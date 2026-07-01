"""
match_friedel_pairs.py
======================
Standalone pipeline script for Friedel pair matching in scanning_3DXRD datasets.
 
Usage — run locally:
    python match_friedel_pairs.py -dsfile /path/to/dataset.h5 \
                                  [-parfile /path/to/params.json] \
                                  [-pairing_options /path/to/options.json] \
                                  [-use2Dpeaks True]
 
Usage — submit to slurm:
    python match_friedel_pairs.py -dsfile /path/to/dataset.h5 \
                                  [-parfile /path/to/params.json] \
                                  [-pairing_options /path/to/options.json] \
                                  [-use2Dpeaks True]
                                  -usecluster True
 
Pairing options JSON keys (all optional — defaults shown):
    # tolerances
    tol_gv             : 0.05
    tol_eta            : 0.2
    tol_logI           : null          (null -> np.inf)
    weights            : {"gx":1,"gy":1,"gz":1,"eta":1,"I":1}
    # pairing strategy
    pair_type          : "omega"       ("omega" | "eta" | "all")
    filter_mode        : "relaxed"
    n_eta_bins         : 360
    drop_unpaired      : false
    extended_bin_search: false
    n_workers          : -1            -1 -> all availables
    timeout            : 120
    n_steps            : 25
    # slurm / cluster options
    slurm_partition    : "nice"
    slurm_mem_G        : 64
    slurm_time         : "02:00:00"
    slurm_cpus         : 16
"""

import os
import sys
import site
import json
import logging
import argparse
import contextlib
import time
import datetime
import subprocess
 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# Python path: needs to be sorted out. I copied the stuff we put at the beginning of notebooks to load ImageD11 from local user folder (cloned from github), but in production this should not be here
# python environment stuff
IMAGED11_PATH = 'ImageD11'  # None means do not use git, otherwise enter the name of the folder to use for the git checkout "ImageD11" or "ImageD11_version_xx", etc
CHECKOUT_PATH = 'ImageD11'  # the name of the git checkout folder within path. None means guess

if IMAGED11_PATH is not None:
    if '/data/id11/nanoscope' not in sys.path:
        sys.path.append('/data/id11/nanoscope')
    import install_ImageD11_from_git
    PYTHONPATH = install_ImageD11_from_git.setup_ImageD11_from_git(IMAGED11_PATH,CHECKOUT_PATH)

import ImageD11.sinograms.dataset
import ImageD11.columnfile
import ImageD11.friedel_pairs as fp

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

 
# =============================================================================
#  Options
# =============================================================================
 
class Options:
    """
    Pairing options container with JSON save/load support.
    All attributes map directly to FriedelPairIndexer constructor arguments
    or match_friedel_pairs / match_friedel_pairs_by_chunks keyword arguments.
    """
 
    def __init__(self):
        # ── tolerances ────────────────────────────────────────────────
        self.tol_gv              = 0.05
        self.tol_eta             = 0.2
        self.tol_logI            = None    # None -> np.inf at runtime
        self.weights             = {'gx': 1., 'gy': 1., 'gz': 1.,
                                    'eta': 1., 'I': 1.}
        # ── pairing strategy ──────────────────────────────────────────
        self.pair_type           = 'omega' # 'omega' | 'eta' | 'all'
        self.filter_mode         = 'relaxed'
        self.n_eta_bins          = 360
        self.drop_unpaired       = False
        self.extended_bin_search = False
        self.n_workers           = -1
        self.timeout             = 120
        self.n_steps             = 25
        # ── slurm / cluster ───────────────────────────────────────────
        self.slurm_partition     = 'nice'
        self.slurm_mem_G         = 64
        self.slurm_time          = '02:00:00'
        self.slurm_cpus          = 16      # should match n_workers
        # y0 for reconstruction
        self.y0                 = None

    @property
    def tol_logI_value(self):
        """Return np.inf when tol_logI is None (JSON cannot store inf)."""
        return np.inf if self.tol_logI is None else float(self.tol_logI)
 
    def to_dict(self):
        return {k: v for k, v in self.__dict__.items()}
 
    def save(self, path='friedel_pairs_opts.json'):
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=4)
        print('[SAVED] Options saved to {}'.format(path))
 
    @classmethod
    def load(cls, path):
        with open(path, 'r') as f:
            data = json.load(f)
        obj = cls()
        for k, v in data.items():
            setattr(obj, k, v)
        return obj
 
    def __repr__(self):
        params = ', '.join('{}={}'.format(k, v) for k, v in self.__dict__.items())
        return '<Options {}>'.format(params)
 
    def __str__(self):
        lines = ['Pairing options:']
        for k, v in self.__dict__.items():
            lines.append('  {:25s} {}'.format(k, v))
        return '\n'.join(lines)
 
 
# =============================================================================
#  Logger context manager — redirects to file
# =============================================================================
 
@contextlib.contextmanager
def log_to_file(log_path, logger_name='ImageD11'):
    """
    Context manager: redirects logging to a file and the terminal.
    Both at INFO level. Noisy third-party loggers capped at WARNING.
    """
    root     = logging.getLogger()

    # silence noisy third-party loggers
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)
 
    # file handler — captures everything at DEBUG level
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        '%(asctime)s  %(levelname)-8s  %(name)s  %(message)s',
        datefmt='%H:%M:%S'))
 
    # console handler — WARNING and above only
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
 
    root.addHandler(fh)
    root.addHandler(ch)
    root.setLevel(logging.INFO)
 
    try:
        yield logging.getLogger(logger_name)
    finally:
        root.removeHandler(fh)
        root.removeHandler(ch)
        fh.close()
        ch.close()
 
 
# =============================================================================
#  Path helpers
# =============================================================================
 
def _basedir(dsfile):
    """Directory containing the dataset file."""
    return os.path.dirname(os.path.abspath(dsfile))
 
 
def _dsname(dsfile):
    """Dataset base name, stripped of '_dataset.h5' / '.h5' suffixes."""
    base = os.path.basename(dsfile)
    for suffix in ('_dataset.h5', '.h5'):
        if base.endswith(suffix):
            base = base[:-len(suffix)]
            break
    return base
 
 
def _slurm_dir(dsfile):
    """Path to the slurm output folder next to the dataset."""
    d = os.path.join(_basedir(dsfile), 'slurm_fpairs')
    os.makedirs(d, exist_ok=True)
    return d
 
 
def _figure_path(dsfile, suffix):
    return os.path.join(_basedir(dsfile),
                        '{}_{}.svg'.format(_dsname(dsfile), suffix))
 
 
def _log_path(dsfile):
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    return os.path.join(_basedir(dsfile),
                        'match_friedel_pairs_{}.log'.format(ts))
 
def _save_figure(fig, path, logger):
    try:
        fig.savefig(path, format='svg', bbox_inches='tight')
        logger.info('Figure saved: %s', path)
    except Exception as e:
        logger.warning('Could not save figure %s: %s', path, e)
    finally:
        plt.close(fig)


# =============================================================================
#  Slurm submission
# =============================================================================
 
def prepare_bash_script(dsfile, parfile, pairing_options_file,
                         use2Dpeaks, opts):
    """ Write a slurm batch script to the slurm/ subfolder """
    
    sdir        = _slurm_dir(dsfile)
    script_path = os.path.join(sdir, '{}_match_friedel_pairs.sh'.format(_dsname(dsfile)))
 
    # build the python argument string — omit optional args when None / empty
    this_script = os.path.abspath(__file__)
    py_args     = [
        '-dsfile        {}'.format(dsfile),
        '-use2Dpeaks    {}'.format(use2Dpeaks),
    ]
    if parfile:
        py_args.append('-parfile        {}'.format(parfile))
    if pairing_options_file:
        py_args.append('-pairing_options {}'.format(pairing_options_file))
 
    py_args_str = ' \\\n    '.join(py_args)
 
    script = (
        '#!/bin/bash\n'
        '#SBATCH --job-name=friedel_{dsname}\n'
        '#SBATCH --nodes=1\n'
        '#SBATCH --ntasks=1\n'
        '#SBATCH --cpus-per-task={cpus}\n'
        '#SBATCH --mem={mem}G\n'
        '#SBATCH --time={time}\n'
        '#SBATCH --partition={partition}\n'
        '#SBATCH --output={sdir}/{dsname}_%j.out\n'
        '#SBATCH --error={sdir}/{dsname}_%j.err\n'
        '#SBATCH --export=ALL,IS_SLURM_PIPELINE=1'
        '\n'
        '# ── job info ─────────────────────────────────────────────────\n'
        'echo "------------------------------------------------------------"\n'
        'echo "Job ID    : $SLURM_JOB_ID"\n'
        'echo "Node      : $SLURM_NODELIST"\n'
        'echo "CPUs      : $SLURM_CPUS_PER_TASK"\n'
        'echo "Started   : $(date)"\n'
        'echo "dsfile    : {dsfile}"\n'
        'echo "------------------------------------------------------------"\n'
        '\n'
        '# ── run pipeline ─────────────────────────────────────────────\n'
        'python {script} \\\n'
        '    {py_args}\n'
        '\n'
        'echo "------------------------------------------------------------"\n'
        'echo "Finished  : $(date)"\n'
        'echo "------------------------------------------------------------"\n'
    ).format(
        dsname    = _dsname(dsfile),
        cpus      = opts.slurm_cpus,
        mem       = opts.slurm_mem_G,
        time      = opts.slurm_time,
        partition = opts.slurm_partition,
        sdir      = sdir,
        dsfile    = dsfile,
        script    = this_script,
        py_args   = py_args_str,
    )
 
    with open(script_path, 'w') as f:
        f.write(script)
 
    os.chmod(script_path, 0o755)
    return script_path
 
 
def submit_to_slurm(dsfile, parfile, pairing_options_file,
                     use2Dpeaks, opts):
    """
    Prepare the batch script and submit it with sbatch.
    Returns job_id : str
    """
    script_path = prepare_bash_script(
        dsfile, parfile, pairing_options_file, use2Dpeaks, opts)
 
    print('Batch script written : {}'.format(script_path))
 
    result = subprocess.run(
        ['sbatch', '--parsable', script_path],
        capture_output=True, text=True)
 
    if result.returncode != 0:
        raise RuntimeError(
            'sbatch failed:\n{}'.format(result.stderr.strip()))
 
    job_id  = result.stdout.strip()
    sdir    = _slurm_dir(dsfile)
    dsname  = _dsname(dsfile)
 
    print('Submitted job        : {}'.format(job_id))
    print('Monitor queue        : squeue -j {}'.format(job_id))
    print('Live stdout          : tail -f {}/{}_{}.out'.format(sdir, dsname, job_id))
    print('Live stderr          : tail -f {}/{}_{}.err'.format(sdir, dsname, job_id))
    print('Cancel               : scancel {}'.format(job_id))
 
    return job_id


# =============================================================================
#  Pairing pipeline
# =============================================================================
 
def _run_pair_type(pair_type, FPIndexer, opts, cf, ds,
                   use_chunks, dsfile, logger):
    """
    Run the full pairing sequence for a single pair_type ('omega' or 'eta').
    """
    logger.info('=' * 60)
    logger.info('PAIRING  [%s]', pair_type)
    logger.info('=' * 60)
 
    chunk_type_map = {'omega': 'scans', 'eta': 'eta_bins'}
 
    if use_chunks:
        # ── chunk-based pairing ───────────────────────────────────────
        chunk_type = chunk_type_map.get(pair_type)
 
        FPIndexer.set_peak_subsets(n_eta_bins=opts.n_eta_bins, y0 = opts.y0)
        FPIndexer.sort_peak_subsets(pair_type=pair_type)
        # symmetry check figure
        try:
            fig_sym, _ = FPIndexer.PeakSubsets.check_symmetry()
            _save_figure(fig_sym,
                         _figure_path(dsfile,
                                      '{}_pair_symmetry'.format(pair_type)),
                         logger)
        except Exception as e:
            logger.warning('check_symmetry failed: %s', e)
 
        cf_paired = FPIndexer.match_friedel_pairs_by_chunks(
            chunk_type=chunk_type,
            extended_bin_search=opts.extended_bin_search,
            filter_mode=opts.filter_mode,
            drop_unpaired=opts.drop_unpaired,
            reset_psub=False,
            doplot=False,
            n_workers=opts.n_workers,
            timeout=opts.timeout)
 
    else:
        # ── global pairing ────────────────────────────────────────────
        cf_paired = FPIndexer.match_friedel_pairs(
            pair_type=pair_type,
            drop_unpaired=opts.drop_unpaired,
            filter_mode=opts.filter_mode,
            doplot=False)
 
    # ── pair distance plot ────────────────────────────────────────────
    try:
        fig_dist, _ = fp.plot_pair_distances(cf=cf_paired, pair_type=pair_type)
        _save_figure(fig_dist,
                     _figure_path(dsfile,
                                  '{}_pair_distances'.format(pair_type)),
                     logger)
    except Exception as e:
        logger.warning('plot_pair_distances failed: %s', e)
 
    # ── sample reconstruction plot ────────────────────────────────────
    try:
        i1, i2 = fp.get_pairs(cf_paired, pair_type)
        if pair_type == 'omega':
            sx, sy = fp.locate_omega_pairs(cf_paired, (i1, i2), ds=ds, y0=ds.y0)
        else:
            sx, sy = fp.locate_eta_pairs(cf_paired, (i1, i2), ds=ds, y0=ds.y0)
 
        valid = np.isfinite(sx) & np.isfinite(sy)
        if valid.sum() > 0 and hasattr(ds, 'ybinedges'):
            hist = np.histogram2d(sx[valid] + ds.y0,
                                  sy[valid] + ds.y0,
                                  bins=ds.ybinedges)[0]
            fig_rec, ax = plt.subplots(1, 1, layout='constrained',
                                       figsize=(6, 6))
            im = ax.pcolormesh(ds.ybinedges, ds.ybinedges, hist,
                               vmax=np.percentile(hist.ravel(), 99), rasterized=True)
            ax.set_aspect(1)
            ax.set(title='{} pairs — sample reconstruction'.format(pair_type),
                   ylabel='Sample Y axis',
                   xlabel='Sample X axis')
            plt.colorbar(im, ax=ax, orientation='vertical',
                         pad=0.04, shrink=0.7)
            _save_figure(fig_rec,
                         _figure_path(dsfile,
                                      '{}_pairs_recon'.format(pair_type)),
                         logger)
    except Exception as e:
        logger.warning('Sample reconstruction failed: %s', e)
 
    return cf_paired
 
# =============================================================================
#  Main pipeline
# =============================================================================
def match_friedel_pairs_pipeline(dsfile,
                                  parfile=None,
                                  pairing_options_file=None,
                                  use2Dpeaks=True):
    """
    Full Friedel pair matching pipeline.
 
    Parameters
    ----------
    dsfile               : str  — path to ImageD11 dataset .h5 file
    parfile              : str  — path to parameters file (optional,
                                  falls back to ds.parfile)
    pairing_options_file : str  — path to JSON options file (optional)
    use2Dpeaks           : bool — if True use 2D peaks, else 4D peaks
    """
    t0 = time.perf_counter()
 
    # ── log file next to dataset ────────────────────────────────────── 
    log_path = _log_path(dsfile)
    
    with log_to_file(log_path) as logger:
        logger.info('match_friedel_pairs_pipeline started')
        logger.info('dsfile      : %s', dsfile)
        logger.info('parfile     : %s', parfile)
        logger.info('options     : %s', pairing_options_file)
        logger.info('use2Dpeaks  : %s', use2Dpeaks)
 
        # ── load options ──────────────────────────────────────────────
        if pairing_options_file is not None:
            opts = Options.load(pairing_options_file)
            logger.info('Loaded pairing options from %s', pairing_options_file)
        else:
            opts = Options()
            logger.info('Using default pairing options')
        logger.info('%s', opts)
 
        # ── load dataset ──────────────────────────────────────────────
        logger.info('Loading dataset...')
        ds = ImageD11.sinograms.dataset.load(dsfile)
 
        # ── load columnfile ───────────────────────────────────────────
        logger.info('Loading columnfile (use2Dpeaks=%s)...', use2Dpeaks)
        if use2Dpeaks:
            cf = ds.get_cf_2d()
        else:
            cf = ds.get_cf_4d()
 
        _parfile = parfile if parfile is not None else ds.parfile
        cf.parameters.loadparameters(_parfile)
        cf.updateGeometry()
        logger.info('Columnfile loaded: %d peaks', cf.nrows)
 
        # ── decide strategy: chunk-based or global ────────────────────
        use_chunks = 'dty' in cf.titles
        if use_chunks:
            logger.info('dty column found — using chunk-based pairing')
        else:
            logger.info('No dty column — using global pairing')
 
        # ── initialise FriedelPairIndexer ─────────────────────────────
        FPIndexer = fp.FriedelPairIndexer(
            cf, ds,
            tol_gv   = opts.tol_gv,
            tol_eta  = opts.tol_eta,
            tol_logI = opts.tol_logI_value,
            weights  = opts.weights,
            n_steps  = opts.n_steps)
 
        # ── run pairing for requested pair_type(s) ────────────────────
        pair_types = ['omega', 'eta'] if opts.pair_type == 'all' \
                     else [opts.pair_type]
 
        cf_paired = cf
        for pt in pair_types:
            cf_paired = _run_pair_type(
                pt, FPIndexer, opts, cf_paired, ds,
                use_chunks, dsfile, logger)
 
        # ── match quadruplets if both pair types were run ─────────────
        if opts.pair_type == 'all':
            logger.info('Matching quadruplets...')
            try:
                quads = FPIndexer.match_quadruplets()
                logger.info('Quadruplets found: %d', len(quads))
            except Exception as e:
                logger.warning('match_quadruplets failed: %s', e)
 
        # ── save columnfile ───────────────────────────────────────────
        logger.info('Saving columnfile...')
        try:
            if use2Dpeaks:
                colfile = ds.col2dfile.replace('peaks_2d','peaks_2d_paired')
            else:
                colfile = ds.col4dfile.replace('peaks_4d','peaks_4d_paired')
            
            if os.path.exists(colfile):
                subprocess.run(['rm', colfile], check=True)
            ImageD11.columnfile.colfile_to_hdf(cf_paired, colfile, name='peaks')
            logger.info('Saved to %s', colfile)
                
        except Exception as e:
            logger.error('Failed to save columnfile: %s', e)
            raise
 
        elapsed = time.perf_counter() - t0
        logger.info('Pipeline completed in %.1f s', elapsed)
        print('Done in {:.1f} s — log: {}'.format(elapsed, log_path))

    # ── remove log file if running inside a slurm job ───────────────── 
    #  the .out file in slurm/ is already a full copy of stdout+stderr
    if os.environ.get('IS_SLURM_PIPELINE') == '1':
        try:
            os.remove(log_path)
        except OSError:
            pass

# =============================================================================
#  CLI entry point
# =============================================================================
 
def _parse_bool(v):
    """Argparse helper: accept 'True'/'False'/'1'/'0' as bool."""
    if isinstance(v, bool):
        return v
    if v.lower() in ('true', '1', 'yes'):
        return True
    if v.lower() in ('false', '0', 'no'):
        return False
    raise argparse.ArgumentTypeError('Boolean value expected, got: {}'.format(v))
 
 
def main():
    parser = argparse.ArgumentParser(
        description='Match Friedel pairs in an ImageD11 dataset.')
    parser.add_argument('-dsfile',
        required=True,
        help='Path to dataset .h5 file')
    parser.add_argument('-parfile',
        default=None,
        help='Path to parameters file (optional)')
    parser.add_argument('-pairing_options',
        default=None,
        help='Path to JSON pairing options file (optional)')
    parser.add_argument('-use2Dpeaks',
        default=True, type=_parse_bool,
        help='Use 2D peaks (True) or 4D peaks (False). Default: True')
    parser.add_argument('-usecluster',
        default=False, type=_parse_bool,
        help='If True, write a slurm script and submit via sbatch. Default: False')
 
    args = parser.parse_args()
 
    # load options now so slurm parameters are available for both paths
    if args.pairing_options is not None:
        opts = Options.load(args.pairing_options)
    else:
        opts = Options()
 
    if args.usecluster:
        # ── slurm path: write script + submit, then exit ──────────────
        submit_to_slurm(
            dsfile               = args.dsfile,
            parfile              = args.parfile,
            pairing_options_file = args.pairing_options,
            use2Dpeaks           = args.use2Dpeaks,
            opts                 = opts)
    else:
        # ── local path: run pipeline in this process ──────────────────
        match_friedel_pairs_pipeline(
            dsfile               = args.dsfile,
            parfile              = args.parfile,
            pairing_options_file = args.pairing_options,
            use2Dpeaks           = args.use2Dpeaks)
 
 
if __name__ == '__main__':
    main()
 

