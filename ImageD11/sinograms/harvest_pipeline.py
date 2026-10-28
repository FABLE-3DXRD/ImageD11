"""
ImageD11.sinograms.harvest_pipeline
===================================

Raw frames -> a cleaned grain set with per-reflection intensities and grain shapes, from a scanning 3DXRD
dataset and its indexed grains. Every policy value is in a JSON config file (see `DEFAULTS`); each step reads
the outputs of the previous ones from cfg['outdir'].

Steps:
  harvest     box integration of every predicted reflection of every grain (Harvester), grains and radii from
              the dataset's grains file and TensorMap labels -> harvest.h5 (or reuse cfg harvest.file)
  duplicates  near-parallel grains sharing their peaks (Harvester.pair_overlaps, duplicate_groups)
              -> duplicates.json
  pass1       harvest of the cleaned set: kept grains plus one union-box grain per duplicate group
              (union_spots), with 1D s, f, omega profiles -> pass1.h5
  sample      sample cross-section: MLEM of the sinogram of all pass-1 box sums (absorption.sample_shape)
              -> sample.h5
  flags       contamination of every box by other grains' peaks (Harvester.deoverlap) -> flags.h5, and the
              grain table (refine or drop) -> grains.csv
  learn       s, f, omega profiles learned from strong clean reflections (learn_profiles) -> learned.h5
  pass2       profile-weighted harvest (Harvester weights) -> pass2.h5
  grains      per grain: shape by MLEM from its strongest reflections (profile_fit.mlem_shape), its projection
              at every reflection, intensities (profile_fit.intensities), absorption (cfg absorption.mu: the
              traced path through the sample cross-section from the grain centre, absorption.path_lengths) and
              the reflections to leave out (drop) -> grains.h5, shapes.h5

Usage:  python -m ImageD11.sinograms.harvest_pipeline config.json [step ...]     (no step: all)
Every output file carries the config, the ImageD11 version (and git commit for a checkout) and the versions of
python, numpy, scipy and numba.
"""

from __future__ import print_function, division

import os
import sys
import csv
import json
import time
import copy
import subprocess
import multiprocessing

import numpy as np
import h5py

STEPS = ['harvest', 'duplicates', 'pass1', 'sample', 'flags', 'learn', 'pass2', 'grains']
COLS = ('grain', 'h', 'k', 'l', 'ds', 'tth', 'eta', 'omega', 's', 'f', 'sc', 'fc', 's0', 'f0', 'hs', 'hf', 'ho', 'bw')

DEFAULTS = {
    "dataset": None,            # ImageD11.sinograms.dataset file
    "phase": None,              # phase name in the parameter file and the grains file
    "raw_data": None,           # folder of the scan files (default: the dataset's)
    "outdir": None,
    "nproc": 8,
    "harvest": {"file": None,   # reuse this first harvest if it exists
                "ds_max": 1.5, "reader": "sparse", "tensormap": None,     # None: TensorMap_<phase>
                "radius_margin": 3.0,                                     # in dty steps
                "box": {}},                                               # ShoeboxModel arguments
    "duplicates": {"symmetry": "cubic", "max_angle": 2.0, "min_shared": 0.3, "in_box": 0.5, "to_partner": 0.5,
                   "min_snr": 20.0},
    "flags": {"deoverlap": "jacobi", "deoverlap_args": {},   # Harvester.deoverlap method ('jacobi', 'solve'), kw
              "contam_frac": 0.05, "contam_nsig": 2.0,
              "min_clean": 500, "min_clean_frac": 0.5},
    "sample": {"omega_bin": 1.0, "niter": 60, "threshold": 0.1},      # threshold: x 99th percentile density
    "absorption": {"mu": None},  # linear absorption coefficient per dty unit (e.g. 1/um); None: no correction
    "exclude": {"min_coverage": 0.98},   # fraction of the grain's projected shape inside the integrated dty rows
    "learn": {"snr_ref": 20.0, "k": 15, "min_ref": 60},
    "intensity": {"t_prof": 3.0, "rej": 4.0},
    "mlem": {"snr_min": 20.0, "nmax": 600, "niter": 40, "outsig": 5.0, "prop": 0.1, "out_at": [5, 15],
             "min_sin_eta": 0.2, "hkl_divisor": 1},
}


def load_config(path):
    """the config file over DEFAULTS (one level of nesting)"""
    cfg = copy.deepcopy(DEFAULTS)
    with open(path) as f:
        user = json.load(f)
    for k, v in user.items():
        if isinstance(v, dict) and isinstance(cfg.get(k), dict):
            cfg[k].update(v)
        else:
            cfg[k] = v
    for k in ('dataset', 'phase', 'outdir'):
        if cfg[k] is None:
            raise ValueError('config: %s is required' % k)
    if not os.path.isdir(cfg['outdir']):
        os.makedirs(cfg['outdir'])
    return cfg


def provenance(cfg):
    import numpy, scipy, numba, ImageD11
    src = os.path.dirname(os.path.dirname(os.path.abspath(ImageD11.__file__)))
    try:
        rev = subprocess.check_output(['git', '-C', src, 'rev-parse', 'HEAD'], stderr=subprocess.STDOUT)
        # content changes only (file mode bits differ on some file systems)
        dirty = subprocess.check_output(['git', '-C', src, '-c', 'core.fileMode=false', 'status', '--porcelain',
                                         '--untracked-files=no', '--',
                                         'ImageD11'], stderr=subprocess.STDOUT)
        rev, dirty = rev.decode().strip(), bool(dirty.decode().strip())
    except (OSError, subprocess.CalledProcessError):
        rev, dirty = '', False
    return dict(config=json.dumps(cfg), imaged11_version=ImageD11.__version__, imaged11_commit=rev,
                imaged11_modified=dirty, python=sys.version.split()[0], numpy=numpy.__version__,
                scipy=scipy.__version__, numba=numba.__version__, date=time.strftime('%Y-%m-%d %H:%M:%S'))


def stamp(fname, cfg):
    with h5py.File(fname, 'a') as f:
        for k, v in provenance(cfg).items():
            f.attrs['provenance_' + k] = v


def out(cfg, name):
    return os.path.join(cfg['outdir'], name)


def log(*a):
    print(*a)
    sys.stdout.flush()


def dataset_and_pars(cfg):
    import ImageD11.sinograms.dataset, ImageD11.parameters
    cwd = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(cfg['dataset'])))     # file names in the dataset are relative to it
    try:
        ds = ImageD11.sinograms.dataset.load(cfg['dataset'])
        P = ImageD11.parameters.parameters()
        P.loadparameters(ds.parfile, phase_name=cfg['phase'])
    finally:
        os.chdir(cwd)
    return ds, P


def raw_data(cfg, ds):
    return cfg['raw_data'] if cfg['raw_data'] else ds.datapath


def first_harvest(cfg):
    return cfg['harvest']['file'] if cfg['harvest']['file'] else out(cfg, 'harvest.h5')


def box_of(fname):
    """the Box or ShoeboxModel (bg_estimator, p_clip) a harvest file was made with, for re-harvests"""
    from ImageD11.sinograms import harvester
    with h5py.File(fname, 'r') as f:
        a = f.attrs
        return (harvester.Box(bg_estimator=str(a.get('box_bg_estimator', 'clipped')), p_clip=float(a['box_p_clip'])),
                float(a['ds_max']))


def harvester_from_file(cfg, fname, **kw):
    """a Harvester with the grains, radii and spot table of a harvest file (no prediction)"""
    import ImageD11.grain
    from ImageD11.sinograms import harvester
    ds, P = dataset_and_pars(cfg)
    box, ds_max = box_of(fname)
    with h5py.File(fname, 'r') as f:
        grains = [ImageD11.grain.grain(u, translation=t) for u, t in zip(f['ubi'][()], f['translation'][()])]
        h = harvester.Harvester(ds, grains, P, ds_max=ds_max, radius=f['radius'][()], datapath=raw_data(cfg, ds),
                                reader=cfg['harvest']['reader'], box=box, **kw)
        h.spots = dict((k, f['spots'][k][()]) for k in f['spots'])
        h.gid = f['gid'][()]
    return h


def harvest_with(cfg, sp, grains, radius, gid, members, fname, **kw):
    """harvest the reflections sp of grains (one pass over the frames)"""
    from ImageD11.sinograms import harvester
    ds, P = dataset_and_pars(cfg)
    box, ds_max = box_of(first_harvest(cfg))
    h = harvester.Harvester(ds, grains, P, ds_max=ds_max, radius=radius, datapath=raw_data(cfg, ds),
                            reader=cfg['harvest']['reader'], box=box, **kw)
    h.spots = sp
    h._static()
    h._overlaps()
    t0 = time.time()
    h.run(nproc=cfg['nproc'], outfile=fname)
    log('harvest %s: %.0f s' % (fname, time.time() - t0))
    with h5py.File(fname, 'a') as o:
        o.create_dataset('gid', data=np.asarray(gid))
        o.create_dataset('members', data=np.array([' '.join(map(str, m)) for m in members], dtype='S'))
        o.create_dataset('ubi', data=np.array([g.ubi for g in grains]))
        o.create_dataset('translation', data=np.array([g.translation for g in grains]))
        o.create_dataset('radius', data=np.asarray(radius))
        o.attrs['datapath'] = raw_data(cfg, ds)
    stamp(fname, cfg)


def labels_recon(grainsfile, group):
    """grain label image of a TensorMap in the reconstruction frame (TensorMap.map_order_to_recon_order)"""
    with h5py.File(grainsfile, 'r') as f:
        m = f[group]['maps']['labels'][0]
    return np.flip(np.swapaxes(m, 0, 1), 1).astype(int)


# ---------------------------------------------------------------------------------------------- steps
def step_harvest(cfg):
    """all grains of the phase in the grains file; radius = the furthest voxel of the grain's largest connected
    region in the TensorMap labels from its position, plus radius_margin dty steps"""
    import scipy.ndimage as ndi
    import ImageD11.grain
    from ImageD11.sinograms import harvester
    from ImageD11.sinograms.sinogram import read_h5
    from ImageD11.sinograms.geometry import recon_bins
    c = cfg['harvest']
    fname = first_harvest(cfg)
    if os.path.exists(fname):
        log('harvest: using', fname)
        return
    ds, P = dataset_and_pars(cfg)
    gss = read_h5(ds.grainsfile, ds, cfg['phase'])
    labels = labels_recon(ds.grainsfile, c['tensormap'] or 'TensorMap_' + cfg['phase'])
    _, bc = recon_bins(ds.ybincens, ds.ymin, ds.ystep, ds.y0)
    gid, grains, radius = [], [], []
    for k, gs in enumerate(gss):
        cc, n = ndi.label(labels == k)
        if n == 0:
            log('harvest: grain %d has no voxels in the labels, left out' % k)
            continue
        main = ndi.binary_fill_holes(cc == (np.argmax(np.bincount(cc.ravel())[1:]) + 1))
        ii, jj = np.nonzero(main)
        t = np.asarray(gs.grain.translation, float)
        radius.append(np.hypot(bc[ii] - t[0], -bc[jj] - t[1]).max() + c['radius_margin'] * ds.ystep)
        grains.append(ImageD11.grain.grain(gs.grain.ubi, translation=t))
        gid.append(k)
    box = harvester.ShoeboxModel(**c['box'])
    h = harvester.Harvester(ds, grains, P, ds_max=c['ds_max'], radius=np.array(radius), datapath=raw_data(cfg, ds),
                            reader=c['reader'], box=box)
    h.predict()
    log('harvest: %d grains, %d reflections, %s' % (len(grains), len(h.spots['h']), box))
    t0 = time.time()
    h.run(nproc=cfg['nproc'], outfile=fname)
    log('harvest %s: %.0f s' % (fname, time.time() - t0))
    with h5py.File(fname, 'a') as f:
        f.create_dataset('gid', data=np.array(gid))
        f.create_dataset('radius', data=np.array(radius))
        f.create_dataset('ubi', data=np.array([g.ubi for g in grains]))
        f.create_dataset('translation', data=np.array([g.translation for g in grains]))
        f.attrs['datapath'] = raw_data(cfg, ds)
    stamp(fname, cfg)


def step_duplicates(cfg):
    from ImageD11.sinograms import harvester
    from ImageD11 import sym_u
    c = cfg['duplicates']
    fn = first_harvest(cfg)
    h = harvester_from_file(cfg, fn)
    with h5py.File(fn, 'r') as f:
        m = dict((k, f['merged'][k][()]) for k in ('I', 'sigI', 's', 'f', 'omega'))
    t0 = time.time()
    M = h.pair_overlaps()
    groups, pairs = harvester.duplicate_groups(h.spots, m, [g.ubi for g in h.grains],
                                               sym_u.getgroup(c['symmetry'])().group, h.omega_step(),
                                               max_angle=c['max_angle'], min_shared=c['min_shared'],
                                               in_box=c['in_box'], to_partner=c['to_partner'], min_snr=c['min_snr'],
                                               pair_overlaps=M)
    gid = h.gid
    res = dict(groups=[[int(gid[x]) for x in g] for g in groups],
               pairs=[dict(gid_i=int(gid[i]), gid_j=int(gid[j]), misorientation=a, in_box=b, to_partner=t)
                      for i, j, a, b, t in pairs], provenance=provenance(cfg))
    with open(out(cfg, 'duplicates.json'), 'w') as fo:
        json.dump(res, fo, indent=1)
    np.save(out(cfg, 'pair_overlaps.npy'), M)
    log('duplicates: %d groups %s (%.0f s)' % (len(groups), res['groups'], time.time() - t0))


def final_set(cfg):
    """spots, grains, radius, gid, members of the cleaned set: the kept grains' rows of the first harvest, then
    one union-box grain per duplicate group (main grain's orientation, radius^2-weighted position)"""
    import ImageD11.grain
    from ImageD11.sinograms import harvester
    with open(out(cfg, 'duplicates.json')) as f:
        groups = json.load(f)['groups']
    h = harvester_from_file(cfg, first_harvest(cfg))
    gid = [int(x) for x in h.gid]
    k_of = dict((g, k) for k, g in enumerate(gid))
    absorbed = set(g for gr in groups for g in gr)
    keep = [k for k, g in enumerate(gid) if g not in absorbed]
    sp0 = h.spots
    kmask = np.isin(sp0['grain'], keep)
    remap = np.full(len(gid), -1)
    remap[keep] = np.arange(len(keep))
    sp = dict((c, sp0[c][kmask]) for c in COLS)
    sp['grain'] = remap[sp['grain']].astype(np.int32)
    ubis = [g.ubi for g in h.grains]
    grains = [ImageD11.grain.grain(ubis[k], translation=h.grains[k].translation) for k in keep]
    radius = list(h.radius[keep])
    out_gid = [gid[k] for k in keep]
    members = [[gid[k]] for k in keep]
    parts = [sp]
    for n, gr in enumerate(groups):
        ks = [k_of[g] for g in gr]
        msp, _ = harvester.union_spots(sp0, ks, ubis, h.omega_step())
        msp['grain'] = np.full(len(msp['h']), len(keep) + n, np.int32)
        parts.append(dict((c, msp[c]) for c in COLS))
        w = np.array([h.radius[k] ** 2 for k in ks])
        T = (w[:, None] * np.array([h.grains[k].translation for k in ks])).sum(0) / w.sum()
        radius.append(max(np.hypot(*(h.grains[k].translation[:2] - T[:2])) + h.radius[k] for k in ks))
        grains.append(ImageD11.grain.grain(ubis[ks[0]], translation=T))
        out_gid.append(gr[0])
        members.append(gr)
    sp = dict((c, np.concatenate([p[c] for p in parts])) for c in COLS)
    return sp, grains, np.array(radius), out_gid, members


def step_pass1(cfg):
    sp, grains, radius, gid, members = final_set(cfg)
    log('pass1: %d grains, %d reflections' % (len(grains), len(sp['h'])))
    harvest_with(cfg, sp, grains, radius, gid, members, out(cfg, 'pass1.h5'), profiles=True)


def step_sample(cfg):
    """sample cross-section from the sinogram of all pass-1 box sums -> sample.h5 (density, mask, sinogram)"""
    from ImageD11.sinograms import absorption
    from ImageD11.sinograms.geometry import recon_bins
    c = cfg['sample']
    ds, _ = dataset_and_pars(cfg)
    ny = len(ds.ybincens)
    om = np.asarray(ds.omega, float)
    edges = np.arange(np.floor(om.min()), np.ceil(om.max()) + c['omega_bin'], c['omega_bin'])
    with h5py.File(out(cfg, 'pass1.h5'), 'r') as f:
        omega = f['spots/omega'][()]
        ps = f['per_scan']
        sino = np.zeros((len(edges) - 1, ny))
        n, chunk = ps['spot'].shape[0], 1 << 25
        for a in range(0, n, chunk):
            sino += absorption.sinogram(omega, ps['spot'][a:a + chunk], ps['scan'][a:a + chunk],
                                        ps['box'][a:a + chunk] - ps['bg'][a:a + chunk], edges, ny)
    nb = len(recon_bins(ds.ybincens, ds.ymin, ds.ystep, ds.y0)[1])
    rho, mask = absorption.sample_shape(sino, 0.5 * (edges[1:] + edges[:-1]), ds.y0, ds.ymin, ds.ystep, (nb, nb),
                                        niter=c['niter'], threshold=c['threshold'])
    with h5py.File(out(cfg, 'sample.h5'), 'w') as o:
        o.create_dataset('density', data=rho.astype(np.float32), compression='gzip')
        o.create_dataset('mask', data=mask, compression='gzip')
        o.create_dataset('sinogram', data=sino.astype(np.float32), compression='gzip')
        o.create_dataset('omega_edges', data=edges)
    stamp(out(cfg, 'sample.h5'), cfg)
    log('sample: %d pixels (equivalent diameter %.1f dty units)' % (mask.sum(), 2 * np.sqrt(mask.sum() / np.pi) * ds.ystep))


def step_flags(cfg):
    from ImageD11.sinograms import harvester
    c = cfg['flags']
    fn = out(cfg, 'pass1.h5')
    h = harvester_from_file(cfg, fn)
    # the box model of the first harvest gives the expected spot widths of the sources
    with h5py.File(first_harvest(cfg), 'r') as f0, h5py.File(fn, 'r') as f:
        I, sigI = f['merged/I'][()], f['merged/sigI'][()]
        members = [x.split() for x in f['members'][()].astype(str)]
        a = f0.attrs
        gid0 = f0['gid'][()].astype(int)
        k0 = dict((g, k) for k, g in enumerate(gid0))
        main = np.array([k0[int(m[0])] for m in members])
        kw = dict((k[4:], a[k]) for k in a if k.startswith('box_') and k[4:] in
                  ('psf', 'bandpass', 'phi0', 'phi0_t', 'phi', 'eps', 'frame', 'sin_eta_min'))
        for k in ('phi', 'eps'):
            kw[k] = np.broadcast_to(np.asarray(kw.get(k, 0.0), float), (len(gid0),))[main]
        model = harvester.ShoeboxModel(**kw)
    sp = h.spots
    sig = list(model.sigmas(sp, h.pars, h.omega_step()))
    merged_grain = np.array([len(m) > 1 for m in members])[sp['grain']]
    # union boxes: their peaks are spread over the whole box
    for i, hk in enumerate(('hs', 'hf', 'ho')):
        sig[i] = np.where(merged_grain, np.maximum(sig[i], (sp[hk] - 1) / 3.0), sig[i])
    t0 = time.time()
    d = h.deoverlap(I, sigI, sigmas=sig, method=c['deoverlap'], verbose=True, **c['deoverlap_args'])
    log('deoverlap %.0f s' % (time.time() - t0))
    flag = d['contamination'] > np.maximum(c['contam_frac'] * np.maximum(I, 0), c['contam_nsig'] * sigI)
    g = sp['grain']
    with h5py.File(out(cfg, 'flags.h5'), 'w') as o:
        for k, v in dict(d, contaminated=flag).items():
            o.create_dataset(k, data=v, compression='gzip')
        o.attrs['rule'] = "contaminated = de-overlapped contamination > max(%g I, %g sigI)" % (
            c['contam_frac'], c['contam_nsig'])
    stamp(out(cfg, 'flags.h5'), cfg)
    sigf = I > 3 * sigI
    rows = []
    for k in range(len(h.grains)):
        m = (g == k) & sigf
        n_clean = int((m & ~flag).sum())
        drop = n_clean < c['min_clean'] or n_clean < c['min_clean_frac'] * m.sum()
        rows.append(dict(gid=int(h.gid[k]), members=' '.join(members[k]), n_significant=int(m.sum()),
                         n_clean=n_clean, frac_contaminated=float(flag[m].mean()) if m.any() else np.nan,
                         refine=not drop))
    with open(out(cfg, 'grains.csv'), 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=['gid', 'members', 'n_significant', 'n_clean', 'frac_contaminated',
                                           'refine'])
        w.writeheader()
        w.writerows(rows)
    log('flags: %.1f%% of significant reflections contaminated; %d of %d grains dropped' % (
        100 * flag[sigf].mean(), sum(not r['refine'] for r in rows), len(rows)))


def step_learn(cfg):
    from ImageD11.sinograms import harvester
    c = cfg['learn']
    with h5py.File(out(cfg, 'pass1.h5'), 'r') as f, h5py.File(out(cfg, 'flags.h5'), 'r') as fl:
        sp = dict((k, f['spots'][k][()]) for k in ('grain', 'hs', 'hf', 'ho', 's', 'f', 's0', 'f0', 'tth', 'eta',
                                                    'edge', 'touch_mask', 'n_masked'))
        I, sigI = f['merged/I'][()], f['merged/sigI'][()]
        prof, off = f['profiles/data'][()].astype(float), f['profiles/offset'][()]
        bad = fl['contaminated'][()] | (fl['coincident_fraction'][()] >= 0.5)
    ref = (I > c['snr_ref'] * sigI) & ~sp['edge'] & ~sp['touch_mask'] & (sp['n_masked'] == 0) & ~bad
    t0 = time.time()
    w, nref, pooled = harvester.learn_profiles(sp, I, sigI, prof, off, ref, k=c['k'], min_ref=c['min_ref'])
    with h5py.File(out(cfg, 'learned.h5'), 'w') as o:
        o.create_dataset('data', data=w.astype(np.float32), compression='gzip')
        o.create_dataset('offset', data=off)
        o.create_dataset('n_ref', data=nref.astype(np.int8))
        o.create_dataset('pooled', data=pooled)
        o.create_dataset('reference', data=ref)
    stamp(out(cfg, 'learned.h5'), cfg)
    log('learn: %d references, %.0f s' % (ref.sum(), time.time() - t0))


def step_pass2(cfg):
    import ImageD11.grain
    with h5py.File(out(cfg, 'learned.h5'), 'r') as w, h5py.File(out(cfg, 'pass1.h5'), 'r') as f:
        weights = (w['data'][()].astype(float), w['offset'][()])
        sp = dict((c, f['spots'][c][()]) for c in COLS)
        gid, members = f['gid'][()], [x.split() for x in f['members'][()].astype(str)]
        radius = f['radius'][()]
        grains = [ImageD11.grain.grain(u, translation=t) for u, t in zip(f['ubi'][()], f['translation'][()])]
    harvest_with(cfg, sp, grains, radius, gid, members, out(cfg, 'pass2.h5'), weights=weights)


_G = {}


def _grain(k):
    """shape (MLEM) and intensities of grain k"""
    from ImageD11.sinograms import profile_fit
    G = _G
    cfg, ds = G['cfg'], G['ds']
    c, ci = cfg['mlem'], cfg['intensity']
    ny = len(ds.ybincens)
    a, b = G['start'][k], G['start'][k + 1]
    ra, rb = G['rstart'][k], G['rstart'][k + 1]
    sp = dict((x, v[a:b]) for x, v in G['sp'].items())
    nr = b - a
    iscans = np.unique(G['scan'][ra:rb])
    si = np.searchsorted(iscans, G['scan'][ra:rb])
    ri = G['spot'][ra:rb] - a
    A = {}
    for x in G['cols']:
        arr = np.zeros((len(iscans), nr))
        arr[si, ri] = G[x][ra:rb]
        A[x] = arr
    seen = A['nframes'] > 0
    w = dict((x, A[x]) for x in ('wA', 'wB', 'wM', 'wC', 'wVbg'))
    # shape from the strongest reflections (box sums; no shape needed for them)
    r0 = profile_fit.row_estimates(A['box'], A['bg'], A['var'], A['nframes'], sp['nbox'])
    I_sum = np.where(seen, r0['net'], 0).sum(0)
    s_sum = np.sqrt(np.where(seen, r0['v_sum'], 0).sum(0))
    H = np.column_stack((sp['h'], sp['k'], sp['l']))
    # boxes cut by the ends of the omega range
    half = (sp['ho'] + 0.5) * G['ostep']
    omega_edge = (sp['omega'] - half < G['omega_range'][0]) | (sp['omega'] + half > G['omega_range'][1])
    use = np.all(H % int(c['hkl_divisor']) == 0, axis=1) & ~G['drop'][a:b] & ~sp['edge'] & ~sp['touch_mask'] & \
        ~omega_edge & (np.abs(np.sin(np.radians(sp['eta']))) >= c['min_sin_eta']) & (I_sum > c['snr_min'] * s_sum)
    idx = np.flatnonzero(use)
    idx = idx[np.argsort(-(I_sum[idx] / s_sum[idx]))][:c['nmax']]
    shape2d = (len(G['bc']), len(G['bc']))
    support = profile_fit.support_disc(G['T'][k], G['radius'][k], shape2d, ds.ystep)
    box, bgs, sn = np.zeros((len(idx), ny)), np.zeros((len(idx), ny)), np.zeros((len(idx), ny))
    box[:, iscans] = A['box'][:, idx].T
    bgs[:, iscans] = A['bg'][:, idx].T
    sn[:, iscans] = seen[:, idx].T
    if len(idx) >= 10:
        rho, nout = profile_fit.mlem_shape(box, bgs, sn > 0, I_sum[idx], sp['omega'][idx], support, ds.y0, ds.ymin,
                                           ds.ystep, niter=c['niter'], outsig=c['outsig'], prop=c['prop'],
                                           out_at=tuple(c['out_at']))
    else:
        rho, nout = support.astype(float), 0       # too few strong reflections: the support disc
    Pall = profile_fit.project(rho, sp['omega'], ds.y0, ds.ymin, ds.ystep, ny)
    P = Pall[:, iscans].T
    # fraction of the grain's projection on the rows that were integrated (the rest is outside the dty window)
    ptot = Pall.sum(1)
    coverage = np.where(ptot > 0, (P * seen).sum(0) / np.where(ptot > 0, ptot, 1), 0.0)
    res = profile_fit.intensities(A['box'], A['bg'], A['var'], A['nframes'], sp['nbox'], P, w=w,
                                  t_prof=ci['t_prof'], rej=ci['rej'])
    net_used = np.where(res['profile_used'][None, :], res['model'], r0['net'] * seen)
    r1 = profile_fit.row_estimates(A['box'], A['bg'], A['var'], A['nframes'], sp['nbox'], w, res['model'])
    if G['mask'] is not None:
        from ImageD11.sinograms import absorption
        L = absorption.path_lengths(G['mask'], ds.ystep, G['T'][k][0], G['T'][k][1], sp['omega'], sp['tth'],
                                    sp['eta'])
        trans = np.exp(-cfg['absorption']['mu'] * L)
    else:
        trans = np.ones(nr)
    overlap = G['drop'][a:b]
    drop = overlap | sp['edge'] | sp['touch_mask'] | omega_edge | (coverage < cfg['exclude']['min_coverage'])
    grp = dict(iscan=iscans, dty=ds.ybincens[iscans], ubi=G['ubi'][k], translation=G['T'][k], refl_hkl=H,
               refl_omega=sp['omega'], refl_s=sp['s'], refl_f=sp['f'], refl_tth=sp['tth'], refl_eta=sp['eta'],
               refl_touch_mask=sp['touch_mask'], refl_edge=sp['edge'], refl_n_overlap=sp['n_overlap'],
               refl_hs=sp['hs'], refl_ho=sp['ho'], box=A['box'].astype(np.float32), bg=A['bg'].astype(np.float32),
               var=A['var'].astype(np.float32), nframes=A['nframes'].astype(np.float32),
               net=(net_used / trans[None, :]).astype(np.float32), I_final=res['I'] / trans,
               sigI_final=res['sigI'] / trans, transmission=trans, coverage=coverage, omega_edge=omega_edge,
               overlap=overlap, I_sum=res['I_sum'],
               sigI_sum=res['sigI_sum'], I_prof=res['I_prof'], sigI_prof=res['sigI_prof'],
               profile_used=res['profile_used'], chi2=res['chi2'], wI=r1['prof'].astype(np.float32),
               wvar=r1['v_prof'].astype(np.float32), drop=drop, members=np.array(G['members'][k], dtype='S'))
    return int(G['gid'][k]), grp, rho.astype(np.float32), len(idx), nout


def step_grains(cfg):
    from ImageD11.sinograms.geometry import recon_bins
    ds, _ = dataset_and_pars(cfg)
    t0 = time.time()
    with h5py.File(out(cfg, 'pass2.h5'), 'r') as f, h5py.File(out(cfg, 'flags.h5'), 'r') as fl:
        ps = f['per_scan']
        spot = ps['spot'][()]
        o = np.argsort(spot, kind='stable')
        _G['spot'], _G['scan'] = spot[o], ps['scan'][()][o]
        _G['cols'] = ('box', 'bg', 'var', 'nframes', 'wA', 'wB', 'wM', 'wC', 'wVbg')
        for x in _G['cols']:
            _G[x] = ps[x][()][o]
        _G['sp'] = dict((x, f['spots'][x][()]) for x in ('grain', 'h', 'k', 'l', 'omega', 's', 'f', 'tth', 'eta',
                                                            'touch_mask', 'edge', 'n_overlap', 'hs', 'ho', 'nbox'))
        _G['gid'], _G['members'] = f['gid'][()], [x.split() for x in f['members'][()].astype(str)]
        _G['ubi'], _G['T'], _G['radius'] = f['ubi'][()], f['translation'][()], f['radius'][()]
        _G['drop'] = fl['contaminated'][()]
    _G['start'] = np.searchsorted(_G['sp']['grain'], np.arange(len(_G['gid']) + 1))
    _G['rstart'] = np.searchsorted(_G['spot'], _G['start'])
    _G['cfg'], _G['ds'] = cfg, ds
    om = np.asarray(ds.omega, float)
    _G['omega_range'] = (om.min(), om.max())
    _G['ostep'] = float(np.median(np.abs(np.diff(om, axis=-1))))
    _G['mask'] = None
    if cfg['absorption']['mu'] is not None:
        with h5py.File(out(cfg, 'sample.h5'), 'r') as f:
            _G['mask'] = f['mask'][()]
    _G['bc'] = recon_bins(ds.ybincens, ds.ymin, ds.ystep, ds.y0)[1]
    log('grains: loaded in %.0f s' % (time.time() - t0))
    # fork: the workers share _G (the default start method may be spawn; python 2 has no get_context and forks)
    mp = multiprocessing.get_context('fork') if hasattr(multiprocessing, 'get_context') else multiprocessing
    pool = mp.Pool(cfg['nproc'])
    try:
        with h5py.File(out(cfg, 'grains.h5'), 'w') as fo, h5py.File(out(cfg, 'shapes.h5'), 'w') as fs:
            for n, (gid, grp, rho, nrows, nout) in enumerate(pool.imap_unordered(_grain, range(len(_G['gid'])))):
                g = fo.create_group(str(gid))
                for x, v in grp.items():
                    g.create_dataset(x, data=v, compression='gzip' if np.ndim(v) == 2 else None)
                gs = fs.create_group(str(gid))
                gs.create_dataset('recon', data=rho, compression='gzip')
                gs.attrs['nrows'], gs.attrs['n_outliers'] = nrows, nout
                if n % 50 == 0:
                    log('%d/%d grain %d: %d MLEM rows (%.0f s)' % (n + 1, len(_G['gid']), gid, nrows,
                                                                    time.time() - t0))
    finally:
        pool.close()
        pool.join()
    stamp(out(cfg, 'grains.h5'), cfg)
    stamp(out(cfg, 'shapes.h5'), cfg)


def run(cfg, steps=None):
    for s in (steps or STEPS):
        t0 = time.time()
        log('=== %s' % s)
        globals()['step_' + s](cfg)
        log('=== %s done in %.0f s' % (s, time.time() - t0))


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or not args[0].endswith('.json'):
        print(__doc__)
        return
    cfg = load_config(args.pop(0))
    bad = [s for s in args if s not in STEPS]
    if bad:
        raise ValueError('unknown steps %s; steps are %s' % (bad, STEPS))
    run(cfg, args)


if __name__ == '__main__':
    main()
