import sys
import numpy as np
from ImageD11.sinograms.point_by_point import PBP, initializer

def load_config(config_path):
    """Load the scalar configuration and forgen info."""
    config = {}
    with open(config_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            key, val = line.split('=', 1)
            if key == 'forgen':
                val = [int(x) for x in val.split(',')]
            elif key in {'y0', 'hkl_tol', 'ds_tol', 'cosine_tol', 'uniqcut'}:
                val = float(val)
            elif key in {'minpks', 'nprocs', 'hmax'}:
                val = int(val)
            else:
                val = val  # parfile, phase_name, symmetry, icolf_filename
            config[key] = val
    return config

def run_chunk(config_path, indices_path, grains_file):
    config = load_config(config_path)

    from ImageD11.sinograms.dataset import load
    dset = load(config['dset_path'])
    ybincens = np.array(dset.ybincens)
    ystep = dset.ystep
    ymin = ybincens.min()
    
    # # Initialize globals needed for indexing
    # initializer(
    #     parfile=config['parfile'],
    #     phase_name=config['phase_name'],
    #     symmetry=config['symmetry'],
    #     colfile=config['icolf_filename'],
    #     loglevel=3
    # )
    
    pbp_args = {
        'parfile': config['parfile'],
        'dset': dset,
        'hkl_tol': config['hkl_tol'],
        'fpks': config['minpks'],
        'ds_tol': config['ds_tol'],
        'cosine_tol': config['cosine_tol'],
        'y0': config['y0'],
        'symmetry':config['symmetry'],
        'forgen': config['forgen'],
        'uniqcut': config['uniqcut'],
        'phase_name': config['phase_name']
    }
    pbp = PBP(**pbp_args)

    # do the stuff that normally happens at the end of setpeaks():
    pbp.hmax = config['hmax']
    pbp.minpks = config['minpks']
    
    # Load the points for this chunk
    points = np.loadtxt(indices_path, dtype=int)
    points = [tuple(row) for row in points]
    
    
    # Run indexing
    pbp.point_by_point(
        grains_filename=grains_file,
        icolf_filename=config['icolf_filename'],
        nprocs=config['nprocs'],
        debugpoints=points
    )

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python -m ImageD11.nbGui.S3DXRD.run_pbp_recon_chunk <config.txt> <indices.txt> <output.txt>")
        sys.exit(1)

    run_chunk(sys.argv[1], sys.argv[2], sys.argv[3])