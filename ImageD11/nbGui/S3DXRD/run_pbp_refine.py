import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import sys
from ImageD11.sinograms.point_by_point import PBPRefine

if __name__ == "__main__":
    pbprefinefile = sys.argv[1]
    print('Loading refinement object from disk')
    refine = PBPRefine.from_h5(pbprefinefile)
    print('Running refinement')
    refine.run_refine(use_cluster=False)
    print('Refinement complete!')
