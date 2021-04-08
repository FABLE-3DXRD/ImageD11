
from __future__ import print_function
import sys
import h5py, timeit
import numpy as np
from ImageD11 import sparseframe

timer= timeit.default_timer

def newpks(hf, g='measurement/frelon3'):
    """ do a peak search in a sparse frame series """
    blbs = []
    with h5py.File(hf,"r") as h:
        # sort as numbers:
        dsu = [ (float(v), v) for v in list(h)]
        dsu.sort()
        for _, scan in dsu:
            print(scan,end=" ")
            sys.stdout.flush()
            if g is None:
                gscan=h[scan]
            else:
                gscan=h[scan][g]
            for fgrp in sorted(list(gscan)):
                f = sparseframe.from_hdf_group(gscan[fgrp])
                sparseframe.sparse_localmax(f)
                blbs.append( sparseframe.sparse_moments(f, "intensity",
                                                        "localmax") )
    return blbs


def main():
    ################### scan specific motorpos ########
    dtyomepos = {}
    oa = np.linspace(0.,360.,3600)
    pks = {}
    start = timer()
    blobs = newpks( sys.argv[1] )
    end = timer()
    print()
    print("Time",end-start, len(blobs))
    return blobs




if __name__=="__main__":
    blobs = main()
