
import os, sys


from ImageD11 import fImageD11, cImageD11

import numpy as np, time

np.random.seed(42)
NPK = 1000000
testvecs = np.random.random( NPK*3 ).reshape((NPK,3))
ubis = (np.random.random(9*20).reshape((20,3,3))-0.5)*10

gvf= np.array( testvecs.T,  order='F')

tim = time.time

for tol in [0.01, 0.1, 0.2, 0.5]:
    int_tmpf = np.zeros( NPK, np.int32 )
    drlv2f = np.zeros( NPK, float)

    start = tim()
    for i, gr in enumerate(ubis):
        fImageD11.assign( gr, gvf , tol, drlv2f, int_tmpf, i )
    end = tim() 
    fortrantime = end -start

    int_tmpc = np.zeros( NPK, np.int32 )
    drlv2c = np.zeros( NPK, float)

    start = tim()


    for i, gr in enumerate(ubis):
        cImageD11.score_and_assign( gr, testvecs , tol, drlv2c, int_tmpc, i )
    end = tim()
    ctime = end - start
    print "tol","fortrantime",fortrantime,"ctime",ctime

    assert  np.allclose( drlv2c, drlv2f )
    assert  (int_tmpc == int_tmpf).all()
