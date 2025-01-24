
from __future__ import print_function



import numpy as np, pylab as pl

def fitstrain( surf, pkcen, pkwid ):
    """
    surf = 2D surface from pyFAI of intensity versus tth/eta
    pkcen = where the peak is on axis 1
    pkwid = how much roi to use
    """
    lo = max(0, int(pkcen - pkwid/2) )
    hi = min(surf.shape[1], int(pkcen + pkwid/2) )
    roi = surf[:,lo:hi]
    bg = roi.min(axis=1)
    roi = (roi.T-bg).T
    sIx = (roi*np.arange(roi.shape[1])).sum(axis=1)
    sI  = roi.sum(axis=1)
    cen = sIx/( sI + (sI==0))
    if 0:
        print(surf.shape, pkcen, pkwid, lo, hi)
        pl.subplot(221)
        pl.imshow(np.log(roi),aspect='auto',interpolation='nearest')
        pl.colorbar()
        pl.subplot(222)
        pl.title(cen)
        pl.plot(cen)
        pl.subplot(223)
        pl.title("tth-I")
        pl.plot(roi.sum(axis=0))
        pl.subplot(224)
        pl.title("sI")
        pl.plot(sI)
        pl.show()
    return cen, sIx, sI, lo, hi

def fitcen( cen, w, azi):
    """
    cen = centroid positions
    w   = weights (intensity or otherwise)
    azi = angle, degrees
    
    Fits to equation :
      cen = cen0 + A*sin(azi) + B*cos(azi) + C*sin(azi*2) + D*cos(azi*2)
    """
    n = len(cen)
    m = 5 # number of variables in equation to fit
    razi = np.radians( azi )
    g = np.zeros( (m, n) )
    g[0] = np.ones(n)
    g[1] = np.sin(razi)
    g[2] = np.cos(razi)
    g[3] = np.sin(2*razi)
    g[4] = np.cos(2*razi)
    rhs = np.zeros(m)
    mat = np.zeros((m,m))
    for i in range(m):
        rhs[i] = (g[i]*cen*w).sum()
        for j in range(i,m):
            mat[i,j] = mat[j,i] = (g[i]*g[j]*w).sum()
    #print mat
    imat = np.linalg.inv(mat)
    sol = np.dot( imat, rhs )
    calc = np.dot( g.T, sol )
    return sol, calc

    

    
