
from __future__ import print_function

import ImageD11.sparseframe
import scipy.sparse
import numpy as np
from timeit import default_timer

def maketestim(N=64, M=2048, SIG=3):
    # tile co-ordinates
    i,j = np.mgrid[0:N,0:N]-N/2.
    def pk(x,y,s):
        # 2D peak shape (with tails)
        dx = i-x
        dy = j-y
        gg= np.exp( -(dx*dx+dy*dy)/s/s )
        ll=s*s/(s*s+dx*dx+dy*dy)
        return ll+gg
    # Image for testing
    im = np.zeros( (M,M), np.float32 )
    # Make a grid of points (tiles)
    for ii in range(0,M,N):
        for jj in range(0,M,N):
            # offset of centroids
            oi = 3*SIG*ii/M
            oj = SIG*ii/M
            # Ratio of heights
            sc = jj*1.0/M
            # put in two peaks, offset increases one direction,
            # scale in the other
            im[ii:ii+N,jj:jj+N] = pk( -oi, -oj, SIG ) + pk( oi, oj, SIG)*sc
    return im

class ticker:
    def __init__(self):
        self.t = default_timer()
    def __call__(self, msg=''):
        now = default_timer()
        print( msg, '%.3f ms'%((now-self.t)*1e3), end=' ')
        self.t = default_timer()

def runthrough( im1, t1, im2, t2):
    t = ticker()
    #    s1 = ImageD11.sparseframe.from_data_mask( im1 > t1, im1.copy(), {} ); t('s1')
    s1 = ImageD11.sparseframe.from_data_cut( im1, t1 ); t('s1c')
    #assert s1 == s1c
    #s2 = ImageD11.sparseframe.from_data_mask( im2 > t2, im2.copy(), {} ); t('s2')
    s2 = ImageD11.sparseframe.from_data_cut( im2, t2 ); t('s2c')
    n1 = ImageD11.sparseframe.sparse_localmax( s1 ); t('n1')
    n2 = ImageD11.sparseframe.sparse_localmax( s2 ); t('n2')
    coo = ImageD11.sparseframe.overlaps(s1, 'localmax', s2, 'localmax'); t('coo')
    olin = ImageD11.sparseframe.overlaps_linear(max(s1.nnz, s2.nnz)+1); 
    novlin, olinlap = olin( s1.row, s1.col, 
                           s1.pixels['localmax'], 
                           s1.meta['localmax']['nlabel'],
                            s2.row, s2.col, 
                           s2.pixels['localmax'], 
                           s2.meta['localmax']['nlabel'] ); t('olinlap')
    assert olinlap.shape[0] == coo.nnz
    assert (coo.row == olinlap[:,0]-1).all()
    assert (coo.col == olinlap[:,1]-1).all()
    assert (coo.data == olinlap[:,2]).all()
    t.t = default_timer()
    omat = ImageD11.sparseframe.overlaps_matrix(max(n1,n2)+1); 
    novmat, olinlap = omat( s1.row, s1.col, 
                           s1.pixels['localmax'], 
                           s1.meta['localmax']['nlabel'],
                            s2.row, s2.col, 
                           s2.pixels['localmax'],
                           s2.meta['localmax']['nlabel'] ); t('omatlap')
    assert olinlap.shape[0] == coo.nnz
    assert (coo.row == olinlap[:,0]-1).all()
    assert (coo.col == olinlap[:,1]-1).all()
    assert (coo.data == olinlap[:,2]).all()
    print('\n nnz1=%d   nnz2=%d  pk1=%d  pk2=%d  nolap=%d'%(s1.nnz, s2.nnz, n1, n2, coo.nnz))

for N in (256, 1024, 2048):
    M = 64
    while M < N:
        im1 = maketestim( M, N, 3)
        im2 = maketestim( M, N, 2.5)
        print(N, M)
        runthrough( im1, 0.1, im2, 0.25)
        runthrough( im1, 0.1, im2, 0.25)
        M  *= 2



