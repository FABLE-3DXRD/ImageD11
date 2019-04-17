
from __future__ import print_function, division

import time
import h5py, scipy.sparse, numpy as np, pylab as pl
from ImageD11 import cImageD11

# see also sandbox/collect_peak_pixels.py


timer = time.time

tic = timer()

def toc(msg):
    global tic
    now = timer()
    print(msg, "%.3f s"%( now - tic))
    tic = now

def read_csr_from_hdf( fname, pname ):
        """
        Reads a sparse matrix from a hdf file(fname)/group(pname)
        """
        h = h5py.File( fname )
        p = h[pname]
        shape = p.attrs['h5sparse_shape']
        assert p.attrs['h5sparse_format'] == 'csr'
        vals = p['data'][:]
        indices = p['indices'][:]
        indptr = p['indptr'][:]
        obj = scipy.sparse.csr_matrix( ( vals, indices, indptr ), shape=shape )
        return obj



    
a40 = read_csr_from_hdf("Au6_s0_040_a.hdf", "Au6_s0_040_a")
toc("read hdf")

# ouch:
frames_a40 = [ a40[i*1024:(i+1)*1024].tocoo() for i in
               range(0,a40.shape[0]//1024) ]

toc("to coo")
# sanity check - frames larger that 2^16=65535 overflow
for f in frames_a40:
    assert cImageD11.sparse_is_sorted( f.row, f.col ) == 0

toc("check sorted")
# allocate labels
labels = [ np.zeros( f.nnz, 'i') for f in frames_a40 ]

toc("alloc labels")
# assign labels at lowest threshold

npx = [ cImageD11.sparse_connectedpixels(  frames_a40[i].data,
                                           frames_a40[i].row,
                                           frames_a40[i].col,
                                           0.,
                                           labels[i] )
        for i in range(len(frames_a40))] 

toc("assign t==0")

p2d = [ cImageD11.sparse_blob2Dproperties( frames_a40[i].data,
                                           frames_a40[i].row,
                                           frames_a40[i].col,
                                           labels[i],
                                           npx[i] )
        for i in range(len(frames_a40))]

toc("propsy t==0")

maxes = []
maxlabels = [ np.zeros( f.nnz, 'i') for f in frames_a40 ]
vmax = [ np.zeros( f.nnz, np.float32) for f in frames_a40 ]
imax = [ np.zeros( f.nnz, 'i') for f in frames_a40 ]
npmx = [ cImageD11.sparse_localmaxlabel(frames_a40[i].data,
                                        frames_a40[i].row,
                                        frames_a40[i].col,
                                        vmax[i], imax[i], maxlabels[i])
         for i in range(len(frames_a40))]

toc("maxlabel t==0")

pm2d = [ cImageD11.sparse_blob2Dproperties( frames_a40[i].data,
                                            frames_a40[i].row,
                                            frames_a40[i].col,
                                            maxlabels[i],
                                            npmx[i] )
         for i in range(len(frames_a40))]

toc("max props")

# Next step : pair max labels with connectedpixels labels
#   con_pk_list : shorter list
#   max_pk_list : longer list
#
# cases :     1 <-> 1 = same
#             1  -> n : split them or not ?
#                   shared boundaries / perimeter / other ?
#
# Parent child relation in this specific case. Localmax is always
# a connected object so that each maxlabel has a parent connected
#
# >>> maxlabels[0]
# array([  1,   1,   2, ..., 668, 668, 668], dtype=int32)
# >>> labels[0]
# array([  1,   1,   2, ..., 347, 347, 347], dtype=int32)
#
#

# need to merge 2D to nD objects ...




if 0:
    pl.plot( npx, "o-")
    pl.show()

