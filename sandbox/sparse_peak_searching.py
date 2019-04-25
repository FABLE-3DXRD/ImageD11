
from __future__ import print_function, division

import time, sys
import h5py, scipy.sparse, numpy as np, pylab as pl
from ImageD11 import cImageD11

# see also sandbox/collect_peak_pixels.py

class sparse_frame( object ):
    """
    Holds the output of an image segmentation
    """
    def __init__(self, row, col, data):
        """
        Sparse matrix, follows scipy.sparse.coo_matrix
        row = slow index (uint16)
        col = fast index (uint16)
        data = intensity
        nnz  = number of entries
        """
        self.row = np.array(row, dtype=np.uint16)
        self.col = np.array(col, dtype=np.uint16)
        self.data = data
        self.nnz = len(self.row)
                
    def check(self):
        """
        Verify data type and that indices are sorted
        """
        assert self.row.dtype == np.uint16
        assert self.col.dtype == np.uint16
        assert cImageD11.sparse_is_sorted( self.row, self.col ) == 0
        for ar in (self.row, self.col, self.data):
            assert len(ar.shape)==1
            assert ar.shape[0] == self.nnz


class cplabels( sparse_frame ):
    """
    The labels for an image segmetation
    """
    def __init__(self, spframe, threshold=0 ):
        self.frame = spframe
        self.threshold = threshold 
        self.data = np.zeros( self.frame.nnz, 'i' )
        self.nlabel = cImageD11.sparse_connectedpixels(
            self.frame.data, self.frame.row, self.frame.col,
            self.threshold,   self.data )
        self.moments = cImageD11.sparse_blob2Dproperties(
            self.frame.data, self.frame.row, self.frame.col,
            self.data, self.nlabel )

    def overlaps(self, other):
        ki = np.empty(  self.frame.nnz, 'i' )
        kj = np.empty( other.frame.nnz, 'i' )
        npx = cImageD11.sparse_overlaps( self.frame.row, self.frame.col, ki,
                                         other.frame.row, other.frame.col, kj)
        row = self.data[ ki[:npx] ]
        col = other.data[ kj[:npx] ]
        ect = np.empty( npx, 'i')
        tj  = np.empty( npx, 'i')
        tmp = np.empty( max(self.nlabel, other.nlabel)+1, 'i')
        
        nedge = cImageD11.compress_duplicates( row, col, ect, tj, tmp )
        crow = row[:nedge]
        ccol = col[:nedge]
        cdata = ect[:nedge]
        cedges = scipy.sparse.coo_matrix( ( cdata, (crow, ccol)) )
        if 0: # debugging
            row = self.data[ ki[:npx] ]
            col = other.data[ kj[:npx] ]
            data = np.ones( npx, 'i')
            if hasattr( scipy.sparse.coo_matrix, "_sum_duplicates" ):
                edges = scipy.sparse.coo_matrix( (data, (row, col)) )
                edges.sum_duplicates()
            else:
                row, col, data = _sum_duplicates( row, col, data )
                edges = scipy.sparse.coo_matrix( (data, (row, col)) )
            assert edges.nnz == nedge
            assert (edges.row == crow).all()
            assert (edges.col == ccol).all()
            assert (edges.data == cdata).all()
        return cedges



class mxlabels( cplabels ):
    """
    The labels for an image segmetation
    """
    def __init__(self, spframe, threshold=0 ):
        self.frame = spframe
        self.threshold = threshold
        nnz = self.frame.nnz 
        self.data = np.zeros(nnz, 'i' )
        vmx = np.zeros( nnz, np.float32 )
        imx = np.zeros( nnz, 'i')
        self.nlabel = cImageD11.sparse_localmaxlabel(
            self.frame.data, self.frame.row, self.frame.col,
            vmx, imx,   self.data )
        self.moments = cImageD11.sparse_blob2Dproperties(
            self.frame.data, self.frame.row, self.frame.col,
            self.data, self.nlabel )
        
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

frames = []
for i in range(0,a40.shape[0]//1024):
    tmp  =  a40[i*1024:(i+1)*1024].tocoo()
    frame = sparse_frame( tmp.row, tmp.col, tmp.data )
    frame.check()
    frames.append( frame )

frames_a40 = frames
toc("to coo")

cl = [ cplabels(f) for f in frames ]

toc("assign and count t==0")

ml = [ mxlabels(f) for f in frames ]

toc("maxlabel")


# Next step : pair up the max labels with connectedpixels labels?
#   con_pk_list : shorter list
#   max_pk_list : longer list
#
# cases :     1 <-> 1 = same
#             1  -> n : split them or not ?
#                   shared boundaries / perimeter / other ?
#
# Parent child relation in this specific case. Localmax is always
# a connected object so that each maxlabel has one parent connected
#
# >>> maxlabels[0]
# array([  1,   1,   2, ..., 668, 668, 668], dtype=int32)
# >>> labels[0]
# array([  1,   1,   2, ..., 347, 347, 347], dtype=int32)
#
# Cases:
#    noisy satellites :
#         identify as s_I(child) or s_1(child) << s_I(parent)
#                     many pixels on the border or close to threshold
#    split peaks:
#         tails are touching : I_shared << I_max[i,j]
#         peaks are overlapping : I_shared ~ I_max * frac
# Classify pixels into:
#    fully inside child (8 neighbors with equal label)
#       9 labels     0, 1, 2,
#                    3, 4, 5,
#                    6, 7, 8
#    can be 0 = border        
#        -> k == neighbor link
#    labels[4] -> labels[j] score v[4] ... when we see j->4 we count v[j]
#    eventually count 8X intensity in links
#
# Graphing :
#   max-labels are nodes in the graph (image_id, node_id)
#     links are then pairs of these.
#

edges = []

for i in range(len(cl)):
    edges.append( cl[i].overlaps( ml[i] ) )
toc("cp mx pairs")
cpedges=[]
for i in range(1, len(cl)):
    cpedges.append( cl[i].overlaps( cl[i-1] ) )
toc("cp cp[-1] pairs")
mxedges=[]
for i in range(1, len(ml)):
    mxedges.append( ml[i].overlaps( ml[i-1] ) )
toc("mx mx[-1] pairs")

#    allpairs.append(
#        pairs(frames_a40[i-1], maxlabels[i-1], frames_a40[i], maxlabels[i] )) 

# foreach frame :
#     row/col/intensity/labels
# foreach pixel frame_i and frame_j
#     overlaps : assign intensity to i,li,j,lj bucket
#     in one only : assign intensity to i,li   bucket




if 0:
    pl.plot( npx, "o-")
    pl.show()

