

import h5py, os
import numpy as np
from ImageD11 import sparseframe, cImageD11
import scipy.sparse
import scipy.sparse.csgraph

# given a file holding pixels in _sparse.h5 format
#  ... we expect to find:
#  ...          shape = image shape
#  ...          nnz = pixels per frame. Shape = number of frames.
#  ...          row/col/intensity


def pairrow( s , row):
    """
    s = SparseScan

    returns SparseScan which adds a dictionary of pairs:
         sourceFrame    destFrame
        [ idty, iomega, idty, iomega ] : nedge, array( (nedge, 3) )
                                                     (src, dest, npixels)
    """
    s.omegaorder = np.argsort( s.motors["omega"] )
    s.sinorow = row
    olap = sparseframe.overlaps_linear( s.nnz.max()+1 )
    pairs = {}
    for i in range(1, len(s.omegaorder)):
        if (s.nnz[s.omegaorder[i]] == 0) or (s.nnz[s.omegaorder[i-1]] == 0):
            continue
        f0 = s.getframe( s.omegaorder[i-1] )
        f1 = s.getframe( s.omegaorder[i] )
        ans = olap( f0.row, f0.col, f0.pixels['labels'], s.nlabels[ s.omegaorder[i-1] ],
                    f1.row, f1.col, f1.pixels['labels'], s.nlabels[ s.omegaorder[i]   ] )
        pairs[ row, s.omegaorder[i-1], row, s.omegaorder[i] ] =  ans
    return pairs



def pairscans( s1, s2, omegatol = 0.051 ):
    olap = sparseframe.overlaps_linear( max(s1.nnz.max(), s2.nnz.max())+1 )
    assert len(s1.nnz) == len(s2.nnz )
    pairs = {}
    for i in range(len(s1.nnz)):
        # check omega angles match
        o1 = s1.motors['omega'][s1.omegaorder[i]]
        o2 = s2.motors['omega'][s2.omegaorder[i]]
        assert abs( o1 - o2 ) < omegatol, 'diff %f, tol %f'%(o1-o2, omegatol)
        if (s1.nnz[s1.omegaorder[i]] == 0) or (s2.nnz[s2.omegaorder[i]] == 0):
            continue
        f0 = s1.getframe( s1.omegaorder[i] )
        f1 = s2.getframe( s2.omegaorder[i] )
        ans = olap( f0.row, f0.col, f0.pixels['labels'], s1.nlabels[ s1.omegaorder[i] ],
                    f1.row, f1.col, f1.pixels['labels'], s2.nlabels[ s2.omegaorder[i] ] )
        pairs[ s1.sinorow, s1.omegaorder[i], s2.sinorow, s2.omegaorder[i] ] =  ans
    return pairs




if 0: # remove unused code in here
  opts = {
        "chunks": (10000,),
        "maxshape": (None,),
        "compression": "lzf",
        "shuffle": True,
       }
    
  def add_connected_peaklabel( h5name, h5path, threshold=0 ):
    # this will load the whole scan into memory.
    scan = sparseframe.SparseScan( h5name, h5path )
    scan.cplabel( threshold, countall=False )
    with h5py.File( h5name, "a" ) as hsp:
        grp = hsp.require_group( h5path )
        if 'cplabel' in grp:
            grp['cplabel'][:] = scan.labels
            grp['Ncplabel'][:] = scan.nlabels
        else:
            grp.create_dataset( name = 'cplabel', data = scan.labels, **opts )
            grp.create_dataset( name = 'Ncplabel', data = scan.nlabels )
    return scan.nlabels, scan.labels

  def add_localmax_peaklabel( h5name, h5path, smooth=True ):
    scan = sparseframe.SparseScan( h5name, h5path )
    scan.lmlabel( smooth=smooth, countall=False )
    with h5py.File( h5name, "a" ) as hsp:
        grp = hsp.require_group( h5path )
        if 'lmlabel' in grp:
            grp['lmlabel'][:] = scan.labels
            grp['Nlmlabel'][:] = scan.nlabels
        else:
            grp.create_dataset( name = 'lmlabel', data = scan.labels, **opts )
            grp.create_dataset( name = 'Nlmlabel', data = scan.nlabels )
    return scan.nlabels, scan.labels





class IncrementalSymmetricCOOMatrix:

    def __init__(self, shape, dtype):
        self.dtype = dtype
        self.shape = shape
        assert self.shape[0] == self.shape[1]
        self.rows = []
        self.cols = []

    def addeye(self):
        r = np.arange(self.shape[0], dtype=np.int32)
        self.append( r, r )

    def append(self, i, j):
        self.rows.append(i)
        self.cols.append(j)

    def concatenate(self):
        print("Concatenating")
        self.rows = np.concatenate(self.rows)
        self.cols = np.concatenate(self.cols)
        self.data = np.ones( len(self.rows), dtype=np.uint8 ) # 1/0

    def tocoo(self):
        return scipy.sparse.coo_matrix((self.data, (self.rows, self.cols)),
                                         shape=self.shape)


def add2mat( ds, r1, r2, pairs ):
    # add entries into a sparse matrix
    #
    # key (500, 2, 501, 2893)
    #  (41, array([[  1,   1,   5],
    #               [  1,   2,   2],
    ct = 0
    p = 0
    for k in pairs.keys():
        ky1, ko1, ky2, ko2 = k
        if ky1 != r1 or ky2 != r2:
            raise Exception("logic error in add2mat")
        assert ky1 == r1
        assert ky2 == r2
        npairs, ijn = pairs[k]
        if npairs == 0:
            continue
        assert ijn.shape[0] == npairs
        ds.coomat.append( ds.firstpk[ky1,ko1] + ijn[:,0] - 1,
                          ds.firstpk[ky2,ko2] + ijn[:,1] - 1 )
        #                          ijn[:,2] )
        ct += npairs



def build_overlap_matrix(ds, hname):
    import tqdm
    names='row','col','intensity','lmlabel'
    lastscan = None
    npk = ds.nlm.sum()
    ds.firstpk = np.zeros( ds.nlm.shape, int )
    s = ds.firstpk.shape
    ds.firstpk.shape=-1
    ds.firstpk[1:] = ds.nlm.ravel().cumsum()[:-1]
    ds.firstpk.shape = s
    ds.coomat = IncrementalSymmetricCOOMatrix( (npk, npk), np.int32 )
    ds.coomat.addeye()

    def rd1( hname, scan, row ):
        spscan = sparseframe.SparseScan(hname, scan, names=names)
        spscan = pairrow( spscan, ds, row )
        return spscan

    def pr2( scan1, scan2 ):
        #assert scan1.sinorow == row1
        #assert scan2.sinorow == row2
        pairs = pairscans( ds, scan1, scan2 )
        return pairs #, row1, row2

    lo = 0
    hi = len(ds.scans)
    s1 = rd1( hname, ds.scans[lo], lo )
    add2mat( ds, lo, lo, s1.pairs )

    for row in tqdm.tqdm( range(lo+1,hi) ):
        s2 = rd1( hname, ds.scans[row], row )
        add2mat( ds, row, row, s2.pairs )
        p = pairscans( ds, s1, s2 )
        add2mat( ds, row-1, row, p )
        s1 = s2
    return ds


def find_global_peaklabel(ds, hname):
    ds = build_overlap_matrix( ds, hname )
    m = ds.coomat
    m.concatenate()
    coo = scipy.sparse.coo_matrix( (np.ones( len(m.rows), dtype=np.uint32 ),
                                    (m.rows, m.cols)), shape=m.shape)
    csr = coo.tocsr()
    csr += csr.T
    cc = scipy.sparse.csgraph.connected_components( csr )
    return cc



