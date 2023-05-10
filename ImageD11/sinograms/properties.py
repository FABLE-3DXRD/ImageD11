
import os
import logging
import sys
import time
import tqdm
from timeit import default_timer
import numpy as np
import h5py
import hdf5plugin
import scipy.sparse
import scipy.sparse.csgraph
import ImageD11.sinograms.dataset
import ImageD11.sparseframe
import numba

### The first part of the code is all to run in parallel on a multicore machine
# 
# The 2D sinogram needs to be processed. Each frame has 4 neighbors.
# 2 in the same rotation and 2 in the previous and next rotation
# Work is split up to be one process doing two rows and joining them
# For each worker added you process a row twice
#
# A sparse matrix is built in shared memory.
# 
# perhaps some refactoring is needed
#
# it is building a sparse matrix of (npks x npks) which tags the overlaps.
# the number of entries in the matrix depends on how many peaks are overlapping.
#  ... so first the overlaps are counted
#  ... then each process shares their data
#  ... for now the matrix is not saved. Should only be about npks * noverlaps.
#

import multiprocessing as mp
from multiprocessing import shared_memory, resource_tracker

#####################################################################################
def remove_shm_from_resource_tracker():
    """Monkey-patch multiprocessing.resource_tracker so SharedMemory won't be tracked
    More details at: https://bugs.python.org/issue38119
    """
    self = None
    def fix_register(name, rtype):
        if rtype == "shared_memory":
            return
        return resource_tracker._resource_tracker.register(name, rtype)
    if resource_tracker.register is not fix_register:
        resource_tracker.register = fix_register
    def fix_unregister(name, rtype):
        if rtype == "shared_memory":
            return
        return resource_tracker._resource_tracker.unregister(name, rtype)
    if resource_tracker.unregister is not fix_register:
        resource_tracker.unregister = fix_register
    resource_tracker.unregister = fix_unregister
    if "shared_memory" in resource_tracker._CLEANUP_FUNCS:
        del resource_tracker._CLEANUP_FUNCS["shared_memory"]
######################################################################################

NPROC = max( 1, int(os.environ['SLURM_CPUS_PER_TASK']) - 1 )

class shared_numpy_array:
    """See: https://bugs.python.org/issue38119
    The multiprocessing pool must stick around until the process exits
    """
    def __init__(self, ary=None, shape=None, dtype=None, shmname=None,
                 fill = None):
        if ary is not None:
            shape = ary.shape
            dtype = ary.dtype
            self.nbytes = ary.nbytes
        else:
            self.nbytes = np.prod(shape) * np.dtype(dtype).itemsize
            
        if shmname is None:
            self.shm = shared_memory.SharedMemory(
                create=True, size = self.nbytes )
            self.creator=True
        else:
            self.shm = shared_memory.SharedMemory(
                create=False,
                name = shmname )
            self.creator=False
        self.array = np.ndarray( shape, dtype, buffer=self.shm.buf)
        if ary is not None:
            self.array[:] = ary
        if fill is not None:
            self.array[:] = fill
            
    def __del__(self):
        del self.array
        self.shm.close()
        if self.creator:
            try:
                self.shm.unlink()
            except Exception as e:
                print('Error: ',e)
                      
    def export(self):
        return { 'shape' : self.array.shape,
                 'dtype' : self.array.dtype,
                 'shmname' : self.shm.name }




class tictoc:
    def __init__(self):
        self.t = default_timer()
    def __call__(self, msg = '' ):
        t = default_timer()
        print("%s : %.6f /s"%( msg, t - self.t ) )
        self.t = default_timer()





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
    olap =  ImageD11.sparseframe.overlaps_linear( s.nnz.max()+1 )
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
    olap =  ImageD11.sparseframe.overlaps_linear( max(s1.nnz.max(), s2.nnz.max())+1 )
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


def props(scan, i, firstframe=None, algorithm='lmlabel' ):
    """
    scan = sparseframe.SparseScan object
    i = sinogram row id : used for tagging pairs
    algorithm = 'lmlabel' | 'cplabel' 

    Labels the peaks with lmlabel
    Assumes a regular scan for labelling frames
    returns ( row, properties[(s1,sI,sRow,sCol,frame),:], pairs, scan )
    """
    scan.sinorow = i
    getattr( scan, algorithm )( countall=False )  
    npks = scan.total_labels
    r = np.empty( (5,npks), np.int64 )
    s = 0
    if firstframe is None:        # not tested or used ...
        j0 = i * scan.shape[0]
    else:
        j0 = firstframe
    for j in range(scan.shape[0]):
        if scan.nnz[j] == 0:
            continue
        f0 = scan.getframe(j)
        e = s + scan.nlabels[j]
        # [1:] means skip the background labels == 0 output
        r[0,s:e] = np.bincount( f0.pixels['labels'] )[1:]
        wt = f0.pixels['intensity'].astype(np.int64)
        signal = np.bincount( f0.pixels['labels'], weights=wt)[1:]
        if signal.min() < 1:
            print("Bad data",scan.hname, scan.scan,i,j,
                  scan.nlabels[j], scan.nnz[j],f0.pixels['intensity'].min())
            raise Exception( 'bad data' )
        r[1,s:e] = signal
        r[2,s:e] = np.bincount( f0.pixels['labels'], weights=f0.row*wt )[1:]
        r[3,s:e] = np.bincount( f0.pixels['labels'], weights=f0.col*wt )[1:]
        r[4,s:e] = j + j0
        s = e
    # Matrix entries for this scan with itself:
    pairs = pairrow( scan, i )
    return r, pairs

###### testing / debug
if 0:
    def checkprops(i):
        hname = 'ds_MA4752_S4_2_XRD_DTL1z60_sparsefull.h5'
        ds = ImageD11.sinograms.dataset.load(hname)
        scan = ImageD11.sparseframe.SparseScan( hname[3:], ds.scans[i] )
        return props( scan, i )
    r, p = checkprops(400)
    for i,(k,(n,v)) in enumerate(p.items()):
        if i < 3 or i > (len(p)-3):
            print(i,k,n)
            print(v.T)
######

def get_start_end( n, p ):
    """ For splitting up the scan row over processes
    n = number of jobs
    p = number of processes
    
    All processes must do at least two rows to get the overlaps between rows.
    The work sharing is based on these overlaps.
    
    rows :  0  1  2  3  4  5   etc.
    overlap  01 12 23 34 45
    """
    overlaps = [ (i-1,i) for i in range(1,n) ]
    assert len(overlaps) >= p
    joins_per_job = np.zeros(p, int)
    for i,o in enumerate(overlaps):
        joins_per_job[i%p] += 1
    assert np.sum(joins_per_job) == len(overlaps)
    start = 0
    slices=[]
    for i in range(p):
        end = start + joins_per_job[i]
        ss = overlaps[start][0]
        se = overlaps[end-1][1]
        slices.append( (ss, se) )
        start += joins_per_job[i]
    return slices

def countse( se ):
    """ counts the points in a block of start/end"""
    k = 0
    for s, e in se:
        for i in range(s, e+1):
            k += 1
    return k

def testit():
    """ sanity check """
    for t in range(95,110):
        pairs = set()
        s = np.arange(t)
        r = set()
        for st, en in get_start_end( t, NPROC):
            last = None
            for i in range(st,en+1):
                x = s[i] # indexerror?
                r.add(i)
                if i > st:
                    pairs.add((i,i-1))
        assert len(pairs) == t-1
        assert len(r) == t
# testit()

def compute_storage( peaks ):
    '''make the cumulative sums of peaks to figure out storage space
    holds the row, npks, nii, nij
    '''
    P = {}
    M = {}
    for row, npks, nii, nij in peaks:
        if row in P:
            # verify duplicates
            assert P[row] == npks
            assert M[row][0] == nii
            if nij > M[row][1]:
                M[row] = nii, nij
        else:
            P[row] = npks
            M[row] = nii, nij
    ks = sorted( P.keys() )
    npk = np.array( [ (P[k],M[k][0],M[k][1]) for k in ks ] )
    return ks, npk

class pks_table:
    def __init__(self, 
                 npk=None,
                 ipk=None,
                 pk_props=None,
                 rc=None,
                 rpk=None,
                 glabel=None,
                 nlabel=0,
                 use_shm=False,
                ):
        """
        Cases: 
           Create from npks counting -> here
           Read from a file          -> classmethod pks_table.load( h5name )
           Read from shared memory   -> classmethod pks_table.fromSHM( h5name )
        """
        self.npk = npk
        self.ipk = ipk
        self.pk_props = pk_props
        self.rc = rc
        self.rpk = rpk
        self.glabel = glabel
        self.nlabel = nlabel
        self.use_shm = use_shm
        self.shared = {}
        # otherwise create
        if self.npk is not None:
            self.create( self.npk )
            
    def share(self, name, *args, **kwds):
        """puts the array into shared memory
        returns the shared copy
        """
        self.shared[name] = shared_numpy_array( *args, **kwds )
        return self.shared[name].array
            
    def create(self, npk):
        # [ nscans, 3 ]
        #       number_of_peaks_in_scan
        #       number_of_pairs_ii
        #       number_of_pairs_ij
        self.npk = npk
        s = npk.sum(axis=0)
        # pointers to peak tables per scan
        ipk = np.empty( npk.shape[0]+1, int )
        ipk[0] = 0
        ipk[1:] = np.cumsum( npk[:,0] )
        self.ipk = self.share( 'ipk', ipk )
        # pointers to r/c positions
        rpk = np.empty( npk.shape[0]+1, int )
        rpk[0] = 0
        rpk[1:] = np.cumsum( npk[:,1]+npk[:,2] )
        self.rpk = self.share( 'rpk', rpk )
        self.pk_props = self.share( 'pk_props', shape=(5, s[0]), dtype=np.int64 )
        self.rc = self.share( 'rc', shape=(3, s[1]+s[2]), dtype=np.int64 )
        
    def export(self):
        return { name : self.shared[name].export() for name in self.shared }
        
    def __del__(self):    
        del self.ipk
        del self.rpk
        del self.pk_props
        del self.rc
        names = list(self.shared.keys())
        for name in names:
            o = self.shared.pop( name )
            del o
            
    def guesschunk(self, ar, m=64*40):
        return None
        """ for parallel / compressed """
        chunk = list(ar.shape)
        nbytes = np.prod(chunk) * ar.dtype.itemsize
        if (nbytes // m) < pow(2,16):
            print('Guessing 1 chunk for', ar.shape, ar.dtype)
            return tuple(chunk)
        axsplit = np.argmax( chunk )
        n = chunk[axsplit]
        chunk[axsplit] = max(1, n//m)
        print('Guessing chunks for', ar.shape, ar.dtype, chunk)
        return tuple(chunk)
            
    def save(self, h5name, group='pks2d'):
        opts = {} # No compression is faster 
        with h5py.File( h5name, 'a' ) as hout:
            grp = hout.require_group( group )
            ds = grp.require_dataset( name = 'ipk', 
                                     shape = self.ipk.shape, 
                                     chunks = self.guesschunk( self.ipk ),
                                     dtype = self.ipk.dtype, **opts )
            ds[:] = self.ipk
            ds.attrs['description']='pointer to start of a scan'
            ds = grp.require_dataset( name = 'pk_props', 
                                     shape = self.pk_props.shape, 
                                     chunks = self.guesschunk( self.pk_props ),
                                     dtype = self.pk_props.dtype, **opts )
            ds[:] = self.pk_props
            ds.attrs['description']='[ ( s1, sI, srI, scI, id ),  Npks ]'
            ds = grp.require_dataset( name = 'npk', 
                                     shape = self.npk.shape, 
                                     chunks = self.guesschunk( self.npk ),
                                     dtype = self.npk.dtype, **opts )
            ds[:] = self.npk
            ds.attrs['description']="[ nscans, (N_peaks_in_scan, N_pairs_ii, N_pairs_ij) ]"
            if hasattr(self,'glabel'):
                ds = grp.require_dataset( name = 'glabel', 
                                         shape = self.glabel.shape, 
                                         chunks = self.guesschunk( self.glabel ),
                                         dtype = self.glabel.dtype, **opts )
                ds[:] = self.glabel
                grp.attrs['nlabel'] = self.nlabel

                
    @classmethod
    def fromSHM(cls, dct):
        names =  'ipk', 'rpk', 'pk_props', 'rc', 'glabel'
        sharrays = { name : shared_numpy_array( **dct[name] ) 
                    for name in names if name in dct }
        arrays = { name : sharrays[name].array 
                  for name in names if name in dct }
        o = cls( **arrays )
        o.shared = sharrays
        return o
    
    @classmethod
    def load(cls, h5name, h5group='pks2d'):
        with h5py.File( h5name, 'r' ) as hin:
            grp = hin[ h5group ]
            ipk = grp[ 'ipk' ][:]
            pk_props = grp[ 'pk_props' ][:]
            if 'glabel' in grp:
                glabel = grp['glabel'][:]
                nlabel = grp.attrs['nlabel']
            else:
                glabel = None
                nlabel = 0
            # rc?
            npk = grp['npk'][:]
            nlabel = grp.attrs['nlabel']
        obj = cls( ipk=ipk, pk_props=pk_props, glabel=glabel, nlabel=nlabel )
        obj.npk = npk # this is ugly. Sending as arg causes allocate.
        return obj                
                
    def find_uniq(self, outputfile=None):
        """ find the unique labels from the rc array """
        t = tictoc()
        n = self.ipk[-1]
        #        print("Row/col sparse array")
        #        for i in range(3):
        #            print(self.rc[i].dtype,self.rc[i].shape)i
        if outputfile is not None:
            with h5py.File(outputfile,"w") as hout:
                hout['data']=self.rc[2]
                hout['i'] = self.rc[0]
                hout['j'] = self.rc[1]
            return None, None
        coo = scipy.sparse.coo_matrix( (self.rc[2], 
                                        (self.rc[0], self.rc[1])), 
                                      shape=(n,n))
        t('coo')
        cc = scipy.sparse.csgraph.connected_components( coo, directed=False, return_labels=True )
        t('find connected components')
        self.cc = cc
        self.nlabel, self.glabel = cc
        return cc
    
    def pk2dmerge(self, omega, dty):
        """
        creates a dictionary of the 3D peaks
        """
        assert omega.shape == dty.shape
        assert omega.size > self.pk_props[4].max()
        
        out = np.zeros( (7, self.nlabel) , float)
        n = numbapkmerge( self.glabel, self.pk_props, omega, dty, out)
        allpks = {
            's_raw' : out[2]/out[1],
            'f_raw' : out[3]/out[1],
            'omega' : out[4]/out[1],
            'Number_of_pixels' : out[0],
            'sum_intensity' : out[1],
            'dty' : out[5]/out[1],
            'spot3d_id' : np.arange(len(out[0])),   # points back to labels in pk2d
            'npk2d': out[6],
        }
        return allpks
    
    def pk2d(self, omega, dty):
        s1, sI, srI, scI, frm = self.pk_props
        allpks = {
            's_raw' : srI / sI,
            'f_raw' : scI / sI,
            'omega' : omega.flat[ frm ],
            'dty' : dty.flat[ frm ],
            'Number_of_pixels' : s1,
            'sum_intensity' : sI,
            'spot3d_id' : self.glabel,
        }
        return allpks
            
if 0:
    def load_and_transpose( hname, itype, vtype ):
        """ Read in a coo file saved by pks_table.find_uniq
        """
        with h5py.File(hname,'r') as hin:
            di = hin['i'] 
            ii = np.empty( len(di)*2, itype )
            jj = np.empty( len(di)*2, itype )
            vv = np.empty( len(di)*2, vtype )
            ii[:len(di)] = hin['i'][:] # this should cast when filling
            ii[len(di):] = hin['j'][:]
            jj[:len(di)] = hin['j'][:]
            jj[len(di):] = hin['i'][:]
            vv[:len(di)] = hin['data'][:]
            vv[len(di):] = hin['data'][:]
        return ii,jj,vv
    

@numba.njit
def numbapkmerge( labels, pks, omega, dty, out):
    for k in range(len(labels)):
        frm = pks[ 4, k ]
        o = omega.flat[ frm ]
        y = dty.flat[ frm ]
        l = labels[ k ]
        out[0,l] += pks[0,k]    # s1 == number of pixels in a peak
        out[1,l] += pks[1,k]    # sI == sum of the intensity
        out[2,l] += pks[2,k]    # srI === sum of intensity * row
        out[3,l] += pks[3,k]    # scI === sum of intensity * column
        out[4,l] += o * pks[1,k]
        out[5,l] += y * pks[1,k]
        out[6,l] += 1           # s0 == number of 2D peaks
    return k
    


def pks_table_from_scan( sparsefilename, ds, row ):
    """
    Labels one rotation scan to a peaks table
    
    sparsefilename = sparse pixels file
    dataset = ImageD11.sinograms.dataset
    row = index for dataset.scan[ row ]
    
    returns a pks_table.
        You might want to call one of "save" or "pk2d" or "pk2dmerge" on the result
    
    This is probably not threadsafe 
    """ 
    sps = ImageD11.sparseframe.SparseScan( sparsefilename, ds.scans[row] )
    peaks, pairs = ImageD11.sinograms.properties.props( sps, row )
    # which frame/peak is which in the peaks array
    # For the 3D merging 
    n1 = sum( pairs[k][0] for k in pairs ) # how many overlaps were found:
    npk = np.array(  [ ( peaks.shape[1], n1, 0), ] )
    pkst = ImageD11.sinograms.properties.pks_table( npk = npk, use_shm=False )
    pkst.pk_props = peaks
    rc =  pkst.rc
    s = 0
    pkid = np.concatenate(([0,], np.cumsum(sps.nlabels)))
    for (row1, frame1, row2, frame2), (npairs, ijn) in pairs.items():
        # add entries into a sparse matrix
        # key:   (500, 2, 501, 2893)
        # value: (41, array([[  1,   1,   5],
        #                    [  1,   2,   2], ...
        if npairs == 0:
            continue
        #assert (row1 == i) and (row2 == i), (row1,row2,i)
        e = s + npairs
        assert e <= rc.shape[1]
        #rc[ 0, s : e ] = ip[row1] + pkid[row1][frame1] + ijn[:,0] - 1       # col
        #rc[ 1, s : e ] = ip[row2] + pkid[row2][frame2] + ijn[:,1] - 1       # row
        rc[ 0, s : e ] = pkid[frame1] + ijn[:,0] - 1       # col
        rc[ 1, s : e ] = pkid[frame2] + ijn[:,1] - 1       # row
        rc[ 2, s : e ] = ijn[:,2]                          # num pixels
        s = e
    uni = pkst.find_uniq()
    
    if 0:   # future TODO : scoring overlaps better. 
        ks = list( pairs.keys() )
        k = ks[100]
        ipt = pkid
        npx1 = peaks[0][ipt[k[1]]:ipt[k[1]+1]] # number of pixels in pk1 (allpeaks)
        npx2 = peaks[0][ipt[k[3]]:ipt[k[3]+1]] # number of pixels in pk2 (allpeaks)
        ijo = pairs[k][1]                      # which peaks overlap each other
        # npx1.shape, npx2.shape, ijo.shape, ijo.max(axis=0)
        n1 = npx1[ijo[:,0]-1]     # number of pixels in pk1 (overlapping)
        n2 = npx2[ijo[:,1]-1]     # number of pixels in pk2 (overlapping)
        no = ijo[:,2]
    
    # pk3d = pkst.pk2dmerge(ds.omega, ds.dty)
    return pkst    
    
    
    
def process(qin, qshm, qout, hname, scans, options):
    remove_shm_from_resource_tracker()
    start, end = qin.get()
    n2 = 0
    # allocate lists to hold results
    nrows = end - start + 1
    pii = {} 
    pij = {}
    mypks = {}
    pkid = {}
    prev = None # suppress flake8 idiocy
    for i in range(start, end+1):
        scan = ImageD11.sparseframe.SparseScan( hname, scans[i] )
        global omega
        scan.motors['omega'] = omega[i] 
        mypks[i], pii[i] = props(scan, i, algorithm = options['algorithm'] )
#        nlabels[i] = scan.nlabels # peaks per frame information
        pkid[i] = np.concatenate(([0,], np.cumsum(scan.nlabels)))
        n1 = sum( pii[i][k][0] for k in pii[i] )
        if i > start:
            pij[i] = pairscans( scan, prev ) 
            n2 = sum( pij[i][k][0] for k in pij[i] )
        # number of pair overlaps required for the big matrix
        qout.put( (i, len(mypks[i][0]), n1, n2) )
        prev = scan
    # Now we are waiting for the shared memory to save the results
    shm = qshm.get()
    pkst = pks_table.fromSHM( shm )
    ip = pkst.ipk
    rc = pkst.rc
    # For each row, save our local results
    for i in range(start, end+1):
        if (i == start) and (i > 0):
            # will be done by the previous worker
            continue
        # now copy our results to shared memory. First the peaks:
        pkst.pk_props[:, ip[i]: ip[i+1] ] = mypks[i]
        #
        # find the unique starting id for each frame:
        # Where to store these in the shared memory
        s = ipstart = pkst.rpk[i]
        ipend = pkst.rpk[i+1]
        #
        if 0: # debugging
            n1 = sum( pii[i][k][0] for k in pii[i] )
            if i > start:
                n2 = sum( pij[i][k][0] for k in pij[i] )
            else:
                n2 = 0
            # # will be done by the previous worker
            # print('Debugging:',i, ipstart, ipend, n1, n2, ipend-ipstart, n1+n2)
            assert (ipend - ipstart)==(n1 + n2)
        #
        for (row1, frame1, row2, frame2), (npairs, ijn) in pii[i].items():
            # add entries into a sparse matrix
            # key:   (500, 2, 501, 2893)
            # value: (41, array([[  1,   1,   5],
            #                    [  1,   2,   2], ...
            if npairs == 0:
                continue
            assert (row1 == i) and (row2 == i)
            e = s + npairs
            assert e <= ipend
            rc[ 0, s : e ] = ip[row1] + pkid[row1][frame1] + ijn[:,0] - 1       # col
            rc[ 1, s : e ] = ip[row2] + pkid[row2][frame2] + ijn[:,1] - 1       # row
            rc[ 2, s : e ] = ijn[:,2]                          # num pixels
            s = e
        if i == 0:
            continue # no pairs to previous row exit
        for (row1, frame1, row2, frame2), (npairs, ijn) in pij[i].items():
            if npairs == 0:
                continue
            assert (row1 == i) and (row2 == i-1)
            e = s + npairs
            assert e <= ipend
            rc[ 0, s : e ] = ip[row1] + pkid[row1][frame1] + ijn[:,0] - 1       # col
            rc[ 1, s : e ] = ip[row2] + pkid[row2][frame2] + ijn[:,1] - 1       # row
            rc[ 2, s : e ] = ijn[:,2]                          # num pixels
            s = e
    qout.put(start)
    
        
def goforit(ds, sparsename, options):
    qin = mp.Queue(maxsize=NPROC)
    qshm = mp.Queue(maxsize=NPROC)
    qresult = mp.Queue(maxsize=len(ds.scans))
    out = []
    with mp.Pool(NPROC, initializer=process, 
               initargs=(qin, qshm, qresult, sparsename, ds.scans, options)) as pool:
        slices = get_start_end( len(ds.scans), NPROC )
        for i,(s,e) in enumerate(slices):
            qin.put((s,e))
        waiting = countse( slices )
        for i in tqdm.tqdm(range(waiting)):
            ans = qresult.get()
            out.append(ans)
        ks, P = compute_storage( out )
        assert len(ks) == len(ds.scans)
        mem = pks_table( P )
        shm = mem.export()
        # workers fill the shared memory
        for i in range(NPROC):
            qshm.put( shm )
        dones = set()
        for i in tqdm.tqdm(range(NPROC)):
            # check done
            f = qresult.get()
            dones.add( f )
            if len(dones) == NPROC:
                break
    return out, ks, P, mem


default_options = { 'algorithm' : 'lmlabel' }
                   
def main( dsname, sparsename, pkname, options = default_options ):
    if os.path.exists(pkname):
        logging.warning("Your output file already exists. May fail. %s"%(pkname))
    global NPROC, omega
    try:
        rmem = None
        t = tictoc()
        ds = ImageD11.sinograms.dataset.load(dsname)
        omega = ds.omega
        t('read ds %s'%(dsname))
        nscans = len(ds.scans)
        if NPROC > nscans-1:
            NPROC = max(1,nscans-1)
        print('Nscans',nscans,'NPROC', NPROC)
        peaks, ks, P, rmem = goforit(ds, sparsename, options)
        t('%d label and pair'%(len(rmem.pk_props[0])))
#        rmem.save( pkname )        
#        t('cache')
        cc = rmem.find_uniq()
        t('%s connected components'%(str(cc[0])))
        rmem.save( pkname )
        t('write hdf5')
    except Exception as e:
        print('Unhandled exception:',e)
        raise
    finally:
        if rmem is not None:
            del rmem
    return



if __name__=="__main__":
    dsname = sys.argv[1]
    sparsename = sys.argv[2]
    pkname = sys.argv[3]
    algorithm = 'lmlabel'
    if len(sys.argv)>=5:
        algorithm = sys.argv[4]
    options = {'algorithm': algorithm }
    main( dsname, sparsename, pkname, options )

    print("Your stuff left in shm:")
    os.system("ls -l /dev/shm | grep %s"%(os.environ['USER']))
    print("... all done")


"""
import numba
@numba.njit
def pkmerge( labels, pks, omega, dty, nlm, out):
    k = 0
    for i in range(nlm.size):
        o = omega[i]
        y = dty[i]
        for j in range(nlm[i]): # peaks on this frame
            l = labels[k]
            out[0,l] += pks[0,k]
            out[1,l] += pks[1,k]
            out[2,l] += pks[2,k]
            out[3,l] += pks[3,k]
            out[4,l] += o * pks[1,k]
            out[5,l] += y * pks[1,k]
            k += 1
    return k    
    
class h5pks_table:
    def __init__(self, h5name, group='pks2d'):
        with h5py.File( h5name, 'r' ) as hin:
            grp = hin[group]
            for name in 'npk ipk pk_props rpk'.split():
                setattr( self, name, grp[name][:] )
                
s1, sI, sr, sc, frm = pks
gpks = np.empty( (6,labels.max()+1), float )
                

pkmerge( labels, pks, ds.omega.ravel(), ds.dty.ravel(), ds.nlm.ravel(), gpks)
sc = gpks[2] / gpks[1]
fc = gpks[3] / gpks[1]
om = gpks[4] / gpks[1]
ty = gpks[5] / gpks[1]
"""
