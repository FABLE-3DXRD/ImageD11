
import os
import sys
import time
import tqdm
from timeit import default_timer
import numpy as np
import h5py
import scipy.sparse
import scipy.sparse.csgraph
import ImageD11.sinograms.dataset
import ImageD11.sinograms.peaklabel
import ImageD11.sparseframe

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
        return resource_tracker._resource_tracker.register(self, name, rtype)
    resource_tracker.register = fix_register
    def fix_unregister(name, rtype):
        if rtype == "shared_memory":
            return
        return resource_tracker._resource_tracker.unregister(self, name, rtype)
    resource_tracker.unregister = fix_unregister
    if "shared_memory" in resource_tracker._CLEANUP_FUNCS:
        del resource_tracker._CLEANUP_FUNCS["shared_memory"]
######################################################################################
remove_shm_from_resource_tracker()
NPROC = max( 1, int(os.environ['SLURM_CPUS_PER_TASK']) - 1 )


class tictoc:
    def __init__(self):
        self.t = default_timer()
    def __call__(self, msg = '' ):
        t = default_timer()
        print("%s : %.6f /s"%( msg, t - self.t ) )
        self.t = default_timer()

def props(scan, i, firstframe=None):
    """
    scan = sparseframe.SparseScan object
    i = sinogram row id : used for tagging pairs

    Labels the peaks with lmlabel
    Assumes a regular scan for labelling frames
    returns ( row, properties[(s1,sI,sRow,sCol,frame),:], pairs, scan )
    """
    scan.sinorow = i
    scan.lmlabel(countall=False)
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
        r[0,s:e] = np.bincount( f0.pixels['labels']-1 )
        r[1,s:e] = np.bincount( f0.pixels['labels']-1, weights=f0.pixels['intensity'] )
        r[2,s:e] = np.bincount( f0.pixels['labels']-1, weights=f0.row*f0.pixels['intensity'] )
        r[3,s:e] = np.bincount( f0.pixels['labels']-1, weights=f0.col*f0.pixels['intensity'] )
        r[4,s:e] = j + j0
        s = e
    # Matrix entries for this scan with itself:
    pairs = ImageD11.sinograms.peaklabel.pairrow( scan, i )
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


class pks_table:

    def __init__(self, npk=None, **kwds):
        if npk is not None:
            self.create( npk )
        else:
            self.ipk = shared_numpy_array( **kwds['ipk'] )
            self.rpk = shared_numpy_array( **kwds['rpk'] )
            self.pk_props = shared_numpy_array( **kwds['pk_props'] )
            self.rc = shared_numpy_array( **kwds['rc'] )

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
        self.ipk = shared_numpy_array( ipk )
        # pointers to r/c positions
        rpk = np.empty( npk.shape[0]+1, int )
        rpk[0] = 0
        rpk[1:] = np.cumsum( npk[:,1]+npk[:,2] )
        self.rpk = shared_numpy_array( rpk )
        self.pk_props = shared_numpy_array( shape=(5, s[0]), dtype=np.int64 )
        self.rc = shared_numpy_array( shape=(3, s[1]+s[2]), dtype=np.int32 )
        
    def export(self):
        return { 'ipk': self.ipk.export(),
                 'rpk': self.rpk.export(),
                 'pk_props' : self.pk_props.export(),
                 'rc' : self.rc.export() }
    
    def __del__(self):
        del self.ipk
        del self.rpk
        del self.pk_props
        del self.rc

    def save(self, h5name, group='pks2d'):
        opts = {'compression':'gzip',
                'compression_opts':2,
                'shuffle':True,
               }
        with h5py.File( h5name, 'a' ) as hout:
            grp = hout.require_group( group )
            for name in 'ipk pk_props'.split():
                # ipk = pointer to start of a scan
                # pk_props = peak properties
                data = getattr( self, name ).array
                ds = grp.require_dataset( name = name, shape = data.shape, dtype = data.dtype, **opts )
                ds[:] = data
            # wtf is npk?
            data = self.npk
            ds = grp.require_dataset( name = 'npk', shape = data.shape, dtype = data.dtype, **opts )
            ds[:] = data
            if hasattr(self,'cc'):
                ds = grp.require_dataset( name = 'glabel', 
                                         shape = self.cc[1].shape, 
                                         dtype = self.cc[1].dtype, **opts )
                ds[:] = self.cc[1]
                grp.attrs['nlabel'] = self.cc[0]
            
    def find_uniq(self):
        """ find the unique labels from the rc array """
        t = tictoc()
        n = self.ipk.array[-1]
        coo = scipy.sparse.coo_matrix( (self.rc.array[2], (self.rc.array[0], self.rc.array[1])), shape=(n,n))
        t('coo')
        csr = coo.tocsr()
        t('tocsr')
        csr += csr.T
        t('transpose')
        cc = scipy.sparse.csgraph.connected_components( csr )
        t('find connected components')
        self.cc = cc
        return cc
    

def process(qin, qshm, qout, hname, scans):
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
        mypks[i], pii[i] = props(scan, i)
#        nlabels[i] = scan.nlabels # peaks per frame information
        pkid[i] = np.concatenate(([0,], np.cumsum(scan.nlabels)))
        n1 = sum( pii[i][k][0] for k in pii[i] )
        if i > start:
            pij[i] = ImageD11.sinograms.peaklabel.pairscans( scan, prev ) 
            n2 = sum( pij[i][k][0] for k in pij[i] )
        # number of pair overlaps required for the big matrix
        qout.put( (i, len(mypks[i][0]), n1, n2) )
        prev = scan
    # Now we are waiting for the shared memory to save the results
    shm = qshm.get()
    pkst = pks_table( **shm )
    ip = pkst.ipk.array
    rc = pkst.rc.array
    # For each row, save our local results
    for i in range(start, end+1):
        if (i == start) and (i > 0):
            # will be done by the previous worker
            continue
        # now copy our results to shared memory. First the peaks:
        pkst.pk_props.array[:, ip[i]: ip[i+1] ] = mypks[i]
        #
        # find the unique starting id for each frame:
        # Where to store these in the shared memory
        s = ipstart = pkst.rpk.array[i]
        ipend = pkst.rpk.array[i+1]
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
    
        
def goforit(ds, sparsename):
    qin = mp.Queue(maxsize=NPROC)
    qshm = mp.Queue(maxsize=NPROC)
    qresult = mp.Queue(maxsize=len(ds.scans))
    out = []
    with mp.Pool(NPROC, initializer=process, 
               initargs=(qin, qshm, qresult, sparsename, ds.scans)) as pool:
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

def main( dsname, sparsename, pkname ):
    global NPROC
    try:
        rmem = None
        t = tictoc()
        ds = ImageD11.sinograms.dataset.load(dsname)
        t('read ds %s'%(dsname))
        nscans = len(ds.scans)
        if NPROC > nscans-1:
            NPROC = nscans-1
        print('Nscans',nscans,'NPROC', NPROC)
        peaks, ks, P, rmem = goforit(ds, sparsename)
        t('%d label and pair'%(len(rmem.pk_props.array[0])))
#        rmem.save( pkname )        
#        t('cache')
        cc = rmem.find_uniq()
        t('%d connected components'%(cc[0]))
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
    main( dsname, sparsename, pkname )

    print("Your stuff left in shm:")
    os.system(f"ls -l /dev/shm | grep {os.environ['USER']}")
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
