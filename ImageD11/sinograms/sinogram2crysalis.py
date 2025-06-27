
"""WARNING: work in in progress"""

from __future__ import print_function

import h5py
import numpy as np
import numba
import sys, os
import ImageD11.cImageD11
import ImageD11.sinograms.dataset
import ImageD11.sinograms.lima_segmenter
import fabio.app.eiger2crysalis, fabio.limaimage
import bslz4_to_sparse
import tqdm
import concurrent.futures
import threading
import hdf5plugin

readlock = threading.Lock()

# Dummy object 
class empty:
    pass
        
# Dummy fabio image for the data
class SinoScan( fabio.limaimage.LimaImage ):
    """ Inherit limaimage to be allowed by converter """
    def __init__(self, dataset, num=0):
        self.header = {}
        self.dataset = dataset
        self.data = self.dataset[num]
        self._nframes = len(self.dataset)
        self.h5 = empty()  # fake this out
        self.h5.attrs = {'default': None}
        
    def close(self):
        pass
    
    def __iter__(self):
        o = empty()
        o.header = {}
        for o.data in self.dataset:
            yield o
            
            
def name_decorator(method):
    def decorate_name(self = None):
        headers = method(self)
        print("Patching header")
        headers["drealpixelsizex"] = 0.075
        headers["drealpixelsizey"] = 0.075
        headers["dexposuretimeinsec"] = 0.1 # fixme
        return headers
    return decorate_name

if not hasattr( fabio.app.eiger2crysalis.Converter, 'patched' ):
    fabio.app.eiger2crysalis.Converter.common_headers = name_decorator(fabio.app.eiger2crysalis.Converter.common_headers)
    
# patch the converter headers
fabio.app.eiger2crysalis.Converter.patched = True

    
@numba.njit(nogil=True)
def accumulate( npixels, values, indices, dest, num ):
    d0 = dest[0].size * num
    for i in range(npixels):
        dest.flat[ d0 + indices[i] ] += values[i]

@numba.njit(parallel=True)
def padd( values, destination ):
    for i in numba.prange( values.size ):
        destination.flat[i] += values.flat[i]
    

def integrate_sino( ds, image_mask, sinomask=None, nomega = 10):
    """
    ds = ImageD11.sinograms.dataset that describes scandata
    sinomask = which frames to use
    nomega = number of omegasteps to combine
    """
    #flake8: global readlock 
    
    if sinomask is not None:
        assert sinomask.shape == ds.shape
    else:
        sinomask = np.ones( ds.shape, bool )
    
    # inverse of omega step size
    istep = (len(ds.obincens)-1)/(ds.obincens[-1]-ds.obincens[0])
    # which omega bin is this frame ?
    obin = np.round( istep * (ds.omega - ds.obincens[0])).astype(int)

    nout = ds.shape[1] // nomega
    destination = np.zeros( (nout, 2162, 2068), np.uint32 )
    address = 0
    npx = 0
    dt = None
    for num_frames, sourcefile in zip( ds.frames_per_file, ds.imagefiles ):
        todo = np.nonzero( sinomask.flat[address : address + num_frames] )[0]
        if len(todo) > 0:
            with readlock:
                with h5py.File(  os.path.join(ds.datapath, sourcefile), 'r' ) as hin:
                    dset = hin[ ds.limapath ]
                    dt = dset.dtype
                    chunks = []
                    for t in todo:
                        num = obin.flat[address + t] // nomega
                        try:
                            chunks.append( (num, dset.id.read_direct_chunk((t,0,0))[1]) )
                        except:
                            print("FAIL",ds.datapath, sourcefile,t,dset)
                            raise
                print('.' ,end='')
            dc = bslz4_to_sparse.chunk2sparse(image_mask, dt)
            for t,(o,b) in zip(todo, chunks):
                npx, (vals, inds ) = dc( b, 0 )
                accumulate( npx, vals, inds, destination, o )
            print(',', end='')
        address += num_frames
    print('/',end='')
    return destination

    
def makecmdline( 
    wvln = 12.3985 / 43.569,
    xbeam = 1024,
    ybeam = 1024,
    omegastart = 90,
    step = 0.1,
    distance = 140,
    savepath = '.',
    run = 0,
    ):
    return "eiger2crysalis /dev/null -w %f --beam %.2f %.2f --omega=-%d+index*%f --distance=%.3f -o %s/esperanto/frame_%d_{index}.esperanto" %(
        wvln, xbeam, ybeam, omegastart, step, distance, savepath, run)


def sum_sinogram( dsname, output_name, nomega = 10, nthreads = None ):
    if nthreads is None:
        nthreads = min( max( ImageD11.cImageD11.cores_available() - 2, 1 ), 20 )
    ds = ImageD11.sinograms.dataset.load( dsname )
    with h5py.File(dsname ,'r') as hin:
        image_mask = fabio.open(hin['lima_segmenter'].attrs['maskfile']).data
    # ds.import_nnz()
    print("Going to sum sinogram of", dsname, 'output to', output_name)
    with concurrent.futures.ThreadPoolExecutor( max_workers=nthreads ) as pool:
        args = []
        for i in range(nthreads):
            rows = np.zeros( ds.shape, bool )
            rows[i::nthreads] = True
            args.append( (ds, image_mask, rows, nomega) )
    
        def fun(args):
            localds, image_mask, sinomask, nomega  = args
            return integrate_sino( localds, image_mask, sinomask=sinomask, nomega=nomega )
    
        spx = None
        for ans in pool.map(fun, args):
            if spx is None:
                spx = ans
            else:
                padd( ans, spx )
    # fill mask?
    mval = pow(2,32)-1
    for i,frm in enumerate(spx):
        spx[i] = np.where( image_mask, spx[i], mval )
     
    with h5py.File( output_name, 'w' ) as hout:
        ds = hout.create_dataset('data', 
                                 shape = spx.shape, 
                                 dtype = spx.dtype,
                                 data = spx,
                                 chunks = (1,spx.shape[1],spx.shape[2]),
                                 **hdf5plugin.Bitshuffle(nelems=0, lz4=True)
                                )
    return spx


def makecmdline( 
    wvln = 12.3985 / 43.569,
    xbeam = 1024,
    ybeam = 1024,
    omegastart = 90,
    step = 0.1,
    distance = 140,
    savepath = '.',
    run = 1,
    ):
    print("cp /data/id11/nanoscope/Eiger/frame_1_.set %s/frame.set"%(savepath))
    return "eiger2crysalis sum_lpcmoX1_slice1.h5 -w %f --flip-lr --calc-mask False --beam %.2f %.2f --omega=-%d+index*%f --distance=%.3f -o %s/frame_%d_{index}.esperanto" %(
        wvln, xbeam, ybeam, omegastart, step, distance, savepath, run)

if __name__=="__main__":
    sample = 'lpcmo_x3'
    dset = 'z50_RT'
    dsname = '../ds_{sample}_{dset}.h5'.format(**locals())
    guessroot = os.getcwd()
    if guessroot.startswith( '/data' ):
        items = guessroot.split('/')
        id11 = items.index('id11')
        dataroot = "/".join( items[:id11+1] ) + "/RAW_DATA"
    print(dataroot)
    if not os.path.exists(dsname):
        ds = ImageD11.sinograms.dataset.DataSet( 
            dataroot = '/data/visitor/hc5185/id11/20230505/RAW_DATA',
            analysisroot = os.getcwd(),
        sample = sample,
        dset = dset, )
        ds.import_all()
        ds.save(dsname)
        ImageD11.sinograms.lima_segmenter.setup(dsname)
        
    spx = sum_sinogram( dsname, "sum_{sample}_{dset}.h5".format(**locals()), nthreads = 20)

    # Dummy file opener
    def myopenimage( filename ):
        #flake8: global spx
        '''flips to reverse order - need to fix this properly ... '''
        return SinoScan( spx[::-1] )

    # replace the file opener by our thing
    fabio.app.eiger2crysalis.fabio_open = myopenimage

    cmd = makecmdline(savepath='lpcmoX1_slice1', step=0.5, )
    sys.argv = cmd.split()
    fabio.app.eiger2crysalis.main()