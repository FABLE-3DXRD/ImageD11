
from __future__ import print_function
import sys, timeit
import numpy as np #, numba
import scipy.ndimage 
import hdf5plugin, h5py, fabio
from ImageD11 import cImageD11, sparseframe, blobcorrector, \
    parameters, transform
import pylab as pl

# Timer (could be a decorator ...)
timer = timeit.default_timer
class ticktock(object):
    def __init__(self):
        self.tick = timer()
    def tock(self, msg = ""):
        tock = timer()
        print(msg, "%.4f /s"%(tock - self.tick), end=" ")
        self.tick = timer()

# Generic filter, keeps input and out around
class imgfilt(object):
    def __init__(self, *args):
        pass
    def __call__(self, input):
        t = ticktock()
        self.input = input
        self.output = self.filt( input )
        t.tock(self.__class__.__name__)
        return self.output
    def filt( self, input ):
        return input.copy()
        

# Dark subtraction and flat normalisation
# faster with omp
#@numba.njit(parallel=True)
#def _darkflat( inp, dark, flmult, out ):
#    for i in numba.prange(inp.shape[0]):
#        for j in range(inp.shape[1]):
#            out[i,j] = ( inp[i,j] - dark[i,j] ) * flmult[i,j]
            
class darkflat(imgfilt):
    def __init__(self, dark, flat):
        self.dark = dark
        self.flmult = 1./flat
        self.out = np.empty( dark.shape, np.float32 )
        
    def filt(self, input):
        cImageD11.uint16_to_float_darkflm( self.out.ravel(),
                                           self.dark.ravel(),
                                           self.flmult.ravel(),
                                           input.ravel())
        return self.out
    
class uniform(imgfilt):
    def __init__(self, n=3, mode='nearest'):
        self.n = n
        self.mode = mode
    def filt( self, input ):
        return scipy.ndimage.uniform_filter( input, self.n, mode = self.mode )

class gauss_bgsub(imgfilt):
    def __init__(self, sigma, mode='nearest' ):
        self.sigma = sigma
        self.mode = mode
    def filt( self, input ):
        bg = scipy.ndimage.gaussian_filter( input.astype(np.float32),
                                            self.sigma, mode = self.mode )
        return input - bg


class gauss(imgfilt):
    def __init__(self, sigma, mode='nearest' ):
        self.sigma = sigma
        self.mode = mode
    def filt( self, input ):
        bg = scipy.ndimage.gaussian_filter( input.astype(np.float32),
                                            self.sigma, mode = self.mode )
        return bg

class laplace(imgfilt):
    def __init__(self, mode='nearest'):
        self.mode = mode
    def filt( self, input ):
        return abs(scipy.ndimage.laplace( input.astype(np.float32),
                                      mode = self.mode ))

class sobel(imgfilt):
    def __init__(self, mode='nearest'):
        self.mode = mode
    def filt( self, input ):
        s0 = scipy.ndimage.sobel( input.astype(np.float32),
                                      mode = self.mode )
        s1 = scipy.ndimage.sobel( input.astype(np.float32),
                                      mode = self.mode )
        return abs(s0) + abs(s1)
    
class laplace_sobel(imgfilt):
    def __init__(self, mode='nearest'):
        self.mode = mode
    def filt( self, input ):
        s0 = scipy.ndimage.sobel( input.astype(np.float32),
                                      mode = self.mode )
        s1 = scipy.ndimage.sobel( input.astype(np.float32),
                                      mode = self.mode )
        s3 = scipy.ndimage.laplace( input.astype(np.float32),
                                   mode = self.mode )
        return abs(s0) + abs(s1) + abs(s3)

class median( imgfilt ):
    def __init__(self, size=3,  mode='nearest'):
        self.mode = mode
        self.size = size
    def filt( self, input ):
        return scipy.ndimage.median_filter( input, size=self.size, mode=self.mode )

class median_bgsub( imgfilt ):
    def __init__(self, size=3,  mode='nearest'):
        self.mode = mode
        self.size = size
    def filt( self, input ):
        bg = scipy.ndimage.median_filter( input, size=self.size, mode=self.mode )
        return input - bg

class median_footprints( imgfilt):
    def __init__(self, footprints, mode='nearest'):
        assert len(footprints)==2
        self.footprints = footprints
        self.mode=mode
    def filt( self, input ):
        bg0 = scipy.ndimage.median_filter( abs(input), footprint = self.footprints[0],
                                          mode=self.mode )
        bg1 = scipy.ndimage.median_filter( abs(input),  footprint = self.footprints[1],
                                          mode=self.mode )
        bg = np.where( bg0 < bg1, bg0, bg1 )
        den = 1/np.sqrt( bg + 1 ) # simd??
        sig = (input-bg)*den
        return np.where( sig > -3 , sig, -3)




class bgsub( imgfilt ):
    def __init__(self, gain=0.5, sigmap=0.2, sigmat=25 ):
        self.sigmap = sigmap
        self.sigmat = sigmat
        self.gain = gain
        self.bg = None
    def filt(self, input):
        if self.bg is None:
            self.msk = np.empty( input.shape, np.uint8 )
            self.bg = np.empty( input.shape, np.float32 )
        cImageD11.bgcalc( input,
                          self.bg,
                          self.msk,
                          self.gain,
                          self.sigmap,
                          self.sigmat )
        #
        #self.bg = bg1d(input.ravel(), self.gain,
        #               self.sigmap, self.sigmat).reshape(input.shape)
        return input - self.bg
    

class variance( imgfilt ):
    def __init__(self, size=5,  mode='nearest'):
        self.mode = mode
        self.size = size
    def filt( self, input ):
        avg = scipy.ndimage.uniform_filter( input.astype(np.float32),
                                            size=self.size,
                                            mode=self.mode )
        d = input - avg
        v = scipy.ndimage.uniform_filter( d*d    ,
                                          size=self.size,
                                          mode=self.mode )
        return np.sqrt(v+1)

    
class radial(imgfilt):
    def __init__(self, splinefile, pars, dims, npx=0.77, mask=None ):
        t, e  = self.compute_tth_eta_lut( splinefile, pars, dims)
        one_px = (t[1:] - t[:-1]).max()
        binsize = one_px * npx
        itth = np.round(t / binsize)
        end = itth.max() + 1
        if mask is not None: # throw to the end
            itth[ mask ] += end
            self.nmask = mask.sum()
        else:
            self.nmask = 0
        self.output = None
        self.lut, self.adrout = self.make_lut( itth, e )

    def filt( self, input ):
        if self.output is None:
            self.output = np.empty_like( input )
        if input.dtype == np.float32:
            cImageD11.reorder_f32_a32( input.ravel(), self.adrout, self.output.ravel())
        elif input.dtype == np.uint16:
            cImageD11.reorder_u16_a32( input.ravel(), self.adrout, self.output.ravel())
        else:
            self.output.flat = input.flat[self.lut]   
        return self.output
         
    def compute_tth_eta_lut( self, splinefile, pars, dims):
        """
        Computes look up values of tth, eta for each pixel
        """
        c = blobcorrector.correctorclass( splinefile )
        p = parameters.read_par_file( pars )
        xp, yp = c.make_pixel_lut( dims )
        t, e = transform.compute_tth_eta( (xp.ravel(), yp.ravel()),
                                      **p.parameters )
        t.shape=dims
        e.shape=dims
        return t, e

    def make_lut( self, itth, eta ):
        ie = 1./361
        sorter = (itth + (ie*eta)).ravel()
        lut = np.argsort( sorter ).astype( np.uint32 ) # read from here
        adrout = np.argsort( lut ).astype( np.uint32 ) # write to here
        return lut, adrout

class undoradial(radial):
    def __init__(self, aradial):
        self.output = None
        self.lut, self.adrout = aradial.adrout, aradial.lut
        self.nmask = aradial.nmask

    def filt(self, input):
        if self.nmask >0:
            input.flat[-self.nmask:] = 0
        return radial.filt( self, input )
        
    
# Masking with a sigma cutoff
class maskcut(imgfilt):
    def __init__(self, sigma=5):
        self.sigma = sigma
        self.output = None
    def filt( self, input ):
        if self.output is None:
            self.output = np.empty( input.shape, np.uint8)
        sigma = self.sigma
        m, s = cImageD11.array_mean_var_msk(input.ravel(), self.output.ravel(), 
                                            cut=self.sigma,
                                            n = 4,
                                            verbose=0)
        return self.output


class clean(imgfilt):
    def filt( self, input ):
        cln  = np.empty( input.shape, 'b' )
        npx = cImageD11.clean_mask( input, cln )
        return cln

class domask( imgfilt ):
    def __init__(self,mask):
        self.mask = mask
    def filt( self, input):
        return np.where( self.mask, 0, input )
        
    

# one in - one out pipeline
class pipeline(imgfilt):
    def __init__(self, filters):
        self.filters = filters
    def filt(self, input):
        x = input.copy()
        for f in self.filters:
            x = f(x)
        return x
    def __getitem__(self, i):
        if i >= 0 and i < len(self.filters):
            return self.filters[i]
    def show( self ):
        n = len(self.filters)
        fig0, ax0 = pl.subplots(1,n+1,sharex=True, sharey = True , figsize=(5*n,4))
        fig1, ax1= pl.subplots(1,n+1,sharex=True, sharey = False, figsize=(5*n,4))
        im1 = ax0[0].imshow( self.input, origin='lower', vmin=0,
                               vmax = self.input.mean() + self.input.std()*3 )
        pl1 = ax1[0].plot( self.input[1024] )
        ax0[0].set_title(self.__class__.__name__+ " input")
        fig0.colorbar( im1, ax=ax0[0] )
        for i in range(n):
            im = self.filters[i].output
            ims = ax0[i+1].imshow( im, origin = 'lower',
                                   vmin = im.mean() - im.std()*3 ,
                                   vmax = im.mean() + im.std()*3 )
            ax0[i+1].set_title( self.filters[i].__class__.__name__ )
            fig0.colorbar( ims, ax=ax0[i+1] )
            ims = ax1[i+1].plot(im[1024])
        fig0.tight_layout()
        fig1.tight_layout()
        fig0.show()


def segment_scan(fname, scan, detector, outname,
                 scanmotors = [],
                 headermotors = []):
    with h5py.File(fname,"r") as hin:
        gin = hin[scan]
        hdr = {}
        with h5py.File( outname, "a") as hout:
            g = hout.create_group( scan )
            for m in scanmotors: # vary : many
                print(m)
                g.create_dataset(m, data = gin['measurement'][m][:] )
            for m in headermotors: # fixed : scalar
                g.create_dataset(m, data = gin['instrument/positioners'][m][()] )
                
            opts = { 'chunks' : (10000,),
                     'maxshape' : (None,),
                     'dtype': np.uint16,
                     'compression':'lzf',
                     'shuffle': True }
            row = g.create_dataset( 'row', (1,), **opts )
            col = g.create_dataset( 'col', (1,), **opts )
            num = g.create_dataset( 'frame', (1,),  **opts )
            opts['dtype'] = np.float32
            sig = g.create_dataset( 'intensity', (1,), **opts )
            frm = gin['measurement'][detector]
            nnz = g.create_dataset( 'nnz', (frm.shape[0],), dtype=np.int32)
            npx = 0
            g.attrs['itype']=np.dtype(np.uint16).name
            g.attrs['nframes']= frm.shape[0]
            g.attrs['shape0'] = frm.shape[1]
            g.attrs['shape1'] = frm.shape[2]            
            for i,frm in enumerate(frm):
                print("%4d"%(i), end=" ")
                msk = mypipeline(frm.astype(np.float32))
                spf = sparseframe.from_data_mask( msk, frm, {} )
                print("%d px"%(spf.nnz))
                if spf.nnz + npx > len(row):
                    row.resize( spf.nnz + npx, axis=0 )
                    col.resize( spf.nnz + npx, axis=0 )
                    sig.resize( spf.nnz + npx, axis=0 )
                    num.resize( spf.nnz + npx, axis=0 )
                row[npx:] = spf.row
                col[npx:] = spf.col
                sig[npx:] = spf.pixels['intensity']
                num[npx:] = i
                nnz[i] = spf.nnz
                npx += spf.nnz
                # spf.to_hdf_group( g.create_group( "%d"%(i) ))

def check( fname, scan, num ):
    with h5py.File( fname, "r") as hin:
        frm = hin[scan][num]
        for i in range(5):
            print("\n",i,end=" ")
            ret = mypipeline( frm )
        sys.stdout.flush()
        mypipeline.show()


#import multiprocessing
#N = multiprocessing.cpu_count() // 2
#cImageD11.cimaged11_omp_set_num_threads( N )
#print("Using",cImageD11.cimaged11_omp_get_max_threads(  ),"threads")

drk = fabio.open("median_image.edf").data
flt = fabio.open("F4M_flood_Ce_aug2010.edf").data
flt = np.where( flt < 0.55 , 1, flt )
mask = fabio.open("mask.edf").data > 0

#m = median_footprints( footprints=[np.ones((9,1)), np.ones((1,9))] )
#df = m( flt ) - flt
#pl.imshow(df*df)
#pl.show()

rad = radial( "frelon4m.spline", "f.par", dims=drk.shape, mask=mask )
unrad = undoradial( rad )


mypipeline = pipeline(
    [ darkflat( dark = drk, flat = flt ),
      bgsub(),
      maskcut( sigma = 10 ),
      clean(),  ]
    )


fname = "fe3o4map01.hdf"
scan = 'entry' #  '%d.1'% iscan # /measurement/frelon3'%(iscan)
outname = "fe3o4_peaks.hdf"
camera = "frelon4m"

if 1:
    segment_scan( fname, scan, camera, outname,
                  scanmotors = ('omega', ),
                  headermotors = () )
else:
    check( fname,
           scan + "/measurement/" + camera,
           0 )
    pl.show()
        
