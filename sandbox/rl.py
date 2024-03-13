
from __future__ import print_function
#!/usr/bin/env fable.python
import numpy as np
from six.moves import input
try:
    import fftw3f
    REAL = np.float32


    class convolver(object):

        plans = {}

        def __init__(self, pointspread, dims):
            """
            dims = target dimensions
            pointspread = pointspread function to use
            """
            self.pointspread = pointspread
            self.dims = dims

            self.N = self.dims[0] * self.dims[1]
        
            self.PSF = centre_psf( pointspread, dims )
    
            self.QPSF = self._do_fft( self.PSF, 'forward' )
            self.eps = 1.0/65535

        def release_plan(self, input):
            """
            release the memory from previous plans
            """
            if input.ctypes.data in self.plans:
                # no idea if it leaks elsewhere
                del self.plans[ input.ctypes.data ]


        def _do_fft( self, input, direct):
            """
            Does the fft using fftw3f
            It caches the plan and result of the fft so that when you 
            do it again it will be faster
            This means you are limited to only using a few arrays before 
            you exhaust the machines memory
            see also release_plan
            """
            key = input.ctypes.data
            if key in self.plans:
                output , plan = self.plans[key]
                plan.execute()
                return output
            # we've never done an fft of this array
            if len(list(self.plans.keys()))>6:
                print("fftw3 is making too much cache")
                print("needs to be programmed to reuse the same arrays")
            if direct == 'forward':
                assert input.dtype == np.float32
                output = np.zeros( input.shape, np.complex64)
            else:
                assert input.dtype == np.complex64
                output = np.zeros( input.shape, np.float32 )
            fp = fftw3f.Plan( 
                    inarray = input, 
                    outarray = output, 
                    direction = direct, 
                    flags = ['estimate'] )
            fp.execute()
            self.plans[ input.ctypes.data ] = output, fp
            return output

        def convolve(self, image):
            """
            Does a convolution
            Currently using wraparound - which should probably be changed to clamp
            """
            assert image.shape == self.dims
            qim = self._do_fft( image, 'forward' )
            np.multiply( qim ,self.QPSF, qim )
            out = self._do_fft( qim, 'backward' )
            np.divide( out, self.N, out)
            return out

        def convolve_transp(self, image):
            """
            Does a convolution
            Currently using wraparound - which should probably be changed to clamp
            """
            assert image.shape == self.dims
            qim = self._do_fft( image, 'forward' )
            np.multiply( qim, self.QPSF.conj(), qim)
            out = self._do_fft( qim, 'backward' )
            np.divide( out, self.N, out)
            return out

  #  raise ImportError()

    print("Using fftw3f, should perform better")

except ImportError:

    print("Using numpy fft for convolution")
    print("You might get better performance from fftw, why not try installing:")
    print("http://prdownload.berlios.de/pyfftw/PyFFTW3-0.2.tar.gz")
    REAL = np.float32

    class convolver(object):

        def __init__(self, pointspread, dims):
            """
            dims = target dimensions
            pointspread = pointspread function to use
            """
            self.pointspread = pointspread
            self.dims = dims
            
            self.PSF = centre_psf( pointspread, dims )
            self.QPSF = np.fft.rfft2( self.PSF )

            self.eps = 1.0/65535

        def convolve(self, image):
            """
            Does a convolution
            Currently using wraparound - which should probably be changed to clamp
            """
            assert image.shape == self.dims
            qim = np.fft.rfft2( image )
            return np.fft.irfft2( qim * self.QPSF )

        def convolve_transp(self, image):
            """
            Does a convolution
            Currently using wraparound - which should probably be changed to clamp
            """
            assert image.shape == self.dims
            qim = np.fft.rfft2( image )
            return np.fft.irfft2( qim * self.QPSF.conj() )





def centre_psf( psf_array, target_image_dims, DEBUG=False ):
    """
    Shifts the centre of mass of the point spread function to (0,0)
    and copies it into a zero padded array of size target_image_dims

    Function is normalised to have integral of 1.0
    """
    assert psf_array.shape[0] < target_image_dims[0]
    assert psf_array.shape[1] < target_image_dims[1]    
    proj_0 = psf_array.sum(axis = 1, dtype=REAL)
    n0 = len(proj_0)/2
    com_0 = n0 + (proj_0*np.arange(-n0, len(proj_0)-n0)).sum() / proj_0.sum()
    com_0 = int(com_0)
    proj_1 = psf_array.sum(axis = 0, dtype=REAL)
    n1 = len(proj_1)/2
    com_1 = n1 + (proj_1*np.arange(-n1, len(proj_1)-n1)).sum() / proj_1.sum()
    com_1 = int(com_1)
    output = np.zeros( target_image_dims, REAL )
    # 
    # 4 corners:
    #
    p0, p1 = psf_array.shape
    output[:p0-com_0, :p1-com_1 ] = psf_array[com_0:, com_1:]
    if com_0 > 0:
        output[:p0-com_0,  -com_1:] = psf_array[com_0:, :com_1]
    if com_1 > 0:
        output[ -com_0:, :p1-com_1] = psf_array[:com_0, com_1:]
    if com_0 > 0 and com_1 > 0:
        output[ -com_0:, -com_1:] = psf_array[:com_0, :com_1]
    return (output/output.sum()).astype(REAL)

def do_scale(image, scale, offset):
    """ scales and offsets an image """
    return (image+offset)*scale

def undo_scale(image, scale, offset):
    """ recovers image undoing scale and offset """
    return image/scale - offset

def prescale(image, tmin = 1./65535, tmax = 1.0):
    """
    Set image to be between tmin and tmax
    """
    mx = image.max()*1.0
    mn = image.min()*1.0
    if (mx-mn) == 0:
        scale = 1
    else:
        scale = 1.0*(tmax - tmin)/(mx - mn)
    offset = tmin - mn
    scale = 1.0*(tmax - tmin)/(mx - mn)
    return do_scale( image, scale, offset), scale, offset




def score( obs, calc):
    """ See how well our model matches the data """
    diff = (obs-calc)
    d2 = (diff*diff).ravel()
    print("score",np.sqrt(d2.mean()), end=' ')

class RichardsonLucy(object):
    def __init__(self, conv, debug=True):
        """
        conv is a convolver object
        """
        self.convolver = conv
        self.debug = debug
        
    def deblur( self, image, niter = 10 ):
        """
        Runs the algorithm to deblur an image
        """
        assert image.shape == self.convolver.dims
        obs, scale, offset = prescale(image) # obs [0->1]
        perfect = obs.copy()
        import time
        calc = obs.copy()
        np.seterr('warn')
        for i in range(niter):
            start = time.time()
            if self.debug: print("RL: iter",i, end=' ')
            calc = self.convolver.convolve( perfect )
            if self.debug: score( obs, calc )
            calc[:] = np.where( calc < 1./65535, 1./65535, calc)
            np.divide( obs, calc, calc ) 
            correction = self.convolver.convolve_transp( calc )
            np.multiply( perfect, correction, perfect)
            print("%.3f"%(time.time()-start))
        rescaled = undo_scale(perfect, scale, offset)
        return rescaled


def gauss2d( dims, sig):
    """ 2D Gaussian """
    x = np.outer( np.ones(dims[0], REAL), np.arange(dims[1]))-dims[1]/2
    y = np.outer( np.arange(dims[0]), np.ones(dims[1], REAL))-dims[0]/2
    arg =  x * x + y * y
    return np.exp( -arg / sig / sig )



def dodeconv( filename_in,
              filename_bg,
              extra_bg,
              filename_out,
              pointspread_image,
              niter = 20):

    from fabio.openimage import openimage

    im_in = openimage( filename_in)
    im_bg = openimage( filename_bg)

    # do background subtraction
    im_in.data = im_in.data.astype(REAL) - im_bg.data - float(extra_bg)
    im_in.data = np.where( im_in.data > 0 , im_in.data, 0 )
    # conv does the image convolution
    conv = convolver( pointspread_image, im_in.data.shape )
    # rl does the deconvolution
    rl = RichardsonLucy( conv )
    # 
    im_in.data = rl.deblur( im_in.data, niter = niter ).astype( np.uint16 )
    im_in.write( filename_out )




############# test stuff

def roughly(x, target):
    return abs(x - target) < 1e-10

def pa(a, name):
    """ debugging routine to print an array, name
    and summary info
    """
    ar = a.ravel()
    print(name,ar.min(), ar.max(), ar.sum()/ar.shape[0])

    
def test_centre_psf( ):
    """ check it does what we want for shifting to 0,0 """
    ps = np.zeros((10,10))
    ps[0,0] = 1
    ret = centre_psf( ps, (100 ,100))
    assert roughly(np.ravel(ret).sum(),1)
    assert roughly(ret[0,0], 1)

    ps = np.zeros((10,10))
    ps[5,5] = 1
    ret = centre_psf( ps, (12 ,15))
    assert roughly(np.ravel(ret).sum(), 1)
    assert roughly(ret[0,0], 1)
    
    ps = np.zeros((7,9))
    ps[4,5] = 8
    ps[5,5] = 10
    ps[6,5] = 8
    ps[5,4] = 4
    ps[5,6] = 4
    ret = centre_psf( ps, (30 ,40))
    assert roughly(np.ravel(ret).sum(), 1.0)
    assert roughly(ret[0,0], 10./(8+10+8+4+4))
    

def test_prescale():
    im = np.zeros((10,10))
    im[5,5]=10
    nim, scale, offset = prescale(im, tmin = 1./65535)
    assert (nim[0,0] - 1./65535)<1./65535 , nim[0]
    assert (nim[5,5] - 1.0) < 1./65535

def testg1():
    from pylab import figure, imshow, colorbar, show
    x = np.arange(20)
    gauss = np.exp( -(x*x+x.T*x.T)/10 )
    dims = (1526,1024)
    ret = centre_psf( gauss, dims )
    from fabio.openimage import openimage
    im = openimage("c_Al_s1_000__quantix_0037.edf")
    
    d = im.data
    figure(1)
    show()
    imshow(d)
    colorbar()
    show()
    ret = centre_psf(  gauss , d.shape)
    rl = RichardsonLucy( gauss, d.shape)
    c = rl.convolve( d )
    def sm(x): return np.ravel(x).sum()
    print(sm(c), sm(d), sm(gauss))
    figure(2)
    imshow( c)
    colorbar()
    input("press a key")


def testg5():

    from fabio.openimage import openimage
    im = openimage("c_Al_s1_000__quantix_0037.edf")
    g = gauss2d
    gaus = g( (50,60), 4.1 )+ g( (50,60), 10.1 )*0.1
    from matplotlib.pylab import imshow, show
    imshow(gaus)
    show()
                                              
    rl = RichardsonLucy( gaus, im.data.shape )
    stuff = rl.deblur( im.data.copy() )
    im.data = stuff
    im.write("deblurred.edf",force_type=np.float32)

def run_tests():
    test_centre_psf()
    test_prescale()
    # testg1()
    testg5()



    
if __name__=="__main__":
    # run_tests()
    
    import sys, os
    from fabio.openimage import openimage

    def usage():
        print(sys.argv[0],"infile bgfilw extra_bg outfile [psfile|Gaussian_sig] [niterations=20]")
        sys.exit()

    filename_in = sys.argv[1]
    filename_bg = sys.argv[2]
    extra_bg = sys.argv[3]
    filename_out = sys.argv[4]
    filename_pointspread = sys.argv[5]
    try:
        niter = int(sys.argv[6])
    except:
        niter = 20
        
    if not os.path.exists(filename_in): usage()
    if not os.path.exists(filename_bg): usage()
    if os.path.exists(filename_pointspread):
        psf =  openimage( filename_pointspread ).data.astype(REAL)
        psfr = psf.ravel()
        # Scale, normalise
        psf = (psf - psfr.min())/psfr.sum()
    else:
        try:
            sig = float( filename_pointspread )
        except:
            usage()
        dims = sig*10, sig*10
        psf = gauss2d( dims, sig )
        
    dodeconv( filename_in,
              filename_bg,
              extra_bg,
              filename_out,
              psf, 
              niter)

    



