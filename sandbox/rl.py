#!/usr/bin/env fable.python
import numpy as np



def centre_psf( psf_array, target_image_dims, DEBUG=False ):
    """
    Shifts the centre of mass of the point spread function to (0,0)
    and copies it into a zero padded array of size target_image_dims

    Function is normalised to have integral of 1.0
    """
    assert psf_array.shape[0] < target_image_dims[0]
    assert psf_array.shape[1] < target_image_dims[1]    
    proj_0 = psf_array.sum(axis = 1, dtype=np.float)
    n0 = len(proj_0)/2
    com_0 = n0 + (proj_0*np.arange(-n0, len(proj_0)-n0)).sum() / proj_0.sum()
    com_0 = int(com_0)
    proj_1 = psf_array.sum(axis = 0, dtype=np.float)
    n1 = len(proj_1)/2
    com_1 = n1 + (proj_1*np.arange(-n1, len(proj_1)-n1)).sum() / proj_1.sum()
    com_1 = int(com_1)
    output = np.zeros( target_image_dims, np.float )
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
    return output/output.sum()

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

def score( obs, calc):
    """ See how well our model matches the data """
    diff = (obs-calc)
    d2 = (diff*diff).ravel()
    print "score",np.sqrt(d2.mean())

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
        for i in range(niter):
            if self.debug: print "RL: iter",i,
            calc = self.convolver.convolve( perfect )
            if self.debug: score( obs, calc )
            factor = obs / calc
            correction = self.convolver.convolve_transp(factor)
            perfect = perfect * correction
        rescaled = undo_scale(perfect, scale, offset)
        return rescaled


def gauss2d( dims, sig):
    """ 2D Gaussian """
    x = np.outer( np.ones(dims[0], np.float), np.arange(dims[1]))-dims[1]/2
    y = np.outer( np.arange(dims[0]), np.ones(dims[1], np.float))-dims[0]/2
    arg =  x * x + y * y
    return np.exp( -arg / sig / sig )



def dodeconv( filename_in,
              filename_out,
              pointspread_image,
              niter = 20):

    from fabio.openimage import openimage
    im_in = openimage( filename_in)
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
    print name,ar.min(), ar.max(), ar.sum()/ar.shape[0]

    
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
    print sm(c), sm(d), sm(gauss)
    figure(2)
    imshow( c)
    colorbar()
    raw_input("press a key")


def testg5():

    from fabio.openimage import openimage
    im = openimage("c_Al_s1_000__quantix_0037.edf")
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
        print sys.argv[0],"infile outfile [psfile|Gaussian_sig] [niterations=20]"
        sys.exit()

    filename_in = sys.argv[1]
    filename_out = sys.argv[2]
    filename_pointspread = sys.argv[3]
    try:
        niter = int(sys.argv[4])
    except:
        niter = 20
        
    if not os.path.exists(filename_in): usage()
    if os.path.exists(filename_pointspread):
        psf =  openimage( filename_pointspread ).data
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
              filename_out,
              psf, 
              niter)

    



