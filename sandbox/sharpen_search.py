
from __future__ import division, print_function
import numpy as np, scipy.ndimage, sys, fabio
from ImageD11 import cImageD11, labelimage

def gauss_1d( N, center, fwhm ):
    """
    Normalised 1D gaussian peak
    """
    sigma = fwhm / 2.35482
    dx = (np.arange(0, N, 1, dtype=np.float32 ) - center)/sigma
    f = np.exp( -dx*dx / 2 ) / sigma / np.sqrt( 2 * np.pi )
    return f

def lorentz_1d( N, center, fwhm ):
    """
    Normalised 1D Lorentz peak
    """
    dx = (np.arange(0, N, 1, dtype=np.float32 ) - center)
    top = fwhm / 2. / np.pi
    bottom = dx*dx + fwhm*fwhm / 4.
    f = top / bottom
    return f




def add_gaussian_2d( img, center, sigma, weight ):
    """ adds a 2d Gaussian into a test image 
    (missing sig_ij for now)
    """
    ns, nf = img.shape
    i, j = np.mgrid[0:ns, 0:nf]
    di = i.astype(np.float32) - center[0]
    dj = j.astype(np.float32) - center[1]
    s2_2 = 2*sigma*sigma
    arg = (di*di+dj*dj)/s2_2
    np.add( img,
            np.exp( -arg )/s2_2/np.pi,
            img )
    return img
    
# gaussian_2d = np.outer( gauss_1d, gauss_1d )
#
# Separation will be interesting for 2theta, azimuth representations...
#

def make_test_image( dims  = (256,256),
                     sigmas = [3.5,],
                     centers = [(128,128),],
                     weights = [1,] ):
    """ produce a test image with some peaks in it """
    imt = np.zeros( dims, np.float32 )
    for (i,j),sigma,weight in zip(centers, sigmas, weights):
        imt = add_gaussian_2d( imt, (i,j), sigma, weight )
    return imt



class convolutepeak(object):
    """
    Gets ready to convolve your image with some functions
       pk - peak width 
       bg - flat
       peak - function to produce 1d peak 
    Uses these convolutions to do a pixel by pixel peak fit (of a local box)
     ... e.g. attempts to strip the background
    """
    def __init__(self, dims, fwhm, peak=gauss_1d ):
        """
        Fits y = pk + bkg    pixel by pixel...
        """
        # Dimension of pk and bg images to convolve
        N = int(3*fwhm+1)*2+1
        p = peak( N, N//2, fwhm )
        o = np.ones(N, np.float32 )
        # 2D peak
        pk = np.outer( p, p )
        bg = np.outer( o, o )
        # Matrix for a least squares problem fitting peak in this box
        m = [ [ (bg*bg).sum(), (bg*pk).sum()],
              [ (bg*pk).sum(), (pk*pk).sum()] ]
        detM = 1./(m[0][0]*m[1][1] - m[0][1]*m[1][0])
        invM = [[ detM*m[1][1],  -detM*m[0][1]],
                [ -detM*m[1][0], detM*m[0][0]]]
        # convolve me to fit peak
        pkc = invM[1][0]*bg + invM[1][1]*pk
        # ...    pk = np.outer( pk1d, pk1d )
        # ...    bg = np.outer( ones, ones )
        # p0 = invM[1][0] * np.ones( N, np.float32 ) + \
        #      invM[1][1] * gauss_1d( N, N//2, sigma )
        #    print( invM )
        #    return p0
        self.p, self.o, self.invM = p, o, invM
        
    def conv( self, im ):
        """ Fits pk + bg """
        p, o, invM = self.p, self.o, self.invM
        # 2D as a pair of 1D
        # pk
        imp = scipy.ndimage.correlate1d( im, invM[1][1] * p,
                                         axis=0, mode= 'nearest' )
        imp = scipy.ndimage.correlate1d( imp, p,
                                         axis=1, mode= 'nearest' )
        # bg
        imb = scipy.ndimage.correlate1d( im, invM[1][0] * o,
                                         axis=0, mode= 'nearest' )
        imb = scipy.ndimage.correlate1d( imb, o,
                                         axis=1, mode= 'nearest' )
        # fitted
        return imb + imp


def process_edf( fname, bkg , sigma_threshold, sigma_peak, sptfile ):
    im = fabio.open( fname ) 
    cpk = convolutepeak( im.data.shape, sigma_peak )
    lio = labelimage.labelimage( im.data.shape, sptfile )
    # for i in range(im.nframes):
    i = 0
    while True:
        try:
            im = im.next()
        except:
            break
        data = im.data.astype( np.float32 )
        np.subtract( data, bkg, data )
        pkdata = cpk.conv( data )
        s0 = pkdata.std()*sigma_threshold
        m0 = abs(pkdata) < s0
        s1 = pkdata[m0].std()*sigma_threshold
        m1 = abs(pkdata) < s1
        s2 = pkdata[m1].std()*sigma_threshold
        lio.peaksearch( pkdata, s2, i )
        print( "%s[%d] %d"%( im.filename, i, lio.npk ) )
        lio.mergelast()
        i+=1
    lio.finalise()
        
        
    
    
if __name__ == "__main__":
    edf      = sys.argv[1]
    bkgname  = sys.argv[2]
    bkg      =  fabio.open( bkgname ).data
    sptfile  = edf.replace(".edf",".flt")
    try:
        sigma_threshold = float(sys.argv[3])
    except:
        sigma_threshold = 5.
    sigma_peak = 2.
    try:
        fwhm_peak = float(sys.argv[4])
    except:
        fwhm_peak = 3.
    process_edf( edf, bkg, sigma_threshold , fwhm_peak , sptfile)

    
