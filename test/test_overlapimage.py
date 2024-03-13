
from __future__ import print_function, division
import sys
import numpy as np, fabio, time, os
from ImageD11 import cImageD11
#print(cImageD11.__file__)

OVIM = os.path.join(
    os.path.dirname(__file__),
    "testoverlaps0000.edf" )

def maketestim(): 
           
    M = 2048   # Overall image size
    N = 64     # tile size
    SIG = 3.   # sigma for peaks (width)
    # tile co-ordinates
    i,j = np.mgrid[0:N,0:N]-N/2.

    def pk(x,y,s):
        # 2D peak shape (with tails)
        dx = i-x
        dy = j-y
        gg= np.exp( -(dx*dx+dy*dy)/s/s )
        ll=s*s/(s*s+dx*dx+dy*dy)
        return ll+gg

    # Image for testing
    im = np.zeros( (M,M) )

    # Make a grid of points (tiles)
    pks = []
    for ii in range(0,M,N):
        for jj in range(0,M,N):
            # offset of centroids
            oi = 3*SIG*ii/M
            oj = SIG*ii/M
            # Ratio of heights
            sc = jj*1.0/M
            # put in two peaks, offset increases one direction,
            # scale in the other
            im[ii:ii+N,jj:jj+N] = pk( -oi, -oj, SIG ) + pk( oi, oj, SIG)*sc
            pks.append( (ii+N//2-oi,jj+N//2-oj) )
            pks.append( (ii+N//2+oi,jj+N//2+oj) )
    if 0:
        import pylab as pl
        pl.imshow( pl.log(im ))
        pl.show()

    s = 30000/im.max()
    int_im = (im*s).astype( np.uint16 )
    if not os.path.exists( OVIM ):
        fabio.edfimage.edfimage( int_im ).write( OVIM )
    return int_im, pks

def bench(int_im, pks):
    im = int_im
    float_im=int_im.astype(np.float32)
    labelc = np.zeros( im.shape, np.int32 )
    labelm = np.zeros( im.shape, np.int32 )
    work   = np.zeros( im.shape, np.int8 )

    end0 = time.time()
    nc    = cImageD11.connectedpixels( float_im, labelc, 1000 )
    end1   = time.time()
    rc    = cImageD11.blobproperties( float_im, labelc, nc, 0 )
    cImageD11.blob_moments(rc)
    end2   = time.time()
    print("cptime",nc,"%.3f %.3f"%((end1-end0)*1000,(end2-end1)*1000))

    wfloatim = np.where( float_im > 1000, float_im, 0)
    #import scipy.ndimage as ndi
    #wfloatim= (float_im - ndi.gaussian_filter( float_im , 32)).clip(0,1e9)

    end0 = time.time()
    nw    = cImageD11.localmaxlabel( wfloatim, labelm, work )
    end1   = time.time()
    rw    = cImageD11.blobproperties( wfloatim, labelm, nw, 0 )
    cImageD11.blob_moments(rw)
    end2   = time.time()
    print("lmtime",nw,"%.3f %.3f"%((end1-end0)*1000,(end2-end1)*1000))
    if "plot" in sys.argv:
        pks=np.array(pks).T
        import pylab as pl
        pl.imshow( pl.log(float_im), origin='lower',
               interpolation='nearest', aspect='auto') 
        pl.plot( pks[1], pks[0],"r+",label="ideal" )
        pl.plot( rc[:,cImageD11.f_raw],rc[:,cImageD11.s_raw],"wx",
        label="connect")
        pl.plot( rw[:,cImageD11.f_raw],rw[:,cImageD11.s_raw],"ko",
             label="localmax",
             markerfacecolor='None')
        pl.legend()
        pl.figure()
        dci= [ np.sqrt(np.min( (pks[1]-x)**2 + (pks[0]-y)**2 ))
           for x,y in zip( rc[:, cImageD11.f_raw],rc[:, cImageD11.s_raw])]
        pl.subplot(121)
        pl.plot( rc[:, cImageD11.s_raw],
                dci , '.', label= "sconnected")
        pl.plot( rc[:, cImageD11.f_raw],
                dci , '.', label= "fconnected")
        pl.ylabel("position error")
        pl.title("Connected pixels = 1 blob")
        pl.legend()
        dwi= [ np.sqrt(np.min( (pks[1]-x)**2 + (pks[0]-y)**2 ))
           for x,y in zip( rw[:, cImageD11.f_raw],rw[:, cImageD11.s_raw])]
        pl.subplot(122)
        pl.plot( rw[:, cImageD11.s_raw],
                dwi , '.', label= "slocalmax")
        pl.plot( rw[:, cImageD11.f_raw],
                dwi , '.', label= "flocalmax")
        pl.title("Localmax = split blobs")
        pl.legend()
        pl.show()


import unittest
class testimexists(unittest.TestCase):
    def setUp(self):
        maketestim()
    def testit(self):
        self.assertTrue(os.path.exists(OVIM))

if __name__=="__main__":
    bench(*maketestim())
    # unittest.main()
    # this is not a unit test !
    
