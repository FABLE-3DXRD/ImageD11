
from __future__ import division, print_function
import numpy as np, time
from ImageD11.cImageD11 import localmaxlabel 
import unittest

def localmax( im, out=None ):
    """ reference implementation 
    ... very slow 
    """
    neighbours = [ (-1,-1),(-1,0), (-1,1),
                   (0,-1),         (0,1),
                   (1,-1), (1,0),  (1,1) ]
    out = np.arange( im.shape[0]*im.shape[1], dtype=np.int )
    out.shape = im.shape
    out0 = out.copy()    # pointers
    print("running python localmax")
    for i in range(1,im.shape[0]-1):
        for j in range(1,im.shape[1]-1):
            mx = im[i,j]
            for k,l in neighbours:
                if im[i+k,j+l] > mx:
                    out[i,j] = out0[i+k,j+l]
                    mx = im[i+k,j+l]
            # now out addresses a neighbour
    # Maximal pixels
    maxpx = (out==out0)
    # Final labels
    npk = 0
    labels = np.zeros( im.shape, np.int32 )-1
    labels[0]=0   # borders
    labels[-1]=0
    labels[:,0]=0
    labels[:,-1]=0
    for i in range(1,im.shape[0]-1):
        for j in range(1,im.shape[1]-1):
            if maxpx[i,j]:
                npk += 1
                labels[i,j]=npk
                out[i,j]=0
    # Now relabel each pixel to point to the max...
    for i in range(1,im.shape[0]-1):
        for j in range(1,im.shape[1]-1):
            adr = out[i,j]
            if adr == 0:
                continue
            while 1:
                adrnew = out.flat[adr]
                if adrnew == 0:
                    break
                adr = adrnew
            l = labels.flat[adr] # take from peak
            labels[i,j] = l
            adr = out[i,j]
            while 1:
                labels.flat[adr] = l
                adrnew = out.flat[adr]
                if adrnew == 0:
                    break
                out.flat[adr] = 0
                adr = adrnew
    return labels


def make_test_image( N=256, w=3 ):
    im = np.zeros( (N, N), np.float32 )
    # patch size to hold peaks
    pN = w*5+1  
    pk = np.zeros( (pN,pN), np.float32 )
    pi,pj = np.mgrid[ 0:pN, 0:pN ]
    # two peaks, separated by i
    ni = N//pN
    s = 5.
    for i in range(0, ni):
        # at max separation we want 0 - pk - pk - pN
        # so:
        dj = (pN/3.0)*(i/(ni-1))
        di =  dj / 5.

        ci1 = pN//2  + di
        ci2 = pN//2  - di
        cj1 = pN//2  + dj
        cj2 = pN//2  - dj
        p1 = np.exp( - (( pi - ci1 )**2 + (pj - cj1)**2)/s)
        p2 = np.exp( - (( pi - ci2 )**2 + (pj - cj2)**2)/s)
        for j in range(0, N//pN):
            im[ i*pN:(i+1)*pN, j*pN:(j+1)*pN ] = p1 + p2*(j+1)/ni
    return im

def demo():
    import pylab as pl
    im = make_test_image()
    l1 = localmax(im)
    l2 = np.zeros(im.shape,np.int32)
    wk = np.zeros(im.shape,np.int8)
    npks = localmaxlabel(im, l2 , wk )
    ax1 = pl.subplot(2,2,1)
    pl.imshow(im)
    pl.xlim( 0,im.shape[0])
    pl.ylim( 0,im.shape[1])
    ax1 = pl.subplot(2,2,2, sharex=ax1, sharey=ax1)
    pl.imshow((l1*13.0)%37)
    pl.xlim( 0,im.shape[0])
    pl.ylim( 0,im.shape[1])
    ax1 = pl.subplot(2,2,3, sharex=ax1, sharey=ax1)
    pl.imshow((l2*13.0)%37)
    pl.xlim( 0,im.shape[0])
    pl.ylim( 0,im.shape[1])
    pl.show()

class test_localmaxlabel( unittest.TestCase ):
    def setUp( self ):
        pass
    def test1( self):
        im = make_test_image()
        start = time.time()
        l1 = localmax(im)
        end = time.time()
        pytime = (end-start)*1000.
        l2 = np.zeros(im.shape,np.int32)
        wk = np.zeros(im.shape,np.int8)
        start = time.time()
        npks = localmaxlabel(im, l2 , wk )
        end=time.time()
        ctime = (end-start)*1000.
        self.assertEqual( (l2==l1).all(), True )
        print("Timing python %.3f ms, c %.3f ms"%(pytime,ctime))

if __name__== "__main__":

    unittest.main()

