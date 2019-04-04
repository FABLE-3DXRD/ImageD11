
from __future__ import division, print_function
import numpy as np, time
from ImageD11.cImageD11 import localmaxlabel
import unittest, sys

def localmax( im, out=None ):
    """ reference implementation 
    ... very slow 
    """
    neighbours = [ (-1,-1),(-1,0), (-1,1),
                   (0,-1), (0,0),   (0,1),
                   (1,-1), (1,0),  (1,1) ]
    out = np.arange( im.shape[0]*im.shape[1], dtype=np.int )
    out.shape = im.shape
    out0 = out.copy()    # pointers
    print("running python localmax", im.shape,end=" " )
    for i in range(1,im.shape[0]-1):
        for j in range(1,im.shape[1]-1):
            mx = -1e10
            for k in range(i-1,i+2):
                for l in range(j-1,j+2):
                    try:
                        if im[k,l] > mx:
                            out[i,j] = out0[k,l]
                            mx = im[k,l]
                    except:
                        print(i,j,k,l,im.shape)
                        raise
            # now out addresses a neighbour
    print("1", end=" " )
    # Maximal pixels
    maxpx = (out==out0)
    # Final labels
    npk = 0
    labels = np.zeros( im.shape, np.int32 )-1
    labels[0]=0   # borders
    labels[-1]=0
    labels[:,0]=0
    labels[:,-1]=0
    print("2", end=" " )
    for i in range(1,im.shape[0]-1):
        for j in range(1,im.shape[1]-1):
            if maxpx[i,j]:
                npk += 1
                labels[i,j]=npk
                out[i,j]=0
    print("3", end=" " )
    sys.stdout.flush()
    # Now relabel each pixel to point to the max...
    pon=False
    for i in range(1,im.shape[0]-1):
        for j in range(1,im.shape[1]-1):
            adr = out[i,j]
            if adr == 0:
                continue
            while 1:
                adrnew = out.flat[adr]
                if adrnew == 0 or labels.flat[adr]==0:
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
    print(".", end=" " )                            
    print ("ok")
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
    npks = localmaxlabel(im, l2 , wk, cpu=0 )
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
        im = make_test_image(256)
        start = time.time()
        l1 = localmax(im)
        end = time.time()
        pytime = (end-start)*1000.
        l2 = np.zeros(im.shape,np.int32)
        l3 = np.zeros(im.shape,np.int32)
        wk = np.zeros(im.shape,np.int8)
        start = time.time()
        npks = localmaxlabel(im.copy(), l2 , wk, cpu=0 )
        end=time.time()
        stime = (end-start)*1000.
        start = time.time()
        npks = localmaxlabel(im, l3 , wk, cpu=1 )
        end=time.time()
        ctime = (end-start)*1000.
        if (l2==l1).all() and (l3==l1).all():
            print("Timing python %.3f ms, c %.3f ms, sse %.3f"%(pytime,ctime,stime))
            self.assertTrue(True)
        else:
            print( "Fail",(l2==l1).all(), (l3==l1).all(), (l2==l3).all())
            self.assertTrue(False)
    def test2( self):
        im = make_test_image(128)[:107,:109].copy()
        start = time.time()
        l1 = localmax(im)
        end = time.time()
        pytime = (end-start)*1000.
        l2 = np.zeros(im.shape,np.int32)
        l3 = np.zeros(im.shape,np.int32)
        wk = np.zeros(im.shape,np.int8)
        start = time.time()
        npks = localmaxlabel(im.copy(), l2 , wk, cpu=0 )
        end=time.time()
        stime = (end-start)*1000.
        start = time.time()
        npks = localmaxlabel(im, l3 , wk, cpu=1 )
        end=time.time()
        ctime = (end-start)*1000.
        if (l2==l1).all() and (l3==l1).all():
            print("Timing python %.3f ms, c %.3f ms, sse %.3f"%(pytime,ctime,stime))
            self.assertTrue(True)
        else:
            print( "Fail",(l2==l1).all(), (l3==l1).all(), (l2==l3).all())
            import pylab as pl
            pl.subplot(221)
            pl.imshow(l1)
            pl.title("reference")
            pl.subplot(222)
            pl.imshow(l2)
            pl.title("ansi")
            pl.subplot(223)
            pl.imshow(l3)
            pl.title("sse")
            pl.subplot(224)
            pl.imshow(im)
            pl.title("problem image")
            pl.show()
            self.assertTrue(False)
    def test3( self):
        print("in test 2")
        sys.stdout.flush()
        im = make_test_image(N=2048)
        l2 = np.zeros(im.shape,np.int32)
        l3 = np.zeros(im.shape,np.int32)
        wk = np.zeros(im.shape,np.int8)
        start = time.time()
        for i in range(10):
            npks = localmaxlabel(im, l2 , wk )
        end=time.time()
        npks = localmaxlabel(im, l2 , wk, 10 )
        ctime = (end-start)*100.
        start = time.time()
        for i in range(10):
            npks = localmaxlabel(im, l3 , wk, 1 )
        end=time.time()
        npks = localmaxlabel(im, l3 , wk, 11 )
        stime = (end-start)*100.
        self.assertEqual( (l2==l3).all(), True )
        print("Timing 2048 c %.3f ms, sse %.3f"%(ctime,stime))



        
if __name__== "__main__":
    import sys
    if len(sys.argv)>1:
        demo()
    unittest.main()

else:
    im = make_test_image(256)
    l2 = np.zeros( im.shape, 'i' )
    wk = np.zeros( im.shape, np.uint8 )
    npks = localmaxlabel(im, l2 , wk )
    

