
from __future__ import division, print_function
import numpy as np
import timeit
from ImageD11.cImageD11 import localmaxlabel, sparse_localmaxlabel
import unittest, sys

timer = timeit.default_timer

def sparselocalmax( im, labels, wk ):
    row, col = np.mgrid[ 0: im.shape[0], 0:im.shape[1] ]
    mask = im > 0
    row = row[mask].astype(np.uint16)
    col = col[mask].astype(np.uint16)
    signal = im[mask].astype(np.float32)
    vmx = np.empty( len(signal), np.float32 )
    imx = np.empty( len(signal), np.int32 )
    slabel = np.empty( len(signal), np.int32 )
    labels.fill(0)
#    print('row',row)
#    print('col',col)
#    print('signal',signal)
    nlabel = sparse_localmaxlabel( signal, row, col, vmx, imx, slabel )
#    print('nlabel',nlabel)
#    print('slabel',slabel)
#    print('vmx   ',vmx)
#    print('imx   ',imx)
    labels[mask] = slabel
    return nlabel


try:
    from numba import jit, prange
    GOTNUMBA = True
    print("Got numba")
except:
    print("No numba")
    GOTNUMBA = False
    def jit(pyfunc=None, **kwds):
        def wrap(func):
            return func
        if pyfunc is not None:
            return wrap(pyfunc)
        else:
            return wrap
    def prange(*a):
        return range(*a)


@jit(nopython=True)
def localmax( im, labels, wk ):
    """ reference implementation 
    ... very slow 
    """
    out = wk
    nf = im.shape[1]
    for i in prange(1,im.shape[0]-1):
        out[i,0] = nf*i
        out[i,nf-1] = nf*(i+1)-1
        for j in range(1,nf-1):
            mx = im[i-1,j-1]
            out[i,j] = (i-1) * nf + j-1
            for k,l in ( (i-1, j  ),
                         (i-1, j+1),
                         (i  , j-1),
                         (i  , j  ),
                         (i  , j+1),
                         (i+1, j-1),
                         (i+1, j  ),
                         (i+1, j+1) ):           
                if im[k,l] > mx:
                    mx = im[k,l]
                    out[i,j] = k * nf + l
            
            # now out addresses a neighbour
    # Final labels
    npk = 0
    labels[0]=0   # borders
    labels[-1]=0
    for i in prange(1,im.shape[0]-1):
        labels[i,0] = 0
        labels[i,nf-1] = 0
        for j in range(1,nf-1):
            if out[i,j] == i*nf + j: # maxpx
                npk += 1
                labels[i,j]=npk
                out[i,j]=0
            else:
                labels[i,j] = -1
    # Now relabel each pixel to point to the max...
    pon=False
    for i in range(1,im.shape[0]-1):
        for j in range(1,nf-1):
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
    l1 = np.zeros(im.shape, np.int32)
    wp = np.zeros(im.shape, np.int32)
    localmax(im,l1,wp)
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
        self.im = make_test_image(128)[:107,:109].copy()
        l1 = np.zeros(self.im.shape, np.int32)
        wp = np.zeros(self.im.shape, np.int32)
        # numba runs
        localmax(self.im, l1, wp)
        self.dopy = True
    def test1( self):
        im = self.im
        if self.dopy:
            l1 = np.zeros(im.shape, np.int32)
            wp = np.zeros(im.shape, np.int32)
            start = timer()
            localmax(im, l1, wp)
            end = timer()
            pytime = (end-start)*1000.
        l = [l1,]
        wk = [0,]
        times = [pytime,]
        names = ["python",]
        for cpu, name in [(0,"ansi")]:#,(1,"sse2"),(2,"avx2")]:
            li = np.zeros(im.shape,np.int32)
            wi = np.zeros(im.shape,np.int8)
            start = timer()
    #        print("cpu",cpu)
            npks = localmaxlabel(im.copy(), li , wi )
            end=timer()
            ctime = (end-start)*1000.
            l.append( li )
            times.append( ctime )
            names.append( name )
        lref = l[0]
        print("Timing",end=" ")
        for i in range(len(names)):
            print( "%s %.3f ms"%( names[i], times[i]), end=" ")
        print()
        for li in l[1:]:
            self.assertTrue((lref == l1).all())
    def test3( self):
     #   print("in test 2")
        im = make_test_image(N=2048)
        l0 = np.zeros(im.shape,np.int32)
        l1 = np.zeros(im.shape,np.int32)
        wp = np.zeros(im.shape,np.int32)
        wk = np.zeros(im.shape,np.int8)
        N = 4
        if GOTNUMBA:
            npks = localmax(im, l0, wp)
            start = timer()
            for i in range(N):
                npks = localmax(im, l0, wp)
            end=timer()
            ptime = (end-start)*1000./N
        else:
            ptime = np.nan
        npks = localmaxlabel(im, l1 , wk )
        start = timer()
        for i in range(N):
            npks = localmaxlabel(im, l1 , wk )
        end=timer()
        ctime = (end-start)*1000./N
        print("Timing 2048 numba %.3f ms c %.3f ms"%(
            ptime, ctime))

    def testbug(self):
        ''' case caught with image from Mariana Mar Lucas,
        caused a segfault'''
        im = np.array([
                [ 0, 0, 0, 0, 0, 0, 0, 0],
                [ 0, 2, 0, 0, 0, 0, 0, 0],
                [ 0, 1, 1, 0, 0, 0, 0, 0],
                [ 0, 0, 0, 0, 0, 0, 0, 0],
                [ 0, 0, 0, 0, 0, 0, 0, 0],
                [ 0, 0, 0, 0, 0, 0, 0, 0]], np.float32)
        la = np.empty( im.shape, 'i' )
        wk = np.empty( im.shape, 'b' )
        localmaxlabel( im, la, wk )
        self.assertEqual(la[1,1],1)

    def testsmall(self):
        ''' sse / avx issue for tiny? '''
#        print("testsmall, fixme, segfault/infinite loop so no CI")
        im = np.array([
                [ 0, 0, 0, 0, 0],
                [ 0, 2, 0, 0, 0],
                [ 0, 1, 1, 0, 0],
                [ 0, 0, 0, 0, 0]], np.float32)
        la = np.empty( im.shape, 'i' )
        wk = np.empty( im.shape, 'b' )
        localmaxlabel( im, la, wk )
        self.assertEqual(la[1,1],1)

        
    def testX(self):
        ''' sse / avx issue for tiny? '''
        im = np.array([
                [ 0, 0, 0, 0, 0],
                [ 0, 1, 0, 1, 0],
                [ 0, 0, 2, 0, 0],
                [ 0, 1, 0, 1, 0],            
                [ 0, 0, 0, 0, 0]], np.float32)
        la = np.empty( im.shape, 'i' )
        wk = np.empty( im.shape, 'b' )
        localmaxlabel( im, la, wk )
        lp = np.empty( im.shape, 'i' )
        wp = np.empty( im.shape, 'b' )
        npks = localmax(im, lp, wp)
        if not (lp == la).all():
            print(lp)
            print(la)
        
        self.assertTrue( (lp == la).all() )

        
    def test_one(self):
        im = np.array([
            [0, 0, 0, 0, 0,],
            [0, 0, 2, 2, 0,],
            [0, 0, 0, 3, 0,],
            [0, 0, 0, 0, 0,],
            [0, 0, 0, 0, 0,],
            [0, 0, 0, 3, 0,],
            [0, 0, 2, 3, 0,],
            [0, 0, 0, 0, 0,],            
        ])
        la = np.empty( im.shape, 'i' )
        wk = np.empty( im.shape, 'b' )
        npks = localmaxlabel( im, la, wk )
        lp = np.empty( im.shape, 'i' )
        wp = np.empty( im.shape, 'b' )
        localmax(im, lp, wp)
        if npks != 2:
            print(la)
            assert(npks == 2), npks
        if not (lp == la).all():
            print(lp)
            print(la)
            
            
    def test_L_sparse(self):
        im = np.array([
            [0, 0, 0, 0, 0,],
            [0, 0, 2, 2, 0,],
            [0, 0, 0, 4, 0,],
            [0, 0, 0, 0, 0,],
            [0, 0, 0, 0, 0,],
            [0, 0, 0, 3, 0,],
            [0, 0, 2, 3, 0,],
            [0, 0, 0, 0, 0,],            
        ])
        la = np.empty( im.shape, 'i' )
        wk = np.empty( im.shape, 'b' )
        npks = localmaxlabel( im, la, wk )
        assert npks == 2
        las = np.empty( im.shape, 'i' )
        npks = sparselocalmax( im, las, wk )
        assert npks == 2
        m = im > 0
        assert (las[m] == la[m]).all()
        
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
    

