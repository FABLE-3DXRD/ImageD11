
from __future__ import print_function
from ImageD11 import cImageD11
import fabio
import numpy as np
import scipy.sparse
import unittest
import time

import platform
if platform.system() == "Windows":
    timer = time.clock
else:
    timer = time.time

class test_array_stats( unittest.TestCase ):
    def test1( self ):
        ar = np.arange( 6*4, dtype=np.float32 ).reshape( 6,4 )
        mini, maxi, mean, var = cImageD11.array_stats( ar.ravel() )
        t = np.allclose( (mini, maxi, mean, var),
                (ar.min(), ar.max(), ar.mean(), ar.var()) )
        if not t:
            print((mini, maxi, mean, var))
            print((ar.min(), ar.max(), ar.mean(), ar.var()) )
        self.assertTrue(t)
        
    def test2( self ):
        ar = np.ones( 6*4, dtype=np.float32 ).reshape( 6,4 )
        mini, maxi, mean, var = cImageD11.array_stats( ar.ravel() )
        self.assertTrue( np.allclose( 
            (mini, maxi, mean, var),
            (ar.min(), ar.max(), ar.mean(), ar.var()) ))
    def test3( self ):
        ar = np.arange( 2048*2048, dtype=np.float32 ).reshape( 2048, 2048 )
        start = timer() 
        mini, maxi, mean, var = cImageD11.array_stats( ar.ravel() )
        endc = timer()
        check = ar.min(), ar.max(), ar.mean(), ar.var()
        endn = timer()
        self.assertTrue( np.allclose( (mini, maxi, mean, var), check ) )
        numpytime = endn - endc
        ctime = endc - start
        print("array_stats: %.3f ms vs %.3f ms, speedup %.1f"%(
            1e3*ctime, 1e3*numpytime, numpytime/ctime))


class test_array_histogram( unittest.TestCase ):
    def test1(self):
        N = 255
        a = (np.random.random( (2048, 2048))*(N-0.1)).astype(np.float32)
        start = timer()
        hc = np.zeros( N, dtype=np.int32 )
        cImageD11.array_histogram( a.ravel(), 0., N, hc )
        endc= timer()
        hn = np.bincount( a.ravel().astype(int) )
        endn = timer()
        if( hc.shape[0] != hn.shape[0]):
            print("shape error")
            print( a.max())
            print(hc)
            print(hn)
        self.assertTrue( (hc == hn).all() )
        numpytime = endn - endc
        ctime = endc - start
        print("array_histogram: %.3f ms vs %.3f ms, speedup %.1f"%(
            1e3*ctime, 1e3*numpytime, numpytime/ctime))


class test_sparse_is_sorted( unittest.TestCase ):
    def setUp( self ):
        self.i = np.arange(10, dtype = np.uint16 )
        self.j = np.arange(10, dtype = np.uint16 )
    def testOK( self ):
        ret = cImageD11.sparse_is_sorted( self.i, self.j)
        self.assertEqual( ret, 0)
    def testrpt( self ):
        self.i[5]=self.i[4]
        self.j[5]=self.j[4]
        ret = cImageD11.sparse_is_sorted( self.i, self.j)
        self.assertEqual( ret, -5)
    def testbadi( self ):
        self.i[5]=self.i[4]-1
        ret = cImageD11.sparse_is_sorted( self.i, self.j)
        self.assertEqual( ret, 5)
    def testbadj( self ):
        self.i[5]=self.i[4]
        self.j[5]=self.j[4]-1
        ret = cImageD11.sparse_is_sorted( self.i, self.j)
        self.assertEqual( ret, 5)


class test_sparse_connected_pixels( unittest.TestCase ):
    def setUp( self ):
        img ="""
00110100000000000000000000000000000000000000000000000000
00001000000000000000000000000000000000000000000000000000
00000000001111000000111110000000000000000000000000000010
00000000011001000001100011000000000000000000000000010010
11000000000000000000000000000000000000000000000000001111
00000000000000000000000000000000020000000000000000000000
00000000000000000000000000000000020000000000000000000000
00000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000001
00000000000000000001100000000000000000000000000000010001
00000000000000000010011000000000000000000000000000001110
00000000000000000100000100000000000000000000000000000000
00000000000000001000200010000000000000000000000000000000
00000000000000001000000100000000000000000000000000000000
00000000000000000100001000000000000000000000000000000000
00000000000000000011110000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000
00001100000000000000000000000000000000000000000000000000
11111100000000000000000000000000000000000000000000000000
00010000000000000000000100000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000
00000000000000000000000010000000000000000000000100000000
"""
        self.a = np.array([ [int(c) for c in line ] for line in img.split() ],
                          dtype=np.float32)
        self.s = scipy.sparse.coo_matrix( self.a )
        self.threshold=0.5
        if 0:
            print(self.s.data)
            print(self.s.row)
            print(self.s.col)
        
    def testOK( self ):
        i = self.s.row.astype( np.uint16 )
        j = self.s.col.astype( np.uint16 )
        v = self.s.data
        self.assertTrue( (self.s.todense() == self.a).all() ) 
        dl = np.zeros( self.a.shape, 'i' )
        nd = cImageD11.connectedpixels(  self.a , dl, self.threshold )
        dstart = timer()
        nd = cImageD11.connectedpixels(  self.a , dl, self.threshold )
        dend = timer()
        sl = np.zeros( v.shape, 'i' )
        ns = cImageD11.sparse_connectedpixels( v, i, j, self.threshold, sl)
        sbegin = timer()
        ns = cImageD11.sparse_connectedpixels( v, i, j, self.threshold, sl)
        send = timer()
        sld = scipy.sparse.coo_matrix( (sl, (i, j)), shape=self.a.shape)
        testcase = (sld.todense() == dl).all()

        densetime = dend - dstart
        sparsetime = send-sbegin
        print("sp_cp: %.3f ms vs %.3f ms, ratio %.1f nnz %d"%(
            1e3*sparsetime, 1e3*densetime, densetime/sparsetime, sld.nnz),
            " t=%.0f"%(self.threshold))
        if ~testcase:
            import pylab as pl
            pl.ioff()
            pl.figure()
            pl.subplot(311)
            pl.title("t=%f"%(self.threshold))
            pl.imshow( dl, interpolation='nearest')
            pl.colorbar()
            pl.subplot(312)
            pl.imshow( sld.todense(), interpolation='nearest')
            pl.colorbar()
            pl.subplot(313)
            pl.imshow( sld.todense() - dl, interpolation='nearest')
            pl.colorbar()
            pl.show()
        self.assertTrue( testcase )
   
class test_sparse_connected_pixels_overlapim1(test_sparse_connected_pixels):
    def setUp(self):
        im = fabio.open("testoverlaps0000.edf").data.astype( np.float32 )
        self.a = np.where( im > 1500, im, 0 )
        self.s = scipy.sparse.coo_matrix( self.a )
        self.threshold = 2000.1

class test_sparse_connected_pixels_overlapim2(test_sparse_connected_pixels):
    def setUp(self):
        im = fabio.open("testoverlaps0000.edf").data.astype( np.float32 )
        self.a = np.where( im > 2000.1, im, 0 )
        self.s = scipy.sparse.coo_matrix( self.a )
        self.threshold = 2000.1

class test_sparse_connected_pixels_overlapim3(test_sparse_connected_pixels):
    def setUp(self):
        im = fabio.open("testoverlaps0000.edf").data.astype( np.float32 )
        self.a = np.where( im > 500.1, im, 0 )
        self.s = scipy.sparse.coo_matrix( self.a )
        self.threshold = 500.1

class test_sparse_connected_pixels_overlapim4(test_sparse_connected_pixels):
    def setUp(self):
        im = fabio.open("testoverlaps0000.edf").data.astype( np.float32 )
        self.a = np.where( im > 5000.1, im, 0 )
        self.s = scipy.sparse.coo_matrix( self.a )
        self.threshold = 5000.1
        

if __name__=="__main__":
    unittest.main()


#    cd .. && python setup.py build && cd test && PYTHONPATH=../build/lib.linux-x86_64-2.7 python test_sparse_image.py
