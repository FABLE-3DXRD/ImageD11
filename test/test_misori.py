
from __future__ import print_function
from ImageD11.cImageD11 import misori_cubic, misori_cubic_pairs, \
    misori_orthorhombic, misori_tetragonal, misori_monoclinic

import xfab.symmetry
import numpy as np
import timeit
import unittest, cProfile, pstats
timer = timeit.default_timer
#    1: Triclinic
#    2: Monoclinic
#    3: Orthorhombic
#    4: Tetragonal
#    5: Trigonal
#    6: Hexagonal
#    7: Cubic
    

def misori_py( u1, u2, sym=7 ):
    # 7 is cubic
    ans = xfab.symmetry.Umis( u1, u2, sym)
    ang = np.min( ans[:,1] )
    trc = np.cos(np.radians(ang))*2+1
    return trc

def make_random_orientations( N ):
    """
    Generate random quaternions and convert to U matrices
    """
    q = np.random.standard_normal( (4, N) )
    s = 1/(q*q).sum( axis=0 )
    U = np.zeros( (N, 3, 3), float )
    r,i,j,k = 0,1,2,3
    U[:,0,0] = 1 - 2*s*(q[j]*q[j]+q[k]*q[k])
    U[:,0,1] = 2*s*(q[i]*q[j]-q[k]*q[r])
    U[:,0,2] = 2*s*(q[i]*q[k]+q[j]*q[r])
    U[:,1,0] = 2*s*(q[i]*q[j]+q[k]*q[r])    
    U[:,1,1] = 1 - 2*s*(q[i]*q[i]+q[k]*q[k])
    U[:,1,2] = 2*s*(q[j]*q[k]-q[i]*q[r])
    U[:,2,0] = 2*s*(q[i]*q[k]-q[j]*q[r])
    U[:,2,1] = 2*s*(q[j]*q[k]+q[i]*q[r])
    U[:,2,2] = 1 - 2*s*(q[i]*q[i]+q[j]*q[j])
    return U




class test_random_orientations( unittest.TestCase ):
    DOPROFILE = False
    DOBENCH   = True
    def setUp(self):
        np.random.seed(42)
        if self.DOPROFILE:
            self.pr = cProfile.Profile()
            self.pr.enable()

    def tearDown(self):
        if self.DOPROFILE: 
            p = pstats.Stats( self.pr )
            p.strip_dirs()
            p.sort_stats ('cumtime')
            p.print_stats ()
            print( "\n--->>>" )

    def test_are_rotations( self ):
        self.U = make_random_orientations( 1000 )
        i = np.eye(3)
        # transpose is inverse
        ids = [ abs( (np.dot( m, m.T )-i ) ).ravel().sum() for m in self.U ]
        self.assertAlmostEqual( np.max( ids ), 0 )
        # determinant is 1
        dts = [ np.linalg.det(m) - 1 for m in self.U ]
        self.assertAlmostEqual( abs(np.array( dts )).max(), 0 )
        
    def test_cubic(self):
        U = make_random_orientations( 20 ) # xfab is a bit slow
        t0 = timer()
        pairs = [ (u1, u2) for u1 in U for u2 in U ]
        c = [    misori_cubic( u1, u2 ) for u1, u2 in pairs ]
        t1 = timer()
        p = [ misori_py( u1, u2, 7 ) for u1, u2 in pairs ]
        t2 = timer()
        if self.DOBENCH:
            print("C time %f s , pytime %f s"%(t1-t0, t2-t1))
        self.assertTrue(np.allclose( np.array(c) ,np.array(c) ))
        
    def test_monoclinic(self):
        U = make_random_orientations( 20 ) # xfab is a bit slow
        t0 = timer()
        pairs = [ (u1, u2) for u1 in U for u2 in U ]
        c = [    misori_monoclinic( u1, u2 ) for u1, u2 in pairs ]
        t1 = timer()
        p = [ misori_py( u1, u2, 2 ) for u1, u2 in pairs ]
        t2 = timer()
        if self.DOBENCH:
            print("C time %f s , pytime %f s"%(t1-t0, t2-t1))
        self.assertTrue(np.allclose( np.array(c) ,np.array(c) ))

    def test_orthorhombic(self):
        U = make_random_orientations( 20 ) # xfab is a bit slow
        t0 = timer()
        pairs = [ (u1, u2) for u1 in U for u2 in U ]
        c = [    misori_orthorhombic( u1, u2 ) for u1, u2 in pairs ]
        t1 = timer()
        p = [ misori_py( u1, u2, 3 ) for u1, u2 in pairs ]
        t2 = timer()
        if self.DOBENCH:
            print("C time %f s , pytime %f s"%(t1-t0, t2-t1))
        self.assertTrue(np.allclose( np.array(c) ,np.array(c) ))

    def test_tetragonal(self):
        U = make_random_orientations( 20 ) # xfab is a bit slow
        t0 = timer()
        pairs = [ (u1, u2) for u1 in U for u2 in U ]
        c = [    misori_tetragonal( u1, u2 ) for u1, u2 in pairs ]
        t1 = timer()
        p = [ misori_py( u1, u2, 4 ) for u1, u2 in pairs ]
        t2 = timer()
        if self.DOBENCH:
            print("C time %f s , pytime %f s"%(t1-t0, t2-t1))
        self.assertTrue(np.allclose( np.array(c) ,np.array(c) ))
        
    def test_cubic_reverse(self):
        N = 250
        U = make_random_orientations( N )
        t0 = timer()
        for u1 in U:
            for u2 in U:
                tr1 = misori_cubic( u1, u2 )
                tr2 = misori_cubic( u2, u1 )
                self.assertAlmostEqual( tr1, tr2 )
        t1 = timer()
        if self.DOBENCH:
            print("C %d pairs in %f s,  %f s per pair"%(
                N*N*2, t1-t0, (t1-t0)/N))

    def test_cubic_pairmat(self):
        N = 500
        U = make_random_orientations( N )
        m0 = []
        t0 = timer()
        for i in range(N):
            for j in range(i+1,N):
                m0.append( misori_cubic( U[i], U[j] ) )
        t1 = timer()
        m1 = np.zeros( len(m0), float )
        misori_cubic_pairs( U, m1 )
        t2 = timer()
        if self.DOBENCH:
            print( "time for distance matrix",t1-t0, t2-t1)
        self.assertTrue( np.allclose( m0, m1 ) )

    
if __name__=="__main__":
    unittest.main()
            
