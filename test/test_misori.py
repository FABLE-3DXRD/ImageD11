
from __future__ import print_function
from ImageD11.cImageD11 import misori_cubic, \
    misori_orthorhombic, misori_tetragonal, misori_monoclinic

import xfab.symmetry
import numpy as np
import numba
import timeit, time
import unittest, cProfile, pstats
timer = timeit.default_timer
#    1: Triclinic
#    2: Monoclinic
#    3: Orthorhombic
#    4: Tetragonal
#    5: Trigonal
#    6: Hexagonal
#    7: Cubic

ROTATIONS = [None] + [np.ascontiguousarray(xfab.symmetry.rotations(i)) for i in range(1,8)]

def isrotation(mat):
    d = np.linalg.det(mat)
    i = np.linalg.inv(mat)
    return np.allclose(d,1) and np.allclose( i, mat.T )

for grp in ROTATIONS[1:]:
    for mat in grp:
        assert isrotation(mat), "Not a rotation"+str(mat)
print("Rotation are all rotation matrices")
#print(ROTATIONS)
def Umis(umat_1, umat_2, crystal_system):
    rot = ROTATIONS[crystal_system]
    return _Umis(umat_1, umat_2, rot)

#@numba.njit(fastmath=True)
def _Umis(umat_1, umat_2, rot):
    misorientations = np.empty( (len(rot),2))
    lengths = 0.5 * (rot * np.dot( umat_1.T, umat_2 )).sum(axis=(2,1)) - 0.5
    misorientations[:,0] = np.arange(len(rot))
    misorientations[:,1] = np.arccos(lengths.clip(-1,1)) * 180./np.pi
    return misorientations

from scipy.spatial.transform import Rotation as R
if not hasattr( R, 'as_matrix' ):
    R.as_matrix = R.as_dcm

from xfab import symmetry
start = time.time()
def test_Umis():
    for crystal_system in range(1,8):
        U1 = R.random(200).as_matrix()
        U2 = R.random(200).as_matrix()
        for u1, u2 in zip(U1, U2):
            m1 = Umis(u1, u2, crystal_system)
            m2 = symmetry.Umis(u1, u2, crystal_system)
            assert np.max(np.abs(m1 - m2)) < 1e-5, str(np.max(np.abs(m1 - m2)))

test_Umis()
print("OK",time.time()-start)

def misori_py( u1, u2, sym=7 ):
    # 7 is cubic
#    ans = xfab.symmetry.Umis( u1, u2, sym)
    ans = Umis( u1, u2, sym )
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


NROT = 200

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
        U = make_random_orientations( NROT ) # xfab is a bit slow
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
        U = make_random_orientations( NROT ) # xfab is a bit slow
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
        U = make_random_orientations( NROT ) # xfab is a bit slow
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
        U = make_random_orientations( NROT ) # xfab is a bit slow
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
        U = make_random_orientations( NROT )
        t0 = timer()
        for u1 in U:
            for u2 in U:
                tr1 = misori_cubic( u1, u2 )
                tr2 = misori_cubic( u2, u1 )
                self.assertAlmostEqual( tr1, tr2 )
        t1 = timer()
        if self.DOBENCH:
            print("C %d pairs in %f s,  %f s per pair"%(
                NROT*NROT*2, t1-t0, (t1-t0)/NROT))


if __name__=="__main__":
    unittest.main()

