
from __future__ import print_function, division

from ImageD11 import cImageD11
import numpy as np, time
import unittest


def norm(v):
    return v/np.linalg.norm(v, axis=1)[:,np.newaxis]

def py_scor_z( ubi,
               ub,
               gv ):
    """
    Decomposing the error into 3 orthogonal directions
      aligned with gvector (length, like two theta)
      normal with gvector and axis (e.g. omega direction)
      normal with that (e.g. azimuth direction)
    """
    nv = gv.shape[0]
    assert gv.shape[1] == 3
    assert ubi.shape == (3,3)
    assert ub.shape == (3,3)
    assert np.allclose( np.dot( ubi, ub), np.eye(3), 6)
    x,y,z = 0,1,2
    # local co-ordinates
    g0 = norm( gv )
    g1 = np.zeros( gv.shape, float )
    g1[:,x] = -g0[:,y]
    g1[:,y] =  g0[:,x]
    g1 = norm( g1 )
    g2 = np.zeros( gv.shape, float )
    g2[:,0] =  g0[:,x]*g0[:,z]
    g2[:,1] =  g0[:,y]*g0[:,z]
    g2[:,2] = -g0[:,x]*g0[:,x]-g0[:,y]*g0[:,y]
    g2 = norm( g2 )
    # error computations
    hkl = np.round( np.dot( ubi, gv.T ).T ) # integer hkl
    d = np.dot(ub, hkl.T).T - gv             # error in g
    e = np.zeros( gv.shape, float )
    e[:,0] = (d*g0).sum(axis=1)
    e[:,1] = (d*g1).sum(axis=1)
    e[:,2] = (d*g2).sum(axis=1)
    # now into local co-ordinates
    return g0, g1, g2, e

        

class test_closest_vec( unittest.TestCase ):

    def setUp( self ):
        
        UBI =np.array( [ [ 9.095614, -4.379317, -4.899265],
                 [ 6.253395 , 6.110302, 2.711559],
                 [ 2.786115, 9.817030 , 14.233333]])
        UB = np.linalg.inv( UBI )

        HKL = np.mgrid[-5:6, -5:6, -5:6].T.copy()
        HKL.shape = HKL.shape[1]*HKL.shape[1]*HKL.shape[1] , 3
        for i in range(len(HKL)):
            if abs(HKL[i]).sum()==0:
                HKL[i] = (1,1,1)
        GVE = np.dot( UB, HKL.T ).T.copy()
        GVE *= ( 0.95 + np.random.random( GVE.shape )/10.  )
        self.UBI = UBI
        self.UB =  UB
        self.GVE = GVE
    

    def test_same_as_python( self ):
        GVE, UBI, UB = self.GVE, self.UBI, self.UB
        pg0,pg1,pg2,pe = py_scor_z( UBI, UB, GVE )
        cg0 = np.zeros( GVE.shape, float )
        cg1 = np.zeros( GVE.shape, float )
        cg2 = np.zeros( GVE.shape, float )
        ce  = np.zeros( GVE.shape, float )
        cImageD11.score_gvec_z( UBI, UB, GVE,
                               cg0, cg1, cg2, ce, 1 )
        self.assertTrue( np.allclose( pg0, cg0, 6 ) )
        self.assertTrue( np.allclose( pg1, cg1, 6 ) )
        self.assertTrue( np.allclose( pg2, cg2, 6 ) )
        self.assertTrue( np.allclose( pe, ce, 6 ) )

if __name__=="__main__":
    unittest.main()
