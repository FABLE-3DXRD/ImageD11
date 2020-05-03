

import unittest
from ImageD11 import finite_strain, indexing
import numpy as np

class setupRandom(object):
    def setUp(self):
        np.random.seed(4567)
        self.ubi = np.random.random((3,3)) + np.eye(3)
        self.ub = np.linalg.inv( self.ubi )
        self.ubi0 = np.random.random((3,3)) + np.eye(3)
        self.ub0 = np.linalg.inv( self.ubi0 )

class test_polar(setupRandom, unittest.TestCase):
    def test_polar(self):
        dgt = finite_strain.DeformationGradientTensor( self.ubi, self.ub0 )
        v,r,s = dgt.VRS
        self.assertTrue( np.allclose( np.dot( v, r), dgt.F ) )
        self.assertTrue( np.allclose( np.dot( r, s), dgt.F ) )

class test_svd(setupRandom, unittest.TestCase):
    def test_svd(self):
        dgt = finite_strain.DeformationGradientTensor( self.ubi, self.ub0 )
        u,s,vt = dgt.SVD
        fcalc = np.dot( np.dot( u, np.diag(s)) , vt )
        self.assertTrue( np.allclose( fcalc, dgt.F ) )

class test_xfab_strain( unittest.TestCase ):
    def setUp(self):
        a0 = 4.04
        b0 = 5.04
        c0 = 6.04
        from scipy.spatial.transform import Rotation
        r =  Rotation.from_euler("XYZ",(123,456,789),degrees=True)
        try:
            self.u = r.as_dcm()
        except:
            self.u = r.as_matrix()
        self.ubi = np.dot(((a0*1.1, 0,0),
                           (0,b0,0),
                           (0,0,c0)),self.u.T)
        self.ubi0 = np.diag((a0,b0,c0))
        self.ub0 = np.diag((1/a0,1/b0,1/c0))
        
    def test_strain_ref(self):
        from xfab.tools import ubi_to_u_and_eps
        dzero_cell = indexing.ubitocellpars( self.ubi0 )
        u, eps_xfab_ref_e6 = ubi_to_u_and_eps( self.ubi, dzero_cell )
        dgt = finite_strain.DeformationGradientTensor( self.ubi, self.ub0 )
        eps_new = dgt.finite_strain_ref( m=0.5 )
        self.assertTrue( np.allclose( dgt.U, u ) )
        self.assertTrue( np.allclose( dgt.U, self.u ) )
        self.assertTrue( np.allclose( finite_strain.symm_to_e6(eps_new), eps_xfab_ref_e6 ) )
        self.assertTrue( np.allclose( finite_strain.e6_to_symm(eps_xfab_ref_e6), eps_new ) )
        

if __name__=="__main__":
    unittest.main()
    



        
