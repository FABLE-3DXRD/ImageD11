

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

class test_m(setupRandom, unittest.TestCase):
    def test_lab(self):
        forwards = finite_strain.DeformationGradientTensor( self.ubi, self.ub0 )
        backward = finite_strain.DeformationGradientTensor( self.ubi0, self.ub )
        for i in range(5):
            m = i * 0.5
            ef = forwards.finite_strain_lab( m )
            # this is a bit misleading - negative m puts you in the other frame
            eb = backward.finite_strain_ref( -m )
            self.assertTrue( np.allclose( ef, -eb ) )
    def test_ref(self):
        forwards = finite_strain.DeformationGradientTensor( self.ubi, self.ub0 )
        backward = finite_strain.DeformationGradientTensor( self.ubi0, self.ub )
        for i in range(5):
            m = i * 0.5
            ef = forwards.finite_strain_ref( m )
            eb = backward.finite_strain_lab( -m )
            self.assertTrue( np.allclose( ef, -eb ) )
            
        
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
        if hasattr(r, "as_matrix"):
            self.u = r.as_matrix()
        else:
            self.u = r.as_dcm()
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

    def test_strain_lab(self):
        from xfab.tools import ubi_to_u_and_eps
        dzero_cell = indexing.ubitocellpars( self.ubi0 )
        u, eps_xfab_ref_e6 = ubi_to_u_and_eps( self.ubi, dzero_cell )
        emat = finite_strain.e6_to_symm( eps_xfab_ref_e6 )
        Exfab = np.dot(np.dot( u, emat), u.T )
        dgt = finite_strain.DeformationGradientTensor( self.ubi, self.ub0 )
        eps_new = dgt.finite_strain_lab( m=0.5 )
        self.assertTrue( np.allclose( dgt.U, u ) )
        self.assertTrue( np.allclose( dgt.U, self.u ) )
        self.assertTrue( np.allclose( eps_new, Exfab ) )

class test_xfab_largestrain( test_xfab_strain ):
    def setUp(self):
        a0 = 4.04
        b0 = 5.04
        c0 = 6.04
        from scipy.spatial.transform import Rotation
        r =  Rotation.from_euler("XYZ",(1230,4560,7890),degrees=True)
        if hasattr(r, "as_matrix"):
            self.u = r.as_matrix()
        else:
            self.u = r.as_dcm()
        self.ubi = np.dot(((a0*12.1, 0,0),
                           (0,b0,0),
                           (0,0,c0)),self.u.T)
        self.ubi0 = np.diag((a0,b0,c0))
        self.ub0 = np.diag((1/a0,1/b0,1/c0))
        

class test_y_strain(unittest.TestCase):
    def setUp(self):
        # longer along the x-axis in the lab
        # this is the a axis in the reference
        self.ubi = np.array([[0,1,0],[-1.1,0,0],[0,0,1]],float) 
        self.ubi0 = np.eye(3)
        self.ub0 = np.linalg.inv(self.ubi0)
        self.dgt = finite_strain.DeformationGradientTensor( self.ubi, self.ub0 )

    def test_strainlab(self):
        elab = self.dgt.finite_strain_lab(m=0.5)
        self.assertAlmostEqual( elab[0,0], 0.1 )
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual( elab[i,j], elab[j,i])
                if i == j == 0:
                    continue
                self.assertAlmostEqual( elab[i,j], 0. )
                
    def test_strainref(self):
        eref = self.dgt.finite_strain_ref(m=0.5)
        self.assertAlmostEqual( eref[1,1], 0.1 )
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual( eref[i,j], eref[j,i])
                if i == j == 1:
                    continue
                self.assertAlmostEqual( eref[i,j], 0. )
        
                
if __name__=="__main__":
    unittest.main()
