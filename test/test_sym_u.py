
import numpy as np
import ImageD11.sym_u
import unittest

class testsyms( unittest.TestCase ):

    def test_cubic(self):
        c = ImageD11.sym_u.cubic()
        assert len(c.group) == 24

    def test_hexagonal(self):
        c = ImageD11.sym_u.hexagonal()
        assert len(c.group) == 12

    def test_tetragonal(self):
        c = ImageD11.sym_u.tetragonal()
        assert len(c.group) == 8
        
    def test_trigonal(self):
        c = ImageD11.sym_u.trigonal()
        assert len(c.group) == 6

    def test_orthorhombic(self):
        c = ImageD11.sym_u.orthorhombic()
        assert len(c.group) == 4

    def test_monoclinic(self):
        c = ImageD11.sym_u.monoclinic_a()
        assert len(c.group) == 2
        c = ImageD11.sym_u.monoclinic_b()
        assert len(c.group) == 2
        c = ImageD11.sym_u.monoclinic_c()
        assert len(c.group) == 2

    def test_triclinic(self):
        c = ImageD11.sym_u.triclinic()
        assert len(c.group) == 1


class trh( object ):
    def test_right_handed(self):
        self.det = np.linalg.det( self.ubi )
        for op in self.c.group:
            ubi_trans = np.dot( op, self.ubi )
            self.assertAlmostEqual( np.linalg.det( ubi_trans ) , self.det )


class test_tetragonal_rh( unittest.TestCase , trh):
    def setUp(self):
        self.ubi = np.array( [ (1.3,0,0), ( 0,1.3,0 ), (0,0,2.1) ] )
        self.c = ImageD11.sym_u.tetragonal()
    
class test_cubic_rh( unittest.TestCase , trh):
    def setUp(self):
        self.ubi = np.array( [ (1.3,0,0), ( 0,1.3,0 ), (0,0,1.3) ] )
        self.c = ImageD11.sym_u.cubic()

class test_hexagonal_rh( unittest.TestCase , trh):
    def setUp(self):
        self.ubi = np.array( [ (1,0,0), ( -0.5, np.sqrt(3)/2, 0 ), (0,0,2.1) ] )
        self.c = ImageD11.sym_u.hexagonal()

class test_trigonal_rh(  unittest.TestCase , trh):
    def setUp(self):
        self.ubi = np.array( [ (1,0,0), ( -0.5, np.sqrt(3)/2, 0 ), (0,0,2.1) ] )
        self.c = ImageD11.sym_u.trigonal()

class test_orthorhombic_rh(  unittest.TestCase , trh):
    def setUp(self):
        self.ubi = np.array( [ (1,0,0), ( 0,1.3,0 ), (0,0,2.1) ] )
        self.c = ImageD11.sym_u.orthorhombic()



if __name__=="__main__":
    unittest.main()

        
        
