
import numpy as np
import ImageD11.sym_u
import ImageD11.unitcell
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


# a cell compatible with each group, in the setting the group expects
CELLS = { 'cubic' : (4.,4.,4.,90.,90.,90.),
          'hexagonal' : (3.,3.,5.,90.,90.,120.),
          'trigonal' : (4.913,4.913,5.405,90.,90.,120.),
          'trigonalH' : (4.913,4.913,5.405,90.,90.,120.),
          'trigonalP' : (5.,5.,5.,70.,70.,70.),
          'rhombohedralP' : (5.,5.,5.,70.,70.,70.),
          'tetragonal' : (3.,3.,5.,90.,90.,90.),
          'orthorhombic' : (3.,4.,5.,90.,90.,90.),
          'monoclinic_a' : (3.,4.,5.,100.,90.,90.),
          'monoclinic_b' : (3.,4.,5.,90.,100.,90.),
          'monoclinic_c' : (3.,4.,5.,90.,90.,100.),
          'triclinic' : (3.,4.,5.,80.,85.,95.) }

class test_ops_keep_cell( unittest.TestCase ):
    """ op . ubi must describe the same lattice : the metric is unchanged """
    def test_metric(self):
        for name, cell in CELLS.items():
            ubi = np.linalg.inv( ImageD11.unitcell.unitcell( cell ).B )
            g = np.dot( ubi, ubi.T )
            for op in ImageD11.sym_u.getgroup( name )().group:
                gop = np.dot( np.dot( op, g ), op.T )
                self.assertTrue( np.allclose( gop, g ), name )

    def test_find_uniq_u_trigonal(self):
        # issue 662 : alpha quartz came back as 4.913, 8.51, 5.405, 90, 90, 30
        rng = np.random.RandomState( 42 )
        B = ImageD11.unitcell.unitcell( CELLS['trigonal'] ).B
        g = np.linalg.inv( np.dot( B.T, B ) )
        grp = ImageD11.sym_u.trigonal()
        for i in range( 100 ):
            U = np.linalg.qr( rng.normal( size = (3,3) ) )[0]
            ubi = np.linalg.inv( np.dot( U, B ) )
            uniq = ImageD11.sym_u.find_uniq_u( ubi, grp )
            self.assertTrue( np.allclose( np.dot( uniq, uniq.T ), g ) )



if __name__=="__main__":
    unittest.main()

        
        
