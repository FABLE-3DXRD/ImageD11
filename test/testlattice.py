
""" Write some testcases for the lattice transforming stuff """
import unittest

import numpy as np
from ImageD11 import refinegrains, unitcell


class testlattices(unittest.TestCase):
    def setUp(self):
        self.uc =  [ 4.31, 4.33, 4.39, 90.1, 90.2, 90.3 ]
    def test_triclinic(self):
        cp = refinegrains.triclinic( self.uc )
        for i,j in zip(cp,self.uc):
            self.assertAlmostEqual( i, j, 6 )
    def test_mono_c(self):
        a,b,c,al,be,ga = refinegrains.monoclinic_c( self.uc )
        self.assertAlmostEqual( al, 90, 5 )
        self.assertAlmostEqual( be, 90, 5 )
    def test_mono_b(self):
        a,b,c,al,be,ga = refinegrains.monoclinic_b( self.uc )
        self.assertAlmostEqual( al, 90, 5 )
        self.assertAlmostEqual( ga, 90, 5 )
    def test_mono_a(self):
        a,b,c,al,be,ga = refinegrains.monoclinic_a( self.uc )
        self.assertAlmostEqual( be, 90, 5 )
        self.assertAlmostEqual( ga, 90, 5 )
    def test_ortho(self):
        a,b,c,al,be,ga = refinegrains.orthorhombic( self.uc )
        self.assertAlmostEqual( al, 90, 5 )
        self.assertAlmostEqual( be, 90, 5 )
        self.assertAlmostEqual( ga, 90, 5 )
    def test_tetra(self):
        a,b,c,al,be,ga = refinegrains.tetragonal( self.uc )
        self.assertAlmostEqual( al, 90, 5 )
        self.assertAlmostEqual( be, 90, 5 )
        self.assertAlmostEqual( ga, 90, 5 )
        self.assertAlmostEqual( a , b,  5 )
    def test_trigP(self):
        a,b,c,al,be,ga = refinegrains.trigonalP( self.uc )
        self.assertAlmostEqual( a , b,  5 )
        self.assertAlmostEqual( a , c,  5 )
        self.assertAlmostEqual( al , be,  5 )
        self.assertAlmostEqual( al , ga,  5 )
    def test_trigH(self):
        """ This one is quite far from the ideal cell ... """
        a,b,c,al,be,ga = refinegrains.trigonalH( self.uc )
        self.assertAlmostEqual( a , b,  5 )
        self.assertAlmostEqual( al , 90,  5 )
        self.assertAlmostEqual( be , 90,  5 )
        self.assertAlmostEqual( ga , 120,  5 )
    def test_hexagonal(self):
        """ This one is quite far from the ideal cell ... """
        a,b,c,al,be,ga = refinegrains.hexagonal( self.uc )
        self.assertAlmostEqual( a , b,  5 )
        self.assertAlmostEqual( al , 90,  5 )
        self.assertAlmostEqual( be , 90,  5 )
        self.assertAlmostEqual( ga , 120,  5 )
        self.assertAlmostEqual( c, self.uc[2] )
    def test_cubic(self):
        """ This one is quite far from the ideal cell ... """
        a,b,c,al,be,ga = refinegrains.cubic( self.uc )
        self.assertAlmostEqual( a, np.mean(self.uc[0:3]), 5)
        self.assertAlmostEqual( a , b,  5 )
        self.assertAlmostEqual( a , c,  5 )
        self.assertAlmostEqual( al , 90,  5 )
        self.assertAlmostEqual( be , 90,  5 )
        self.assertAlmostEqual( ga , 90,  5 )
        
        

        
if __name__ == "__main__":
    unittest.main()
