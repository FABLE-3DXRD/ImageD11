

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

        

if __name__=="__main__":
    unittest.main()

        
        
