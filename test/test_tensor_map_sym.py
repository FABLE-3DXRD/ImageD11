
import numpy as np
import unittest
from ImageD11 import unitcell
from ImageD11.sinograms.tensor_map import sym_ops_from_phase


class test_sym_ops_from_phase( unittest.TestCase ):

    def test_laue_classes(self):
        # 4/m, -3, -31m, 6/m and m-3 used to raise ValueError
        for cell, sg, n in [ ((3,3,5,90,90,90), 88, 4),
                             ((4.9,4.9,5.4,90,90,120), 148, 3),
                             ((4.9,4.9,5.4,90,90,120), 162, 6),
                             ((3,3,5,90,90,120), 176, 6),
                             ((4,4,4,90,90,90), 205, 12),
                             ((4,4,4,90,90,90), 225, 24),
                             ((5,5,5,70,70,70), 166, 6) ]:
            ops = sym_ops_from_phase( unitcell.unitcell( cell, sg ) )
            self.assertEqual( len(ops), n, sg )
            for s in ops:
                self.assertTrue( np.allclose( np.dot( s, s.T ), np.eye(3) ) )
                self.assertAlmostEqual( np.linalg.det( s ), 1 )

    def test_nonstandard_setting(self):
        # P21/c is unique axis b: a cell with gamma != 90 is another setting
        with self.assertRaises( ValueError ):
            sym_ops_from_phase( unitcell.unitcell( (3,4,5,90,90,100), 14 ) )


if __name__=="__main__":
    unittest.main()
