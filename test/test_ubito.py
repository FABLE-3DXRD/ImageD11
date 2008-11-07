
import unittest
from ImageD11.indexing import ubitoU,ubitoB,ubitoRod
import numpy as n


U_GS =  n.array([[-0.78834909,  0.61401788,  0.03857147],
                 [ 0.55263593,  0.67919891,  0.48299314],
                 [ 0.27036873,  0.40208318, -0.87477418]])


UBI_GS =  n.array([[-61.964236942, 43.437183163, 21.250981607],
                   [48.261804673, 53.385032994, 31.603737470],
                   [1.425215782, 17.846596993, -32.322906616]])

Rod_GS = n.array([5.033078533, 14.419161711, 3.818320914])


B = n.array([[  1.27226463e-02,   0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   1.27226463e-02,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   2.70635995e-02]])



class test1(unittest.TestCase):
    def test_ubitoB(self):

        b = ubitoB(UBI_GS)
        assert(n.abs(B-b) <10**-7).all(), 'some value of B is too far off'

    def test_ubitoU(self):

        u = ubitoU(UBI_GS)
        assert(n.abs(U_GS-u) <10**-7).all(), 'some value of U is too far from original U'
    def test_ubitoRod(self):

        Rod = ubitoRod(UBI_GS)
        assert(n.abs(Rod_GS-Rod) <10**-7).all(), 'some value of the Rodrigues vector  is too far from the original'
        
if __name__ == '__main__':
    unittest.main()
