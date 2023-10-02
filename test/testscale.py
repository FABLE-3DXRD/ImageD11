""" Write some testcases for the image scaling stuff """
import unittest

import numpy
from ImageD11.scale import scale

class testscale(unittest.TestCase):
    """old testcase """
    def testscaleimage(self):
        """check the object works """
        im1 = numpy.ones((10, 10), float)
        im1[2:8, 2:8] *= 2
        im1[4:6, 4:6] *= 2
        im2 = im1 * 2. + 3.
        o = scale(im1)
        a , b = o.scale(im2)
        self.assertEqual( abs(a-2.)< 1e-6 and abs(b-3.) < 1e-6 , True)
        scaled = o.scaleimage(im2)
        diff = numpy.ravel(scaled - im1)
        self.assertEqual( numpy.sum(diff*diff) < 1e-6 , True)

if __name__ == "__main__":
    unittest.main()
