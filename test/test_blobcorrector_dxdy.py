import unittest
import numpy as np
from ImageD11 import blobcorrector


class test_dxdy(unittest.TestCase):
    def setUp(self):
        s, f = np.mgrid[0:200, 0:300].astype(float)
        self.dx = 0.5 * np.sin(f / 40.) + 0.001 * s      # fast
        self.dy = 0.3 * np.cos(s / 30.) - 0.002 * f      # slow

    def test_invert(self):
        rng = np.random.default_rng(1)
        s = rng.uniform(5, 190, 500)
        f = rng.uniform(5, 290, 500)
        sc, fc = blobcorrector.correct_dxdy_bilinear(s, f, self.dx, self.dy)
        s2, f2 = blobcorrector.invert_dxdy(sc, fc, self.dx, self.dy, iterations=8)
        # the inverse uses nearest pixel look ups, so is good to a fraction of the map gradient
        self.assertLess(np.abs(s2 - s).max(), 0.05)
        self.assertLess(np.abs(f2 - f).max(), 0.05)

    def test_bilinear_vs_nearest(self):
        ds = blobcorrector.detector_spatial.__new__(blobcorrector.detector_spatial)
        ds.dx, ds.dy = self.dx, self.dy
        s, f = np.arange(10.0, 100.0), np.arange(20.0, 110.0)     # integer pixels: identical
        pk = ds({'s_raw': s, 'f_raw': f})
        sc, fc = ds.correct_bilinear(s, f)
        self.assertTrue(np.allclose(sc, pk['sc']) and np.allclose(fc, pk['fc']))

    def test_pixel_lut(self):
        """pixel_lut agrees with __call__ at every pixel, and returns (slow, fast)"""
        ds = blobcorrector.detector_spatial.__new__(blobcorrector.detector_spatial)
        ds.dx, ds.dy = self.dx, self.dy
        sl, fa = ds.pixel_lut()
        s, f = np.mgrid[0:200, 0:300]
        pk = ds({'s_raw': s.ravel().astype(float), 'f_raw': f.ravel().astype(float)})
        self.assertTrue(np.allclose(sl.ravel(), pk['sc']))
        self.assertTrue(np.allclose(fa.ravel(), pk['fc']))
        self.assertTrue(np.allclose(sl, s + self.dy) and np.allclose(fa, f + self.dx))

    def test_none(self):
        s, f = np.array([1.5]), np.array([2.5])
        for fn in (lambda: blobcorrector.invert_dxdy(s, f, None, None),
                   lambda: blobcorrector.correct_dxdy_bilinear(s, f, None, None)):
            a, b = fn()
            self.assertTrue(a[0] == 1.5 and b[0] == 2.5)


if __name__ == "__main__":
    unittest.main()
