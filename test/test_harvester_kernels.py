"""Regression tests for harvester.py kernels (_poisson_upper, _weight_sums_kernel)"""
from __future__ import print_function, division

import unittest

import numpy as np

from ImageD11.sinograms import harvester


class TestPoissonUpper(unittest.TestCase):
    def test_against_scipy(self):
        from scipy.stats import poisson
        for p in (1e-2, 1e-3, 1e-6):
            for mean in (0.5, 5.0, 100.0, 699.0, 700.0, 745.0, 800.0, 1000.0, 5000.0, 1e5):
                t = harvester._poisson_upper(mean, p)
                self.assertEqual(t, poisson.isf(p, mean), (mean, p))
                self.assertLessEqual(poisson.sf(t, mean), p)
                self.assertGreater(poisson.sf(t - 1, mean), p)

    def test_zero(self):
        self.assertEqual(harvester._poisson_upper(0.0, 1e-3), 0.0)


class TestWeightSums(unittest.TestCase):
    def test_against_numpy(self):
        rs = np.random.RandomState(0)
        H, W, n = 60, 50, 300
        hs, hf, ho = rs.randint(1, 6, n), rs.randint(1, 6, n), rs.randint(1, 4, n)
        s0, f0 = rs.randint(-6, H + 6, n), rs.randint(-6, W + 6, n)     # some boxes outside the detector
        good = rs.rand(H, W) > 0.1
        off = np.concatenate([[0], np.cumsum(2 * hs + 1 + 2 * hf + 1 + 2 * ho + 1)]).astype(np.int64)
        flat = rs.rand(off[-1])
        SB, SW, SC = harvester._weight_sums_kernel(flat, off, good, hs, hf, s0, f0)
        for k in range(n):
            r0, c0 = s0[k] - hs[k], f0[k] - hf[k]
            a0, a1, b0, b1 = max(r0, 0), min(s0[k] + hs[k] + 1, H), max(c0, 0), min(f0[k] + hf[k] + 1, W)
            if a1 <= a0 or b1 <= b0:
                self.assertEqual((SB[k], SW[k], SC[k]), (0, 0, 0))
                continue
            S = flat[off[k]:off[k] + 2 * hs[k] + 1]
            F = flat[off[k] + 2 * hs[k] + 1:off[k] + 2 * hs[k] + 2 * hf[k] + 2]
            w = np.outer(S[a0 - r0:a1 - r0], F[b0 - c0:b1 - c0]) * good[a0:a1, b0:b1]
            np.testing.assert_allclose([SB[k], SW[k], SC[k]], [(w * w).sum(), w.sum(), (w * w * w).sum()],
                                       rtol=1e-12, atol=1e-300)


if __name__ == '__main__':
    unittest.main()
