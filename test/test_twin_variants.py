import os
import unittest

import numpy as np

import ImageD11.sym_u
import ImageD11.unitcell
from ImageD11.refinegrains import best_twin_variant, lf, polarization
from ImageD11.sinograms.sinogram import StructureFactorCorrection

HERE = os.path.dirname(os.path.abspath(__file__))


def tetragonal_4_m_sf(hmax=6):
    """|F|^2 with the 4/m symmetry of I4_1/a (88) but not 4/mmm:
    h k (h^2 - k^2) is kept by (h,k,l)->(-k,h,l) and by the inversion, not by h<->k."""
    r = np.arange(-hmax, hmax + 1)
    h, k, l = np.meshgrid(r, r, r, indexing="ij")
    t = h * k * (h * h - k * k)
    sf = 1.5 + t / (1.0 + np.abs(t))
    valid = (h + k + l) % 2 == 0  # I centring
    return StructureFactorCorrection(sf, valid, (-hmax, -hmax, -hmax))


class test_best_twin_variant(unittest.TestCase):

    def setUp(self):
        rng = np.random.RandomState(42)
        U = np.linalg.qr(rng.normal(size=(3, 3)))[0]
        if np.linalg.det(U) < 0:
            U = -U
        self.ubi = np.linalg.inv(np.dot(U, ImageD11.unitcell.unitcell([5.2, 5.2, 11.4, 90, 90, 90]).B))
        self.sf = tetragonal_4_m_sf()
        # peaks of the grain: every tabulated reflection that is not absent
        r = np.arange(-6, 7)
        hkl = np.array(np.meshgrid(r, r, r, indexing="ij")).reshape(3, -1)
        hkl = hkl[:, self.sf.valid_mask(hkl) & (np.abs(hkl).sum(axis=0) > 0)]
        self.gv = np.dot(np.linalg.inv(self.ubi), hkl).T
        n = hkl.shape[1]
        self.tth = rng.uniform(2, 20, n)
        self.eta = rng.uniform(-170, 170, n)
        # observed = |F|^2 / (lf / polarization), with 2 % noise
        self.intensity = (self.sf.sf_values(hkl) * polarization(self.tth, self.eta)
                          / lf(self.tth, self.eta) * rng.normal(1, 0.02, n))

    def test_finds_the_true_variant(self):
        laue = ImageD11.sym_u.laue_rotations(88)
        for start in ImageD11.sym_u.twin_variants(self.ubi, 88):
            best, scores = best_twin_variant(start, self.gv, self.intensity,
                                             self.tth, self.eta, 88, self.sf)
            self.assertEqual(len(scores), 2)
            self.assertGreater(np.nanmax(scores), 0.95)
            self.assertTrue(laue.isMember(np.dot(best, np.linalg.inv(self.ubi))))

    def test_no_ambiguity_keeps_ubi(self):
        best, scores = best_twin_variant(self.ubi, self.gv, self.intensity,
                                         self.tth, self.eta, 139, self.sf)
        self.assertEqual(len(scores), 1)
        self.assertTrue(np.allclose(best, self.ubi))


class test_structure_factors_from_cif(unittest.TestCase):

    def test_fcc_aluminium(self):
        try:
            import Dans_Diffraction  # noqa: F401
        except ImportError:
            self.skipTest("Dans_Diffraction is not installed")
        sf = StructureFactorCorrection.from_cif(
            os.path.join(HERE, "twinprob", "Al.cif"), wavelength=0.3, tth_max=20)
        hkl = np.array([[1, 1, 1], [2, 0, 0], [1, 0, 0], [1, 1, 0]]).T
        # F centring: mixed parity reflections are absent
        self.assertEqual(list(sf.valid_mask(hkl)), [True, True, False, False])
        f111, f200 = sf.sf_values(hkl[:, :2])
        self.assertGreater(f111, f200)


if __name__ == "__main__":
    unittest.main()
