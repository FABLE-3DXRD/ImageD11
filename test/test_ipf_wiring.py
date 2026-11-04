"""IPF colours through unitcell, grain, TensorMap and the IPF plot (ImageD11.rgb_colors)."""

import sys
import unittest

import numpy as np

from ImageD11 import grain, rgb_colors, sym_u, unitcell
from ImageD11.ipf_keys import GROUPS

HEX = (3.0, 3.0, 5.0, 90.0, 90.0, 120.0)


def rotations(rng, n):
    q = np.array([np.linalg.qr(m)[0] for m in rng.normal(size=(n, 3, 3))])
    return q * np.sign(np.linalg.det(q))[:, None, None]


def strained_ubis(ucell, rng, n, strain=1e-3):
    """ubi = inv(U.B) with a small strain of the cell"""
    eps = rng.uniform(-strain, strain, (n, 3, 3))
    UB = np.einsum("nij,jk,nkl->nil", rotations(rng, n), ucell.B, np.eye(3) + (eps + eps.transpose(0, 2, 1)) / 2)
    return np.linalg.inv(UB)


class TestSpacegroup(unittest.TestCase):

    def test_laue_class(self):
        """the point group of each space group has the Laue class sym_u gives"""
        for sgno in range(1, 231):
            name = rgb_colors.spacegroup_group(sgno)
            self.assertEqual(GROUPS[name]["laue"].replace("12/m1", "2/m"), sym_u.laue_class(sgno), sgno)
            r = rgb_colors.spacegroup_group(sgno, "rhombohedral")
            self.assertEqual(GROUPS[r]["laue"].replace("12/m1", "2/m"), sym_u.laue_class(sgno, "rhombohedral"), sgno)

    def test_unitcell(self):
        self.assertEqual(unitcell.unitcell(HEX, 194).ipf_group, "6/mmm")
        self.assertEqual(unitcell.unitcell((4.9, 4.9, 13.6, 90, 90, 120), 166).ipf_group, "-3m1")
        self.assertEqual(unitcell.unitcell((5, 5, 5, 70, 70, 70), 166).ipf_group, "-3mR")
        self.assertEqual(unitcell.unitcell((4, 5, 6, 90, 100, 90), 14).ipf_group, "12/m1")
        with self.assertRaises(AttributeError):
            unitcell.unitcell(HEX, "P").ipf_group


class TestColours(unittest.TestCase):

    def test_unitcell(self):
        rng = np.random.RandomState(1)
        uc = unitcell.unitcell(HEX, 194)
        ubis = strained_ubis(uc, rng, 20)
        U = rotations(rng, 20)
        for key in ("tsl", "hsv"):
            rgb = rgb_colors.ubi_to_rgb(ubis, "6/mmm", key=key, axis=(1, 0, 0))
            self.assertLess(np.abs(uc.get_ipf_colour_from_ubi(ubis, (1, 0, 0), key) - rgb).max(), 1e-15)
            self.assertLess(np.abs(uc.get_ipf_colour(np.linalg.inv(ubis), (1, 0, 0), key) - rgb).max(), 1e-12)
            self.assertLess(np.abs(uc.get_ipf_colour_fast(U, (1, 0, 0), key) -
                                   uc.get_ipf_colour(np.dot(U, uc.B), (1, 0, 0), key)).max(), 1e-12)

    def test_grain(self):
        rng = np.random.RandomState(2)
        uc = unitcell.unitcell((4.0, 4.0, 4.0, 90.0, 90.0, 90.0), 225)
        for ubi in strained_ubis(uc, rng, 5):
            g = grain.grain(ubi)
            g.ref_unitcell = uc
            for key in ("tsl", "hsv"):
                rgb = rgb_colors.ubi_to_rgb([ubi], "m-3m", key=key, axis=(0, 1, 0))[0]
                self.assertLess(np.abs(g.get_ipf_colour((0, 1, 0), key) - rgb).max(), 1e-15)

    def test_orix_orien(self):
        """orix orientations go back to the ubi they came from"""
        try:
            import orix  # noqa: F401
        except ImportError:
            self.skipTest("no orix")
        rng = np.random.RandomState(3)
        for cell, sg in [(HEX, 194), ((4, 5, 6, 90, 100, 90), 14), ((5, 5, 5, 70, 70, 70), 166)]:
            uc = unitcell.unitcell(cell, sg)
            U = rotations(rng, 10)
            ubis = np.linalg.inv(np.dot(U, uc.B))
            # through orix's B matrix and quaternions: measured 2.7e-12 (R axes)
            back = uc.get_ubi_from_orix_orien(uc.get_orix_orien_fast(U))
            self.assertLess(np.abs(back - ubis).max() / np.abs(ubis).max(), 1e-11, sg)
            back = uc.get_ubi_from_orix_orien(uc.get_orix_orien(np.dot(U, uc.B)))
            self.assertLess(np.abs(back - ubis).max() / np.abs(ubis).max(), 1e-11, sg)


class TestTensorMap(unittest.TestCase):

    def test_get_ipf_maps(self):
        if sys.version_info[0] == 2:
            self.skipTest("tensor_map does not import on Python 2")
        from ImageD11.sinograms.tensor_map import TensorMap
        rng = np.random.RandomState(4)
        uc = unitcell.unitcell(HEX, 194)
        ubi = strained_ubis(uc, rng, 6).reshape(1, 2, 3, 3, 3)
        ubi[0, 1, 2] = np.nan
        phase = np.zeros((1, 2, 3), int)
        phase[0, 0, 1] = -1
        tm = TensorMap(maps={"UBI": ubi, "phase_ids": phase}, phases={0: uc})
        for key in ("tsl", "hsv"):
            tm.get_ipf_maps(key)
            ok = (phase == 0) & np.isfinite(ubi).all(axis=(-2, -1))
            for i, letter in enumerate("xyz"):
                rgb = rgb_colors.ubi_to_rgb(ubi[ok], "6/mmm", key=key, axis=np.eye(3)[i])
                self.assertLess(np.abs(tm["ipf_" + letter][ok] - rgb).max(), 1e-15)
                self.assertTrue((tm["ipf_" + letter][~ok] == 0).all())


class TestIpfPoints(unittest.TestCase):

    def test_on_the_key(self):
        """a grain's point lands on the pixel of key_image that has its colour"""
        rng = np.random.RandomState(5)
        n = 401
        for cell, name, key in [((4, 4, 4, 90, 90, 90), "m-3m", "tsl"), (HEX, "6/mmm", "hsv"),
                                ((4, 5, 6, 80, 95, 110), "-1", "tsl"), ((4, 5, 6, 100, 90, 90), "2/m11", "hsv"),
                                ((5, 5, 5, 70, 70, 70), "-3mR", "tsl"), ((4, 5, 6, 80, 95, 110), "1", "hsv")]:
            ubi0 = np.linalg.inv(unitcell.unitcell(cell, "P").B)
            ubis = np.einsum("ij,njk->nik", ubi0, rotations(rng, 200))
            images, _ = rgb_colors.key_image(ubis[0], name, key=key, n=n)
            rgb = rgb_colors.ubi_to_rgb(ubis, name, key=key, axis=(0.3, -0.5, 0.8))
            x, y, panel = rgb_colors.ipf_points(ubis, name, key=key, axis=(0.3, -0.5, 0.8))
            col = np.round((x + 1) / 2 * (n - 1)).astype(int)
            row = np.round((y + 1) / 2 * (n - 1)).astype(int)
            pix = np.array([images[min(p, len(images) - 1)][r, c] for p, r, c in zip(panel, row, col)])
            inside = pix[:, 3] == 1
            self.assertGreater(inside.mean(), 0.95, name)
            # the colour changes by about 0.03 across one pixel at n = 401
            self.assertLess(np.abs(pix[inside, :3] - rgb[inside]).max(), 0.05, (name, key))

    def test_plot(self):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            self.skipTest("no matplotlib")
        rng = np.random.RandomState(6)
        ubis = np.einsum("ij,njk->nik", np.linalg.inv(unitcell.unitcell(HEX, "P").B), rotations(rng, 50))
        for key in ("tsl", "hsv"):
            axes = rgb_colors.plot_ipf(ubis, "6mm", key=key, n=51)
            self.assertEqual(sum(len(ax.collections[0].get_offsets()) for ax in axes), 50)
            plt.close("all")


if __name__ == "__main__":
    unittest.main()
