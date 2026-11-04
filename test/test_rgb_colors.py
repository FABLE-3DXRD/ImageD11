"""IPF colours: check colours (ipf_oracle.json, from test/ipf_keys/make_shipped.py): exact, and MTEX 7.2,
one colour for every equivalent ubi, rhombohedral == hexagonal axes, the strain warning."""

import json
import os
import unittest

import numpy as np

from ImageD11 import rgb_colors, sym_u
from ImageD11.ipf_keys import GROUPS

HERE = os.path.dirname(os.path.abspath(__file__))
T = np.array([[1, -1, 0], [0, 1, -1], [1, 1, 1]])   # ubi_H = T.ubi_R (obverse)
R_OF_H = {"3": "3R", "-3": "-3R", "321": "32R", "3m1": "3mR", "-3m1": "-3mR"}


def cell_ubi(cell):
    """a ubi (rows a, b, c) of cell a b c alpha beta gamma (deg)"""
    a, b, c = cell[:3]
    al, be, ga = np.radians(cell[3:])
    g = np.array([[a * a, a * b * np.cos(ga), a * c * np.cos(be)],
                  [a * b * np.cos(ga), b * b, b * c * np.cos(al)],
                  [a * c * np.cos(be), b * c * np.cos(al), c * c]])
    return np.linalg.cholesky(g)


def rotations(rng, n):
    q = np.array([np.linalg.qr(m)[0] for m in rng.normal(size=(n, 3, 3))])
    return q * np.sign(np.linalg.det(q))[:, None, None]


def turn_to_z(v, spin):
    """a rotation taking v onto +z, then turning by spin about z"""
    z = v / np.linalg.norm(v)
    x = np.cross(z, [0.3, 0.5, 0.7])
    x /= np.linalg.norm(x)
    Q = np.array([x, np.cross(z, x), z])
    c, s = np.cos(spin), np.sin(spin)
    return np.dot([[c, -s, 0], [s, c, 0], [0, 0, 1]], Q)


def proper(ops):
    return [o for o in ops if np.linalg.det(o) > 0]


class TestOracle(unittest.TestCase):

    def test_mtex_check_colours(self):
        """each check direction of each key, put along z by turning the grain"""
        with open(os.path.join(HERE, "ipf_oracle.json")) as f:
            oracle = json.load(f)
        rng = np.random.RandomState(1)
        for name, o in oracle.items():
            ubi0 = cell_ubi(o["test_cell"])
            ubis = []
            for u in o["check_uvw"]:
                Q = turn_to_z(np.dot(ubi0.T, u), rng.uniform(0, 2 * np.pi))
                ubis.append(np.dot(ubi0, Q.T))
            for kind in ("tsl", "hsv"):
                for k in o[kind]:
                    for view in k["views"]:
                        rgb = rgb_colors.ubi_to_rgb(ubis, name, view, kind)
                        self.assertLess(np.abs(rgb - k["rgb"]).max(), 1e-7, (name, kind, view))
                        # MTEX sums the hue integral in 1000 steps (3e-5), with 4 corners
                        # takes the mean of a 1 deg grid as the centre (5e-3), and its HSV
                        # hue correction starts at F(0) = 1.7e-3 (up to 6 x 1.7e-3 in rgb)
                        mtex = 1e-2 if (kind == "tsl" and name in ("m-3", "23")) else 1e-4
                        mtex = 1.1e-2 if kind == "hsv" else mtex
                        self.assertLess(np.abs(rgb - k["mtex_rgb"]).max(), mtex, (name, kind, view))
        self.assertEqual(len(oracle), len(GROUPS))

    def test_sym_u_groups(self):
        for f, name in [(sym_u.cubic, "432"), (sym_u.hexagonal, "622"), (sym_u.trigonal, "321"),
                        (sym_u.rhombohedralP, "32R"), (sym_u.tetragonal, "422"),
                        (sym_u.orthorhombic, "222"), (sym_u.monoclinic_b, "121"),
                        (sym_u.triclinic, "1")]:
            self.assertEqual(rgb_colors.group_name(f()), name)
        for name in GROUPS:
            self.assertEqual(rgb_colors.group_name(rgb_colors.point_group(name)), name)


class TestInvariance(unittest.TestCase):

    def test_equivalent_ubis(self):
        """a strained grain, indexed every equivalent way, has one colour"""
        rng = np.random.RandomState(2)
        for name in GROUPS:
            with open(os.path.join(HERE, "ipf_oracle.json")) as f:
                cell = json.load(f)[name]["test_cell"]
            eps = rng.uniform(-3e-3, 3e-3, (3, 3))
            ubi = np.dot(np.dot(cell_ubi(cell), np.eye(3) + (eps + eps.T) / 2), rotations(rng, 1)[0])
            ops = rgb_colors.point_group(name).group
            for kind in ("tsl", "hsv"):
                eq = proper(list(ops) + [-o for o in ops]) if kind == "tsl" else proper(ops)
                ubis = [np.dot(o, ubi) for o in eq]
                for view in set(v for k in GROUPS[name][kind] for v in k["views"]):
                    rgb = rgb_colors.ubi_to_rgb(ubis, name, view, kind, axis=rng.normal(size=3))
                    self.assertLess(np.abs(rgb - rgb[0]).max(), 1e-9, (name, kind, view))

    def test_rhombohedral_axes(self):
        """an R crystal on rhombohedral axes has the colours of its hexagonal axes"""
        rng = np.random.RandomState(3)
        ubi0 = np.dot(cell_ubi([4.9, 4.9, 13.6, 90, 90, 120]), np.eye(3) + 1e-3 * rng.normal(size=(3, 3)))
        ubi_h = np.einsum("ij,nkj->nik", ubi0, rotations(rng, 50))
        ubi_r = np.einsum("ij,njk->nik", np.linalg.inv(T), ubi_h)
        for h, r in R_OF_H.items():
            for kind in ("tsl", "hsv"):
                for view in set(v for k in GROUPS[h][kind] for v in k["views"]):
                    a = rgb_colors.ubi_to_rgb(ubi_h, h, view, kind)
                    b = rgb_colors.ubi_to_rgb(ubi_r, r, view, kind)
                    self.assertLess(np.abs(a - b).max(), 1e-9, (h, kind, view))


class TestHueTable(unittest.TestCase):
    """the even-hue table is made once per key, from the first grain"""

    def figure(self, k, gs):
        """dot products of the unit vectors the colour uses: planes, corners,
        centre and hue reference, (n,m,m)"""
        L = np.linalg.cholesky(gs)
        Li = np.linalg.inv(L)
        N, V, edges = rgb_colors._region(L, Li, k.planes if k.kind == "hsv" else k.normals)
        c = rgb_colors._centre(N, V, edges)
        for s in N[:, len(k.normals):].transpose(1, 0, 2):
            c = np.where((rgb_colors._dot(c, s) < 1e-5)[:, None], c - 2 * rgb_colors._dot(c, s)[:, None] * s, c)
        rx = rgb_colors._vec(L, Li, k.oos) - c
        flat = np.sqrt(rgb_colors._dot(rx, rx) - rgb_colors._dot(rx, c) ** 2) < 1e-10
        rx[flat] = (rgb_colors._vec(L, Li, k.hue_ref) - c)[flat]
        X = np.concatenate([N, c[:, None], rgb_colors._unit(rx - rgb_colors._dot(rx, c)[:, None] * c)[:, None]] +
                           ([] if V is None else [V]), axis=1)
        return np.einsum("nai,nbi->nab", X, X)

    def test_same_figure(self):
        """the sector is the same figure in every metric its group allows"""
        with open(os.path.join(HERE, "ipf_oracle.json")) as f:
            oracle = json.load(f)
        rng = np.random.RandomState(4)
        for name in GROUPS:
            eps = rng.uniform(-0.2, 0.2, (50, 3, 3))
            ubis = np.einsum("ij,njk->nik", cell_ubi(oracle[name]["test_cell"]),
                             np.eye(3) + (eps + eps.transpose(0, 2, 1)) / 2)
            gs, _ = rgb_colors.symmetric_metric(ubis, rgb_colors.point_group(name))
            for kind in ("tsl", "hsv"):
                for k in GROUPS[name][kind]:
                    G = self.figure(rgb_colors._Key(name, kind, k["views"][0]), gs)
                    self.assertLess(np.abs(G - G[:1]).max(), 1e-12, (name, kind, k["views"]))

    def test_first_grain(self):
        """the colours do not depend on which grain made the table"""
        rng = np.random.RandomState(6)
        for name, cell in [("6/mmm", [3, 3, 5, 90, 90, 120]), ("112/m", [4, 5, 6, 90, 90, 100]),
                           ("m-3", [4, 4, 4, 90, 90, 90])]:
            eps = rng.uniform(-3e-3, 3e-3, (100, 3, 3))
            ubis = np.einsum("ij,njk,nkl->nil", cell_ubi(cell), np.eye(3) + (eps + eps.transpose(0, 2, 1)) / 2,
                             rotations(rng, 100))
            for kind in ("tsl", "hsv"):
                a = rgb_colors.ubi_to_rgb(ubis, name, key=kind)
                b = rgb_colors.ubi_to_rgb(ubis[::-1], name, key=kind)[::-1]
                self.assertLess(np.abs(a - b).max(), 1e-12, (name, kind))


class TestKeyImage(unittest.TestCase):

    def test_area(self):
        """the sector drawn in equal area covers 4 pi / (operations folded with)"""
        with open(os.path.join(HERE, "ipf_oracle.json")) as f:
            oracle = json.load(f)
        n = 101
        for name in GROUPS:
            ubi = cell_ubi(oracle[name]["test_cell"])
            for kind in ("tsl", "hsv"):
                for view in [k["views"][0] for k in GROUPS[name][kind]]:
                    imgs, _ = rgb_colors.key_image(ubi, name, view, kind, "equal area", n)
                    area = sum(im[..., 3].sum() for im in imgs) * (2. / n) ** 2 * 2
                    ops = rgb_colors._Key(name, kind, view).ops
                    nops = len(set(tuple(o.ravel()) for o in ops))
                    # pixels on the boundary: 5.4 % at worst for n = 101
                    self.assertLess(abs(area * nops / (4 * np.pi) - 1), 0.07, (name, kind, view))

    def test_rhombohedral_axes(self):
        ubi_h = cell_ubi([4.9, 4.9, 13.6, 90, 90, 120])
        ubi_r = np.dot(np.linalg.inv(T), ubi_h)
        for h, r in R_OF_H.items():
            for kind in ("tsl", "hsv"):
                for view in set(v for k in GROUPS[h][kind] for v in k["views"]):
                    a, la = rgb_colors.key_image(ubi_h, h, view, kind, n=100)
                    b, lb = rgb_colors.key_image(ubi_r, r, view, kind, n=100)
                    self.assertEqual(len(a), len(b))
                    for x, y in zip(a, b):
                        self.assertLess(np.abs(x - y).max(), 1e-9, (h, kind, view))
                    self.assertEqual(len(la), len(lb))

    def test_labels(self):
        _, labels = rgb_colors.key_image(cell_ubi([4, 4, 4, 90, 90, 90]), "m-3m", n=11)
        self.assertEqual(sorted(l[3] for l in labels), ["[0 0 1]", "[1 0 1]", "[1 1 1]"])


class TestStrain(unittest.TestCase):

    def test_symmetric_metric(self):
        cubic = rgb_colors.point_group("m-3m")
        ubi = cell_ubi([4, 4, 4, 90, 90, 90])
        gs, e = rgb_colors.symmetric_metric(ubi, cubic)
        self.assertLess(np.abs(e).max(), 1e-15)
        # c/a is allowed in tetragonal, not in cubic: e = (eps, eps, -2 eps) about the mean
        ubi = np.dot(ubi, np.diag([1 - 2e-3, 1 - 2e-3, 1 + 4e-3]))
        gs, e = rgb_colors.symmetric_metric(ubi, cubic)
        self.assertAlmostEqual(np.abs(e).max(), 4e-3, 4)
        self.assertLess(abs(e.sum()), 1e-12)
        gs, e = rgb_colors.symmetric_metric(ubi, rgb_colors.point_group("4/mmm"))
        self.assertLess(np.abs(e).max(), 1e-15)

    def test_warning(self):
        if not hasattr(self, "assertLogs"):
            self.skipTest("no assertLogs")
        ubi = np.dot(cell_ubi([4, 4, 4, 90, 90, 90]), np.diag([1, 1, 1.02]))
        with self.assertLogs(level="WARNING"):
            rgb_colors.ubi_to_rgb([ubi], "m-3m")


if __name__ == "__main__":
    unittest.main()
