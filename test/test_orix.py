import numpy as np
from orix.vector import Miller
from scipy.spatial.transform import Rotation as R

from ImageD11.grain import grain
from ImageD11.unitcell import unitcell
import unittest


class TestOrixSG(unittest.TestCase):
    def test_no_spacegroup_fails(self):
        # set up B matrix
        cubic_ucell_array = np.array([3., 3., 3., 90., 90., 90.])
        cubic_ucell = unitcell(cubic_ucell_array, symmetry="F")

        B = cubic_ucell.B

        # set up U matrix

        # we want:
        # a* to point in +Z lab
        # b* to point in +X lab
        # c* to point in +Y lab

        U = np.squeeze(R.from_euler('zyx', [-90, 0, -90], degrees=True).as_matrix())

        UBI = np.linalg.inv(np.dot(U, B))

        g = grain(UBI)

        with self.assertRaises(NameError):
            sg = g.spacegroup


class TestOrixPhase(unittest.TestCase):
    def test_phase_recomputed(self):
        # set up B matrix
        cubic_ucell_array = np.array([3., 3., 3., 90., 90., 90.])
        cubic_ucell = unitcell(cubic_ucell_array, symmetry="F")

        B = cubic_ucell.B

        # set up U matrix

        # we want:
        # a* to point in +Z lab
        # b* to point in +X lab
        # c* to point in +Y lab

        U = np.squeeze(R.from_euler('zyx', [-90, 0, -90], degrees=True).as_matrix())

        UBI = np.linalg.inv(np.dot(U, B))

        g = grain(UBI)

        # set a spacegroup
        g.spacegroup = 1

        # set an orix phase
        phase_1 = g.orix_phase

        # change the spacegroup
        g.spacegroup = 225

        # orix phase should be different
        phase_2 = g.orix_phase

        self.assertNotEqual(phase_1, phase_2)


class TestOrixOrien(unittest.TestCase):
    def test_hexagonal(self):
        # set up hexagonal B matrix
        hex_ucell_array = np.array([3, 3, 4, 90., 90., 120.])
        hex_ucell = unitcell(hex_ucell_array, symmetry="P")

        B = hex_ucell.B

        # set up U matrix

        # we want:
        # a* to point in +Z lab
        # b* to point in +X lab
        # c* to point in +Y lab

        U = np.squeeze(R.from_euler('zyx', [-90, 0, -90], degrees=True).as_matrix())

        UBI = np.linalg.inv(np.dot(U, B))

        g = grain(UBI)

        # set a spacegroup
        g.spacegroup = 194

        # a* should be pointing in lab Z

        # get direction of (100) normal in cartesian lab frame as Numpy array
        ra_direc_lab = np.array((~g.orix_orien * Miller(hkl=(1,0,0), phase=g.orix_phase)).xyz).T
        lab_z = np.array([0, 0, 1])

        ra_direc_lab_unit = ra_direc_lab / np.linalg.norm(ra_direc_lab)

        self.assertTrue(np.allclose(ra_direc_lab_unit, lab_z))


class TestIPFCol(unittest.TestCase):
    def test_cubic(self):
        # set up B matrix
        cubic_ucell_array = np.array([3., 3., 3., 90., 90., 90.])
        cubic_ucell = unitcell(cubic_ucell_array, symmetry="F")

        B = cubic_ucell.B

        # set up U matrix

        # we want:
        # a* to point in +Z lab
        # b* to point in +X lab
        # c* to point in +Y lab

        U = np.squeeze(R.from_euler('zyx', [-90, 0, -90], degrees=True).as_matrix())

        UBI = np.linalg.inv(np.dot(U, B))

        g = grain(UBI)

        # set the spacegroup
        g.spacegroup = 225

        # each of the 100, 010, 001, etc directions should have red colours (cubic symmetry)

        axes = np.array([
            [1., 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [-1, 0, 0],
            [0, -1, 0],
            [0, 0, -1]
        ])

        red = np.array([[1., 0., 0.]])

        for axis in axes:
            rgb = g.get_ipf_colour(axis=axis)
            self.assertTrue(np.allclose(rgb, red))

    def test_cubic_110(self):
        cubic_ucell_array = np.array([3., 3., 3., 90., 90., 90.])
        cubic_ucell = unitcell(cubic_ucell_array, symmetry="F")

        B = cubic_ucell.B

        # set up U matrix

        # we want:
        # some (110) normal to point along X
        # rotate 45 degrees in Z
        # should mean 0.5 e1 + 0.5 e2 + 0.0 e3

        U = np.squeeze(R.from_euler('z', [44], degrees=True).as_matrix())

        UBI = np.linalg.inv(np.dot(U, B))

        g = grain(UBI)

        # set the spacegroup
        g.spacegroup = 225

        # 110 should be pointing along X
        # so in an IPF-X plot, the colour should be green

        rgb = g.get_ipf_colour(axis=np.array([1, 0, 0]))
        green = np.array([0., 1., 0.])

        self.assertTrue(np.allclose(rgb, green, atol=0.05))

    def test_cubic_111(self):
        # set up B matrix
        cubic_ucell_array = np.array([3., 3., 3., 90., 90., 90.])
        cubic_ucell = unitcell(cubic_ucell_array, symmetry="F")

        B = cubic_ucell.B

        # set up U matrix

        # we want:
        # some (111) normal to point along X
        # rotate 45 degrees in Z
        # should mean 0.5 e1 + 0.5 e2 + 0.0 e3

        U = np.squeeze(R.from_euler('zy', [45, 45], degrees=True).as_matrix())

        UBI = np.linalg.inv(np.dot(U, B))

        g = grain(UBI)

        # set the spacegroup
        g.spacegroup = 225

        # 111 should be pointing along X
        # so in an IPF-X plot, the colour should be blue-ish

        rgb = g.get_ipf_colour(axis=np.array([1, 0, 0]))
        blue = np.array([0., 0., 1.])

        self.assertTrue(np.allclose(rgb, blue, atol=0.3))

    def test_hexagonal(self):
        # set up hexagonal B matrix
        hex_ucell_array = np.array([3, 3, 4, 90., 90., 120.])
        hex_ucell = unitcell(hex_ucell_array, symmetry="P")

        B = hex_ucell.B

        # set up U matrix

        # we want:
        # a* to point in +Z lab
        # b* to point in +X lab
        # c* to point in +Y lab

        U = np.squeeze(R.from_euler('zyx', [-90, 0, -90], degrees=True).as_matrix())

        UBI = np.linalg.inv(np.dot(U, B))

        g = grain(UBI)

        # set a spacegroup
        g.spacegroup = 194

        # on an IPF-Z:
        # a* should be parallel to Z
        # (100) (10-10) should be along Z
        # on TSL, (10-10) is blue

        rgb = g.get_ipf_colour(axis=np.array([0, 0, 1]))
        blue = np.array([0., 0., 1.])

        self.assertTrue(np.allclose(rgb, blue, atol=0.3))

        # on an IPF-Y, c* should be parallel to Y
        # c* is (0001) normal

        rgb = g.get_ipf_colour(axis=np.array([0, 1, 0]))
        red = np.array([1., 0., 0.])

        self.assertTrue(np.allclose(rgb, red, atol=0.3))
