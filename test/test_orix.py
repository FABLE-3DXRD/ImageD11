import numpy as np
from scipy.spatial.transform import Rotation as R

from ImageD11.grain import grain
from ImageD11.unitcell import unitcell
import unittest


class TestIPFCol(unittest.TestCase):
    def test_cubic(self):
        # set up B matrix
        cubic_ucell_array = np.array([3., 3., 3., 90., 90., 90.])
        cubic_ucell = unitcell(cubic_ucell_array, symmetry="F")

        B = cubic_ucell.B

        # set up U matrix

        # we want:
        # a to point in +Z
        # b to point in +X
        # c to point in +Y

        U = R.from_euler('zyx', [-90, 0, -90], degrees=True).as_matrix()

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
        # set up B matrix
        cubic_ucell_array = np.array([3., 3., 3., 90., 90., 90.])
        cubic_ucell = unitcell(cubic_ucell_array, symmetry="F")

        B = cubic_ucell.B

        # set up U matrix

        # we want:
        # some 110 direction to point along X
        # rotate 45 degrees in Z
        # should mean 0.5 e1 + 0.5 e2 + 0.0 e3

        U = R.from_euler('z', [45], degrees=True).as_matrix()

        UBI = np.linalg.inv(np.dot(U, B))

        g = grain(UBI)

        # set the spacegroup
        g.spacegroup = 225

        # 110 should be pointing along X
        # so in an IPF-X plot, the colour should be green

        rgb = g.get_ipf_colour(axis=np.array([1, 0, 0]))
        green = np.array([0., 1., 0.])

        self.assertTrue(np.allclose(rgb, green, atol=0.01))

    def test_cubic_111(self):
        # set up B matrix
        cubic_ucell_array = np.array([3., 3., 3., 90., 90., 90.])
        cubic_ucell = unitcell(cubic_ucell_array, symmetry="F")

        B = cubic_ucell.B

        # set up U matrix

        # we want:
        # some 111 direction to point along X
        # rotate 45 degrees in Z
        # should mean 0.5 e1 + 0.5 e2 + 0.0 e3

        U = R.from_euler('zy', [45, 45], degrees=True).as_matrix()

        UBI = np.linalg.inv(np.dot(U, B))

        g = grain(UBI)

        # set the spacegroup
        g.spacegroup = 225

        # 111 should be pointing along X
        # so in an IPF-X plot, the colour should be blue-ish

        rgb = g.get_ipf_colour(axis=np.array([1, 0, 0]))
        blue = np.array([0., 0., 1.])

        self.assertTrue(np.allclose(rgb, blue, atol=0.3))