import sys
import unittest

if int(sys.version_info.major) == 2:
    raise unittest.SkipTest("Skipping tensor_map tests on Python 2")

import numpy as np
from ImageD11 import unitcell
from ImageD11.sinograms.tensor_map import TensorMap
from scipy.spatial.transform import Rotation

try:
    import orix
except ImportError:
    orix = None


def make_map(nan_voxel=None):
    cell = unitcell.unitcell((4.05, 4.05, 4.05, 90, 90, 90), 225)
    U = Rotation.random(6, random_state=42).as_matrix()
    ubi = np.array([np.linalg.inv(np.dot(u, cell.B)) for u in U])
    UBI = ubi.reshape(1, 2, 3, 3, 3)
    if nan_voxel is not None:
        UBI[nan_voxel] = np.nan
    maps = {"UBI": UBI, "phase_ids": np.zeros((1, 2, 3), int)}
    return TensorMap(maps=maps, phases={0: cell})


class test_ipf_nan_voxels(unittest.TestCase):
    """A voxel with a nan UBI must not shift the colours of the others (#654)"""

    def check(self, backend):
        ref = make_map()
        ref.get_ipf_maps(backend=backend)
        tm = make_map(nan_voxel=(0, 0, 1))
        tm.get_ipf_maps(backend=backend)
        for letter in "xyz":
            a = ref.maps["ipf_" + letter]
            b = tm.maps["ipf_" + letter]
            self.assertTrue(np.all(b[0, 0, 1] == 0))
            good = np.ones((1, 2, 3), bool)
            good[0, 0, 1] = False
            self.assertTrue(np.allclose(a[good], b[good]), letter)

    def test_rgb_colors(self):
        self.check("rgb_colors")

    @unittest.skipIf(orix is None, "orix not installed")
    def test_orix(self):
        self.check("orix")


if __name__ == "__main__":
    unittest.main()
