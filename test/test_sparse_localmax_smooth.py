"""
Checks that sparseframe.sparse_localmax(frame, smooth=...) reproduces the
per-frame labelling that SparseScan.lmlabel(smooth=...) already does
internally, for both the smoothed and unsmoothed paths.
"""
from __future__ import print_function

import os
import unittest

import numpy as np

from ImageD11 import sparseframe

SPARSEFILE = os.path.join(
    os.path.dirname(__file__), "pixelmapper", "silicon_fullscan_strong_sparse.h5"
)


class TestSparseLocalmaxSmooth(unittest.TestCase):
    def setUp(self):
        self.scan = sparseframe.SparseScan(SPARSEFILE, "1.1")
        # a handful of frames with plenty of pixels, so localmax has real work to do
        self.frame_ids = np.argsort(self.scan.nnz)[-5:]
        for i in self.frame_ids:
            self.assertGreater(self.scan.nnz[i], 0)

    def check_smooth_option(self, smooth):
        # path A: label the whole scan with SparseScan.lmlabel (countall=False
        # so labels restart from 1 on each frame, matching a single sparse_frame)
        scan = sparseframe.SparseScan(SPARSEFILE, "1.1")
        scan.lmlabel(smooth=smooth, countall=False)

        for i in self.frame_ids:
            # path B: pull out a single frame and label it with sparse_localmax
            frame = self.scan.getframe(int(i))
            nlabel = sparseframe.sparse_localmax(frame, smooth=smooth)

            labels_a = scan.labels[scan.ipt[i] : scan.ipt[i + 1]]
            labels_b = frame.pixels["localmax"]

            self.assertGreater(nlabel, 0)
            self.assertEqual(scan.nlabels[i], nlabel)
            np.testing.assert_array_equal(labels_a, labels_b)

    def test_smooth_true(self):
        self.check_smooth_option(True)

    def test_smooth_false(self):
        self.check_smooth_option(False)


if __name__ == "__main__":
    unittest.main()
