from __future__ import print_function, division

"""Issue 664: with cut=0 the top_pixels thresholds were all 0, so a busy
frame could not be cut down to howmany pixels."""

import sys
import unittest

if sys.version_info[0] == 2 and "h5py" in sys.modules:
    # hdf5plugin 1.4.1 (setup.py pin for python 2) will not import after h5py
    raise unittest.SkipTest("hdf5plugin 1.4.1 cannot import after h5py")

import numpy as np

from ImageD11.sinograms import lima_segmenter


class TestCutZero(unittest.TestCase):
    def test_thresholds(self):
        opts = lima_segmenter.SegmenterOptions(cut=0)
        opts.setup()
        self.assertEqual(opts.thresholds, (1, 2, 4, 8, 16, 32))
        opts = lima_segmenter.SegmenterOptions(cut=3)
        opts.setup()
        self.assertEqual(opts.thresholds, (3, 6, 12, 24, 48, 96))

    def test_top_pixels(self):
        opts = lima_segmenter.SegmenterOptions(cut=0, howmany=10)
        opts.setup()
        img = np.ones((8, 8), np.uint16)
        img[0, :5] = 50  # 5 strong pixels in a frame of 64 above cut
        n = img.size
        row = np.empty(n, np.uint16)
        col = np.empty(n, np.uint16)
        val = np.empty(n, img.dtype)
        nnz = lima_segmenter.select(img, 1, row, col, val, opts.cut)
        self.assertEqual(nnz, 64)
        nnz = lima_segmenter.top_pixels(
            nnz, row, col, val, opts.howmany, opts.thresholds
        )
        self.assertEqual(nnz, 5)
        self.assertTrue((val[:nnz] == 50).all())


if __name__ == "__main__":
    unittest.main()
