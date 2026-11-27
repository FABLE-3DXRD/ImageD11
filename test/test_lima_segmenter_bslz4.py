from __future__ import print_function, division

"""lima_segmenter.reader gives the same sparse frames with bslz4_to_sparse
(compressed chunks, or frames with a background) as with the numpy code,
when top_pixels has to cut frames."""

import os
import shutil
import tempfile
import unittest

import numpy as np
import h5py
import hdf5plugin

from ImageD11.sinograms import lima_segmenter


def run_reader(frms, options, use_bslz4):
    lima_segmenter.OPTIONS = options
    saved = lima_segmenter.chunk2sparse, lima_segmenter.CHUNK2SPARSE_HOWMANY
    if not use_bslz4:
        lima_segmenter.chunk2sparse = None
        lima_segmenter.CHUNK2SPARSE_HOWMANY = False
    try:
        return list(lima_segmenter.reader(frms, options.mask, options.cut))
    finally:
        lima_segmenter.chunk2sparse, lima_segmenter.CHUNK2SPARSE_HOWMANY = saved


@unittest.skipIf(lima_segmenter.chunk2sparse is None, "needs bslz4_to_sparse")
@unittest.skipIf(not hasattr(hdf5plugin, "Bitshuffle"), "needs hdf5plugin.Bitshuffle")
class TestReaderBslz4(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="id11_lima_segmenter_")
        rng = np.random.RandomState(42)
        nf, ni, nj = 10, 64, 128
        frames = rng.poisson(0.5, (nf, ni, nj)).astype(np.uint16)
        for f in frames:  # some bright 3x3 spots
            for i, j in zip(rng.randint(2, ni - 2, 20), rng.randint(2, nj - 2, 20)):
                spot = rng.randint(5, 500, (3, 3)).astype(np.uint16)
                f[i - 1 : i + 2, j - 1 : j + 2] += spot
        self.h5name = os.path.join(self.tmpdir, "frames.h5")
        with h5py.File(self.h5name, "w") as h:
            h.create_dataset(
                "data", data=frames, chunks=(1, ni, nj), **hdf5plugin.Bitshuffle()
            )
        self.options = lima_segmenter.SegmenterOptions(cut=1, howmany=300)
        self.options.setup()
        self.options.mask = np.ones((ni, nj), np.uint8)
        self.options.mask[10:20, 30:40] = 0
        # poisson(0.5) puts about 9% of 8192 pixels above 1: top_pixels must cut
        self.assertGreater((frames * self.options.mask > 1).sum(axis=(1, 2)).min(), 300)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def compare(self):
        with h5py.File(self.h5name, "r") as h:
            got = run_reader(h["data"], self.options, use_bslz4=True)
            ref = run_reader(h["data"], self.options, use_bslz4=False)
        self.assertEqual(len(got), len(ref))
        for a, b in zip(got, ref):
            self.assertTrue(0 < a.nnz <= self.options.howmany)
            self.assertEqual(a.nnz, b.nnz)
            self.assertTrue((a.row == b.row).all())
            self.assertTrue((a.col == b.col).all())
            self.assertEqual(a.pixels["intensity"].dtype, b.pixels["intensity"].dtype)
            self.assertTrue((a.pixels["intensity"] == b.pixels["intensity"]).all())

    def test_chunks(self):
        self.compare()

    def test_frames_with_background(self):
        self.options.bg = np.random.RandomState(1).uniform(0, 2, (64, 128))
        self.compare()

    def test_frames_with_3d_background(self):
        self.options.bg = np.random.RandomState(1).uniform(0, 2, (10, 64, 128))
        self.compare()


if __name__ == "__main__":
    unittest.main()
