from __future__ import print_function, division

import os
import tempfile
import unittest

import h5py
import numpy as np

from ImageD11.sparseframe import SparseScan
from ImageD11.sinograms.properties import props, pks_table


def write_sparse(hname, frames, shape=(16, 16)):
    """Minimal sparse file: one group per frame set, pixels row-major sorted.

    frames = list of [(row, col, intensity), ...], one list per frame
    """
    nnz = np.array([len(f) for f in frames], np.uint32)
    row = np.concatenate([[p[0] for p in f] for f in frames]).astype(np.uint16)
    col = np.concatenate([[p[1] for p in f] for f in frames]).astype(np.uint16)
    val = np.concatenate([[p[2] for p in f] for f in frames]).astype(np.uint32)
    with h5py.File(hname, "w") as h:
        g = h.create_group("1.1")
        g.attrs["nframes"] = len(frames)
        g.attrs["shape0"] = shape[0]
        g.attrs["shape1"] = shape[1]
        g.attrs["itype"] = val.dtype.name
        g["nnz"] = nnz
        g["row"] = row
        g["col"] = col
        g["intensity"] = val
        g["measurement/rot"] = np.arange(len(frames), dtype=float)


class TestPropsMax(unittest.TestCase):
    """The peak maximum must be recorded before wtmax clips the pixel values."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.hname = os.path.join(self.tmp, "sparse.h5")
        # two frames, one three-pixel blob each; the second contains a pixel
        # that wtmax will clip
        write_sparse(self.hname,
                     [[(5, 5, 10), (5, 6, 20), (5, 7, 10)],
                      [(5, 5, 10), (5, 6, 100), (5, 7, 10)]])

    def test_max_is_recorded_before_clipping(self):
        scan = SparseScan(self.hname, "1.1")
        r, pairs = props(scan, 0, algorithm="cplabel", wtmax=50)

        self.assertEqual(r.shape[0], 6, "pk_props should have six rows")
        # true maxima, unaffected by wtmax
        self.assertEqual(r[5, 0], 20)
        self.assertEqual(r[5, 1], 100)
        # sum_intensity does see the clipping: 10 + 50 + 10
        self.assertEqual(r[1, 0], 40)
        self.assertEqual(r[1, 1], 70)

    def test_max_without_wtmax(self):
        scan = SparseScan(self.hname, "1.1")
        r, pairs = props(scan, 0, algorithm="cplabel", wtmax=None)
        self.assertEqual(r[5, 1], 100)
        self.assertEqual(r[1, 1], 120)


class TestPkTableIMax(unittest.TestCase):
    """IMax_int comes out of pk2d/pk2dmerge, and old five-row tables still work."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        # three 2D peaks on two frames, the first two merging into one spot
        self.pk_props6 = np.array([[3, 3, 4],            # s1
                                   [40, 70, 50],         # sI
                                   [200, 350, 250],      # srI
                                   [240, 420, 300],      # scI
                                   [0, 1, 1],            # frm
                                   [20, 100, 30]],       # mxI
                                  dtype=np.int64)
        self.glabel = np.array([0, 0, 1], np.int64)
        self.omega = np.array([[0.0, 1.0]])
        self.dty = np.array([[0.0, 0.0]])

    def table(self, pk_props):
        t = pks_table(ipk=np.array([0, pk_props.shape[1]], np.int64),
                      pk_props=pk_props, glabel=self.glabel, nlabel=2)
        t.npk = np.array([[pk_props.shape[1], 0, 0]], np.int64)
        return t

    def test_pk2d_reports_imax(self):
        p = self.table(self.pk_props6).pk2d(self.omega, self.dty)
        self.assertIn("IMax_int", p)
        self.assertTrue((p["IMax_int"] == [20, 100, 30]).all())

    def test_pk2dmerge_takes_the_largest(self):
        p = self.table(self.pk_props6).pk2dmerge(self.omega, self.dty)
        self.assertIn("IMax_int", p)
        # spot 0 merges the 20 and the 100
        self.assertEqual(p["IMax_int"][0], 100)
        self.assertEqual(p["IMax_int"][1], 30)

    def test_five_row_table_still_works(self):
        """Peaks tables written before this change have no sixth row."""
        t = self.table(self.pk_props6[:5])
        p = t.pk2d(self.omega, self.dty)
        self.assertNotIn("IMax_int", p)
        self.assertEqual(len(p["s_raw"]), 3)
        self.assertNotIn("IMax_int", t.pk2dmerge(self.omega, self.dty))

    def test_roundtrip_through_h5(self):
        hname = os.path.join(self.tmp, "pks.h5")
        self.table(self.pk_props6).save(hname)
        back = pks_table.load(hname)
        self.assertEqual(back.pk_props.shape[0], 6)
        self.assertTrue((back.pk2d(self.omega, self.dty)["IMax_int"]
                         == [20, 100, 30]).all())


if __name__ == "__main__":
    unittest.main()
