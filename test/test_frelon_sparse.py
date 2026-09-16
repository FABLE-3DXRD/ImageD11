"""
Tests for frelon_peaksearch.segment_dataset_to_sparse: the frelon segmenter
writing a sparse pixel file that the existing sparse readers understand.

The synthetic bliss dataset is adapted from the work of @jadball in PR #633.
"""
from __future__ import print_function

import os
import shutil
import sys
import tempfile
import unittest

import h5py
import numpy as np

# python2 lacks functools.lru_cache, which frelon_peaksearch uses
IMPORT_ERROR = None
try:
    import ImageD11.sinograms.dataset
    import ImageD11.sparseframe
    import ImageD11.frelon_peaksearch as frelon_peaksearch
except Exception as e:  # pragma: no cover
    IMPORT_ERROR = e
SKIP = (sys.version_info[0] < 3) or (IMPORT_ERROR is not None)
REASON = "needs python3 (import error: %s)" % (IMPORT_ERROR,)

NFRAMES = 12
SHAPE = (96, 96)
OSTEP = 1.0
DTYS = (-0.5, 0.0)
# row, col, first frame, last frame inclusive
SPOTS = [
    (30.0, 30.0, 2, 6),
    (70.0, 70.0, 1, 3),
    (60.0, 20.0, 3, 5),
    (60.0, 26.0, 6, 8),
]
WORKER_ARGS = {
    "bgfile": None,
    "maskfile": None,
    "darkfile": None,
    "flatfile": None,
    "threshold": 70,
    "smoothsigma": 1.0,
    "bgc": 0.9,
    "minpx": 3,
    "m_offset_thresh": 100,
    "m_ratio_thresh": 150,
}


def make_frames(iscan):
    r, c = np.mgrid[: SHAPE[0], : SHAPE[1]]
    frames = np.full((NFRAMES,) + SHAPE, 100.0)
    for k, (r0, c0, f0, f1) in enumerate(SPOTS):
        for f in range(f0, f1 + 1):
            dr = 0.1 * np.sin(f + k + iscan)
            amp = 2000.0 * (1 + 0.2 * np.cos(f + iscan))
            frames[f] += amp * np.exp(
                -0.5 * ((r - r0 - dr) ** 2 + (c - c0) ** 2) / 1.5 ** 2)
    return frames.round().astype(np.uint16)


def make_bliss(dataroot, sample, dset, nscans):
    """Minimal bliss layout: masterfile with VDS frames pointing at lima files"""
    dsname = sample + "_" + dset
    path = os.path.join(dataroot, sample, dsname)
    os.makedirs(path)
    master = os.path.join(path, dsname + ".h5")
    with h5py.File(master, "w") as hm:
        for i in range(nscans):
            scan = "%d.1" % (i + 1)
            limaname = "scan%04d/frelon3_0000.h5" % (i + 1)
            os.makedirs(os.path.join(path, "scan%04d" % (i + 1)))
            with h5py.File(os.path.join(path, limaname), "w") as hl:
                hl["/entry_0000/measurement/data"] = make_frames(i)
            g = hm.create_group(scan)
            g["title"] = "fscan diffrz 0 1 %d 0.1" % NFRAMES
            m = g.create_group("measurement")
            m["diffrz"] = np.arange(NFRAMES) * OSTEP + 0.5 * OSTEP
            g["instrument/positioners/diffty"] = DTYS[i]
            layout = h5py.VirtualLayout(shape=(NFRAMES,) + SHAPE, dtype=np.uint16)
            vsrc = h5py.VirtualSource(
                limaname, "/entry_0000/measurement/data", shape=(NFRAMES,) + SHAPE)
            layout[0:NFRAMES] = vsrc[0:NFRAMES]
            m.create_virtual_dataset("frelon3", layout)
    return master


def make_dataset(tmp, dset, nscans):
    dataroot = os.path.join(tmp, "RAW")
    analysisroot = os.path.join(tmp, "PROCESSED")
    make_bliss(dataroot, "S", dset, nscans)
    ds = ImageD11.sinograms.dataset.DataSet(
        dataroot=dataroot,
        analysisroot=analysisroot,
        sample="S",
        dset=dset,
        detector="frelon3",
        omegamotor="diffrz",
        dtymotor="diffty",
    )
    ds.import_all()
    ds.save()
    return ds


@unittest.skipIf(SKIP, REASON)
class TestSegmentToSparse(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="id11_frelon_sparse_")
        cls.ds = make_dataset(cls.tmp, "d", len(DTYS))
        cls.sparsefile = frelon_peaksearch.segment_dataset_to_sparse(
            cls.ds, WORKER_ARGS, num_cpus=1)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_file_was_written(self):
        self.assertTrue(os.path.exists(self.sparsefile))
        self.assertEqual(self.sparsefile, self.ds.sparsefile)

    def test_layout_matches_the_sparse_readers(self):
        with h5py.File(self.sparsefile, "r") as h:
            self.assertIn("h5input", h.attrs)
            for scan in ("1.1", "2.1"):
                g = h[scan]
                for name in ("row", "col", "intensity", "labels",
                             "nnz", "nlabel", "title"):
                    self.assertIn(name, g, "%s missing from %s" % (name, scan))
                for a in ("itype", "nframes", "shape0", "shape1", "npx"):
                    self.assertIn(a, g.attrs)
                self.assertEqual(g.attrs["nframes"], NFRAMES)
                self.assertEqual((g.attrs["shape0"], g.attrs["shape1"]), SHAPE)
                self.assertIn("measurement/diffrz", g)
                self.assertIn("instrument/positioners/diffty", g)

    def test_title_is_kept(self):
        """the scan title says fscan / fscan2d / f2scan and is needed later"""
        with h5py.File(self.sparsefile, "r") as h:
            title = h["1.1"]["title"].asstr()[()]
        self.assertTrue(title.startswith("fscan"), title)

    def test_dtypes(self):
        with h5py.File(self.sparsefile, "r") as h:
            g = h["1.1"]
            self.assertEqual(g["row"].dtype, np.uint16)
            self.assertEqual(g["col"].dtype, np.uint16)
            # follows worker.cor, which is float32 after background subtraction
            self.assertEqual(g["intensity"].dtype, np.float32)
            self.assertEqual(g["labels"].dtype, np.int32)
            # itype describes the detector frames, as in assemble_label. It
            # must not be set to the intensity dtype: from_hdf_group reads it
            # as the row/col index type.
            self.assertEqual(np.dtype(g.attrs["itype"]), np.uint16)

    def test_some_peaks_were_found(self):
        with h5py.File(self.sparsefile, "r") as h:
            for scan in ("1.1", "2.1"):
                self.assertGreater(h[scan]["nlabel"][:].sum(), 0)
                self.assertEqual(h[scan]["nnz"][:].sum(), h[scan].attrs["npx"])

    def test_pixels_are_sorted_within_each_frame(self):
        """the C overlap code is a merge and requires (row, col) order"""
        with h5py.File(self.sparsefile, "r") as h:
            for scan in ("1.1", "2.1"):
                g = h[scan]
                nnz = g["nnz"][:]
                row, col = g["row"][:], g["col"][:]
                ipt = ImageD11.sparseframe.nnz_to_pointer(nnz)
                for i in range(len(nnz)):
                    r = row[ipt[i]:ipt[i + 1]].astype(np.int64)
                    c = col[ipt[i]:ipt[i + 1]].astype(np.int64)
                    key = r * (SHAPE[1] + 1) + c
                    self.assertTrue((np.diff(key) > 0).all(),
                                    "frame %d of %s is not sorted" % (i, scan))

    def test_labels_run_one_to_nlabel_per_frame(self):
        with h5py.File(self.sparsefile, "r") as h:
            for scan in ("1.1", "2.1"):
                g = h[scan]
                nnz, nlab = g["nnz"][:], g["nlabel"][:]
                labels = g["labels"][:]
                ipt = ImageD11.sparseframe.nnz_to_pointer(nnz)
                for i in range(len(nnz)):
                    lab = labels[ipt[i]:ipt[i + 1]]
                    if len(lab) == 0:
                        self.assertEqual(nlab[i], 0)
                        continue
                    self.assertEqual(lab.min(), 1)
                    self.assertEqual(lab.max(), nlab[i])
                    self.assertEqual(len(np.unique(lab)), nlab[i])

    def test_sparsescan_can_read_it(self):
        """the load half: the existing reader must understand what we wrote"""
        scan = ImageD11.sparseframe.SparseScan(self.sparsefile, "1.1")
        self.assertEqual(len(scan.nnz), NFRAMES)
        self.assertEqual(scan.shape, (NFRAMES,) + SHAPE)
        frame = None
        for i in range(NFRAMES):
            frame = scan.getframe(i)
            if frame is not None:
                break
        self.assertIsNotNone(frame, "no non-empty frame to read back")
        self.assertEqual(frame.row.dtype, np.uint16)

    def test_refuses_to_overwrite(self):
        self.assertRaises(
            ValueError,
            frelon_peaksearch.segment_dataset_to_sparse,
            self.ds, WORKER_ARGS, self.sparsefile, 1)


if __name__ == "__main__":
    unittest.main()
