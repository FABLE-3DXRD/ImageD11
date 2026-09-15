"""
Tests for the unified sparse labelling / merging:
  frelon_peaksearch.segment_dataset_to_sparse  -> sparse pixels + labels
  sinograms.properties.main(algorithm='stored', max_centroid_dist=...)
"""
from __future__ import print_function
import os
import shutil
import tempfile
import unittest
import warnings

import h5py
import numpy as np

import ImageD11.sinograms.dataset
import ImageD11.sinograms.properties as properties
import ImageD11.frelon_peaksearch as frelon_peaksearch

NFRAMES = 12
SHAPE = (96, 96)
OSTEP = 1.0
DTYS = (-0.5, 0.0)
# row, col, first frame, last frame (inclusive)
SPOTS = [
    (30.0, 30.0, 2, 6),  # isolated
    (70.0, 70.0, 1, 3),  # isolated
    (60.0, 20.0, 3, 5),  # C : fades out on frame 5 ...
    (60.0, 26.0, 6, 8),  # D : ... appears on frame 6, pixels touch C, centroid 6px away
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
        # sub-pixel wobble between frames and scans so centroids differ a bit
        for f in range(f0, f1 + 1):
            dr = 0.1 * np.sin(f + k + iscan)
            amp = 2000.0 * (1 + 0.2 * np.cos(f + iscan))
            frames[f] += amp * np.exp(-0.5 * ((r - r0 - dr) ** 2 + (c - c0) ** 2) / 1.5**2)
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
            frame = np.arange(NFRAMES)
            m["fpico4"] = (1e5 * (1 + 0.1 * np.sin(frame + i))).astype(np.float32)  # bliss counters are float32
            g["instrument/positioners/diffty"] = DTYS[i]
            layout = h5py.VirtualLayout(shape=(NFRAMES,) + SHAPE, dtype=np.uint16)
            vsrc = h5py.VirtualSource(
                limaname, "/entry_0000/measurement/data", shape=(NFRAMES,) + SHAPE
            )
            layout[0:NFRAMES] = vsrc[0:NFRAMES]  # explicit selection, like bliss
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


def match(old, new, keys=("s_raw", "f_raw", "omega", "sum_intensity", "Number_of_pixels")):
    """sort both peak dicts by (s_raw, f_raw, omega) and compare"""
    o = np.lexsort((np.round(old["omega"], 3), np.round(old["f_raw"], 1), np.round(old["s_raw"], 1)))
    n = np.lexsort((np.round(new["omega"], 3), np.round(new["f_raw"], 1), np.round(new["s_raw"], 1)))
    return {k: (np.asarray(old[k])[o], np.asarray(new[k])[n]) for k in keys}


class TestCentroidGate(unittest.TestCase):
    def test_gate(self):
        # three peaks, sI = 10, centroids at (0,0), (1,1), (5,0)
        s1 = [3, 3, 3]
        sI = [10, 10, 10]
        srI = [0, 10, 50]
        scI = [0, 10, 0]
        frm = [0, 1, 1]
        pk_props = np.array([s1, sI, srI, scI, frm], np.int64)
        rc = np.array([[0, 0], [1, 2], [4, 1]], np.int64)
        keep = properties.centroid_gate(rc, pk_props, 1.6)
        self.assertEqual(list(keep), [True, False])
        keep = properties.centroid_gate(rc, pk_props, 5.0)
        self.assertEqual(list(keep), [True, True])

    def test_find_uniq_keep(self):
        npk = np.array([[3, 2, 0]])
        pkst = properties.pks_table(npk=npk)
        pkst.pk_props[:] = np.array(
            [[3, 3, 3], [10, 10, 10], [0, 10, 50], [0, 10, 0], [0, 1, 1]]
        )
        pkst.rc[:] = np.array([[0, 0], [1, 2], [4, 1]])
        n, labels = pkst.find_uniq()
        self.assertEqual(n, 1)
        keep = properties.merge_mask(pkst, 1.6, verbose=0)
        n, labels = pkst.find_uniq(keep=keep)
        self.assertEqual(n, 2)
        self.assertEqual(labels[0], labels[1])
        self.assertNotEqual(labels[0], labels[2])
        self.assertIsNone(properties.merge_mask(pkst, None))
        del pkst

    def test_bad_options(self):
        with self.assertRaises(Exception):
            properties.main("nofile", options={"max_centroid_dist_typo": 1})
        with self.assertRaises(ValueError):
            properties.main("nofile", options={"algorithm": "nope"})
        with self.assertRaises(ValueError):
            properties.main("nofile", options={"max_centroid_dist": -1})


class TestSparseFrelon(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def run_sparse(self, ds, gate, monitor=None, name=None):
        """segment once, merge into a separate pksfile per gate"""
        if not os.path.exists(ds.sparsefile):
            frelon_peaksearch.segment_dataset_to_sparse(
                ds, WORKER_ARGS, num_cpus=2, monitor_name=monitor
            )
        ds.pksfile = os.path.join(ds.analysispath, "pks_%s.h5" % name)
        ds.save()
        properties.main(
            ds.dsfile, options={"algorithm": "stored", "max_centroid_dist": gate}
        )
        ds2 = ImageD11.sinograms.dataset.load(ds.dsfile)
        return ds2

    def test_4d(self):
        ds = make_dataset(self.tmp, "fourd", 2)
        self.assertEqual(ds.shape, (2, NFRAMES))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cf2_old, cf_old = frelon_peaksearch.segment_dataset(
                ds, WORKER_ARGS, num_cpus=2, scan_number=range(2),
                monitor_name="fpico4", disable=True,
            )
        # any overlap: C and D chain together
        ds_any = self.run_sparse(ds, None, monitor="fpico4", name="any")
        self.assertEqual(ds_any.peaks_table.nlabel, len(SPOTS) - 1)
        # the sparse file has the pixels and labels
        with h5py.File(ds.sparsefile, "r") as h:
            for scan in ("1.1", "2.1"):
                g = h[scan]
                for name in ("row", "col", "intensity", "labels", "nnz", "nlabel"):
                    self.assertIn(name, g)
                self.assertEqual(g["nnz"][:].sum(), len(g["row"]))
                self.assertIn("fpico4", g["measurement"])
                labels = g["labels"][:]
                nlabel = g["nlabel"][:]
                ip = np.concatenate(([0], np.cumsum(g["nnz"][:]))).astype(np.int64)
                for f in range(NFRAMES):
                    lf = labels[ip[f]: ip[f + 1]]
                    if nlabel[f]:
                        self.assertEqual(sorted(set(lf)), list(range(1, nlabel[f] + 1)))

        with h5py.File(ds_any.pksfile, "r") as h:
            self.assertEqual(h["pks2d"].attrs["algorithm"], "stored")
            self.assertEqual(h["pks2d"].attrs["max_centroid_dist"], "None")

        # gated: matches the old frelon 4D merge
        ds_gate = self.run_sparse(ds, 1.6, name="gate")
        with h5py.File(ds_gate.pksfile, "r") as h:
            self.assertEqual(h["pks2d"].attrs["max_centroid_dist"], 1.6)
        self.assertIsNotNone(ds_gate.monitor)
        pk4d = ds_gate.pk4d
        self.assertEqual(len(pk4d["s_raw"]), len(SPOTS))
        self.assertEqual(cf_old.nrows, len(SPOTS))
        old = {k: cf_old.getcolumn(k) for k in cf_old.titles}
        for k, (a, b) in match(old, pk4d).items():
            if k == "sum_intensity":
                # new: sums truncated per peak. old: float sums
                np.testing.assert_allclose(a, b, rtol=1e-4, err_msg=k)
            elif k == "Number_of_pixels":
                np.testing.assert_array_equal(a, b, err_msg=k)
            else:
                np.testing.assert_allclose(a, b, atol=0.01, err_msg=k)
        # both dty rows were merged into each spot
        np.testing.assert_allclose(sorted(pk4d["dty"]), [np.mean(DTYS)] * len(SPOTS), atol=0.1)
        # the 2D peaks agree too
        pk2d = ds_gate.pk2d
        self.assertEqual(len(pk2d["s_raw"]), cf2_old.nrows)
        np.testing.assert_allclose(
            np.sort(pk2d["sum_intensity"]), np.sort(cf2_old.sum_intensity), rtol=1e-3
        )

    def test_3d_single_scan(self):
        ds = make_dataset(self.tmp, "threed", 1)
        self.assertEqual(len(ds.scans), 1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cf2_old, cf3_old = frelon_peaksearch.segment_dataset(
                ds, WORKER_ARGS, num_cpus=2, scan_number=0, disable=True
            )
        ds_gate = self.run_sparse(ds, 1.6, name="gate3d")
        pk3d = ds_gate.pk4d
        self.assertEqual(len(pk3d["s_raw"]), len(SPOTS))
        old = {k: cf3_old.getcolumn(k) for k in cf3_old.titles}
        for k, (a, b) in match(old, pk3d).items():
            if k in ("sum_intensity",):
                np.testing.assert_allclose(a, b, rtol=1e-4, err_msg=k)
            else:
                np.testing.assert_allclose(a, b, atol=0.01, err_msg=k)
        pk2d = ds_gate.pk2d
        self.assertEqual(len(pk2d["s_raw"]), cf2_old.nrows)
        ds_any = self.run_sparse(ds, None, name="any3d")
        self.assertEqual(ds_any.peaks_table.nlabel, len(SPOTS) - 1)

    def test_refuse_overwrite(self):
        ds = make_dataset(self.tmp, "again", 1)
        frelon_peaksearch.segment_dataset_to_sparse(ds, WORKER_ARGS, num_cpus=1)
        with self.assertRaises(ValueError):
            frelon_peaksearch.segment_dataset_to_sparse(ds, WORKER_ARGS, num_cpus=1)


if __name__ == "__main__":
    unittest.main()
