import unittest

import os
import tempfile
import h5py
import numpy as np

import sys
if int(sys.version_info.major) == 2:
    raise unittest.SkipTest('Skipping voxel_mask tests on Python 2')
else:
    from ImageD11.sinograms.voxel_mask import (
        VoxelSinoMasker, fill_voxel_idx, max_candidates,
        save_peak_index_cache, load_peak_index_cache, _PEAK_H5GROUP, _CACHE_VERSION)


YSTEP = 2.0
Y0 = 0.0
YBIN_MIN = -274.0          # stands in for ds.ybincens.min()
NDTY = 275
NFRAMES = 400


def brute_force(xi0, yi0, sinomega, cosomega, dty, y0=Y0, ystep=YSTEP):
    """The predicate the partition must reproduce, tested against every peak."""
    d = np.abs(y0 - xi0 * sinomega - yi0 * cosomega - dty)
    idx = np.nonzero(d <= ystep)[0]
    return idx, d[idx]


def make_peaks(npk=200000, first_row=0, seed=0):
    """Synthetic sinogram.

    first_row > 0 leaves the outer dty rows empty, which is what real data
    looks like and what a test must exercise: the partition origin then
    differs from dty.min().
    """
    rng = np.random.default_rng(seed)
    iframe = rng.integers(0, NFRAMES, npk)
    frame_omega = 0.0626 + np.arange(NFRAMES) * 0.45
    omega = frame_omega[iframe]
    rows = rng.integers(first_row, NDTY - first_row, npk)
    dty = YBIN_MIN + rows * YSTEP
    return omega, dty, iframe, frame_omega


class TestPartitionSelection(unittest.TestCase):
    """The partition is a pre-filter; the selection must be exactly the
    brute-force predicate, at every omega binning."""

    def _check(self, masker, nvoxel=20, seed=1):
        rng = np.random.default_rng(seed)
        sinomega = np.sin(np.radians(masker.omega))
        cosomega = np.cos(np.radians(masker.omega))
        for _ in range(nvoxel):
            xi0 = rng.uniform(-280, 280)
            yi0 = rng.uniform(-280, 280)
            idx, ydist = masker.mask(xi0, yi0, YSTEP, Y0)
            eidx, eydist = brute_force(xi0, yi0, sinomega, cosomega,
                                       masker.dty)
            self.assertTrue(np.array_equal(np.sort(idx), np.sort(eidx)))
            self.assertTrue(np.allclose(np.sort(ydist), np.sort(eydist)))

    def test_heuristic_binsize(self):
        omega, dty, _, _ = make_peaks()
        m = VoxelSinoMasker(omega.copy(), dty.copy(), YSTEP, ymin=YBIN_MIN)
        m.partition()
        self._check(m)

    def test_supplied_frame_bins(self):
        omega, dty, iframe, frame_omega = make_peaks()
        m = VoxelSinoMasker(omega.copy(), dty.copy(), YSTEP, ymin=YBIN_MIN)
        m.partition(omega_bin_indices=iframe, omega_bins=frame_omega)
        self._check(m)

    def test_coarse_binsize(self):
        omega, dty, _, _ = make_peaks()
        m = VoxelSinoMasker(omega.copy(), dty.copy(), YSTEP, ymin=YBIN_MIN)
        m.partition(omega_binsize=4.0)
        self._check(m)


class TestOmegaBinCentres(unittest.TestCase):
    """Supplied bin indices must keep their correspondence with the supplied
    centres. Renumbering the indices without slicing the centres to match
    shifts every bin, which silently drops peaks."""

    def _max_dev(self, m):
        ob = np.asarray(m.omega_bins)
        worst = 0.0
        for i in range(len(ob)):
            lo, hi = m.omega_partitions[i], m.omega_partitions[i + 1]
            if hi > lo:
                worst = max(worst, np.abs(m.omega[lo:hi] - ob[i]).max())
        return worst

    def test_first_frames_populated(self):
        omega, dty, iframe, frame_omega = make_peaks()
        m = VoxelSinoMasker(omega.copy(), dty.copy(), YSTEP, ymin=YBIN_MIN)
        m.partition(omega_bin_indices=iframe, omega_bins=frame_omega)
        self.assertLess(self._max_dev(m), 1e-9)

    def test_leading_frames_empty(self):
        # no peaks in the first 7 frames -- min(bin index) is 7, not 0
        omega, dty, iframe, frame_omega = make_peaks()
        keep = iframe >= 7
        omega, dty, iframe = omega[keep], dty[keep], iframe[keep]
        m = VoxelSinoMasker(omega.copy(), dty.copy(), YSTEP, ymin=YBIN_MIN)
        m.partition(omega_bin_indices=iframe, omega_bins=frame_omega)
        self.assertLess(self._max_dev(m), 1e-9)


class TestYminOrigin(unittest.TestCase):
    """The dty bin grid must start where the caller says, not at the lowest
    observed dty. A caller comparing bin indices against a dtyi column built
    from ds.ybincens.min() is otherwise offset by the number of empty outer
    rows."""

    def test_origin_is_respected(self):
        omega, dty, _, _ = make_peaks(first_row=4)
        self.assertGreater(dty.min(), YBIN_MIN)
        m = VoxelSinoMasker(omega.copy(), dty.copy(), YSTEP, ymin=YBIN_MIN)
        m.partition()
        self.assertAlmostEqual(m.ymin, YBIN_MIN)
        expect = np.round((m.dty - YBIN_MIN) / YSTEP).astype(int)
        got = np.round((m.dty - m.ymin) / YSTEP).astype(int)
        self.assertTrue(np.array_equal(expect, got))

    def test_default_origin_is_dty_min(self):
        omega, dty, _, _ = make_peaks(first_row=4)
        m = VoxelSinoMasker(omega.copy(), dty.copy(), YSTEP)
        self.assertAlmostEqual(m.ymin, dty.min())


class TestCallerArraysUntouched(unittest.TestCase):
    """partition() must not reorder the arrays it was given: they are usually
    columns of a columnfile, and every other column would be misaligned."""

    def test_not_modified(self):
        omega, dty, _, _ = make_peaks()
        om_ref, dty_ref = omega.copy(), dty.copy()
        m = VoxelSinoMasker(omega, dty, YSTEP, ymin=YBIN_MIN)
        m.partition()
        self.assertTrue(np.array_equal(omega, om_ref))
        self.assertTrue(np.array_equal(dty, dty_ref))

    def test_inplace_does_modify(self):
        omega, dty, _, _ = make_peaks()
        om_ref = omega.copy()
        m = VoxelSinoMasker(omega, dty, YSTEP, ymin=YBIN_MIN)
        m.partition(inplace=True)
        self.assertFalse(np.array_equal(omega, om_ref))


class TestBufferOverflow(unittest.TestCase):
    """fill_voxel_idx returns -1 rather than raising: it is called from
    inside a prange, where an exception is not reliably propagated."""

    def test_returns_sentinel(self):
        omega, dty, _, _ = make_peaks()
        m = VoxelSinoMasker(omega.copy(), dty.copy(), YSTEP, ymin=YBIN_MIN)
        m.partition()
        n = fill_voxel_idx(
            0.0, 0.0, Y0, YSTEP, m.ymin,
            m.omega_partitions, m.dty_partitions, m.dty,
            m.sinomega, m.cosomega, m.sinomega_bins, m.cosomega_bins,
            np.empty(4, np.int64), np.empty(4, np.float64))
        self.assertEqual(n, -1)


class TestMaxCandidates(unittest.TestCase):
    """max_candidates must bound what fill_voxel_idx actually returns, or the
    buffers sized from it can overflow."""

    def test_is_an_upper_bound(self):
        omega, dty, _, _ = make_peaks()
        m = VoxelSinoMasker(omega.copy(), dty.copy(), YSTEP, ymin=YBIN_MIN)
        m.partition()
        rng = np.random.default_rng(2)
        xs = rng.uniform(-280, 280, 40)
        ys = rng.uniform(-280, 280, 40)
        bound = int(max_candidates(
            np.ascontiguousarray(xs), np.ascontiguousarray(ys),
            Y0, YSTEP, m.ymin, m.dty_partitions,
            m.sinomega_bins, m.cosomega_bins))
        worst = max(len(m.mask(xs[k], ys[k], YSTEP, Y0)[0])
                    for k in range(len(xs)))
        self.assertGreaterEqual(bound, worst)


class TestIndexCache(unittest.TestCase):
    """The cache must round-trip, and must refuse a file it cannot trust."""
    def _idx(self):
        rng = np.random.default_rng(0)
        return dict(order=rng.permutation(500),
                    omega_partitions=np.arange(0, 501, 50),
                    dty_partitions=rng.integers(0, 50, (10, 12)),
                    usin=np.sin(np.linspace(0, 6, 10)),
                    ucos=np.cos(np.linspace(0, 6, 10)),
                    nom=10, ndty=11, ymin=YBIN_MIN,
                    ray_margin=2.0, dev=0.0, maxlocal=1234)
 
    def setUp(self):
        self.idx = self._idx()
        self.fn = os.path.join(tempfile.mkdtemp(), "idx.h5")
        save_peak_index_cache(self.idx, self.fn)
 
    def test_round_trip(self):
        got = load_peak_index_cache(self.fn)
        self.assertIsNotNone(got)
        for k in ("order", "omega_partitions", "dty_partitions", "usin", "ucos"):
            self.assertTrue(np.array_equal(np.asarray(got[k]), self.idx[k]), k)
        self.assertEqual(got["maxlocal"], 1234)
        self.assertAlmostEqual(got["ymin"], YBIN_MIN)
        self.assertEqual(got["nom"], 10)
 
    def test_mmap_matches(self):
        mm = load_peak_index_cache(self.fn, mmap=True)
        self.assertIsNotNone(mm)
        for k in ("order", "omega_partitions", "dty_partitions", "usin", "ucos"):
            self.assertTrue(np.array_equal(np.asarray(mm[k]), self.idx[k]), k)
 
    def test_missing_file(self):
        self.assertIsNone(load_peak_index_cache(self.fn + ".nope"))
 
    def test_wrong_version(self):
        # a file from an older ImageD11 must be refused, not half-read
        with h5py.File(self.fn, "a") as h:
            h[_PEAK_H5GROUP].attrs["version"] = _CACHE_VERSION - 1
        self.assertIsNone(load_peak_index_cache(self.fn))
 
    def test_no_group(self):
        with h5py.File(self.fn, "a") as h:
            del h[_PEAK_H5GROUP]
        self.assertIsNone(load_peak_index_cache(self.fn))
 
    def test_maxlocal_absent_is_none(self):
        # initializer must be able to tell "no maxlocal" from a real value:
        # np.empty(None) raises an unhelpful TypeError
        with h5py.File(self.fn, "a") as h:
            del h[_PEAK_H5GROUP].attrs["maxlocal"]
        self.assertIsNone(load_peak_index_cache(self.fn)["maxlocal"])


if __name__ == "__main__":
    unittest.main()