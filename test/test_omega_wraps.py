"""
DataSet.omega_wraps: when the scan turns more than once, omega has been folded
into 0-360, so a 4D peak can have 2D peaks at 359 and at 1 degrees. The merge
then has to average omega on the circle rather than linearly. The flag defaults
to False and is set True for f2scan or by hand.

The Welford running mean is what does the averaging: the update needs the
difference o - mean, so that difference is taken the short way round the circle,
and the running mean itself is the reference the unwrapping is measured against.
"""
from __future__ import print_function

import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

IMPORT_ERROR = None
try:
    from ImageD11.sinograms.dataset import DataSet
    from ImageD11.sinograms.properties import pks_table, numbapkmerge
except Exception as e:  # pragma: no cover
    IMPORT_ERROR = e
SKIP = (sys.version_info[0] < 3) or (IMPORT_ERROR is not None)
REASON = "needs python3 (import error: %s)" % (IMPORT_ERROR,)


def merge_one(hand_built, omega, dty, weights, scale_factor=None,
              omega_wraps=False):
    """One 4D peak made of len(weights) 2D peaks. pk_props has float frm on
    purpose, as it can be on the float route."""
    n = len(weights)
    pk_props = np.zeros((5, n))
    pk_props[0] = 1                     # s1, pixels per 2D peak
    pk_props[1] = weights               # sI, the weight
    pk_props[2] = weights * 3.0         # srI
    pk_props[3] = weights * 4.0         # scI
    pk_props[4] = np.arange(n, dtype=float)  # frm, float on purpose
    glabel = np.zeros(n, np.int64)      # all one 4D peak
    pkst = pks_table(pk_props=pk_props, glabel=glabel, nlabel=1)
    return pkst, np.asarray(omega, float), np.asarray(dty, float)


@unittest.skipIf(SKIP, REASON)
class TestCircularMerge(unittest.TestCase):
    """numbapkmerge / pk2dmerge averaging omega the short way round"""

    def test_default_is_linear(self):
        """omega_wraps=False must not perturb the existing arithmetic"""
        o = np.array([100.0, 101.0, 102.0, 103.0])
        w = np.array([1.0, 5.0, 4.0, 2.0])
        pkst, omega, dty = merge_one(None, o, np.zeros_like(o), w)
        got = pkst.pk2dmerge(omega, dty)["omega"][0]
        self.assertAlmostEqual(got, (o * w).sum() / w.sum(), places=12)

    def test_circular_mean_crosses_the_seam(self):
        """the same peak seen at 359 and 1 degrees averages to 0, not 180 --
        which is exactly what the linear mean gets wrong"""
        o = np.array([358.0, 359.0, 0.0, 1.0])
        w = np.array([1.0, 5.0, 4.0, 2.0])
        pkst, omega, dty = merge_one(None, o, np.zeros_like(o), w)
        unwrapped = np.array([-2.0, -1.0, 0.0, 1.0])
        want = (unwrapped * w).sum() / w.sum() % 360.0
        got = pkst.pk2dmerge(omega, dty, omega_wraps=True)["omega"][0]
        self.assertAlmostEqual(got, want, places=10)
        linear = pkst.pk2dmerge(omega, dty)["omega"][0]
        self.assertGreater(abs(((linear - want + 180.0) % 360.0) - 180.0), 90.0)

    def test_result_is_in_range(self):
        o = np.array([359.5, 0.5])
        w = np.array([1.0, 1.0])
        pkst, omega, dty = merge_one(None, o, np.zeros_like(o), w)
        got = pkst.pk2dmerge(omega, dty, omega_wraps=True)["omega"][0]
        self.assertTrue(0.0 <= got < 360.0, got)

    def test_circular_mean_is_order_independent(self):
        rng = np.random.RandomState(7)
        o = (rng.random_sample(30) * 40.0 - 20.0) % 360.0
        w = rng.random_sample(30) * 1000.0 + 1.0
        vals = []
        for i in range(5):
            perm = rng.permutation(len(o))
            pkst, omega, dty = merge_one(None, o[perm], np.zeros_like(o[perm]),
                                         w[perm])
            vals.append(pkst.pk2dmerge(omega, dty, omega_wraps=True)["omega"][0])
        self.assertLess(max(vals) - min(vals), 1e-9)

    def test_dty_is_not_treated_as_an_angle(self):
        """dty is a translation: it keeps the plain weighted mean"""
        o = np.array([358.0, 2.0])
        dty = np.array([10.0, -10.0])
        w = np.array([1.0, 1.0])
        pkst, omega, dty = merge_one(None, o, dty, w)
        got = pkst.pk2dmerge(omega, dty, omega_wraps=True)["dty"][0]
        self.assertAlmostEqual(got, (dty * w).sum() / w.sum(), places=10)

    def test_scale_factor_reaches_the_centroid(self):
        """scale_factor normalises against the incident beam monitor, so
        the centroid is weighted by the normalised intensity"""
        o = np.array([10.0, 20.0])
        w = np.array([1.0, 1.0])
        sf = np.array([1.0, 3.0])
        pkst, omega, dty = merge_one(None, o, np.zeros_like(o), w)
        want = (o * w * sf).sum() / (w * sf).sum()
        got = pkst.pk2dmerge(omega, dty, scale_factor=sf,
                             omega_wraps=True)["omega"][0]
        self.assertAlmostEqual(got, want, places=10)

    def test_numbapkmerge_float_frm_with_scale_factor(self):
        """out[4] is a mean, but frm is float: Welford must still work"""
        n = 4
        pk_props = np.zeros((5, n))
        pk_props[1] = np.array([1.0, 2.0, 3.0, 4.0])
        pk_props[4] = np.arange(n, dtype=float)
        out = np.zeros((7, 1), float)
        o = np.array([0.0, 5.0, 355.0, 10.0])
        dty = np.zeros(n)
        sf = np.ones(n)
        labels = np.zeros(n, np.int64)
        numbapkmerge(labels, pk_props, o, dty, out,
                     scale_factor=sf, omega_wraps=True)
        self.assertEqual(out[4, 0], out[4, 0])  # no exception, finite result
        self.assertTrue(np.isfinite(out[4, 0]))


@unittest.skipIf(SKIP, REASON)
class TestDataSetOmegaWraps(unittest.TestCase):
    """the flag on DataSet: default False, settable, saved and loaded"""

    def make_ds(self, nturn, nomega=36, ny=5):
        ds = DataSet(sample="S", dset="d")
        ds.shape = (ny, nomega * nturn)
        o = np.arange(nomega * nturn) * (360.0 / nomega)
        ds.omega = np.tile(o, (ny, 1))
        ds.dty = np.repeat(np.linspace(-1, 1, ny), nomega * nturn).reshape(ds.shape)
        return ds

    def test_defaults_to_false(self):
        ds = self.make_ds(3)
        ds.guessbins()
        self.assertFalse(ds.omega_wraps)

    def test_settable_via_init(self):
        ds = DataSet(sample="S", dset="d", omega_wraps=True)
        self.assertTrue(ds.omega_wraps)

    def test_wraps_false_does_not_fold_even_if_span_exceeds_360(self):
        """an fscan2d that overruns 360 (say -1 to 361) is not periodic: its
        ends are used for sample alignment, not merged, so guessbins must not
        fold it just because the raw span exceeds 360"""
        ds = DataSet(sample="S", dset="d")
        ds.shape = (1, 1440)
        ds.omega = np.linspace(-1, 361, 1440).reshape(1, 1440)
        ds.dty = np.zeros((1, 1440))
        ds.guessbins()
        self.assertFalse(ds.omega_wraps)
        self.assertTrue((ds.omega_for_bins == ds.omega).all())
        self.assertEqual(len(ds.obincens), 1440)

    def test_wraps_true_folds_and_bins_match_the_shape(self):
        """when omega does wrap, one bin per frame closes the circle: there is
        no extra 0/360 bin (the off-by-one) and the folded omega stays in range"""
        ds = self.make_ds(1, nomega=1440)
        ds.omega = (np.arange(1440) * (360.0 / 1440)).reshape(1, 1440)
        ds.dty = np.zeros((1, 1440))
        ds.omega_wraps = True
        ds.guessbins()
        self.assertEqual(len(ds.obincens), 1440)
        self.assertEqual(len(ds.obinedges), 1441)
        self.assertAlmostEqual(len(ds.obincens) * ds.ostep, 360.0, places=9)
        # all frames land in a bin
        io = np.digitize(ds.omega_for_bins, ds.obinedges) - 1
        self.assertTrue(((io >= 0) & (io < len(ds.obincens))).all())

    def test_pk4d_passes_omega_wraps(self):
        """the flag reaches pk2dmerge"""
        ds = self.make_ds(1)
        ds.monitor = None
        ds.omega_for_bins = np.zeros(3)
        ds.dty = np.zeros(3)
        ds.obinedges = np.array([-0.125, 0.125])
        ds.omega_wraps = True
        calls = {}

        class Fake(object):
            def pk2dmerge(self, omega, dty, scale_factor=None, omega_wraps=False,
                          omega0=0.0):
                calls["omega_wraps"] = omega_wraps
                calls["omega0"] = omega0
                return {}

        ds._peaks_table = Fake()
        ds.pk4d
        self.assertTrue(calls["omega_wraps"])
        self.assertEqual(calls["omega0"], ds.obinedges[0])

    def test_flag_survives_save_and_load(self):
        tmp = tempfile.mkdtemp(prefix="id11_omega_wraps_")
        try:
            for flag in (False, True):
                ds = self.make_ds(1 if flag else 3)
                ds.omega_wraps = flag
                h5name = os.path.join(tmp, "ds_%s.h5" % flag)
                ds.save(h5name)
                self.assertEqual(DataSet(filename=h5name).omega_wraps, flag)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_old_file_without_flag_loads_as_false(self):
        """a file saved before omega_wraps existed has no attribute: it must
        load, and default to False"""
        import h5py
        tmp = tempfile.mkdtemp(prefix="id11_omega_wraps_old_")
        try:
            ds = self.make_ds(1)
            h5name = os.path.join(tmp, "old.h5")
            ds.save(h5name)
            with h5py.File(h5name, "a") as h:
                del h["/"].attrs["omega_wraps"]
            self.assertFalse(DataSet(filename=h5name).omega_wraps)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
