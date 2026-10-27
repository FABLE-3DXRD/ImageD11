"""Synthetic test for ImageD11.sinograms.harvester: Gaussian spots with known intensity and position are
painted on Poisson background frames at the predicted positions of a grain, written with bitshuffle
compression, and harvested back."""
import os
import shutil
import tempfile
import unittest

import numpy as np

try:
    import hdf5plugin      # before h5py: old versions refuse otherwise
except ImportError:
    hdf5plugin = None
import h5py

import ImageD11.grain
from ImageD11.sinograms import harvester

PARS = dict(wavelength=0.2, distance=100e3, y_center=128.0, z_center=128.0, y_size=75.0, z_size=75.0,
            tilt_x=0.0, tilt_y=0.0, tilt_z=0.0, o11=1, o12=0, o21=0, o22=-1, wedge=0.0, chi=0.0,
            t_x=0.0, t_y=0.0, t_z=0.0, omegasign=1.0)
SHAPE = (256, 256)
NFRAMES = 180


class FakeDataset(object):
    """the DataSet attributes the harvester uses"""

    def __init__(self, folder):
        self.datapath = folder
        self.imagefiles = ['scan0001/frames.h5']
        self.limapath = '/data'
        self.scans = ['1.1']
        self.omega = (np.arange(NFRAMES) * 1.0 - 89.5)[None, :]
        self.dty = np.zeros((1, NFRAMES))
        self.y0 = 0.0
        self.ystep = 1.0
        self.shape = self.omega.shape
        self.monitor = None
        self.monitor_ref = None
        self.maskfile = None
        self.e2dxfile = None


def write_frames(ds, spots, intensity, sig_px, sig_om, bg, mask, seed=0):
    rng = np.random.RandomState(seed)
    ss, ff = np.mgrid[0:SHAPE[0], 0:SHAPE[1]]
    frames = np.zeros((NFRAMES,) + SHAPE, float)
    for k in range(len(spots['s'])):
        wf = np.exp(-0.5 * ((ds.omega[0] - spots['omega'][k]) / sig_om) ** 2)
        wf /= wf.sum()
        img = np.exp(-0.5 * (((ss - spots['s'][k]) / sig_px) ** 2 + ((ff - spots['f'][k]) / sig_px) ** 2))
        img /= img.sum()
        for j in np.flatnonzero(wf > 1e-6):
            frames[j] += intensity[k] * wf[j] * img
    data = rng.poisson(frames + bg).astype(np.uint16)
    data[:, mask] = 0
    os.makedirs(os.path.join(ds.datapath, 'scan0001'))
    with h5py.File(os.path.join(ds.datapath, ds.imagefiles[0]), 'w') as h:
        try:
            # bitshuffle (filter 32008): 8192 byte blocks (4096 uint16) + lz4; the same on all hdf5plugin versions
            h.create_dataset('data', data=data, chunks=(1,) + SHAPE, compression=32008, compression_opts=(4096, 2))
            return True
        except (ValueError, RuntimeError):     # no usable bitshuffle plugin (old hdf5plugin / HDF5)
            h.create_dataset('data', data=data, chunks=(1,) + SHAPE)
            return False


class TestHarvester(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.ds = FakeDataset(cls.tmp)
        a = 4.1
        cls.grain = ImageD11.grain.grain(np.eye(3) * a, translation=np.zeros(3))
        cls.mask = np.zeros(SHAPE, bool)
        cls.mask[:, 120:123] = True          # a masked gap through the middle
        # same prediction settings as Harvester.predict
        orange = (cls.ds.omega.min() - 1, cls.ds.omega.max() + 1)
        sp = harvester.predict_spots(cls.grain.ubi, cls.grain.translation, PARS, 0.75, SHAPE, omega_range=orange)
        rng = np.random.RandomState(1)
        cls.truth = dict(sp)
        # only spots well inside the detector and the omega range are painted; the others must come out ~ 0
        cls.painted = ((sp['omega'] > -85) & (sp['omega'] < 85) & (sp['s'] > 20) & (sp['s'] < SHAPE[0] - 20) &
                       (sp['f'] > 20) & (sp['f'] < SHAPE[1] - 20))
        cls.truth['I'] = np.where(cls.painted, rng.uniform(200, 2000, len(sp['s'])), 0.0)
        cls.truth['s'] = sp['s'] + rng.uniform(-0.4, 0.4, len(sp['s']))     # real spots off the prediction
        cls.truth['f'] = sp['f'] + rng.uniform(-0.4, 0.4, len(sp['s']))
        cls.compressed = write_frames(cls.ds, cls.truth, cls.truth['I'], 1.0, 0.7, 0.2, cls.mask)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def harvest(self, reader):
        h = harvester.Harvester(self.ds, [self.grain], PARS, ds_max=0.75, radius=5.0, mask=self.mask, reader=reader,
                                box=harvester.Box(half_s=5, half_f=5, half_omega=4, bg_width=5))
        h.predict()
        return h, h.run(nproc=1, verbose=False)

    def check(self, h, m):
        sp = h.spots
        self.assertGreater(len(sp['s']), 10)
        # spots and truth are the same list (same prediction call)
        self.assertEqual(len(sp['s']), len(self.truth['s']))
        clean = ~sp['touch_mask'] & ~sp['edge'] & (sp['n_overlap'] == 0) & self.painted
        self.assertGreater(clean.sum(), 5)
        empty = ~self.painted & ~sp['edge'] & (m['nframes'] > 0)
        if empty.any():
            self.assertLess(np.abs(m['I'][empty] / m['sigI'][empty]).max(), 5)
        ratio = m['I'][clean] / self.truth['I'][clean]
        z = (m['I'][clean] - self.truth['I'][clean]) / m['sigI'][clean]
        self.assertLess(abs(np.median(ratio) - 1), 0.03)
        self.assertLess(np.std(z), 2.0)
        ds_ = m['s'][clean] - self.truth['s'][clean]
        df_ = m['f'][clean] - self.truth['f'][clean]
        self.assertLess(np.sqrt(np.mean(ds_ ** 2)), 0.1)
        self.assertLess(np.sqrt(np.mean(df_ ** 2)), 0.1)
        # rms widths: typical value 1 px, most within 0.3 (a few noisy spots are allowed)
        self.assertLess(abs(np.median(m['sig_s'][clean]) - 1.0), 0.1)
        self.assertGreater(np.mean(np.abs(m['sig_s'][clean] - 1.0) < 0.3), 0.75)
        # boxes on the masked gap are flagged
        gap = (np.abs(sp['f0'] - 121) <= h.box.half_f + 1) & ~sp['edge']
        self.assertTrue(np.all(sp['touch_mask'][gap]))

    def test_h5py(self):
        if hdf5plugin is None:
            self.skipTest('hdf5plugin is needed to read bitshuffle frames')
        self.check(*self.harvest('h5py'))

    def test_shoebox_model(self):
        if hdf5plugin is None:
            self.skipTest('hdf5plugin is needed to read bitshuffle frames')
        model = harvester.ShoeboxModel(psf=1.0, bandpass=1e-4, phi0=0.3, nsigma=3.5, margin_px=0.5, margin_frames=1.0)
        h = harvester.Harvester(self.ds, [self.grain], PARS, ds_max=0.75, radius=5.0, mask=self.mask, reader='h5py', box=model)
        h.predict()
        sp = h.spots
        self.assertTrue(np.all(sp['hs'] >= 4) and np.all(sp['hs'] <= 6))
        # the rocking width grows as 1/|sin eta|
        se = np.abs(np.sin(np.radians(sp['eta'])))
        self.assertGreater(np.corrcoef(1 / np.maximum(se, 0.05), sp['ho'])[0, 1], 0.5)
        m = h.run(nproc=1, verbose=False)
        clean = ~sp['touch_mask'] & ~sp['edge'] & (sp['n_overlap'] == 0) & self.painted
        self.assertLess(abs(np.median(m['I'][clean] / self.truth['I'][clean]) - 1), 0.03)

    def test_profiles_and_weights(self):
        if hdf5plugin is None:
            self.skipTest('hdf5plugin is needed to read bitshuffle frames')
        box = harvester.Box(half_s=5, half_f=5, half_omega=4, bg_width=5)
        h = harvester.Harvester(self.ds, [self.grain], PARS, ds_max=0.75, radius=5.0, mask=self.mask, reader='h5py',
                                box=box, profiles=True)
        h.predict()
        m = h.run(nproc=1, verbose=False)
        sp = h.spots
        prof, off = h.profile_data
        seen = m['nframes'] > 0
        for k in np.flatnonzero(seen):
            for sl in harvester.profile_slices(sp['hs'][k], sp['hf'][k], sp['ho'][k], off[k]):
                self.assertAlmostEqual(prof[sl].sum(), m['I'][k], delta=1e-6 * max(1.0, abs(m['I'][k])))
        # weights = the true profile of each spot (normalised over its box)
        om = self.ds.omega[0]
        step = om[1] - om[0]
        flat = np.zeros(off[-1])
        for k in range(len(sp['s'])):
            a, b, c = harvester.profile_slices(sp['hs'][k], sp['hf'][k], sp['ho'][k], off[k])
            js = sp['s0'][k] + np.arange(-sp['hs'][k], sp['hs'][k] + 1)
            jf = sp['f0'][k] + np.arange(-sp['hf'][k], sp['hf'][k] + 1)
            t = np.arange(-sp['ho'][k], sp['ho'][k] + 1)
            dom = om[np.clip(np.searchsorted(om, sp['omega'][k] + t * step - step / 2), 0, len(om) - 1)] - sp['omega'][k]
            for sl, w in ((a, np.exp(-0.5 * (js - self.truth['s'][k]) ** 2)),
                          (b, np.exp(-0.5 * (jf - self.truth['f'][k]) ** 2)),
                          (c, np.exp(-0.5 * (dom / 0.7) ** 2))):
                flat[sl] = w / w.sum()
        hw = harvester.Harvester(self.ds, [self.grain], PARS, ds_max=0.75, radius=5.0, mask=self.mask,
                                 reader='h5py', box=box, weights=(flat, off))
        hw.predict()
        out = os.path.join(self.tmp, 'weights.h5')
        hw.run(nproc=1, verbose=False, outfile=out)
        with h5py.File(out, 'r') as f:
            ps = f['per_scan']
            spot = ps['spot'][()]
            wA, wB, wM, wC, wV = (ps[k][()].astype(float) for k in ('wA', 'wB', 'wM', 'wC', 'wVbg'))
        clean = (~sp['touch_mask'] & ~sp['edge'] & (sp['n_overlap'] == 0) & self.painted)[spot] & (wB > 0)
        I = self.truth['I'][spot][clean]
        Iw = wA[clean] / wB[clean]
        sig = np.sqrt((wM[clean] + I * wC[clean] + wV[clean]) / wB[clean] ** 2)
        z = (Iw - I) / sig
        self.assertLess(abs(np.median(Iw / I) - 1), 0.03)
        # the variance model is calibrated, without scale factors
        self.assertGreater(np.std(z), 0.8)
        self.assertLess(np.std(z), 1.25)
        # and the profile-fitted intensities are more precise than the box sums
        Is = m['I'][spot][clean]
        self.assertLess(np.std(Iw - I), np.std(Is - I))

    def test_duplicates_and_learning(self):
        if hdf5plugin is None:
            self.skipTest('hdf5plugin is needed to read bitshuffle frames')
        from ImageD11 import sym_u
        a = np.radians(0.3)
        R = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
        twin = ImageD11.grain.grain(np.dot(self.grain.ubi, R.T), translation=np.zeros(3))   # the same crystal
        box = harvester.Box(half_s=5, half_f=5, half_omega=4, bg_width=5)
        h = harvester.Harvester(self.ds, [self.grain, twin], PARS, ds_max=0.75, radius=5.0, mask=self.mask,
                                reader='h5py', box=box, profiles=True)
        h.predict()
        m = h.run(nproc=1, verbose=False)
        sp = h.spots
        n = np.bincount(sp['grain'])
        M = h.pair_overlaps()
        self.assertGreater(M[0, 1] / float(n[0]), 0.5)
        groups, pairs = harvester.duplicate_groups(sp, m, [g.ubi for g in h.grains], sym_u.getgroup('cubic')().group,
                                                   h.omega_step(), pair_overlaps=M)
        self.assertEqual(len(groups), 1)
        self.assertEqual(sorted(groups[0]), [0, 1])
        merged, match = harvester.union_spots(sp, groups[0], [g.ubi for g in h.grains], h.omega_step())
        k0 = groups[0][0]
        rows0 = np.flatnonzero(sp['grain'] == k0)
        self.assertTrue(np.all(merged['hs'] >= sp['hs'][rows0]))
        mm = match[groups[0][1]]
        ok = mm >= 0
        self.assertGreater(ok.mean(), 0.9)
        # every matched box lies inside the union box
        for c, hc in (('s0', 'hs'), ('f0', 'hf')):
            self.assertTrue(np.all(np.abs(sp[c][mm[ok]] - merged[c][ok]) + sp[hc][mm[ok]] <= merged[hc][ok]))
        d = h.deoverlap(m['I'])
        painted = np.concatenate([self.painted, self.painted])
        self.assertGreater(np.median(d['coincident_fraction'][painted & ~sp['edge']]), 0.5)
        # profiles learned from the strong, clean spots match the true (Gaussian, 1 px) s profile
        prof, off = h.profile_data
        ref = (m['I'] > 20 * m['sigI']) & ~sp['touch_mask'] & ~sp['edge'] & (sp['grain'] == 0)
        learned, nref, pooled = harvester.learn_profiles(sp, m['I'], m['sigI'], prof, off, ref, k=5, min_ref=5)
        corr = []
        for k in np.flatnonzero(ref):
            a_, _, _ = harvester.profile_slices(sp['hs'][k], sp['hf'][k], sp['ho'][k], off[k])
            js = sp['s0'][k] + np.arange(-sp['hs'][k], sp['hs'][k] + 1)
            corr.append(np.corrcoef(learned[a_], np.exp(-0.5 * (js - sp['s'][k]) ** 2))[0, 1])
        self.assertGreater(np.median(corr), 0.9)

    def test_sparse(self):
        if harvester.bslz4_to_sparse is None:
            self.skipTest('bslz4_to_sparse is not installed')
        if not self.compressed:
            self.skipTest('no bitshuffle plugin for h5py: the test file is not compressed')
        self.check(*self.harvest('sparse'))


class TestDeoverlap(unittest.TestCase):

    def test_solve(self):
        import scipy.sparse
        # 0, 1, 2 overlap partially; 3, 4 hold each other's peaks completely; 5 alone
        F = np.zeros((6, 6))
        F[0, 1], F[1, 0], F[2, 0] = 0.3, 0.1, 0.05
        F[3, 4] = F[4, 3] = 1.0
        T = np.array([100.0, 50.0, 10.0, 80.0, 20.0, 7.0])
        I = T + np.dot(F, T)
        d = harvester.deoverlap_solve(scipy.sparse.csr_matrix(F), I, np.sqrt(I))
        self.assertTrue(np.allclose(d['I_deoverlapped'][:3], T[:3], rtol=1e-4))
        self.assertTrue(np.allclose(d['contamination'][:3], np.dot(F, T)[:3], rtol=1e-4))
        # the split of 3 and 4 is not measured: even, and their sum is kept
        self.assertTrue(np.allclose(d['I_deoverlapped'][3:5], 50.0, rtol=1e-4))
        self.assertEqual(d['I_deoverlapped'][5], 7.0)
        self.assertEqual(list(d['largest_row']), [1, 0, 0, 4, 3, -1])
        self.assertEqual(d['cluster'][5], -1)
        self.assertEqual(d['cluster'][0], d['cluster'][2])
        self.assertNotEqual(d['cluster'][0], d['cluster'][3])

    def test_jacobi(self):
        import scipy.sparse
        F = np.zeros((6, 6))
        F[0, 1], F[1, 0], F[2, 0] = 0.3, 0.1, 0.05
        F[3, 4] = F[4, 3] = 1.0
        T = np.array([100.0, 50.0, 10.0, 80.0, 20.0, 7.0])
        I = T + np.dot(F, T)
        d = harvester.deoverlap_jacobi(scipy.sparse.csr_matrix(F), I, niter=40)
        self.assertTrue(np.allclose(d['I_deoverlapped'], [100, 50, 10, 50, 50, 7], rtol=1e-6))
        self.assertTrue(np.allclose(d['contamination'][:3], np.dot(F, T)[:3], rtol=1e-6))
        self.assertEqual(list(d['largest_row']), [1, 0, 0, 4, 3, -1])


if __name__ == '__main__':
    unittest.main()
