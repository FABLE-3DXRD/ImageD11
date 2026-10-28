"""End-to-end run of ImageD11.sinograms.harvest_pipeline (steps duplicates .. grains) on the synthetic frames of
test_harvester: one grain, and a copy of it rotated by 0.3 degrees that must be merged with it."""
import os
import sys
import json
import shutil
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_harvester import FakeDataset, write_frames, hdf5plugin, PARS, SHAPE     # noqa: E402

import h5py                                     # noqa: E402
import ImageD11.grain                           # noqa: E402
from ImageD11.sinograms import harvester, harvest_pipeline as hp      # noqa: E402


class TestHarvestPipeline(unittest.TestCase):

    def setUp(self):
        if hdf5plugin is None:
            self.skipTest('hdf5plugin is needed to read bitshuffle frames')
        self.tmp = tempfile.mkdtemp()
        ds = FakeDataset(self.tmp)
        ds.ybincens, ds.ymin = np.array([0.0]), 0.0
        g0 = ImageD11.grain.grain(np.eye(3) * 4.1, translation=np.zeros(3))
        a = np.radians(0.3)
        R = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
        g1 = ImageD11.grain.grain(np.dot(g0.ubi, R.T), translation=np.zeros(3))
        orange = (ds.omega.min() - 1, ds.omega.max() + 1)
        sp = harvester.predict_spots(g0.ubi, g0.translation, PARS, 0.75, SHAPE, omega_range=orange)
        painted = ((sp['omega'] > -85) & (sp['omega'] < 85) & (sp['s'] > 20) & (sp['s'] < SHAPE[0] - 20) &
                   (sp['f'] > 20) & (sp['f'] < SHAPE[1] - 20))
        self.truth = dict(sp, I=np.where(painted, np.random.RandomState(1).uniform(200, 2000, len(sp['s'])), 0.0))
        write_frames(ds, self.truth, self.truth['I'], 1.0, 0.7, 0.2, np.zeros(SHAPE, bool))
        self.ds = ds
        self._dp = hp.dataset_and_pars
        hp.dataset_and_pars = lambda cfg: (ds, PARS)
        first = os.path.join(self.tmp, 'harvest.h5')
        h = harvester.Harvester(ds, [g0, g1], PARS, ds_max=0.75, radius=5.0, reader='h5py',
                                box=harvester.Box(half_s=5, half_f=5, half_omega=4, bg_width=5))
        h.predict()
        h.run(nproc=1, outfile=first, verbose=False)
        with h5py.File(first, 'a') as f:
            f['gid'], f['radius'] = np.array([0, 1]), np.array([5.0, 5.0])
            f['ubi'], f['translation'] = np.array([g0.ubi, g1.ubi]), np.zeros((2, 3))
        cfgfile = os.path.join(self.tmp, 'config.json')
        with open(cfgfile, 'w') as f:
            json.dump({"dataset": "fake", "phase": "fake", "outdir": os.path.join(self.tmp, 'out'), "nproc": 1,
                       "harvest": {"file": first, "reader": "h5py"},
                       "duplicates": {"min_snr": 5.0},
                       "flags": {"min_clean": 5},
                       "learn": {"k": 5, "min_ref": 5},
                       "mlem": {"snr_min": 5.0}}, f)
        self.cfg = hp.load_config(cfgfile)

    def tearDown(self):
        if hasattr(self, '_dp'):
            hp.dataset_and_pars = self._dp
        shutil.rmtree(self.tmp)

    def test_run(self):
        hp.run(self.cfg, hp.STEPS[1:])
        with open(hp.out(self.cfg, 'duplicates.json')) as f:
            self.assertEqual(json.load(f)['groups'], [[0, 1]])
        with h5py.File(hp.out(self.cfg, 'pair_overlaps.h5'), 'r') as f:
            self.assertEqual(f['pair_overlaps'].shape, (2, 2))
            self.assertGreater(f['pair_overlaps'][0, 1], 0.5 * f['nrefl'][0])
        with h5py.File(hp.out(self.cfg, 'grains.h5'), 'r') as f:
            self.assertEqual(list(f.keys()), ['0'])
            g = f['0']
            self.assertIn('provenance_config', f.attrs)
            ok = np.isfinite(g['I_final'][()]) & ~g['refl_edge'][()] & (g['I_sum'][()] > 20 * g['sigI_sum'][()])
            self.assertGreater(ok.sum(), 10)
            # the union boxes hold the painted peaks: their intensities come back
            t = dict(((h, k, l), I) for h, k, l, I in zip(self.truth['h'], self.truth['k'], self.truth['l'],
                                                           self.truth['I']))
            H = g['refl_hkl'][()][ok]
            ratio = g['I_final'][()][ok] / np.array([t.get(tuple(x), np.nan) for x in H])
            self.assertLess(abs(np.nanmedian(ratio) - 1), 0.05)
            I0 = g['I_final'][()]
            self.assertTrue(np.all(g['transmission'][()] == 1))
        with h5py.File(hp.out(self.cfg, 'sample.h5'), 'r') as f:
            self.assertGreater(f['mask'][()].sum(), 0)
        # with absorption: the same intensities divided by the traced transmission
        try:
            import xraylib       # noqa: F401
        except ImportError:
            return
        self.cfg['absorption'].update(formula='Si', density=2.33, dty_unit='mm')
        hp.run(self.cfg, ['grains'])
        with h5py.File(hp.out(self.cfg, 'grains.h5'), 'r') as f:
            T = f['0/transmission'][()]
            self.assertTrue(np.all((T > 0) & (T <= 1)) and np.any(T < 1))
            self.assertTrue(np.allclose(f['0/I_final'][()] * T, I0, equal_nan=True))



class TestHarvestStep(unittest.TestCase):
    """the harvest step on the synthetic frames: the grains file and its labels replaced by stand-ins; the
    ShoeboxModel instrumental terms are fitted to a pilot harvest (the painted spots are 1 px wide)"""

    def setUp(self):
        if hdf5plugin is None:
            self.skipTest('hdf5plugin is needed to read bitshuffle frames')
        from ImageD11.sinograms import sinogram
        self.tmp = tempfile.mkdtemp()
        ds = FakeDataset(self.tmp)
        ds.ybincens, ds.ymin, ds.grainsfile = np.array([0.0]), 0.0, 'none'
        g0 = ImageD11.grain.grain(np.eye(3) * 4.1, translation=np.zeros(3))
        orange = (ds.omega.min() - 1, ds.omega.max() + 1)
        sp = harvester.predict_spots(g0.ubi, g0.translation, PARS, 0.75, SHAPE, omega_range=orange)
        painted = ((sp['omega'] > -85) & (sp['omega'] < 85) & (sp['s'] > 20) & (sp['s'] < SHAPE[0] - 20) &
                   (sp['f'] > 20) & (sp['f'] < SHAPE[1] - 20))
        I = np.where(painted, np.random.RandomState(1).uniform(2000, 20000, len(sp['s'])), 0.0)
        write_frames(ds, sp, I, 1.0, 0.7, 0.2, np.zeros(SHAPE, bool))

        class GS(object):
            grain = g0
        self.saved = (hp.dataset_and_pars, hp.labels_recon, sinogram.read_h5)
        hp.dataset_and_pars = lambda cfg: (ds, PARS)
        hp.labels_recon = lambda f, g: np.zeros((3, 3), int)
        sinogram.read_h5 = lambda f, d, p: [GS()]
        cfgfile = os.path.join(self.tmp, 'config.json')
        with open(cfgfile, 'w') as f:
            json.dump({"dataset": "fake", "phase": "fake", "outdir": os.path.join(self.tmp, 'out'), "nproc": 1,
                       "harvest": {"reader": "h5py", "ds_max": 0.75, "fit_every": 1,
                                   "box": {"psf": 2.0, "phi0": 0.2, "bandpass": 5e-3}}}, f)
        self.cfg = hp.load_config(cfgfile)

    def tearDown(self):
        from ImageD11.sinograms import sinogram
        if hasattr(self, 'saved'):
            hp.dataset_and_pars, hp.labels_recon, sinogram.read_h5 = self.saved
        shutil.rmtree(self.tmp)

    def test_harvest_fit(self):
        hp.step_harvest(self.cfg)
        with h5py.File(hp.out(self.cfg, 'harvest.h5'), 'r') as f:
            self.assertEqual(list(f['gid'][()]), [0])
            self.assertIn('shoebox_fit_psf', f.attrs)
            psf = f.attrs['shoebox_fit_psf']
            self.assertLess(abs(psf - 1.0), 0.3)          # started at 2 px; the spots are 1 px wide
            self.assertGreater(f['merged/I'][()].max(), 1000)


class TestGrainSpreads(unittest.TestCase):

    def test_spreads(self):
        rng = np.random.RandomState(0)
        n = 400
        ubi0 = np.eye(3) * 4.0
        ubis, epss = [], []
        for k in range(n):
            ax = rng.normal(size=3)
            ax /= np.linalg.norm(ax)
            t = np.radians(rng.normal(0, 0.1))
            K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
            R = np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * np.dot(K, K)
            ubis.append(np.dot(ubi0, R.T))
            e = rng.normal(0, 1e-3, (3, 3))
            epss.append((e + e.T) / 2)
        u, phi, eps = hp.grain_spreads(np.array(ubis), np.array(epss), np.ones(n, bool))
        self.assertTrue(np.allclose(u, ubi0, atol=0.01))
        self.assertLess(abs(phi - 0.1) / 0.1, 0.1)          # rotations of N(0, 0.1 deg) about random axes
        self.assertGreater(eps, 0)


if __name__ == '__main__':
    unittest.main()
