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


if __name__ == '__main__':
    unittest.main()
