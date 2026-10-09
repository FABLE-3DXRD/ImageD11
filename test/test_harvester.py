"""Synthetic test for ImageD11.sinograms.harvester: Gaussian spots with known intensity and position are
painted on Poisson background frames at the predicted positions of a grain, written with bitshuffle
compression, and harvested back."""
import os
import shutil
import tempfile
import unittest

import numpy as np
import h5py

try:
    import hdf5plugin
except ImportError:
    hdf5plugin = None

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
    rng = np.random.default_rng(seed)
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
    kw = {}
    if hdf5plugin is not None:
        kw = dict(hdf5plugin.Bitshuffle(nelems=0, cname='lz4'))
    with h5py.File(os.path.join(ds.datapath, ds.imagefiles[0]), 'w') as h:
        h.create_dataset('data', data=data, chunks=(1,) + SHAPE, **kw)


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
        rng = np.random.default_rng(1)
        cls.truth = dict(sp)
        # only spots well inside the detector and the omega range are painted; the others must come out ~ 0
        cls.painted = ((sp['omega'] > -85) & (sp['omega'] < 85) & (sp['s'] > 20) & (sp['s'] < SHAPE[0] - 20) &
                       (sp['f'] > 20) & (sp['f'] < SHAPE[1] - 20))
        cls.truth['I'] = np.where(cls.painted, rng.uniform(200, 2000, len(sp['s'])), 0.0)
        cls.truth['s'] = sp['s'] + rng.uniform(-0.4, 0.4, len(sp['s']))     # real spots off the prediction
        cls.truth['f'] = sp['f'] + rng.uniform(-0.4, 0.4, len(sp['s']))
        write_frames(cls.ds, cls.truth, cls.truth['I'], 1.0, 0.7, 0.2, cls.mask)

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
        self.assertTrue(np.all(np.abs(m['sig_s'][clean] - 1.0) < 0.3))
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

    def test_sparse(self):
        if harvester.bslz4_to_sparse is None:
            self.skipTest('bslz4_to_sparse is not installed')
        self.check(*self.harvest('sparse'))


if __name__ == '__main__':
    unittest.main()
