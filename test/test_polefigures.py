
from __future__ import print_function

import unittest
import sys

if int(sys.version_info.major) == 2:
    raise unittest.SkipTest('Skipping polefigures tests on Python 2')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pl
from ImageD11.sinograms import polefigures as pf
import ImageD11.columnfile


class Pars(object):
    def get(self, name):
        return {'wavelength': 1.0}[name]


class PMap(object):
    def __init__(self):
        ubi = np.full((3, 3, 3, 3), np.nan)
        ubi[0, 0] = np.eye(3) * 4.0
        ubi[1, 2] = np.eye(3) * 4.0
        self.best_ubi = ubi


class FakeDS(object):
    obinedges = np.arange(0, 361, 1.0)
    ybinedges = np.linspace(-5, 5, 11)
    ystep = 1.0


class TestWrap(unittest.TestCase):
    def test_wrap(self):
        w = pf.wrap_omega(np.array([-170., 0., 10., 359.]), 0.)
        self.assertTrue(np.allclose(w, [190., 0., 10., 359.]))
        w = pf.wrap_omega(np.array([-90., 0., 89., 190., 270.]), -90.)
        self.assertTrue(np.allclose(w, [-90., 0., 89., 190., -90.]))


class TestCalc(unittest.TestCase):
    def test_unreachable_is_nan(self):
        hkls = np.array([[1., 1., 0.], [10., 0., 0.]])
        mask, tth, e1, e2, o1, o2 = pf.calc_gve_angles(PMap(), hkls, Pars())
        self.assertEqual(mask.sum(), 2)
        for a in (tth, e1, e2, o1, o2):
            self.assertEqual(a.shape, (2, 2))
            self.assertTrue(np.isfinite(a[:, 0]).all())
            self.assertTrue(np.isnan(a[:, 1]).all())


class TestPlots(unittest.TestCase):
    def setUp(self):
        hkls = np.array([[1., 1., 0.], [10., 0., 0.]])
        self.ds = FakeDS()
        self.mask, tth, self.e1, self.e2, self.o1, self.o2 = \
            pf.calc_gve_angles(PMap(), hkls, Pars())
        self.dty1, self.dty2 = pf.calc_dty(self.ds, self.mask, self.o1, self.o2, 0.)
        rs = np.random.RandomState(0)
        n = 50
        cf = ImageD11.columnfile.newcolumnfile(
            ['ds', 'eta', 'omega', 'dty', 'sum_intensity'])
        cf.bigarray = np.array([np.full(n, 0.25 * np.sqrt(2)),
                                rs.uniform(-180, 180, n),
                                rs.uniform(0, 360, n),
                                rs.uniform(-4, 4, n),
                                rs.uniform(1, 10, n)])
        cf.nrows = n
        cf.ncols = 5
        cf.set_attributes()
        self.cf = cf

    class UC(object):
        ringds = [0.25 * np.sqrt(2)]
        ringhkls = {0.25 * np.sqrt(2): [(1, 1, 0)]}

    def tearDown(self):
        pl.close('all')

    def test_polefigures(self):
        pf.plot_ring_polefigures(self.UC, [1], self.e1[:, :1], self.e2[:, :1],
                                 self.o1[:, :1], self.o2[:, :1], self.cf,
                                 ds=self.ds)
        ax = pl.gcf().axes
        self.assertTrue(ax[0].get_shared_x_axes().joined(ax[0], ax[1]))
        self.assertTrue(ax[0].get_shared_y_axes().joined(ax[0], ax[1]))
        self.assertEqual(ax[0].get_ylim(), ax[1].get_ylim())
        # every reachable calc peak lands in a bin (none lost below 0 deg)
        self.assertEqual(ax[0].collections[0].get_array().sum(),
                         np.isfinite(self.o1[:, :1]).sum() * 2)
        pf.plot_ring_polefigures(self.UC, [1], self.e1[:, :1], self.e2[:, :1],
                                 self.o1[:, :1], self.o2[:, :1], self.cf)

    def test_sinograms(self):
        pf.plot_ring_sinograms(self.ds, self.cf, self.UC, [1],
                               self.o1[:, :1], self.o2[:, :1],
                               self.dty1[:, :1], self.dty2[:, :1])
        ax = pl.gcf().axes
        self.assertEqual(ax[0].collections[0].get_array().sum(),
                         np.isfinite(self.o1[:, :1]).sum() * 2)
        self.assertTrue(ax[0].get_shared_x_axes().joined(ax[0], ax[1]))
        self.assertTrue(ax[0].get_shared_y_axes().joined(ax[0], ax[1]))


if __name__ == '__main__':
    unittest.main()
