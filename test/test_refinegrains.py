from __future__ import print_function, division

import os
import sys
import unittest
import numpy as np
from ImageD11 import refinegrains

HERE = os.path.dirname(os.path.abspath(__file__))


def total_intensity_full_columns(colfile, indices, tth_range, ntrim=2):
    """ compute_total_intensity as it was: products of full columns, then take """
    raw = colfile.npixels * colfile.avg_intensity
    intensities = np.take(raw * colfile.Lorentz_per_grain, indices)
    sigma_i = np.sum(intensities)
    tth_vals = np.take(colfile.tth_per_grain, indices)
    intensities = np.compress((tth_vals > min(tth_range)) & (tth_vals < max(tth_range)), intensities)
    if len(intensities) < 1:
        return "no peaks"
    intensities.sort()
    if len(intensities) > ntrim * 2 + 1:
        intensities = intensities[ntrim:-ntrim]
    return "sum_of_all = %f , middle %d from %f to %f in tth: median = %f , min = %f , max = %f , mean = %f , std = %f , n = %d" % (
        sigma_i, len(intensities), min(tth_range), max(tth_range), intensities[len(intensities) // 2],
        intensities.min(), intensities.max(), intensities.mean(), intensities.std(), intensities.shape[0])


class TestAssignLabels(unittest.TestCase):
    def setUp(self):
        d = os.path.join(HERE, 'makemap')
        self.o = refinegrains.refinegrains(tolerance=0.05, intensity_tth_range=(0, 180))
        ss = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        try:
            self.o.loadparameters(os.path.join(d, 'test.prm'))
            self.o.readubis(os.path.join(d, 'map.ubi'))
            self.o.loadfiltered(os.path.join(d, 'test.flt'))
            self.o.generate_grains()
            self.o.assignlabels(quiet=True)
        finally:
            sys.stdout.close()
            sys.stdout = ss

    def test_grain_peaks_and_intensity(self):
        # each grain's peaks are the rows labelled with it, in order, and its
        # intensity summary is the one from the full columns
        s = self.o.scannames[0]
        cf = self.o.scandata[s]
        assigned = 0
        for g in self.o.grainnames:
            gr = self.o.grains[(g, s)]
            ind = np.flatnonzero(cf.labels == g)
            self.assertTrue(np.array_equal(gr.ind, ind))
            self.assertEqual(gr.intensity_info, total_intensity_full_columns(cf, ind, (0, 180)))
            assigned += len(ind)
        self.assertGreater(assigned, 0)


if __name__ == '__main__':
    unittest.main()
