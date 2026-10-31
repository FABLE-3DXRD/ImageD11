"""Synthetic tests for ImageD11.sinograms.profile_fit: a grain of known shape, reflections of known intensity seen
versus dty with Poisson counts on a background, a few rows carrying another grain's peak."""
import unittest

import numpy as np

from ImageD11.sinograms import profile_fit

NY, YSTEP, Y0, YMIN = 61, 1.0, 0.0, -30.0
SHAPE = (61, 61)


def grain(shape=SHAPE):
    ii, jj = np.mgrid[0:shape[0], 0:shape[1]]
    rho = (((ii - 33) / 6.0) ** 2 + ((jj - 27) / 4.0) ** 2 <= 1).astype(float)     # an off-centre ellipse
    return rho


def rows(rho, I, omegas, bg, rng, outliers=0):
    """(nscan=NY, nrefl) box / bg / var / nframes arrays for reflections of intensity I, all dty bins measured,
    the normalised profiles, and the columns given an outlier"""
    P = profile_fit.project(rho, omegas, Y0, YMIN, YSTEP, NY)          # (nrefl, NY)
    P /= P.sum(1, keepdims=True)
    mean = (I[:, None] * P).T                                          # (NY, nrefl)
    box = rng.poisson(mean + bg).astype(float)
    bgest = np.full_like(box, float(bg))
    c = rng.randint(0, len(I), outliers)
    if outliers:
        box[rng.randint(0, NY, outliers), c] += 50 * bg + 200          # another grain's peak in the box
    return box, bgest, box + bg / 100.0, np.ones_like(box), P, c


class TestProfileFit(unittest.TestCase):

    def test_projector_matches_project(self):
        rho = grain() * np.random.RandomState(0).uniform(0.5, 1.5, SHAPE)
        om = np.linspace(-90, 90, 13)
        A, (ri, rj) = profile_fit.shape_projector(rho > 0, om, Y0, YMIN, YSTEP, NY)
        p1 = A.dot(rho[ri, rj]).reshape(len(om), NY)
        p2 = profile_fit.project(rho, om, Y0, YMIN, YSTEP, NY)
        self.assertLess(np.abs(p1 - p2).max(), 1e-9 * p2.max())
        # the beam footprint conserves the density
        self.assertTrue(np.allclose(p2.sum(1), rho.sum()))

    def test_fit_dty_unbiased_and_calibrated(self):
        rng = np.random.RandomState(1)
        rho = grain()
        n = 400
        om = rng.uniform(-90, 90, n)
        I = rng.uniform(0, 60, n)                                       # weak, on a background of 20 per row
        box, bg, var, nfr, P, hitc = rows(rho, I, om, 20.0, rng)
        r = profile_fit.row_estimates(box, bg, var, nfr, np.ones(n) * 100)
        f = profile_fit.fit_dty(r['net'], r['v_sum'], nfr > 0, P.T)
        z = (f['I'] - I) / f['sigI']
        self.assertLess(abs(np.mean(z)), 0.15)
        self.assertGreater(np.std(z), 0.85)
        self.assertLess(np.std(z), 1.15)
        # more precise than the box sum over all rows
        s_sum = np.sqrt(np.sum(r['v_sum'], 0))
        self.assertLess(np.median(f['sigI']), 0.8 * np.median(s_sum))

    def test_outlier_rows_rejected(self):
        rng = np.random.RandomState(2)
        rho = grain()
        n = 200
        om = rng.uniform(-90, 90, n)
        I = rng.uniform(0, 60, n)
        box, bg, var, nfr, P, hitc = rows(rho, I, om, 20.0, rng, outliers=60)
        r = profile_fit.row_estimates(box, bg, var, nfr, np.ones(n) * 100)
        f = profile_fit.fit_dty(r['net'], r['v_sum'], nfr > 0, P.T)
        Is = r['net'].sum(0)
        hit = np.zeros(n, bool)
        hit[hitc] = True
        # the contaminated reflections keep their intensity in the fit, not in the box sum
        self.assertLess(np.median(np.abs(f['I'][hit] - I[hit])), 0.3 * np.median(np.abs(Is[hit] - I[hit])))

    def test_intensities_hybrid(self):
        rng = np.random.RandomState(3)
        rho = grain()
        n = 300
        om = rng.uniform(-90, 90, n)
        I = np.concatenate([rng.uniform(0, 30, n // 2), rng.uniform(5000, 20000, n // 2)])
        box, bg, var, nfr, P, hitc = rows(rho, I, om, 20.0, rng)
        out = profile_fit.intensities(box, bg, var, nfr, np.ones(n) * 100, P.T)
        weak = I < 30
        self.assertTrue(np.all(~out['profile_used'][~weak]))          # strong: plain box sums
        self.assertGreater(out['profile_used'][weak].mean(), 0.8)
        self.assertTrue(np.allclose(out['I'][~weak], out['I_sum'][~weak]))
        z = (out['I'] - I) / out['sigI']
        self.assertLess(abs(np.mean(z[weak])), 0.2)
        self.assertLess(np.std(z[weak]), 1.2)

    def test_mlem_shape(self):
        rng = np.random.RandomState(4)
        rho = grain()
        n = 120
        om = rng.uniform(-90, 90, n)
        I = rng.uniform(2000, 8000, n)
        box, bg, var, nfr, P, hitc = rows(rho, I, om, 20.0, rng, outliers=40)
        support = profile_fit.support_disc((3.0, 3.0), 14.0, SHAPE, YSTEP)
        self.assertTrue(np.all(support[rho > 0]))
        rec, nout = profile_fit.mlem_shape(box.T, bg.T, nfr.T > 0, I, om, support, Y0, YMIN, YSTEP, niter=60)
        self.assertGreater(nout, 10)
        found = rec > 0.5 * rec[rho > 0].mean()
        iou = (found & (rho > 0)).sum() / float((found | (rho > 0)).sum())
        self.assertGreater(iou, 0.8)
        self.assertLess(rec[~(rho > 0)].sum() / rec.sum(), 0.1)


if __name__ == '__main__':
    unittest.main()
