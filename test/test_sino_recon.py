"""Reconstruction-level tests for the scanning 3DXRD geometry.

Derived from ImageD11/sandbox/scanning_geometry.ipynb. test_geometry.py covers the
coordinate transforms on their own; this file drives them through an actual
reconstruction, because several of the conventions (where the rotation axis sits in
the output array, how `shift` is applied, how wide the detector has to be) only bite
once you reconstruct.

Two acquisition types are covered:

  full  - 180 degrees of omega, dty scanned across the whole sample
  half  - 360 degrees of omega, dty scanned across half the sample

crossed with both parities of the reconstruction size. The parity matters:
roi_iradon puts the rotation axis on integer index N//2, ASTRA centres the volume on
the origin at (N-1)/2, and those only agree when N is odd.

sino_shift_and_pad now always returns a pad that makes N = ny + pad odd, so the even
cases are reached by adding one to its pad ("extra_pad" below). That is not a
contrived situation: the notebooks let the user set pad by hand through
update_recon_parameters, so the even path still has to reconstruct correctly.
"""

from __future__ import print_function

import unittest
import warnings

import numpy as np

from ImageD11.sinograms import geometry
from ImageD11.sinograms.roi_iradon import run_iradon

try:
    import astra
except ImportError:
    astra = None

if astra is not None:
    from ImageD11.sinograms.sinogram import run_astra

HAVE_ASTRA = astra is not None
HAVE_CUDA = bool(astra.use_cuda()) if HAVE_ASTRA else False

if not HAVE_ASTRA:
    warnings.warn("astra not installed - skipping all ASTRA reconstruction tests")
elif not HAVE_CUDA:
    warnings.warn("no GPU available to astra - skipping the *_CUDA reconstruction tests")

YSTEP = 10.0            # um per dty step
SX, SY = 300.0, -400.0  # grain centre of mass in the sample frame (radius 500 um = 50 px)

# how close the reconstructed grain has to land to geometry.sample_to_recon.
# well inside a pixel: a half-pixel convention slip has to fail this.
TOL_PX = 0.35

# Each case is a synthetic acquisition. y0 is the true rotation axis position; the
# scan is not centred on it, which is the y0 error. The two y0 values in each pair
# differ by half a dty step, so the pair also covers two different sub-pixel shifts.
# extra_pad is added to the pad from sino_shift_and_pad to force an even output size.
SCANS = {
    "full_odd":  dict(omega=np.arange(0.0, 180.0, 1.0),
                      ymin=13400.0, ymax=14600.0, y0=13955.0, extra_pad=0),
    "full_even": dict(omega=np.arange(0.0, 180.0, 1.0),
                      ymin=13400.0, ymax=14600.0, y0=13950.0, extra_pad=1),
    "half_odd":  dict(omega=np.arange(0.0, 360.0, 1.0),
                      ymin=13400.0, ymax=14060.0, y0=13995.0, extra_pad=0),
    "half_even": dict(omega=np.arange(0.0, 360.0, 1.0),
                      ymin=13400.0, ymax=14060.0, y0=13990.0, extra_pad=1),
}


def build_scan(name, sx=SX, sy=SY, ystep=YSTEP):
    """Synthesise the sinogram of a single grain for one of the SCANS above.

    Returns a dict with everything the reconstruction routines need, plus the
    (ri, rj) position the grain is expected to reconstruct at.
    """
    p = SCANS[name]
    ymin, ymax, y0, omega = p["ymin"], p["ymax"], p["y0"], p["omega"]

    ny = int((ymax - ymin) // ystep) + 1
    ybincens = np.linspace(ymin, ymax, ny)

    dty = geometry.dty_values_grain_in_beam(sx, sy, y0, omega)
    dtyi = geometry.dty_to_dtyi(dty, ystep, ybincens[0])

    # a sharp ridge following the grain's sinusoid, as in the notebook
    rows = np.arange(ny)[:, None]
    sino = 1.0 / (50.0 * np.cbrt(np.abs(rows - dtyi[None, :])) + 0.01)
    sino /= sino.max()   # SIRT/BP clamp the reconstruction to [0, 1]

    shift, pad = geometry.sino_shift_and_pad(y0, ny, ymin, ystep)
    pad = pad + p["extra_pad"]
    recon_shape = (ny + pad, ny + pad)

    return dict(name=name, sino=sino, omega=omega, ny=ny, ystep=ystep,
                y0=y0, ymin=ymin, dty=dty, dtyi=dtyi, ybincens=ybincens,
                shift=shift, pad=pad, recon_shape=recon_shape,
                expected=geometry.sample_to_recon(sx, sy, recon_shape, ystep))


def peak_position(recon):
    """Sub-pixel position of the brightest voxel, by parabolic fit on each axis."""
    i, j = np.unravel_index(np.argmax(recon), recon.shape)
    if i < 1 or j < 1 or i >= recon.shape[0] - 1 or j >= recon.shape[1] - 1:
        return float(i), float(j)

    def offset(lo, mid, hi):
        denom = lo - 2.0 * mid + hi
        return 0.0 if denom == 0 else 0.5 * (lo - hi) / denom

    return (i + offset(recon[i - 1, j], recon[i, j], recon[i + 1, j]),
            j + offset(recon[i, j - 1], recon[i, j], recon[i, j + 1]))


class TestScanMatrix(unittest.TestCase):
    """The cases have to actually cover what they claim to cover."""


    def test_both_parities_present(self):
        parities = {}
        for name in SCANS:
            n = build_scan(name)["recon_shape"][0]
            parities[name] = n % 2
        self.assertEqual(parities["full_odd"], 1)
        self.assertEqual(parities["full_even"], 0)
        self.assertEqual(parities["half_odd"], 1)
        self.assertEqual(parities["half_even"], 0)

    def test_y0_error_is_real(self):
        # a shift of zero would make these tests vacuous
        for name in SCANS:
            s = build_scan(name)
            self.assertGreater(abs(s["shift"]), 0.5, name)

    def test_half_scans_are_actually_half(self):
        # the grain's sinusoid should run off the end of a half acquisition
        # and stay inside a full one
        for name in ("full_odd", "full_even"):
            s = build_scan(name)
            self.assertLessEqual(s["dtyi"].max(), s["ny"] - 1, name)
            self.assertGreaterEqual(s["dtyi"].min(), 0, name)
        for name in ("half_odd", "half_even"):
            s = build_scan(name)
            self.assertGreater(s["dtyi"].max(), s["ny"] - 1, name)

    # we now enforce odd sinograms in sino_shift_and_pad
    def test_sino_shift_and_pad_is_always_odd(self):
        # the odd-size invariant is what lets recon_to_step, recon_bins and ASTRA
        # agree without a half-pixel correction, so pin it down directly rather
        # than only through the scans above
        ystep = YSTEP
        for ny in range(40, 140):
            for offset_steps in np.arange(-3.0, 3.01, 0.125):
                y0 = (ny // 2 + offset_steps) * ystep
                shift, pad = geometry.sino_shift_and_pad(y0, ny, 0.0, ystep)
                self.assertEqual((ny + pad) % 2, 1,
                                 "ny=%d offset=%.3f gave even size %d"
                                 % (ny, offset_steps, ny + pad))
                # must still be wide enough to hold the whole sample
                self.assertGreaterEqual(pad, 2 * abs(shift),
                                        "ny=%d offset=%.3f pad too small"
                                        % (ny, offset_steps))

class TestGeometryRoundTrips(unittest.TestCase):

    def test_sample_recon_round_trip(self):
        for name in SCANS:
            s = build_scan(name)
            ri, rj = geometry.sample_to_recon(SX, SY, s["recon_shape"], YSTEP)
            sx, sy = geometry.recon_to_sample(ri, rj, s["recon_shape"], YSTEP)
            self.assertAlmostEqual(sx, SX, places=6, msg=name)
            self.assertAlmostEqual(sy, SY, places=6, msg=name)

    def test_sx_sy_y0_recovered_from_dty(self):
        for name in SCANS:
            s = build_scan(name)
            sx, sy, y0 = geometry.sx_sy_y0_from_dty_omega(s["dty"], s["omega"])
            self.assertAlmostEqual(sx, SX, places=4, msg=name)
            self.assertAlmostEqual(sy, SY, places=4, msg=name)
            self.assertAlmostEqual(y0, s["y0"], places=4, msg=name)


class ReconstructionMixin(object):
    """Shared checks: whatever the backend, the grain has to land on the marker."""

    def check_locates_grain(self, recon, scan, label):
        self.assertEqual(recon.shape, scan["recon_shape"],
                         "%s %s: reconstruction changed shape" % (scan["name"], label))
        ri, rj = peak_position(recon)
        ei, ej = scan["expected"]
        self.assertLess(abs(ri - ei), TOL_PX,
                        "%s %s: i off by %.3f px" % (scan["name"], label, ri - ei))
        self.assertLess(abs(rj - ej), TOL_PX,
                        "%s %s: j off by %.3f px" % (scan["name"], label, rj - ej))


class TestIradon(unittest.TestCase, ReconstructionMixin):

    def test_locates_grain(self):
        for name in SCANS:
            s = build_scan(name)
            recon = run_iradon(s["sino"], s["omega"], pad=s["pad"],
                               shift=s["shift"], workers=1)
            self.check_locates_grain(recon, s, "iradon")

    def test_recon_position_reproduces_dty(self):
        # walk the whole loop: sinogram -> reconstruction -> back to a dty curve
        for name in SCANS:
            s = build_scan(name)
            recon = run_iradon(s["sino"], s["omega"], pad=s["pad"],
                               shift=s["shift"], workers=1)
            ri, rj = np.unravel_index(np.argmax(recon), recon.shape)
            dty_calc = geometry.recon_omega_to_dty(ri, rj, s["omega"], s["y0"],
                                                  s["recon_shape"], s["ystep"])
            self.assertLess(np.abs(s["dty"] - dty_calc).max(), s["ystep"], name)


@unittest.skipUnless(HAVE_ASTRA, "astra is not installed")
class TestAstraCPU(unittest.TestCase, ReconstructionMixin):

    METHODS = (("FBP", 1), ("BP", 1), ("SIRT", 150))

    def test_locates_grain(self):
        for name in SCANS:
            s = build_scan(name)
            for method, niter in self.METHODS:
                recon = run_astra(s["sino"], s["omega"], shift=s["shift"],
                                  pad=s["pad"], niter=niter, astra_method=method)
                self.check_locates_grain(recon, s, method)

    def test_agrees_with_iradon(self):
        for name in SCANS:
            s = build_scan(name)
            ref = peak_position(run_iradon(s["sino"], s["omega"], pad=s["pad"],
                                           shift=s["shift"], workers=1))
            for method, niter in self.METHODS:
                got = peak_position(run_astra(s["sino"], s["omega"], shift=s["shift"],
                                              pad=s["pad"], niter=niter,
                                              astra_method=method))
                self.assertLess(abs(got[0] - ref[0]), TOL_PX,
                                "%s %s vs iradon: i" % (name, method))
                self.assertLess(abs(got[1] - ref[1]), TOL_PX,
                                "%s %s vs iradon: j" % (name, method))

    def test_mask_is_honoured(self):
        # FBP silently ignores astra's ReconstructionMaskId, so run_astra has to
        # apply the mask itself. This catches it if that handling is dropped.
        s = build_scan("half_odd")
        n = s["recon_shape"][0]
        ii, jj = np.mgrid[0:n, 0:n]
        mask = ((ii - n // 2) ** 2 + (jj - n // 2) ** 2) < (n // 4) ** 2
        for method, niter in self.METHODS:
            recon = run_astra(s["sino"], s["omega"], shift=s["shift"], pad=s["pad"],
                              mask=mask, niter=niter, astra_method=method)
            self.assertEqual(recon.shape, mask.shape, method)
            self.assertEqual(np.abs(recon[~mask]).max(), 0.0,
                             "%s: reconstructed outside the mask" % method)


@unittest.skipUnless(HAVE_CUDA, "no GPU available to astra")
class TestAstraCUDA(unittest.TestCase, ReconstructionMixin):

    METHODS = (("FBP_CUDA", 1), ("BP_CUDA", 1), ("SIRT_CUDA", 150), ("EM_CUDA", 150))

    def test_locates_grain(self):
        for name in SCANS:
            s = build_scan(name)
            for method, niter in self.METHODS:
                recon = run_astra(s["sino"], s["omega"], shift=s["shift"],
                                  pad=s["pad"], niter=niter, astra_method=method)
                self.check_locates_grain(recon, s, method)

    def test_mask_is_honoured(self):
        # FBP_CUDA accepts ReconstructionMaskId and never applies it, without even
        # a warning; EM_CUDA rejects it outright. Both have to be masked by hand.
        s = build_scan("half_odd")
        n = s["recon_shape"][0]
        ii, jj = np.mgrid[0:n, 0:n]
        mask = ((ii - n // 2) ** 2 + (jj - n // 2) ** 2) < (n // 4) ** 2
        for method, niter in self.METHODS:
            recon = run_astra(s["sino"], s["omega"], shift=s["shift"], pad=s["pad"],
                              mask=mask, niter=niter, astra_method=method)
            self.assertEqual(recon.shape, mask.shape, method)
            self.assertEqual(np.abs(recon[~mask]).max(), 0.0,
                             "%s: reconstructed outside the mask" % method)


if __name__ == "__main__":
    unittest.main()