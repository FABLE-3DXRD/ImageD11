"""
ImageD11.sinograms.harvester
============================

Integrate the predicted Bragg reflections of a list of grains directly from the raw detector frames of
an ImageD11 dataset (scanning 3DXRD, or 3DXRD with a single dty), reading every frame only once.

Segmentation (a threshold, then connected pixels) loses weak reflections and the weak tails of strong
ones. Here every predicted reflection of every grain is summed in a fixed box ("shoebox") on the raw
frames, with a local background, so weak and absent reflections are measured too. One pass over the
data serves all grains.

Typical use::

    import ImageD11.sinograms.dataset, ImageD11.grain
    from ImageD11.sinograms import harvester
    ds = ImageD11.sinograms.dataset.load(dsfile)
    grains = [ImageD11.grain.grain(ubi, translation=t) for ubi, t in ...]
    h = harvester.Harvester(ds, grains, pars, ds_max=1.5, radius=radii)
    box = h.auto_box()                   # or h.box = harvester.Box(half_s=6, half_f=6, half_omega=3)
    h.run(nproc=8, outfile='harvest.h5') # per (reflection, scan) rows written to outfile
    spots = h.merged()                   # per reflection: intensity, sigma, centroid, flags

The pieces
----------
* Prediction (`predict_spots`): every hkl with ``|g| <= ds_max`` and both omega solutions. The detector
  position is computed for the grain at its lab position when it is in the beam: the diffraction origin
  is (lab_x, 0, 0), a per-peak origin along x (a per-peak distance), as in point_by_point. The spatial
  distortion is inverted with the e2dx/e2dy maps.
* Which frames: a reflection is integrated on the frames within ``half_omega`` frames of its predicted
  omega, in the scans where its grain is in the beam (``|dty - dty_pred(omega)| <= radius``).
* Reading: each frame is decoded to its non-zero, unmasked pixels, either with ``bslz4_to_sparse``
  (bitshuffle lz4 / zstd chunks decoded straight to a sparse list; no threshold, cut = 0) or with
  h5py + hdf5plugin (dense read, then the non-zero pixels). Masked pixels never enter the sums.
* Box size: either fixed (`Box`), or computed for each reflection from the expected spot shape
  (`ShoeboxModel`): rocking width in omega from the grain orientation spread, the strain spread and the
  energy bandpass, scaled by 1/|sin eta|; radial / tangential extent on the detector from the strain and
  bandpass (2 tan(theta) eps) and the orientation spread, plus the detector point spread; then
  nsigma x sigma plus a margin for the prediction error. The background + peak box is bg_factor times the
  peak box. `fit_shoebox_model` fits the instrumental terms to the moments of strong reflections.
* Box and background: ``(2 half_s + 1) x (2 half_f + 1)`` pixels on each of the
  ``2 half_omega + 1`` frames. The background is estimated per frame from an annulus of width
  ``bg_width`` pixels around the box, with estimator ``'mean'`` or ``'clipped'``. ``'clipped'`` drops
  annulus pixels above the Poisson upper limit (probability ``p_clip``) of the annulus mean, so a
  neighbouring reflection does not raise the background. `auto_box` measures spot sizes and
  prediction offsets on strong reflections with a large pilot box and sets the box from them.
* Masks: each reflection records the number of masked pixels in its box (``n_masked``), whether its
  box touches the mask or the detector edge (``touch_mask``, ``edge``), and the net intensity on
  pixels adjacent to a masked pixel (``I_adj_mask``). ``n_overlap`` counts the boxes of other grains
  that overlap in (s, f, omega) while both grains can be in the beam.
* Positions: background-subtracted centroids (s, f, omega) and second moments of each box, merged
  over scans per reflection and corrected for distortion to (sc, fc). These can be used to refine
  the grain (they include weak reflections and are not cut by a threshold).
* Viewing: `spot_frames` returns the pixel stack of one reflection in one scan, and `show_spot`
  plots it with the box and the annulus marked.
"""

from __future__ import print_function, division

import math
import os
import sys
import time
import multiprocessing

import numpy as np
try:
    import hdf5plugin  # noqa: F401  (registers the bitshuffle / zstd filters for h5py; old versions need it before h5py)
except ImportError:
    hdf5plugin = None

import h5py
import numba

import ImageD11.transform
import ImageD11.blobcorrector
import ImageD11.unitcell
import ImageD11.columnfile

try:
    import bslz4_to_sparse
except ImportError:
    bslz4_to_sparse = None


# ------------------------------------------------------------------------------------- options
class Box(object):
    """Integration box ("shoebox") and background options.

    half_s, half_f : half sizes of the box in detector pixels (slow, fast): (2 n + 1) pixels.
    half_omega     : half size of the box in frames: (2 n + 1) frames.
    bg_width       : width in pixels of the background annulus around the box on each frame.
    bg_estimator   : 'mean' (all annulus pixels) or 'clipped' (annulus pixels above the Poisson upper
                     limit of the annulus mean, at probability p_clip, are left out).
    p_clip         : tail probability for 'clipped'.
    """

    def __init__(self, half_s=6, half_f=6, half_omega=3, bg_width=6, bg_estimator='clipped', p_clip=1e-3):
        self.half_s = int(half_s)
        self.half_f = int(half_f)
        self.half_omega = int(half_omega)
        self.bg_width = int(bg_width)
        if bg_estimator not in ('mean', 'clipped'):
            raise ValueError("bg_estimator must be 'mean' or 'clipped'")
        self.bg_estimator = bg_estimator
        self.p_clip = float(p_clip)

    def __repr__(self):
        return ('Box(half_s=%d, half_f=%d, half_omega=%d, bg_width=%d, bg_estimator=%r, p_clip=%g)' %
                (self.half_s, self.half_f, self.half_omega, self.bg_width, self.bg_estimator, self.p_clip))

    def sizes(self, spots, pars, omega_step):
        """per-reflection (half_s, half_f, half_omega, bg_width): the same for all"""
        n = len(spots['s'])
        return (np.full(n, self.half_s, np.int64), np.full(n, self.half_f, np.int64),
                np.full(n, self.half_omega, np.int64), np.full(n, self.bg_width, np.int64))


class ShoeboxModel(object):
    """Box size for each reflection from the expected spot shape (all spreads are rms, 1 sigma).

    Rocking width (rotation about z, beam along x), in omega:
        sigma_omega^2 = (frame * step)^2 + omega0^2 + (phi_tot^2 + (2 tan(theta) eps_tot)^2) / sin(eta)^2
    On the detector, radial and tangential spreads (converted to pixels with the local gradients of tth
    and eta), plus the point spread:
        radial:      d(2 theta) = 2 tan(theta) eps_tot
        tangential:  d(eta)     = phi_t / cos(theta)
    with phi_tot^2 = phi0^2 + phi_grain^2 (instrument + grain orientation spread) and
         eps_tot^2 = eps_grain^2 + (bandpass / 2.355)^2 (strain spread + energy bandpass FWHM).

    psf            : detector point spread (pixels, rms)
    bandpass       : relative energy bandpass dE/E (FWHM)
    phi0           : instrumental rocking spread (deg, rms): beam divergence and unresolved mosaicity
    phi0_t         : instrumental spread perpendicular to the scattering plane (deg); None = phi0
    omega0         : rocking spread that does not depend on eta (deg, rms): wobble or jitter of the rotation
    phi, eps       : grain orientation spread (deg) and strain spread, rms: scalars or arrays over grains
                     (e.g. the voxel misorientation and strain spreads of a point-by-point refinement)
    frame          : rms of the omega integration of one frame, in frames (1/sqrt(12))
    nsigma         : half size = nsigma * sigma + margin
    margin_px, margin_frames : allowance for the error of the predicted position
    bg_factor      : background + peak box = bg_factor x the peak box (half sizes); annulus >= bg_min px
    max_half_frames, sin_eta_min : limits near eta = 0, 180 where the rocking width diverges
    bg_estimator, p_clip : as for Box
    """

    def __init__(self, psf=0.9, bandpass=2e-3, phi0=0.085, phi0_t=None, phi=0.0, eps=0.0, frame=1 / np.sqrt(12),
                 nsigma=3.0, margin_px=1.0, margin_frames=2.0, bg_factor=2.0, bg_min=4, max_half_frames=40,
                 sin_eta_min=0.05, bg_estimator='clipped', p_clip=1e-3, omega0=0.0):
        self.psf, self.bandpass, self.phi0, self.phi0_t, self.omega0 = psf, bandpass, phi0, phi0_t, omega0
        self.phi, self.eps, self.frame = phi, eps, frame
        self.nsigma, self.margin_px, self.margin_frames = nsigma, margin_px, margin_frames
        self.bg_factor, self.bg_min = bg_factor, int(bg_min)
        self.max_half_frames, self.sin_eta_min = int(max_half_frames), sin_eta_min
        if bg_estimator not in ('mean', 'clipped'):
            raise ValueError("bg_estimator must be 'mean' or 'clipped'")
        self.bg_estimator, self.p_clip = bg_estimator, float(p_clip)

    def __repr__(self):
        return ('ShoeboxModel(psf=%g px, bandpass=%g, phi0=%g deg, omega0=%g deg, phi=%s deg, eps=%s, nsigma=%g, '
                'bg_factor=%g)' %
                (self.psf, self.bandpass, self.phi0, self.omega0, np.round(np.median(self.phi), 4),
                 np.round(np.median(self.eps), 6),
                 self.nsigma, self.bg_factor))

    def sigmas(self, spots, pars, omega_step):
        """rms widths (sigma_s, sigma_f in pixels, sigma_omega in frames) of each reflection"""
        g = spots['grain']
        phi_g = np.broadcast_to(np.asarray(self.phi, float), (g.max() + 1,))[g] if np.ndim(self.phi) else self.phi
        eps_g = np.broadcast_to(np.asarray(self.eps, float), (g.max() + 1,))[g] if np.ndim(self.eps) else self.eps
        phi_tot = np.radians(np.sqrt(self.phi0 ** 2 + np.square(phi_g)))
        phi_t = np.radians(np.sqrt((self.phi0 if self.phi0_t is None else self.phi0_t) ** 2 + np.square(phi_g)))
        eps_tot = np.sqrt(np.square(eps_g) + (self.bandpass / 2.355) ** 2)
        th = np.radians(spots['tth'] / 2)
        d2th = np.degrees(2 * np.tan(th) * eps_tot)            # radial, deg of 2 theta
        deta = np.degrees(phi_t / np.cos(th))                  # tangential, deg of eta
        # local gradients of tth and eta on the detector (deg per pixel), for a grain at the rotation axis
        sc, fc = spots['sc'], spots['fc']
        def te(ds_, df_):
            return tth_eta_from_pixels(sc + ds_, fc + df_, 0.0, pars)
        t_sp, e_sp = te(0.5, 0); t_sm, e_sm = te(-0.5, 0); t_fp, e_fp = te(0, 0.5); t_fm, e_fm = te(0, -0.5)
        gt = np.array([t_sp - t_sm, t_fp - t_fm])
        ge = np.array([(e_sp - e_sm + 180) % 360 - 180, (e_fp - e_fm + 180) % 360 - 180])
        R = d2th / np.linalg.norm(gt, axis=0)                   # radial spread in pixels
        T = deta / np.linalg.norm(ge, axis=0)                   # tangential spread in pixels
        r = gt / np.linalg.norm(gt, axis=0); t = ge / np.linalg.norm(ge, axis=0)
        sig_s = np.sqrt(self.psf ** 2 + (R * r[0]) ** 2 + (T * t[0]) ** 2)
        sig_f = np.sqrt(self.psf ** 2 + (R * r[1]) ** 2 + (T * t[1]) ** 2)
        se = np.maximum(np.abs(np.sin(np.radians(spots['eta']))), self.sin_eta_min)
        rock = np.sqrt(phi_tot ** 2 + (2 * np.tan(th) * eps_tot) ** 2) / se
        sig_w = np.sqrt(self.frame ** 2 + (self.omega0 / omega_step) ** 2 + (np.degrees(rock) / omega_step) ** 2)
        return sig_s, sig_f, sig_w

    def sizes(self, spots, pars, omega_step):
        """per-reflection (half_s, half_f, half_omega, bg_width)"""
        sig_s, sig_f, sig_w = self.sigmas(spots, pars, omega_step)
        hs = np.ceil(self.nsigma * sig_s + self.margin_px).astype(np.int64)
        hf = np.ceil(self.nsigma * sig_f + self.margin_px).astype(np.int64)
        ho = np.minimum(np.ceil(self.nsigma * sig_w + self.margin_frames), self.max_half_frames).astype(np.int64)
        bw = np.maximum(np.ceil((self.bg_factor - 1) * np.maximum(hs, hf)), self.bg_min).astype(np.int64)
        return hs, hf, ho, bw


def fit_shoebox_model(outfile, pars, model=None, min_snr=30.0, sin_eta_min=0.3, nmax=200000, seed=0, omega_step=None):
    """Fit the instrumental terms of a ShoeboxModel (psf, phi0 (rocking), phi0_t (perpendicular to the scattering
    plane), omega0 (rocking, independent of eta) and the instrumental d-spacing spread, reported as `bandpass` =
    2.355 x eps0) to the per-scan box moments of strong, clean reflections in a harvest file written by
    Harvester.run(outfile=...). The grain terms (model.phi, model.eps) are kept fixed. The residuals are
    log(model / measured) widths: symmetric, so the broad scatter of the measured rocking widths does not bias
    the fit (a residual relative to the measured width favours the narrow ones).
    Rows with |sin eta| < sin_eta_min are left out: their rocking curves are long and may be truncated by
    the box used for the harvest. Returns (fitted ShoeboxModel, dict of fit details)."""
    import copy
    import scipy.optimize
    model = ShoeboxModel() if model is None else model
    with h5py.File(outfile, 'r') as h:
        sp = dict((k, h['spots'][k][()]) for k in ('grain', 'tth', 'eta', 'sc', 'fc', 'n_overlap', 'touch_mask', 'edge'))
        ps = h['per_scan']
        spot = ps['spot'][()]
        I = ps['box'][()] - ps['bg'][()]
        good = I > min_snr * np.sqrt(np.maximum(ps['var'][()], 1))
        good &= (sp['n_overlap'][spot] == 0) & ~sp['touch_mask'][spot] & ~sp['edge'][spot]
        good &= np.abs(np.sin(np.radians(sp['eta'][spot]))) >= sin_eta_min
        rows = np.flatnonzero(good)
        if len(rows) > nmax:
            rows = np.sort(np.random.RandomState(seed).choice(rows, nmax, replace=False))
        mom = dict((k, ps[k][()][rows]) for k in ('ms', 'mss', 'mf', 'mff', 'mw', 'mww'))
        step = omega_step
        if step is None and 'box_omega_step' in h.attrs:
            step = float(h.attrs['box_omega_step'])
    I = I[rows]
    sub = dict((k, v[spot[rows]]) for k, v in sp.items())
    meas_s = np.sqrt(np.maximum(mom['mss'] / I - (mom['ms'] / I) ** 2, 0))
    meas_f = np.sqrt(np.maximum(mom['mff'] / I - (mom['mf'] / I) ** 2, 0))
    meas_w = np.sqrt(np.maximum(mom['mww'] / I - (mom['mw'] / I) ** 2, 0))          # degrees
    if step is None:
        raise ValueError('the harvest file has no omega step (box_omega_step attribute): pass omega_step')
    meas_w = meas_w / step
    trial = copy.copy(model)

    ok = (meas_s > 0) & (meas_f > 0) & (meas_w > 0)
    sub = dict((k, v[ok]) for k, v in sub.items())
    meas_s, meas_f, meas_w = meas_s[ok], meas_f[ok], meas_w[ok]

    def setp(p):
        trial.psf, trial.phi0, trial.phi0_t, trial.bandpass, trial.omega0 = p[0], p[1], p[2], p[3] * 2.355, p[4]

    def resid(p):
        setp(p)
        a, b, c = trial.sigmas(sub, pars, step)
        return np.log(np.concatenate([a / meas_s, b / meas_f, c / meas_w]))
    p0 = [model.psf, model.phi0, model.phi0 if model.phi0_t is None else model.phi0_t, model.bandpass / 2.355,
          model.omega0]
    r = scipy.optimize.least_squares(resid, p0, bounds=([0, 0, 0, 0, 0], [10, 5, 5, 0.05, 5]), loss='soft_l1',
                                     f_scale=0.3)
    setp(r.x)
    a, b, c = trial.sigmas(sub, pars, step)
    info = dict(nrows=int(ok.sum()), psf=r.x[0], phi0=r.x[1], phi0_t=r.x[2], eps0=r.x[3], omega0=r.x[4],
                median_ratio_s=float(np.median(meas_s / a)), median_ratio_f=float(np.median(meas_f / b)),
                median_ratio_omega=float(np.median(meas_w / c)))
    return trial, info


# ------------------------------------------------------------------------------------- frame readers
class H5FrameReader(object):
    """Dense h5py read (any filter hdf5plugin can decode), reduced to the non-zero unmasked pixels."""

    method = 'h5py'

    def __init__(self, filename, path, good):
        self.h = h5py.File(filename, 'r')
        self.d = self.h[path]
        self.good = good.ravel().astype(bool)
        self.buf = np.empty(self.d.shape[1:], self.d.dtype)

    def __call__(self, i):
        """Frame i -> (flat indices, ascending; values as float32) of non-zero unmasked pixels"""
        self.d.read_direct(self.buf, np.s_[i])
        b = self.buf.ravel()
        b[~self.good] = 0
        idx = np.flatnonzero(b)
        return idx, b[idx].astype(np.float32)

    def close(self):
        self.h.close()


class SparseFrameReader(object):
    """bslz4_to_sparse: bitshuffle-lz4 or -zstd chunks decoded straight to the non-zero unmasked pixels"""

    method = 'sparse'

    def __init__(self, filename, path, good):
        if bslz4_to_sparse is None:
            raise ImportError('bslz4_to_sparse is not available')
        self.h = h5py.File(filename, 'r')
        self.d = self.h[path]
        kw = {}
        if hasattr(bslz4_to_sparse, 'detect_codec'):     # releases up to 0.0.16 are lz4 only, with no codec argument
            kw['codec'] = bslz4_to_sparse.detect_codec(self.d)
        self.dec = bslz4_to_sparse.chunk2sparse(np.ascontiguousarray(good, dtype=np.uint8), dtype=self.d.dtype, **kw)

    def __call__(self, i):
        _, chunk = self.d.id.read_direct_chunk((i, 0, 0))
        npx, (vals, idx) = self.dec(chunk, 0)
        idx = idx[:npx].astype(np.int64)
        vals = vals[:npx].astype(np.float32)
        if npx > 1 and np.any(idx[1:] < idx[:-1]):      # the kernel needs ascending indices
            o = np.argsort(idx, kind='stable')
            idx, vals = idx[o], vals[o]
        return idx, vals

    def close(self):
        self.h.close()


def open_reader(filename, path, good, method='auto'):
    """Frame reader for one image file. method: 'sparse', 'h5py' or 'auto' (sparse if it can decode the
    file, otherwise h5py)."""
    if method == 'h5py':
        return H5FrameReader(filename, path, good)
    if method == 'sparse':
        return SparseFrameReader(filename, path, good)
    if method != 'auto':
        raise ValueError("method must be 'auto', 'sparse' or 'h5py'")
    if bslz4_to_sparse is not None:
        try:
            r = SparseFrameReader(filename, path, good)
            r(0)
            return r
        except Exception:
            pass
    return H5FrameReader(filename, path, good)


# ------------------------------------------------------------------------------------- prediction
def lab_x(translation, omega):
    """Lab x (along the beam) of a grain at sample position `translation` at rotation `omega` (deg), when it
    is in the beam: its diffraction origin is (lab_x, 0, 0). As in point_by_point / pbp_refine, this is
    used as a per-peak origin along x, i.e. a per-peak sample to detector distance (distance - lab_x)."""
    w = np.radians(omega)
    return translation[0] * np.cos(w) - translation[1] * np.sin(w)


def _no_translation(pars):
    """pars without the sample translation, which is handled by the per-ray origin"""
    return {k: v for k, v in pars.items() if k not in ('t_x', 't_y', 't_z')}


def detector_position(tth, eta, origin_x, pars):
    """Corrected detector pixel (sc, fc) of rays (tth, eta) leaving the lab point (origin_x, 0, 0);
    origin_x may be an array (one per ray)."""
    tth = np.asarray(tth, float)
    origin = np.zeros((3, tth.size))
    origin[0] = origin_x
    fc, sc = ImageD11.transform.compute_xyz_from_tth_eta(tth, eta, None, origin=origin, **_no_translation(pars))
    return sc, fc


def tth_eta_from_pixels(sc, fc, origin_x, pars):
    """tth, eta of corrected pixels (sc, fc) for diffraction origins (origin_x, 0, 0): compute_xyz_lab with
    a per-peak distance, as in point_by_point"""
    p = _no_translation(pars)
    p['distance'] = pars['distance'] - np.asarray(origin_x, float)
    xyz = ImageD11.transform.compute_xyz_lab(np.array([sc, fc], float), **p)
    return ImageD11.transform.compute_tth_eta_from_xyz(xyz, None)


def predict_spots(ubi, translation, pars, ds_max, shape, e2dx=None, e2dy=None, omega_range=None, margin=10):
    """All reflections of one grain up to ds_max that fall on the detector.

    ubi         : (3, 3) real-space cell of the grain (rows a, b, c), as in ImageD11.grain.ubi
    translation : (x, y, z) sample-frame position of the grain (same units as the distance)
    pars        : ImageD11 parameters dict (wavelength, wedge, chi, detector geometry)
    shape       : detector shape (slow, fast)
    omega_range : (min, max) in degrees, or None for all
    Returns a dict of arrays: h, k, l, ds, tth, eta, omega, sc, fc (corrected) and s, f (raw pixels).
    """
    ubi = np.asarray(ubi, float)
    B = np.linalg.inv(ubi)
    n = int(ds_max * np.linalg.norm(ubi, axis=1).max()) + 1
    r = np.arange(-n, n + 1)
    hkl = np.array(np.meshgrid(r, r, r, indexing='ij')).reshape(3, -1).T
    g = hkl.dot(B.T)
    dsv = np.linalg.norm(g, axis=1)
    keep = (dsv > 0) & (dsv <= ds_max)
    hkl, g, dsv = hkl[keep], g[keep], dsv[keep]
    tth, eta2, omega2 = ImageD11.transform.uncompute_g_vectors(g.T, pars['wavelength'],
                                                               wedge=pars.get('wedge', 0.0),
                                                               chi=pars.get('chi', 0.0))
    out = []
    for sol in (0, 1):
        ok = np.isfinite(omega2[sol]) & (tth > 0)
        out.append((hkl[ok], dsv[ok], tth[ok], eta2[sol][ok], omega2[sol][ok]))
    hkl = np.concatenate([o[0] for o in out])
    dsv = np.concatenate([o[1] for o in out])
    tth = np.concatenate([o[2] for o in out])
    eta = np.concatenate([o[3] for o in out])
    om = np.concatenate([o[4] for o in out])
    om = (om + 180.0) % 360.0 - 180.0
    if omega_range is not None:
        lo, hi = omega_range
        k = ((om >= lo) & (om <= hi)) | ((om + 360 >= lo) & (om + 360 <= hi)) | ((om - 360 >= lo) & (om - 360 <= hi))
        om = np.where((om < lo), om + 360, np.where(om > hi, om - 360, om))
        hkl, dsv, tth, eta, om = hkl[k], dsv[k], tth[k], eta[k], om[k]
    sc, fc = detector_position(tth, eta, lab_x(translation, om), pars)
    s, f = ImageD11.blobcorrector.invert_dxdy(sc, fc, e2dx, e2dy)
    on = (s > -margin) & (s < shape[0] + margin) & (f > -margin) & (f < shape[1] + margin)
    return dict(h=hkl[on, 0], k=hkl[on, 1], l=hkl[on, 2], ds=dsv[on], tth=tth[on], eta=eta[on],
                omega=om[on], sc=sc[on], fc=fc[on], s=s[on], f=f[on])


# ------------------------------------------------------------------------------------- integration kernel
# accumulator columns for each (reflection, scan)
ACC = ['box', 'bg', 'var', 'nframes', 'ms', 'mf', 'mw', 'mss', 'mff', 'mww', 'nclip', 'I_adj_mask', 'maxval',
       'nbgpix', 'nnz']
NACC = len(ACC)
_I = dict((k, i) for i, k in enumerate(ACC))


@numba.njit(cache=True)
def _poisson_upper(mean, p):
    """smallest integer t with P(X > t) <= p for X ~ Poisson(mean)"""
    if mean <= 0:
        return 0.0
    term = np.exp(-mean)
    cdf = term
    t = 0
    while 1.0 - cdf > p and t < 100000:
        t += 1
        term *= mean / t
        cdf += term
    return float(t)


@numba.njit(cache=True)
def _integrate_frame(idx, vals, nfast, adj, spots, s0, f0, hs, hf, bw, nbox, nann, S1s, S1f, S2s, S2f,
                     domega, scale, clipped, p_clip, acc):
    """Accumulate one frame into acc[spots] (see ACC); no profiles. See _integrate_frame_core."""
    z1 = np.zeros(1)
    zi = np.zeros(1, np.int64)
    _integrate_frame_core(idx, vals, nfast, adj, spots, s0, f0, hs, hf, bw, nbox, nann, S1s, S1f, S2s, S2f,
                          domega, scale, clipped, p_clip, acc, 0, 1.0, zi, np.zeros((1, 1)), z1, zi, z1, zi, z1, z1,
                          z1, np.zeros((1, NACCW)))


NACCW = 7   # profile-weighted sums per (reflection, scan) (see Harvester, weights)
ACCW = ['wA', 'wB', 'wQ', 'wW', 'wVbg', 'wM', 'wC']


@numba.njit(cache=True)
def _integrate_frame_core(idx, vals, nfast, adj, spots, s0, f0, hs, hf, bw, nbox, nann, S1s, S1f, S2s, S2f,
                          domega, scale, clipped, p_clip, acc, mode, ostep, ho, gint, prof, poff, wflat, woff, SB, SW,
                          SC, accw):
    """Accumulate one frame into acc[spots] (see ACC). idx ascending flat pixel indices, vals counts.
    mode & 1: 1D profiles. prof[poff[k]:] holds, for reflection k, the background-subtracted sums over the box
    projected on s (2 hs + 1), on f (2 hf + 1) and on omega (2 ho + 1 frames, by round(domega / ostep)).
    gint is the integral image of the good-pixel mask, for the background per box row / column.
    mode & 2: profile-weighted sums with the per-pixel weight S(s) F(f) W(omega) taken from wflat[woff[k]:]
    (same layout as prof). SB[k], SW[k], SC[k] = sums over the good box pixels of (S F)^2, S F, (S F)^3. Adds to
    accw[k]: A = sum w (x - bg), B = sum w^2, Q = sum w^2 x, Wt = sum w, Vbg = sum (W SW)^2 var(bg per pixel),
    M = sum w^2 bg, C = sum w^3 (x in monitor-scaled counts). With w the profile of the peak, the variance of the
    profile-fitted intensity A / B of a peak of intensity I is (M + I C + Vbg) / B^2.
    hs, hf, bw: per-reflection half sizes and annulus width.
    S1s.. are sums over the good box pixels of (s - s0), (f - f0) and their squares (for the background
    correction of the moments); domega[k] = omega_frame - omega_pred of spot k."""
    for kk in range(spots.shape[0]):
        k = spots[kk]
        hsk = hs[k]
        hfk = hf[k]
        r0 = s0[k] - hsk - bw[k]
        r1 = s0[k] + hsk + bw[k]
        c0 = f0[k] - hfk - bw[k]
        c1 = f0[k] + hfk + bw[k]
        lo = np.searchsorted(idx, max(r0, 0) * nfast)
        hi = np.searchsorted(idx, (r1 + 1) * nfast)
        # pass 1: annulus sum, box sums
        asum = 0.0
        bsum = 0.0
        bs = 0.0
        bf = 0.0
        bss = 0.0
        bff = 0.0
        badj = 0.0
        vmax = 0.0
        nnz = 0.0
        wa = 0.0
        wq = 0.0
        for p in range(lo, hi):
            r = idx[p] // nfast
            c = idx[p] - r * nfast
            if c < c0 or c > c1:
                continue
            nnz += 1.0
            v = vals[p]
            ds_ = r - s0[k]
            df_ = c - f0[k]
            if abs(ds_) <= hsk and abs(df_) <= hfk:
                if mode & 1:
                    prof[poff[k] + ds_ + hsk] += v * scale
                    prof[poff[k] + 2 * hsk + 1 + df_ + hfk] += v * scale
                if mode & 2:
                    wsf = wflat[woff[k] + ds_ + hsk] * wflat[woff[k] + 2 * hsk + 1 + df_ + hfk]
                    wa += wsf * v
                    wq += wsf * wsf * v
                bsum += v
                bs += v * ds_
                bf += v * df_
                bss += v * ds_ * ds_
                bff += v * df_ * df_
                if adj[idx[p]]:
                    badj += v
                if v > vmax:
                    vmax = v
            else:
                asum += v
        # background per pixel
        na = nann[k]
        if na <= 0:
            m = 0.0
            nused = 0.0
            nclip = 0.0
        else:
            m = asum / na
            nused = na
            nclip = 0.0
            if clipped:
                t = _poisson_upper(m, p_clip)
                csum = 0.0
                for p in range(lo, hi):
                    r = idx[p] // nfast
                    c = idx[p] - r * nfast
                    if c < c0 or c > c1:
                        continue
                    if abs(r - s0[k]) <= hsk and abs(c - f0[k]) <= hfk:
                        continue
                    if vals[p] > t:
                        csum += vals[p]
                        nclip += 1.0
                if nclip > 0:
                    nused = na - nclip
                    m = (asum - csum) / nused if nused > 0 else 0.0
        bg = m * nbox[k]
        net = bsum - bg
        sc2 = scale * scale
        if mode & 3:
            iw = int(np.floor(domega[kk] / ostep + 0.5)) + ho[k]
            nwk = 2 * ho[k] + 1
            if mode & 1:
                nrow = gint.shape[0] - 1
                ncol = gint.shape[1] - 1
                ca = max(f0[k] - hfk, 0)
                cb = min(f0[k] + hfk, ncol - 1)
                ra = max(s0[k] - hsk, 0)
                rb = min(s0[k] + hsk, nrow - 1)
                for r in range(ra, rb + 1):
                    if cb >= ca:
                        ng = gint[r + 1, cb + 1] - gint[r, cb + 1] - gint[r + 1, ca] + gint[r, ca]
                        prof[poff[k] + r - s0[k] + hsk] -= m * ng * scale
                for c in range(ca, cb + 1):
                    if rb >= ra:
                        ng = gint[rb + 1, c + 1] - gint[ra, c + 1] - gint[rb + 1, c] + gint[ra, c]
                        prof[poff[k] + 2 * hsk + 1 + c - f0[k] + hfk] -= m * ng * scale
                if 0 <= iw < nwk:
                    prof[poff[k] + 2 * hsk + 1 + 2 * hfk + 1 + iw] += net * scale
            if mode & 2:
                ww = wflat[woff[k] + 2 * hsk + 1 + 2 * hfk + 1 + iw] if 0 <= iw < nwk else 0.0
                vm = m / nused if nused > 0 else 0.0
                accw[k, 0] += ww * (wa - m * SW[k]) * scale
                accw[k, 1] += ww * ww * SB[k]
                accw[k, 2] += ww * ww * wq * sc2
                accw[k, 3] += ww * SW[k]
                accw[k, 4] += (ww * SW[k]) ** 2 * vm * sc2
                accw[k, 5] += ww * ww * SB[k] * m * sc2
                accw[k, 6] += ww * ww * ww * SC[k] * scale
        acc[k, 0] += bsum * scale
        acc[k, 1] += bg * scale
        acc[k, 2] += (bsum + (m / nused if nused > 0 else 0.0) * nbox[k] * nbox[k]) * sc2
        acc[k, 3] += 1.0
        acc[k, 4] += (bs - m * S1s[k]) * scale
        acc[k, 5] += (bf - m * S1f[k]) * scale
        acc[k, 6] += net * domega[kk] * scale
        acc[k, 7] += (bss - m * S2s[k]) * scale
        acc[k, 8] += (bff - m * S2f[k]) * scale
        acc[k, 9] += net * domega[kk] * domega[kk] * scale
        acc[k, 10] += nclip
        acc[k, 11] += badj * scale
        if vmax > acc[k, 12]:
            acc[k, 12] = vmax
        acc[k, 13] += nused
        acc[k, 14] += nnz


@numba.njit(cache=True)
def _count_overlaps(order, key, wf, s0, f0, hs, hf, ho, grain, dty, rad, ncf, wmax):
    """for reflections sorted by (grid cell, omega): count boxes of other grains overlapping in s, f, omega"""
    n = key.shape[0]
    out = np.zeros(n, np.int32)
    for i in range(n):
        for dc in (-ncf - 1, -ncf, -ncf + 1, -1, 0, 1, ncf - 1, ncf, ncf + 1):
            kk = key[i] + dc
            lo = np.searchsorted(key, kk)
            hi = np.searchsorted(key, kk + 1)
            if hi <= lo:
                continue
            a = lo + np.searchsorted(wf[lo:hi], wf[i] - wmax)
            b = lo + np.searchsorted(wf[lo:hi], wf[i] + wmax, side='right')
            for j in range(a, b):
                if j == i or grain[j] == grain[i]:
                    continue
                if abs(wf[j] - wf[i]) > ho[i] + ho[j]:
                    continue
                if abs(s0[j] - s0[i]) > hs[i] + hs[j] or abs(f0[j] - f0[i]) > hf[i] + hf[j]:
                    continue
                if abs(dty[j] - dty[i]) >= rad[i] + rad[j]:
                    continue
                out[i] += 1
    return out


@numba.njit(cache=True)
def _overlap_contamination(key, wf, s, f, s0, f0, hs, hf, ho, sig_s, sig_f, sig_w, grain, dty, rad, I, ncf, wmax):
    """for reflections sorted by (grid cell, omega): intensity of other grains' reflections expected inside each box,
    sum_j max(I_j, 0) * (fraction of j's Gaussian profile, widths sig_*, inside i's box)"""
    n = key.shape[0]
    out = np.zeros(n)
    nmax = np.zeros(n)
    for i in range(n):
        for dc in (-ncf - 1, -ncf, -ncf + 1, -1, 0, 1, ncf - 1, ncf, ncf + 1):
            kk = key[i] + dc
            lo = np.searchsorted(key, kk)
            hi = np.searchsorted(key, kk + 1)
            if hi <= lo:
                continue
            a = lo + np.searchsorted(wf[lo:hi], wf[i] - wmax)
            b = lo + np.searchsorted(wf[lo:hi], wf[i] + wmax, side='right')
            for j in range(a, b):
                if j == i or grain[j] == grain[i] or I[j] <= 0:
                    continue
                if abs(dty[j] - dty[i]) >= rad[i] + rad[j]:
                    continue
                fr = 1.0
                for c, w, sg, h in ((s[j] - s0[i], 0.0, sig_s[j], hs[i]), (f[j] - f0[i], 0.0, sig_f[j], hf[i]),
                                    (wf[j] - wf[i], 0.0, sig_w[j], ho[i])):
                    lo_ = (-h - 0.5 - c) / (np.sqrt(2.0) * sg)
                    hi_ = (h + 0.5 - c) / (np.sqrt(2.0) * sg)
                    fr *= 0.5 * (math.erf(hi_) - math.erf(lo_))
                if fr > 1e-4:
                    out[i] += fr * I[j]
                    if fr * I[j] > nmax[i]:
                        nmax[i] = fr * I[j]
    return out, nmax


@numba.njit(cache=True)
def _pair_overlaps(key, wf, s0, f0, hs, hf, ho, grain, dty, rad, ncf, wmax, ng):
    """for reflections sorted by (grid cell, omega): M[gi, gj] = number of reflections of grain gi whose box is
    overlapped by at least one box of grain gj while both grains can be in the beam"""
    n = key.shape[0]
    M = np.zeros((ng, ng), np.int32)
    last = np.full(ng, -1, np.int64)
    for i in range(n):
        for dc in (-ncf - 1, -ncf, -ncf + 1, -1, 0, 1, ncf - 1, ncf, ncf + 1):
            kk = key[i] + dc
            lo = np.searchsorted(key, kk)
            hi = np.searchsorted(key, kk + 1)
            if hi <= lo:
                continue
            a = lo + np.searchsorted(wf[lo:hi], wf[i] - wmax)
            b = lo + np.searchsorted(wf[lo:hi], wf[i] + wmax, side='right')
            for j in range(a, b):
                gj = grain[j]
                if gj == grain[i] or last[gj] == i:
                    continue
                if abs(wf[j] - wf[i]) > ho[i] + ho[j]:
                    continue
                if abs(s0[j] - s0[i]) > hs[i] + hs[j] or abs(f0[j] - f0[i]) > hf[i] + hf[j]:
                    continue
                if abs(dty[j] - dty[i]) >= rad[i] + rad[j]:
                    continue
                last[gj] = i
                M[grain[i], gj] += 1
    return M


@numba.njit(cache=True)
def _profile_in_box(i, j, s, f, wf, s0, f0, hs, hf, ho, sig_s, sig_f, sig_w):
    """fraction of reflection j's Gaussian profile (widths sig_*) inside the box of reflection i"""
    fr = 1.0
    for c, sg, h in ((s[j] - s0[i], sig_s[j], hs[i]), (f[j] - f0[i], sig_f[j], hf[i]),
                     (wf[j] - wf[i], sig_w[j], ho[i])):
        lo_ = (-h - 0.5 - c) / (np.sqrt(2.0) * sg)
        hi_ = (h + 0.5 - c) / (np.sqrt(2.0) * sg)
        fr *= 0.5 * (math.erf(hi_) - math.erf(lo_))
    return fr


@numba.njit(cache=True, parallel=True)
def _overlap_fractions(key, wf, s, f, s0, f0, hs, hf, ho, sig_s, sig_f, sig_w, grain, dty, rad, ncf, wmax, frmin,
                       ptr, col, frac):
    """for reflections sorted by (grid cell, omega): CSR matrix F[i, j] = fraction of reflection j's profile inside
    the box of i, other grains only, while both can be in the beam, F >= frmin. With col empty, ptr[i + 1] gets the
    number of entries of row i (call again with the cumulated ptr and col, frac allocated to fill them)."""
    n = key.shape[0]
    fill = col.shape[0] > 0
    for i in numba.prange(n):
        k = ptr[i] if fill else 0
        for dc in (-ncf - 1, -ncf, -ncf + 1, -1, 0, 1, ncf - 1, ncf, ncf + 1):
            kk = key[i] + dc
            lo = np.searchsorted(key, kk)
            hi = np.searchsorted(key, kk + 1)
            if hi <= lo:
                continue
            a = lo + np.searchsorted(wf[lo:hi], wf[i] - wmax)
            b = lo + np.searchsorted(wf[lo:hi], wf[i] + wmax, side='right')
            for j in range(a, b):
                if j == i or grain[j] == grain[i]:
                    continue
                if abs(dty[j] - dty[i]) >= rad[i] + rad[j]:
                    continue
                fr = _profile_in_box(i, j, s, f, wf, s0, f0, hs, hf, ho, sig_s, sig_f, sig_w)
                if fr < frmin:
                    continue
                if fill:
                    col[k] = j
                    frac[k] = fr
                k += 1
        if not fill:
            ptr[i + 1] = k


@numba.njit(cache=True)
def _average_profiles(nbr, hs, hf, ho, ds_, df_, off, data, I, out):
    """out[off[t]:] = mean over the neighbours r = nbr[t] (-1: none) of the profiles of r normalised by I[r],
    resampled on t's box grid relative to its predicted position; then clipped at 0 and normalised to sum 1 per
    dimension (flat if nothing was learned)"""
    for t in range(nbr.shape[0]):
        ns, nf, nw = 2 * hs[t] + 1, 2 * hf[t] + 1, 2 * ho[t] + 1
        o = off[t]
        nn = 0
        for q in range(nbr.shape[1]):
            r = nbr[t, q]
            if r < 0:
                continue
            nn += 1
            a = off[r]
            nsr, nfr, nwr = 2 * hs[r] + 1, 2 * hf[r] + 1, 2 * ho[r] + 1
            for dim in range(2):
                n_t = ns if dim == 0 else nf
                n_r = nsr if dim == 0 else nfr
                h_t = hs[t] if dim == 0 else hf[t]
                h_r = hs[r] if dim == 0 else hf[r]
                d_t = ds_[t] if dim == 0 else df_[t]
                d_r = ds_[r] if dim == 0 else df_[r]
                base_r = a if dim == 0 else a + nsr
                base_t = o if dim == 0 else o + ns
                for j in range(n_t):
                    x = j - h_t - d_t + d_r + h_r      # target pixel, relative to the prediction, on r's grid
                    i0 = int(np.floor(x))
                    w1 = x - i0
                    v = 0.0
                    if 0 <= i0 < n_r:
                        v += (1 - w1) * data[base_r + i0]
                    if 0 <= i0 + 1 < n_r:
                        v += w1 * data[base_r + i0 + 1]
                    out[base_t + j] += v / I[r]
            for j in range(nw):
                xw = j - ho[t] + ho[r]
                if 0 <= xw < nwr:
                    out[o + ns + nf + j] += data[a + nsr + nfr + xw] / I[r]
        if nn > 0:
            for j in range(ns + nf + nw):
                out[o + j] /= nn
        for b0, n in ((o, ns), (o + ns, nf), (o + ns + nf, nw)):
            tot = 0.0
            for j in range(n):
                if out[b0 + j] < 0:
                    out[b0 + j] = 0.0
                tot += out[b0 + j]
            for j in range(n):
                out[b0 + j] = out[b0 + j] / tot if tot > 0 else 1.0 / n


def _box_sums(integral, s0, f0, hs, hf):
    """sums of an image over the boxes [s0-hs, s0+hs] x [f0-hf, f0+hf] from its integral image"""
    H, W = integral.shape[0] - 1, integral.shape[1] - 1
    a0, a1 = np.clip(s0 - hs, 0, H), np.clip(s0 + hs + 1, 0, H)
    b0, b1 = np.clip(f0 - hf, 0, W), np.clip(f0 + hf + 1, 0, W)
    return integral[a1, b1] - integral[a0, b1] - integral[a1, b0] + integral[a0, b0]


def _integral(img):
    out = np.zeros((img.shape[0] + 1, img.shape[1] + 1), float)
    out[1:, 1:] = img.cumsum(0).cumsum(1)
    return out


def profile_slices(hs, hf, ho, offset):
    """(s, f, omega) slices of one reflection's profiles in the flat profile array (see Harvester profiles)"""
    a = int(offset)
    b = a + 2 * int(hs) + 1
    c = b + 2 * int(hf) + 1
    return slice(a, b), slice(b, c), slice(c, c + 2 * int(ho) + 1)


def learn_profiles(spots, I, sigI, profiles, offset, reference, k=15, min_ref=60):
    """Expected s, f, omega profile of every reflection, learned from strong ones (the Harvester profiles option).

    spots     : dict of arrays (grain, hs, hf, ho, s, f, s0, f0, tth, eta), one entry per reflection, grains sorted
    I, sigI   : merged intensities (from the run that collected the profiles)
    profiles, offset : the flat profiles and their offsets (Harvester.profile_data or the outfile group)
    reference : bool, the reflections to learn from (e.g. I / sigI > 20, not touching the mask or another grain)
    For each reflection: the mean of the normalised profiles of the k nearest references of the same grain in
    (2theta, eta) (omega widths go as 1 / |sin eta|), relative to the predicted position, resampled onto its own box,
    clipped at 0 and normalised to 1 per dimension. Grains with fewer than min_ref references use all grains.
    A reference does not learn from itself. Returns (weights, n_ref, pooled): weights in the profile layout, to
    use as Harvester(weights=(weights, offset))."""
    from scipy.spatial import cKDTree
    g = np.asarray(spots['grain'])
    feat = np.column_stack((spots['tth'] / 2.0, 3 * np.cos(np.radians(spots['eta'])),
                            3 * np.sin(np.radians(spots['eta']))))
    ref = np.asarray(reference, bool) & (I > 0)
    nbr = np.full((len(g), k), -1, np.int64)
    pooled = np.zeros(len(g), bool)
    iall = np.flatnonzero(ref)
    tall = cKDTree(feat[iall])
    start = np.searchsorted(g, np.arange(g.max() + 2))
    for gi in range(g.max() + 1):
        rows = np.arange(start[gi], start[gi + 1])
        if len(rows) == 0:
            continue
        r = rows[ref[rows]]
        if len(r) >= min_ref:
            cand = r[cKDTree(feat[r]).query(feat[rows], k + 1)[1]]
        else:
            cand = iall[tall.query(feat[rows], k + 1)[1]]
            pooled[rows] = True
        cand = np.where(cand == rows[:, None], -1, cand)
        o = np.argsort(cand < 0, axis=1, kind='stable')
        nbr[rows] = np.take_along_axis(cand, o, 1)[:, :k]
    out = np.zeros(int(offset[-1]))
    _average_profiles(nbr, np.asarray(spots['hs'], np.int64), np.asarray(spots['hf'], np.int64),
                      np.asarray(spots['ho'], np.int64), np.asarray(spots['s'] - spots['s0'], float),
                      np.asarray(spots['f'] - spots['f0'], float), np.asarray(offset, np.int64),
                      np.asarray(profiles, float), np.where(I != 0, I, 1.0).astype(float), out)
    return out, (nbr >= 0).sum(1), pooled


def misorientation(U1, U2, group):
    """smallest rotation angle (deg) between the orientations U1 and U2 (crystal -> lab, orthonormal) under the
    rotations of group (e.g. ImageD11.sym_u.getgroup('cubic')().group)"""
    d = np.dot(U1.T, U2)
    tr = max(np.trace(np.dot(np.asarray(op, float), d)) for op in group)
    return np.degrees(np.arccos(np.clip((tr - 1) / 2, -1, 1)))


def misorientation_axis(U1, U2, group):
    """(angle in deg, unit rotation axis in the crystal frame of U1) of the smallest rotation between the
    orientations U1 and U2 under the rotations of group; the axis is (0, 0, 0) for no rotation"""
    d = np.dot(U1.T, U2)
    R = max((np.dot(np.asarray(op, float), d) for op in group), key=np.trace)
    ang = np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))
    ax = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    n = np.linalg.norm(ax)
    if n < 1e-9:
        if ang < 1e-6:
            return 0.0, np.zeros(3)
        # 180 degrees: the axis is the eigenvector of R with eigenvalue 1
        w, v = np.linalg.eigh((R + R.T) / 2)
        ax, n = v[:, np.argmax(w)], 1.0
    return np.degrees(ang), ax / n


def twin_law_distance(U1, U2, group, angle, axis):
    """distance (deg) of the orientation U2 from U1 rotated by the twin law (angle in deg about the crystal axis,
    e.g. 60, (1, 1, 1)), over the symmetry-equivalent axes of the law and of U2: 0 for an ideal twin"""
    k = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    t = np.radians(angle)
    R = np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * np.dot(K, K)
    ops = [np.asarray(op, float) for op in group]
    return min(misorientation(np.dot(U1, np.dot(op, np.dot(R, op.T))), U2, group) for op in ops)


def _match(H_a, ubi_a, ubi_b, H_b, om_a, om_b):
    """for each reflection of grain a: the row of grain b with the same reflection (h_b = round(ubi_b inv(ubi_a))
    h_a, the omega solution closest to a's), or -1"""
    M = np.round(np.dot(ubi_b, np.linalg.inv(ubi_a))).astype(int)
    hb = np.dot(H_a, M.T)
    lut = {}
    for q, key in enumerate(map(tuple, H_b)):
        lut.setdefault(key, []).append(q)
    out = np.full(len(H_a), -1, np.int64)
    for p, key in enumerate(map(tuple, hb)):
        c = lut.get(key)
        if c:
            out[p] = min(c, key=lambda q: abs((om_b[q] - om_a[p] + 180) % 360 - 180))
    return out


def union_spots(spots, members, ubis, omega_step, max_sep=4.0):
    """One box per reflection for a group of duplicate grains (members: grain indices, the first is the main one):
    each reflection of the main grain is matched to the same reflection of every other member (Miller indices
    through round(ubi_b inv(ubi_a)), nearest omega solution, predicted positions within max_sep box half sizes)
    and its box becomes the bounding box (s, f, omega) of all the matched boxes, so every pixel is counted once.
    Returns (spots dict for the merged grain, {member: matched row of each main reflection, -1 if none})."""
    start = np.searchsorted(spots['grain'], np.arange(np.max(spots['grain']) + 2))
    blk = [dict((c, np.asarray(spots[c])[start[m]:start[m + 1]]) for c in spots) for m in members]
    A = blk[0]
    H = np.column_stack((A['h'], A['k'], A['l']))
    lo_s, hi_s = A['s0'] - A['hs'], A['s0'] + A['hs']
    lo_f, hi_f = A['f0'] - A['hf'], A['f0'] + A['hf']
    lo_w, hi_w = A['omega'] - A['ho'] * omega_step, A['omega'] + A['ho'] * omega_step
    bw = A['bw'].copy()
    match = {members[0]: np.arange(start[members[0]], start[members[0] + 1])}
    for m, B in zip(members[1:], blk[1:]):
        j = _match(H, ubis[members[0]], ubis[m], np.column_stack((B['h'], B['k'], B['l'])), A['omega'], B['omega'])
        a = np.flatnonzero(j >= 0)
        b = j[a]
        sep = np.max([np.abs(B['s'][b] - A['s'][a]) / A['hs'][a], np.abs(B['f'][b] - A['f'][a]) / A['hf'][a],
                      np.abs(B['omega'][b] - A['omega'][a]) / omega_step / A['ho'][a]], axis=0)
        a, b = a[sep < max_sep], b[sep < max_sep]
        lo_s[a] = np.minimum(lo_s[a], B['s0'][b] - B['hs'][b]); hi_s[a] = np.maximum(hi_s[a], B['s0'][b] + B['hs'][b])
        lo_f[a] = np.minimum(lo_f[a], B['f0'][b] - B['hf'][b]); hi_f[a] = np.maximum(hi_f[a], B['f0'][b] + B['hf'][b])
        lo_w[a] = np.minimum(lo_w[a], B['omega'][b] - B['ho'][b] * omega_step)
        hi_w[a] = np.maximum(hi_w[a], B['omega'][b] + B['ho'][b] * omega_step)
        bw[a] = np.maximum(bw[a], B['bw'][b])
        mm = np.full(len(H), -1, np.int64)
        mm[a] = start[m] + b
        match[m] = mm
    sp = dict((c, A[c].copy()) for c in ('h', 'k', 'l', 'ds', 'tth', 'eta', 's', 'f', 'sc', 'fc'))
    sp['s0'] = np.round((lo_s + hi_s) / 2.0).astype(np.int64)
    sp['f0'] = np.round((lo_f + hi_f) / 2.0).astype(np.int64)
    sp['hs'] = np.ceil(np.maximum(hi_s - sp['s0'], sp['s0'] - lo_s)).astype(np.int64)
    sp['hf'] = np.ceil(np.maximum(hi_f - sp['f0'], sp['f0'] - lo_f)).astype(np.int64)
    sp['omega'] = (lo_w + hi_w) / 2.0
    sp['ho'] = np.ceil((hi_w - lo_w) / 2.0 / omega_step - 1e-6).astype(np.int64)
    sp['bw'] = bw
    return sp, match


def duplicate_groups(spots, merged, ubis, group, omega_step, max_angle=2.0, min_shared=0.3, in_box=0.5,
                     to_partner=0.5, min_snr=20.0, pair_overlaps=None):
    """Groups of grains that integrate the same peaks (duplicate indexings of one crystal, or domains that the boxes
    cannot resolve). A pair is a duplicate when the misorientation (under the rotations of group) is below max_angle,
    they share boxes for more than min_shared of their reflections (pair_overlaps, from Harvester.pair_overlaps), and
    either the partner's predicted peak lies inside the box for more than in_box of the reflections, or the measured
    centroids (s, f, omega, in box half sizes) of strong reflections (I / sigI > min_snr) lie closer to the partner's
    prediction than to their own for more than to_partner of them. merged: dict with I, sigI, s, f, omega. Pairs are joined into groups (union-find). Returns a list of lists of grain
    indices, the main grain (most strong reflections) first, and the list of (gi, gj, angle, in_box, to_partner)."""
    ng = len(ubis)
    U = []
    for u in ubis:
        a, _, vt = np.linalg.svd(np.linalg.inv(u))
        U.append(np.dot(a, vt))
    F = pair_overlaps / np.maximum(np.bincount(spots['grain'], minlength=ng)[:, None], 1).astype(float)
    start = np.searchsorted(spots['grain'], np.arange(ng + 1))
    parent = list(range(ng))

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    pairs = []
    cand = np.argwhere(np.triu(np.maximum(F, F.T) > min_shared, 1))
    for i, j in cand:
        ang = misorientation(U[i], U[j], group)
        if ang >= max_angle:
            continue
        best = (0.0, 0.0)
        for a, b in ((i, j), (j, i)):
            ra = np.arange(start[a], start[a + 1])
            rb = np.arange(start[b], start[b + 1])
            H = np.column_stack((spots['h'][ra], spots['k'][ra], spots['l'][ra]))
            m = _match(H, ubis[a], ubis[b], np.column_stack((spots['h'][rb], spots['k'][rb], spots['l'][rb])),
                       spots['omega'][ra], spots['omega'][rb])
            ok = m >= 0
            pa, pb = ra[ok], rb[m[ok]]
            hs, hf, ho = spots['hs'][pa], spots['hf'][pa], spots['ho'][pa]
            inb = np.mean((np.abs(spots['s'][pb] - spots['s'][pa]) <= hs) & (np.abs(spots['f'][pb] - spots['f'][pa]) <= hf)
                          & (np.abs(spots['omega'][pb] - spots['omega'][pa]) <= ho * omega_step)) if ok.any() else 0.0
            st = merged['I'][pa] > min_snr * merged['sigI'][pa]
            wst = ho * omega_step
            own = np.sqrt(((merged['s'][pa] - spots['s'][pa]) / hs) ** 2 + ((merged['f'][pa] - spots['f'][pa]) / hf) ** 2
                          + ((merged['omega'][pa] - spots['omega'][pa]) / wst) ** 2)
            par = np.sqrt(((merged['s'][pa] - spots['s'][pb]) / hs) ** 2 + ((merged['f'][pa] - spots['f'][pb]) / hf) ** 2
                          + ((merged['omega'][pa] - spots['omega'][pb]) / wst) ** 2)
            sel = st & np.isfinite(own) & np.isfinite(par)
            tp = np.mean(par[sel] < own[sel]) if sel.sum() > 10 else 0.0
            best = (max(best[0], inb), max(best[1], tp))
        if best[0] > in_box or best[1] > to_partner:
            pairs.append((int(i), int(j), ang, best[0], best[1]))
            parent[find(i)] = find(j)
    groups = {}
    for i, j, _, _, _ in pairs:
        for x in (i, j):
            groups.setdefault(find(x), set()).add(x)
    nstrong = np.bincount(spots['grain'], weights=merged['I'] > min_snr * merged['sigI'], minlength=ng)
    out = [sorted(v, key=lambda x: -nstrong[x]) for v in groups.values()]
    return out, pairs


# ------------------------------------------------------------------------------------- the harvester
_WORKER = {}


def deoverlap_solve(F, I, sigI=None, ridge=1e-3, link=0.5, niter=50, tol=1e-3, verbose=False):
    """The measured I of a box is its own reflection plus the parts of its neighbours' peaks inside it:
    I = (1 + F) T, F[i, j] = fraction of reflection j's profile inside the box of i (Harvester.overlap_fractions).
    T >= 0 is the weighted (1 / sigI^2) least-squares solution with a ridge (ridge T / sigI): where two boxes
    hold each other's peaks (F[i, j] ~ F[j, i] ~ 1) the split of their sum is not measured and comes out even.
    Solved by blocks: reflections linked by contributions that can matter, F[i, j] (max(I[j], 0) + sigI[j]) >=
    link sigI[i] (either way), form a block (so boxes holding each other's peaks are always in one block), solved exactly (non-negative least squares) with the contributions from outside the block fixed at
    the current T; repeated until T changes by less than tol sigI (at most niter times).
    Returns a dict of arrays per reflection: contamination (F T), I_deoverlapped (T), largest (largest
    F[i, j] T[j]), largest_row (its j, -1: none), coincident_fraction (largest F[i, j], geometry alone),
    coincident_row, cluster (block label; -1: no overlap)."""
    import scipy.sparse, scipy.sparse.csgraph, scipy.optimize
    t0 = time.time()
    F = scipy.sparse.csr_matrix(F)
    I = np.asarray(I, float)
    n = len(I)
    sig = np.maximum(np.sqrt(np.maximum(np.abs(I), 1.0)) if sigI is None else np.asarray(sigI, float), 1.0)
    T = np.maximum(I, 0)
    Fc = F.tocoo()
    r_, c_, v_ = Fc.row, Fc.col, Fc.data
    strong = v_ * (T[c_] + sig[c_]) >= link * sig[r_]
    G = scipy.sparse.csr_matrix((np.ones(strong.sum()), (r_[strong], c_[strong])), shape=(n, n))
    G = (G + G.T).tocsr()
    _, lab = scipy.sparse.csgraph.connected_components(G, directed=False)
    nb = np.bincount(lab)
    single = nb[lab] == 1
    lab = np.where(single & (np.diff((F + F.T).tocsr().indptr) == 0), -1, lab)
    inside = lab[r_] == lab[c_]
    Fout = scipy.sparse.csr_matrix((v_[~inside], (r_[~inside], c_[~inside])), shape=(n, n))
    # blocks of more than one reflection: members (sorted by block) and the entries inside each block
    multi = np.flatnonzero(~single)
    mo = multi[np.argsort(lab[multi], kind='stable')]
    bstart = np.flatnonzero(np.r_[True, lab[mo][1:] != lab[mo][:-1], True])
    loc = np.zeros(n, np.int64)
    loc[mo] = np.arange(len(mo)) - np.repeat(bstart[:-1], np.diff(bstart))
    bid = np.full(n, -1)
    bid[mo] = np.repeat(np.arange(len(bstart) - 1), np.diff(bstart))
    ein = np.flatnonzero(inside & (bid[r_] >= 0))
    ein = ein[np.argsort(bid[r_[ein]], kind='stable')]
    estart = np.searchsorted(bid[r_[ein]], np.arange(len(bstart)))
    blocks = []
    for k in range(len(bstart) - 1):
        m = mo[bstart[k]:bstart[k + 1]]
        e = ein[estart[k]:estart[k + 1]]
        A = np.eye(len(m))
        A[loc[r_[e]], loc[c_[e]]] += v_[e]
        w = 1.0 / sig[m]
        blocks.append((m, np.vstack([w[:, None] * A, ridge * np.diag(w)]), w))
    if verbose:
        print('deoverlap: %d single, %d blocks (largest %d), %d links outside blocks (%.0f s)' % (
            single.sum(), len(blocks), np.diff(bstart).max() if len(blocks) else 0, Fout.nnz, time.time() - t0))
    s1 = np.flatnonzero(single & (lab >= 0))        # no overlap at all (lab -1): T = max(I, 0)
    for it in range(niter):
        r = I - Fout.dot(T)
        new = T.copy()
        new[s1] = np.maximum(r[s1], 0) / (1 + ridge ** 2)
        for m, Aw, w in blocks:
            new[m] = scipy.optimize.nnls(Aw, np.concatenate([w * r[m], np.zeros(len(m))]))[0]
        change = np.max(np.abs(new - T) / sig)
        T = new
        if verbose:
            print('deoverlap: iteration %d, max change %.3g sigma (%.0f s)' % (it, change, time.time() - t0))
        if change < tol:
            break

    return _deoverlap_result(F, T, cluster=lab)


def _rowmax(M):
    n = M.shape[0]
    mx, arg = np.zeros(n), np.full(n, -1)
    r_ = np.repeat(np.arange(n), np.diff(M.indptr))
    k = np.lexsort((-M.data, r_))
    if len(k):
        first = k[np.r_[True, r_[k][1:] != r_[k][:-1]]]
        mx[r_[first]], arg[r_[first]] = M.data[first], M.indices[first]
    return mx, arg


def _deoverlap_result(F, T, **extra):
    largest, jmax = _rowmax(F.multiply(T[None, :]).tocsr())
    frmax, jfr = _rowmax(F)
    d = dict(contamination=F.dot(T), I_deoverlapped=T, largest=largest, largest_row=np.where(largest > 0, jmax, -1),
             coincident_fraction=frmax, coincident_row=jfr)
    d.update(extra)
    return d


def deoverlap_jacobi(F, I, niter=8, damping=0.5, verbose=False):
    """I = (1 + F) T solved by the damped fixed point T <- damping T + (1 - damping) max(I - F T, 0) from
    T = max(I, 0), niter times (not iterated to convergence: where boxes hold each other's peaks it oscillates).
    Same F and returned dict as deoverlap_solve (no cluster); contamination (and largest) is F T with the T of
    the last update's input."""
    import scipy.sparse
    F = scipy.sparse.csr_matrix(F)
    I = np.asarray(I, float)
    T = prev = np.maximum(I, 0)
    for it in range(niter):
        prev = T
        T = damping * prev + (1 - damping) * np.maximum(I - F.dot(prev), 0)
        if verbose:
            print('deoverlap iteration %d: max change %.3g' % (it, np.max(np.abs(T - prev))))
    d = _deoverlap_result(F, prev)
    d['I_deoverlapped'] = T
    return d


def _work_scan(iscan):
    try:
        return _WORKER['h']._harvest_scan(iscan)
    except Exception as e:          # e.g. an image file that is unreadable or still being restored
        return iscan, None, '%s: %r' % (_WORKER['h'].scan_file(iscan), e), 0, 0.0, None


class Harvester(object):
    """Box integration of the predicted reflections of many grains in one pass over the frames.

    ds        : ImageD11.sinograms.dataset.DataSet (scans, omega, dty, imagefiles, limapath, monitor,
                maskfile, e2dxfile, e2dyfile)
    grains    : list of objects with .ubi and .translation (e.g. ImageD11.grain.grain)
    pars      : ImageD11 parameters dict (or ImageD11.parameters.parameters)
    ds_max    : maximum |g| (1/d, 1/Angstrom) of the predicted reflections
    radius    : per-grain in-beam half width in dty units (array or scalar); a grain is in the beam on a
                frame when |dty - dty_pred(omega)| <= radius. None: 3 ystep.
    mask      : detector mask (True / 1 = masked); None reads ds.maskfile (ImageD11 convention 1 = masked)
    datapath  : folder holding ds.imagefiles (default ds.datapath), e.g. a copy of the raw data
    reader    : 'auto', 'sparse' (bslz4_to_sparse) or 'h5py'
    monitor   : normalise each frame by ds.monitor_ref / ds.monitor (if the dataset has a monitor)
    box       : Box (see auto_box)
    profiles  : also collect, for every reflection, its background-subtracted 1D profiles on s, f and omega summed
                over scans (the dty profile is the per-scan rows). Written to group 'profiles' of the outfile:
                data[offset[k]:offset[k + 1]] = s profile (2 hs + 1), f profile (2 hf + 1), omega profile
                (2 ho + 1 frames, relative to the predicted omega). See profile_slices.
    weights   : None, or (flat, offset) per-pixel profile weights in the same layout as the profiles (S(s), F(f),
                W(omega) for each reflection). Adds the profile-weighted sums of ACCW to each (reflection, scan)
                row: the profile-fitted intensity of the row is wA / wB, its variance (wM + I wC + wVbg) / wB^2 for
                a peak of intensity I (use a model I, not the noisy row: count-based weights bias weak means).
    """

    def __init__(self, ds, grains, pars, ds_max, radius=None, mask=None, datapath=None, reader='auto',
                 monitor=True, box=None, e2dx=None, e2dy=None, profiles=False, weights=None):
        self.ds = ds
        self.grains = list(grains)
        if hasattr(pars, 'get_parameters'):
            pars = pars.get_parameters()
        self.pars = dict(pars)
        self.ds_max = float(ds_max)
        ng = len(self.grains)
        if radius is None:
            radius = 3 * ds.ystep
        self.radius = np.broadcast_to(np.asarray(radius, float), (ng,)).copy()
        self.datapath = ds.datapath if datapath is None else datapath
        self.reader = reader
        self.box = Box() if box is None else box
        if mask is None:
            if getattr(ds, 'maskfile', None):
                import fabio
                mask = fabio.open(ds.maskfile).data != 0
            else:
                mask = None
        self.mask = None if mask is None else np.asarray(mask).astype(bool)
        if e2dx is None and getattr(ds, 'e2dxfile', None):
            import fabio
            e2dx = fabio.open(ds.e2dxfile).data.astype(float)
            e2dy = fabio.open(ds.e2dyfile).data.astype(float)
        self.e2dx, self.e2dy = e2dx, e2dy
        self.scale = np.ones(ds.shape)
        if monitor and getattr(ds, 'monitor', None) is not None:
            ref = ds.monitor_ref if ds.monitor_ref is not None else np.median(ds.monitor)
            self.scale = float(ref) / np.asarray(ds.monitor, float)
        self.spots = None
        self.results = None
        self.profiles = bool(profiles)
        self.weights = weights

    # ---------------------------------------------------------------- files
    def scan_file(self, iscan):
        return os.path.join(self.datapath, self.ds.imagefiles[iscan])

    def present_scans(self):
        """scan indices whose image file exists"""
        return [i for i in range(len(self.ds.scans)) if os.path.exists(self.scan_file(i))]

    def detector_shape(self):
        if self.mask is not None:
            return self.mask.shape
        for i in self.present_scans():
            with h5py.File(self.scan_file(i), 'r') as h:
                return h[self.ds.limapath].shape[1:]
        raise RuntimeError('no image files found')

    # ---------------------------------------------------------------- prediction
    def predict(self):
        """Predict the reflections of all grains (sets self.spots, a dict of arrays, one entry per reflection)"""
        shape = self.detector_shape()
        orange = (np.min(self.ds.omega) - 1, np.max(self.ds.omega) + 1)
        cols = {}
        for gid, g in enumerate(self.grains):
            p = predict_spots(g.ubi, g.translation, self.pars, self.ds_max, shape, self.e2dx, self.e2dy, orange)
            p['grain'] = np.full(len(p['h']), gid, np.int32)
            for k, v in p.items():
                cols.setdefault(k, []).append(v)
        sp = dict((k, np.concatenate(v)) for k, v in cols.items())
        sp['s0'] = np.round(sp['s']).astype(np.int64)
        sp['f0'] = np.round(sp['f']).astype(np.int64)
        self.spots = sp
        sp['hs'], sp['hf'], sp['ho'], sp['bw'] = self.box.sizes(sp, self.pars, self.omega_step())
        self._static()
        self._overlaps()
        return sp

    def profile_offsets(self):
        """offsets of each reflection's [s | f | omega] profile in the flat profile array (length n + 1)"""
        sp = self.spots
        n = (2 * sp['hs'] + 1) + (2 * sp['hf'] + 1) + (2 * sp['ho'] + 1)
        return np.concatenate([[0], np.cumsum(n)]).astype(np.int64)

    def _weight_sums(self):
        """SB, SW, SC = sums over the good box pixels of (S F)^2, S F and (S F)^3, per reflection, for the weights"""
        sp = self.spots
        flat, off = self.weights
        good = self._good
        H, W = good.shape
        n = len(sp['hs'])
        SB, SW, SC = np.zeros(n), np.zeros(n), np.zeros(n)
        for k in range(n):
            hs, hf, s0, f0 = int(sp['hs'][k]), int(sp['hf'][k]), int(sp['s0'][k]), int(sp['f0'][k])
            S = flat[off[k]:off[k] + 2 * hs + 1]
            F = flat[off[k] + 2 * hs + 1:off[k] + 2 * hs + 2 * hf + 2]
            r0, c0 = s0 - hs, f0 - hf
            a0, a1, b0, b1 = max(r0, 0), min(s0 + hs + 1, H), max(c0, 0), min(f0 + hf + 1, W)
            if a1 <= a0 or b1 <= b0:
                continue
            g = good[a0:a1, b0:b1]
            w = np.outer(S[a0 - r0:a1 - r0], F[b0 - c0:b1 - c0]) * g
            SB[k], SW[k], SC[k] = (w * w).sum(), w.sum(), (w * w * w).sum()
        return SB, SW, SC

    def omega_step(self):
        om = np.asarray(self.ds.omega[0], float)
        return float(np.median(np.abs(np.diff(om)))) if om.size > 1 else 1.0

    def _overlaps(self):
        """n_overlap: number of other predicted reflections whose box overlaps this one (in s, f and omega)
        and whose grain can be in the beam in the same scan (|dty_pred difference| < sum of radii)"""
        sp = self.spots
        step = self.omega_step()
        T = np.array([g.translation for g in self.grains], float)
        g = sp['grain']
        w = np.radians(sp['omega'])
        dty = self.ds.y0 - T[g, 0] * np.sin(w) - T[g, 1] * np.cos(w)
        cell = int(2 * max(sp['hs'].max(), sp['hf'].max()) + 1)
        cs, cf = sp['s0'] // cell, sp['f0'] // cell
        ncf = int(cf.max()) + 3
        key = (cs + 1) * ncf + (cf + 1)
        wf = sp['omega'] / step
        order = np.lexsort((wf, key))
        sp['n_overlap'] = _count_overlaps(order, key[order], wf[order], sp['s0'][order], sp['f0'][order], sp['hs'][order],
                                          sp['hf'][order], sp['ho'][order], g[order], dty[order],
                                          self.radius[g][order], ncf, int(2 * sp['ho'].max()))[np.argsort(order)]

    def _static(self):
        """Per-reflection box geometry from the mask: good pixels in box and annulus, coordinate sums for
        the background correction of the moments, masked pixels, edge flag."""
        sp = self.spots
        hs, hf, bw = sp['hs'], sp['hf'], sp['bw']
        shape = self.detector_shape()
        good = np.ones(shape, bool) if self.mask is None else ~self.mask
        ss, ff = np.mgrid[0:shape[0], 0:shape[1]].astype(float)
        s0, f0 = sp['s0'], sp['f0']
        I_g = _integral(good.astype(float))
        n_box = _box_sums(I_g, s0, f0, hs, hf)
        n_out = _box_sums(I_g, s0, f0, hs + bw, hf + bw)
        Ss = _box_sums(_integral(good * ss), s0, f0, hs, hf)
        Sf = _box_sums(_integral(good * ff), s0, f0, hs, hf)
        Sss = _box_sums(_integral(good * ss * ss), s0, f0, hs, hf)
        Sff = _box_sums(_integral(good * ff * ff), s0, f0, hs, hf)
        sp['nbox'] = n_box
        sp['nann'] = n_out - n_box
        sp['S1s'] = Ss - s0 * n_box
        sp['S1f'] = Sf - f0 * n_box
        sp['S2s'] = Sss - 2 * s0 * Ss + s0 * s0 * n_box
        sp['S2f'] = Sff - 2 * f0 * Sf + f0 * f0 * n_box
        full = (2 * hs + 1) * (2 * hf + 1)
        inside = ((s0 - hs >= 0) & (s0 + hs < shape[0]) & (f0 - hf >= 0) & (f0 + hf < shape[1]))
        sp['edge'] = ~inside
        sp['n_masked'] = np.where(inside, full - n_box, 0).astype(np.int32)
        # box including a one-pixel rim: "in contact" with a masked pixel
        n_rim = _box_sums(I_g, s0, f0, hs + 1, hf + 1)
        sp['touch_mask'] = inside & (n_rim < (2 * hs + 3) * (2 * hf + 3))
        if self.mask is not None:
            import scipy.ndimage
            self._adj = scipy.ndimage.binary_dilation(self.mask, np.ones((3, 3), bool)).ravel().astype(np.uint8)
        else:
            self._adj = np.zeros(int(np.prod(shape)), np.uint8)
        self._good = good

    def overlap_contamination(self, I, sigmas=None):
        """Intensity from other grains' reflections expected inside each box, from the harvested intensities I
        (one per reflection, e.g. merged()['I']): sum over overlapping reflections j of max(I_j, 0) times the
        fraction of j's Gaussian profile inside this box (widths from the ShoeboxModel, or box / (2 * 3)).
        Returns (contamination, largest single contribution); compare with I and sigI."""
        sp = self.spots
        step = self.omega_step()
        if sigmas is None:
            if hasattr(self.box, 'sigmas'):
                sigmas = self.box.sigmas(sp, self.pars, step)
            else:
                sigmas = (sp['hs'] / 3.0, sp['hf'] / 3.0, sp['ho'] / 3.0)
        sig_s, sig_f, sig_w = [np.maximum(np.asarray(x, float), 0.3) for x in sigmas]
        T = np.array([g.translation for g in self.grains], float)
        g = sp['grain']
        w = np.radians(sp['omega'])
        dty = self.ds.y0 - T[g, 0] * np.sin(w) - T[g, 1] * np.cos(w)
        cell = int(2 * max(sp['hs'].max(), sp['hf'].max()) + 1)
        cs, cf = sp['s0'] // cell, sp['f0'] // cell
        ncf = int(cf.max()) + 3
        key = (cs + 1) * ncf + (cf + 1)
        wf = sp['omega'] / step
        o = np.lexsort((wf, key))
        c, cmax = _overlap_contamination(key[o], wf[o], sp['s'][o], sp['f'][o], sp['s0'][o], sp['f0'][o], sp['hs'][o],
                                         sp['hf'][o], sp['ho'][o], sig_s[o], sig_f[o], sig_w[o], g[o], dty[o],
                                         self.radius[g][o], np.asarray(I, float)[o], ncf, int(2 * sp['ho'].max()))
        inv = np.argsort(o)
        return c[inv], cmax[inv]

    def _grid_order(self):
        """(key, omega in frames, sort order, ncf, dty of the grain at each reflection) for the neighbour kernels"""
        sp = self.spots
        T = np.array([g.translation for g in self.grains], float)
        g = sp['grain']
        w = np.radians(sp['omega'])
        dty = self.ds.y0 - T[g, 0] * np.sin(w) - T[g, 1] * np.cos(w)
        cell = int(2 * max(sp['hs'].max(), sp['hf'].max()) + 1)
        cs, cf = sp['s0'] // cell, sp['f0'] // cell
        ncf = int(cf.max()) + 3
        key = (cs + 1) * ncf + (cf + 1)
        wf = sp['omega'] / self.omega_step()
        return key, wf, np.lexsort((wf, key)), ncf, dty

    def pair_overlaps(self):
        """(ngrains, ngrains): M[i, j] = number of reflections of grain i whose box overlaps a box of grain j (in s,
        f, omega, while both can be in the beam). M[i, j] / (reflections of i) is the shared fraction."""
        sp = self.spots
        key, wf, o, ncf, dty = self._grid_order()
        g = sp['grain']
        return _pair_overlaps(key[o], wf[o], sp['s0'][o], sp['f0'][o], sp['hs'][o], sp['hf'][o], sp['ho'][o], g[o],
                              dty[o], self.radius[g][o], ncf, int(2 * sp['ho'].max()), len(self.grains))

    def overlap_fractions(self, sigmas=None, frmin=1e-4):
        """sparse (nrefl, nrefl) CSR matrix F[i, j] = fraction of reflection j's profile (Gaussian, widths sigmas
        = (s, f, omega in frames) per reflection; default ShoeboxModel.sigmas, else box half widths / 3) inside the
        box of reflection i, for reflections of other grains while both can be in the beam (entries >= frmin)"""
        import scipy.sparse
        sp = self.spots
        if sigmas is None:
            if hasattr(self.box, 'sigmas'):
                sigmas = self.box.sigmas(sp, self.pars, self.omega_step())
            else:
                sigmas = (sp['hs'] / 3.0, sp['hf'] / 3.0, sp['ho'] / 3.0)
        sig_s, sig_f, sig_w = [np.maximum(np.asarray(x, float), 0.3) for x in sigmas]
        key, wf, o, ncf, dty = self._grid_order()
        g = sp['grain']
        n = len(o)
        args = (key[o], wf[o], sp['s'][o], sp['f'][o], sp['s0'][o], sp['f0'][o], sp['hs'][o], sp['hf'][o],
                sp['ho'][o], sig_s[o], sig_f[o], sig_w[o], g[o], dty[o], self.radius[g][o], ncf,
                int(2 * sp['ho'].max()), float(frmin))
        ptr = np.zeros(n + 1, np.int64)
        _overlap_fractions(*(args + (ptr, np.zeros(0, np.int64), np.zeros(0))))
        ptr = np.cumsum(ptr)
        col, frac = np.zeros(ptr[-1], np.int64), np.zeros(ptr[-1])
        _overlap_fractions(*(args + (ptr, col, frac)))
        rows = np.repeat(o, np.diff(ptr))     # back to the spot order
        return scipy.sparse.csr_matrix((frac, (rows, o[col])), shape=(n, n))

    def deoverlap(self, I, sigI=None, sigmas=None, frmin=1e-4, F=None, method='solve', verbose=False, **kw):
        """Contamination of every box by other grains' peaks, with the intensities corrected for it, on the
        overlap_fractions matrix F: method 'solve' (deoverlap_solve) or 'jacobi' (deoverlap_jacobi); kw go to
        that function. Returns its dict."""
        t0 = time.time()
        if F is None:
            F = self.overlap_fractions(sigmas, frmin)
        if verbose:
            print('deoverlap: %d overlaps (%.0f s)' % (F.nnz, time.time() - t0))
        if method == 'jacobi':
            return deoverlap_jacobi(F, I, verbose=verbose, **kw)
        if method == 'solve':
            return deoverlap_solve(F, I, sigI, verbose=verbose, **kw)
        raise ValueError('deoverlap method %r' % method)

    # ---------------------------------------------------------------- planning
    def plan_scan(self, iscan):
        """(frame indices, CSR pointers, reflection ids, omega offsets) for one scan: which reflections
        are integrated on which frames"""
        ds, sp = self.ds, self.spots
        om_f = np.asarray(ds.omega[iscan], float)
        dty = float(np.mean(ds.dty[iscan]))
        g = sp['grain']
        T = np.array([gr.translation for gr in self.grains], float)
        w = np.radians(sp['omega'])
        dty_pred = ds.y0 - T[g, 0] * np.sin(w) - T[g, 1] * np.cos(w)
        cand = np.flatnonzero(np.abs(dty - dty_pred) <= self.radius[g])
        if len(cand) == 0:
            return None
        order = np.argsort(om_f)
        oms = om_f[order]
        step = np.median(np.abs(np.diff(oms))) if len(oms) > 1 else 1.0
        half = (sp['ho'][cand] + 0.5) * step
        lo = np.searchsorted(oms, sp['omega'][cand] - half)
        hi = np.searchsorted(oms, sp['omega'][cand] + half)
        n = hi - lo
        ok = n > 0
        cand, lo, n = cand[ok], lo[ok], n[ok]
        if len(cand) == 0:
            return None
        rep = np.repeat(np.arange(len(cand)), n)
        pos = np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n) + lo[rep]
        frame = order[pos]
        spot = cand[rep]
        dom = om_f[frame] - sp['omega'][spot]
        o = np.argsort(frame, kind='stable')
        frame, spot, dom = frame[o], spot[o], dom[o]
        frames, ptr = np.unique(frame, return_index=True)
        ptr = np.append(ptr, len(frame))
        return frames, ptr, spot, dom

    # ---------------------------------------------------------------- one scan
    def _harvest_scan(self, iscan):
        plan = self.plan_scan(iscan)
        if plan is None:
            return iscan, None, None, 0, 0.0, None
        frames, ptr, spot, dom = plan
        sp, b = self.spots, self.box
        ids = np.unique(spot)
        local = np.searchsorted(ids, spot).astype(np.int64)
        acc = np.zeros((len(ids), NACC))
        arr = dict((k, np.ascontiguousarray(sp[k][ids])) for k in ('s0', 'f0', 'hs', 'hf', 'bw', 'nbox', 'nann', 'S1s', 'S1f',
                                                                      'S2s', 'S2f'))
        t0 = time.time()
        mode = (1 if self.profiles else 0) | (2 if self.weights is not None else 0)
        ho = np.ascontiguousarray(sp['ho'][ids])
        nprof = (2 * arr['hs'] + 1) + (2 * arr['hf'] + 1) + (2 * ho + 1)
        poff = np.concatenate([[0], np.cumsum(nprof)[:-1]]).astype(np.int64)
        prof = np.zeros(int(nprof.sum()) if mode & 1 else 1)
        if mode & 2:
            wflat, woff = self.weights[0], np.ascontiguousarray(self.weights[1][ids]).astype(np.int64)
            SB, SW = np.ascontiguousarray(self._SB[ids]), np.ascontiguousarray(self._SW[ids])
            SC = np.ascontiguousarray(self._SC[ids])
        else:
            wflat, woff, SB, SW, SC = np.zeros(1), np.zeros(1, np.int64), np.zeros(1), np.zeros(1), np.zeros(1)
        accw = np.zeros((len(ids) if mode & 2 else 1, NACCW))
        gint = self._gint if mode & 1 else np.zeros((1, 1))
        rd = open_reader(self.scan_file(iscan), self.ds.limapath, self._good, self.reader)
        nfast = self._good.shape[1]
        step = self.omega_step()
        try:
            for j, fr in enumerate(frames):
                idx, vals = rd(int(fr))
                a, z = ptr[j], ptr[j + 1]
                _integrate_frame_core(idx, vals, nfast, self._adj, local[a:z], arr['s0'], arr['f0'], arr['hs'],
                                      arr['hf'], arr['bw'], arr['nbox'], arr['nann'], arr['S1s'], arr['S1f'],
                                      arr['S2s'], arr['S2f'], dom[a:z], float(self.scale[iscan, fr]),
                                      b.bg_estimator == 'clipped', b.p_clip, acc, mode, step, ho, gint, prof, poff,
                                      wflat, woff, SB, SW, SC, accw)
        finally:
            rd.close()
        extra = {}
        if mode & 1:
            extra['prof'] = prof
        if mode & 2:
            extra['accw'] = accw
        return iscan, ids, acc, len(frames), time.time() - t0, extra

    # ---------------------------------------------------------------- run
    def run(self, scans=None, nproc=None, outfile=None, verbose=True):
        """Integrate all predicted reflections over the scans (default: all whose files exist).

        Per-reflection totals are accumulated in memory (see merged()). With outfile, the per
        (reflection, scan) rows are written to an HDF5 file as they arrive (group 'per_scan', with the
        reflection table in 'spots'), for sinogram / tomographic use. Scans whose file cannot be read are
        skipped and listed in self.failed_scans (and 'failed_scans' in outfile)."""
        if self.spots is None:
            self.predict()
        if scans is None:
            scans = self.present_scans()
        nsp = len(self.spots['h'])
        tot = np.zeros((nsp, NACC))
        nscans = np.zeros(nsp, np.int32)
        if self.profiles:
            gpoff = self.profile_offsets()
            gprof = np.zeros(int(gpoff[-1]))
            good = self._good.astype(float)
            self._gint = np.zeros((good.shape[0] + 1, good.shape[1] + 1))
            self._gint[1:, 1:] = good.cumsum(0).cumsum(1)
        if self.weights is not None:
            self._SB, self._SW, self._SC = self._weight_sums()
        _WORKER['h'] = self
        h5 = None
        if outfile is not None:
            h5 = h5py.File(outfile, 'w')
            gsp = h5.create_group('spots')
            for k, v in self.spots.items():
                gsp.create_dataset(k, data=v)
            gps = h5.create_group('per_scan')
            dsets = dict(spot=gps.create_dataset('spot', (0,), np.int32, maxshape=(None,), chunks=(65536,)),
                         scan=gps.create_dataset('scan', (0,), np.int32, maxshape=(None,), chunks=(65536,)))
            for k in ACC + (ACCW if self.weights is not None else []):
                dsets[k] = gps.create_dataset(k, (0,), np.float32, maxshape=(None,), chunks=(65536,))
            h5.attrs['box'] = repr(self.box)
            h5.attrs['box_omega_step'] = self.omega_step()
            for k, v in vars(self.box).items():
                if v is not None:
                    h5.attrs['box_' + k] = v
            h5.attrs['ds_max'] = self.ds_max
        t0 = time.time()
        nframes = 0
        failed = []
        if nproc is None:
            nproc = max(1, len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else
                        multiprocessing.cpu_count())
        if nproc > 1:
            # fork: the workers share the Harvester (python 2 has no get_context and forks)
            mp = multiprocessing.get_context('fork') if hasattr(multiprocessing, 'get_context') else multiprocessing
            pool = mp.Pool(nproc)
            it = pool.imap_unordered(_work_scan, scans)
        else:
            pool = None
            it = map(_work_scan, scans)
        try:
            for k, (iscan, ids, acc, nf, dt, extra) in enumerate(it):
                if ids is None:
                    if acc is not None:
                        failed.append((iscan, acc))
                        if verbose:
                            print('skipped scan %d: %s' % (iscan, acc))
                            sys.stdout.flush()
                    continue
                nframes += nf
                if self.profiles:
                    n_ = np.diff(gpoff)[ids]
                    src = np.arange(int(n_.sum()))
                    dst = np.repeat(gpoff[ids], n_) + src - np.repeat(np.cumsum(n_) - n_, n_)
                    gprof[dst] += extra['prof']
                for c in range(NACC):
                    if ACC[c] == 'maxval':
                        tot[ids, c] = np.maximum(tot[ids, c], acc[:, c])
                    else:
                        tot[ids, c] += acc[:, c]
                nscans[ids] += 1
                if h5 is not None:
                    n0 = dsets['spot'].shape[0]
                    for name in dsets:
                        dsets[name].resize((n0 + len(ids),))
                    dsets['spot'][n0:] = ids
                    dsets['scan'][n0:] = iscan
                    for c, name in enumerate(ACC):
                        dsets[name][n0:] = acc[:, c]
                    if self.weights is not None:
                        for c, name in enumerate(ACCW):
                            dsets[name][n0:] = extra['accw'][:, c]
                if verbose and (k % 20 == 0 or k == len(scans) - 1):
                    el = time.time() - t0
                    print('%d/%d scans, %d frames, %.0f s, %.0f frames/s' % (k + 1, len(scans), nframes, el, nframes / max(el, 1e-9)))
                    sys.stdout.flush()
        finally:
            if pool is not None:
                pool.close()
                pool.join()
        self.results = dict(total=tot, nscans=nscans)
        self.failed_scans = failed
        if failed and verbose:
            print('%d scans could not be read and were skipped: %s' % (len(failed), [f[0] for f in failed]))
        m = self.merged()
        if self.profiles:
            self.profile_data = (gprof, gpoff)
        if h5 is not None:
            if self.profiles:
                gp = h5.create_group('profiles')
                gp.create_dataset('offset', data=gpoff)
                gp.create_dataset('data', data=gprof.astype(np.float32), chunks=True, compression='gzip')
                gp.attrs['layout'] = ('data[offset[k]:offset[k+1]] = [s profile (2 hs + 1) | f profile (2 hf + 1) | '
                                      'omega profile (2 ho + 1 frames, centred on the predicted omega)], net counts '
                                      'summed over scans; dty profile = per_scan rows')
            h5.create_dataset('failed_scans', data=np.array([f[0] for f in failed], np.int32))
            gm = h5.create_group('merged')
            for k, v in m.items():
                gm.create_dataset(k, data=v)
            h5.close()
        return m

    # ---------------------------------------------------------------- results
    def merged(self):
        """Per-reflection results (dict of arrays): I (net), sigI, bg, box, nscans, nframes, centroid
        (s, f, omega) and its distortion-corrected (sc, fc), rms widths (sig_s, sig_f, sig_omega), the
        prediction, the flags n_masked, touch_mask, edge, I_adj_mask, nclip, and nnz (non-zero pixels read in
        the box and annulus over all frames: the size of a sparse copy of the spots)."""
        sp, t = self.spots, self.results['total']
        out = dict((k, sp[k]) for k in ('grain', 'h', 'k', 'l', 'ds', 'tth', 'eta', 'omega', 's', 'f', 'sc', 'fc',
                                        'n_masked', 'touch_mask', 'edge', 'n_overlap'))
        out = dict(('pred_' + k if k in ('tth', 'eta', 'omega', 's', 'f', 'sc', 'fc') else k, v) for k, v in out.items())
        I = t[:, _I['box']] - t[:, _I['bg']]
        out['I'] = I
        out['sigI'] = np.sqrt(np.maximum(t[:, _I['var']], 0))
        out['box'] = t[:, _I['box']]
        out['bg'] = t[:, _I['bg']]
        out['nframes'] = t[:, _I['nframes']].astype(np.int32)
        out['nscans'] = self.results['nscans']
        out['nclip'] = t[:, _I['nclip']]
        out['I_adj_mask'] = t[:, _I['I_adj_mask']]
        out['maxval'] = t[:, _I['maxval']]
        out['nnz'] = t[:, _I['nnz']]
        with np.errstate(divide='ignore', invalid='ignore'):
            ds_ = t[:, _I['ms']] / I
            df_ = t[:, _I['mf']] / I
            dw_ = t[:, _I['mw']] / I
            out['s'] = sp['s0'] + ds_
            out['f'] = sp['f0'] + df_
            out['omega'] = sp['omega'] + dw_
            out['sig_s'] = np.sqrt(np.maximum(t[:, _I['mss']] / I - ds_ ** 2, 0))
            out['sig_f'] = np.sqrt(np.maximum(t[:, _I['mff']] / I - df_ ** 2, 0))
            out['sig_omega'] = np.sqrt(np.maximum(t[:, _I['mww']] / I - dw_ ** 2, 0))
        ok = np.isfinite(out['s']) & np.isfinite(out['f'])
        sc, fc = ImageD11.blobcorrector.correct_dxdy_bilinear(np.where(ok, out['s'], 0), np.where(ok, out['f'], 0), self.e2dx, self.e2dy)
        out['sc'] = np.where(ok, sc, np.nan)
        out['fc'] = np.where(ok, fc, np.nan)
        return out

    def gvectors(self, m, use='measured'):
        """g-vectors (n, 3) of merged reflections from their measured centroids (use='measured') or the
        predicted positions ('predicted'), for each reflection's own grain position"""
        if use == 'measured':
            sc, fc, om = m['sc'], m['fc'], m['omega']
        else:
            sc, fc, om = m['pred_sc'], m['pred_fc'], m['pred_omega']
        T = np.array([g.translation for g in self.grains], float)[m['grain']]
        p = self.pars
        tth, eta = tth_eta_from_pixels(sc, fc, lab_x(T.T, om), p)
        gv = ImageD11.transform.compute_g_vectors(tth, eta, om, p['wavelength'], wedge=p.get('wedge', 0.0),
                                                  chi=p.get('chi', 0.0))
        return gv.T

    def to_columnfile(self, m=None):
        """merged results as an ImageD11 columnfile"""
        if m is None:
            m = self.merged()
        return ImageD11.columnfile.colfile_from_dict(dict((k, np.asarray(v, float)) for k, v in m.items()))

    # ---------------------------------------------------------------- box sizing
    def auto_box(self, nscans=24, pilot=None, min_snr=20.0, nsigma=3.0, quantile=0.95, trim=4.0, nproc=None, verbose=True):
        """Choose the box from the data: integrate a sample of scans with a large pilot box, take the
        reflections with I / sigI > min_snr, and set each half size to
            ceil(nsigma * median rms width + `quantile` of |centroid - prediction|)
        in s, f and omega (frames). With many grains a large pilot box often catches a neighbouring
        reflection too, so only the core of the strong reflections is used: offsets below `trim` times the
        median offset (iterated). The annulus width is kept from the pilot. Sets and returns self.box; the
        measurements are kept in self.box_stats."""
        if pilot is None:
            pilot = Box(half_s=12, half_f=12, half_omega=8, bg_width=6, bg_estimator='clipped')
        keep_box = self.box
        self.box = pilot
        self.predict()
        present = self.present_scans()
        sel = [present[i] for i in np.linspace(0, len(present) - 1, min(nscans, len(present))).astype(int)]
        m = self.run(scans=sel, nproc=nproc, verbose=False)
        step = float(np.median(np.abs(np.diff(self.ds.omega[sel[0]]))))
        strong = (m['I'] > min_snr * m['sigI']) & ~m['touch_mask'] & ~m['edge']
        if strong.sum() < 20:
            self.box = keep_box
            raise RuntimeError('auto_box: only %d strong reflections in the pilot scans' % strong.sum())
        # With many grains, a large pilot box often also catches a neighbouring reflection, which drags the
        # centroid and inflates the moments. Keep the core: offsets within `trim` robust widths of the median.
        d_s = np.abs(m['s'] - m['pred_s'])
        d_f = np.abs(m['f'] - m['pred_f'])
        d_w = np.abs(m['omega'] - m['pred_omega']) / step
        core = strong.copy()
        for _ in range(3):
            lim = [np.median(d[core]) * trim + 0.5 for d in (d_s, d_f, d_w)]
            core = strong & (d_s < lim[0]) & (d_f < lim[1]) & (d_w < lim[2])
        ws = np.median(m['sig_s'][core]); wf = np.median(m['sig_f'][core]); ww = np.median(m['sig_omega'][core]) / step
        os_ = np.quantile(d_s[core], quantile)
        of_ = np.quantile(d_f[core], quantile)
        ow_ = np.quantile(d_w[core], quantile)
        box = Box(half_s=int(np.ceil(nsigma * ws + os_)), half_f=int(np.ceil(nsigma * wf + of_)),
                  half_omega=int(np.ceil(nsigma * ww + ow_)), bg_width=pilot.bg_width, bg_estimator=pilot.bg_estimator,
                  p_clip=pilot.p_clip)
        self.box_stats = dict(n_strong=int(strong.sum()), n_core=int(core.sum()), rms_s=ws, rms_f=wf, rms_omega_frames=ww,
                              offset_s=os_, offset_f=of_, offset_omega_frames=ow_, omega_step=step)
        if verbose:
            print('auto_box: %d strong reflections, %d in the core (%.0f%% trimmed as contaminated); rms width s %.2f '
                  'f %.2f px, omega %.2f frames; %.0f%% |offset| s %.2f f %.2f px, omega %.2f frames -> %r'
                  % (strong.sum(), core.sum(), 100 * (1 - core.mean() / strong.mean()), ws, wf, ww, quantile * 100,
                     os_, of_, ow_, box))
        self.box = box
        self.spots = None
        self.results = None
        self.predict()
        return box

    # ---------------------------------------------------------------- viewing
    def spot_frames(self, ispot, iscan, extra=4):
        """Pixel stack (frames, rows, cols) around reflection ispot in scan iscan, covering the box and
        annulus (plus `extra` frames either side), the omega of each frame, and the box corner."""
        sp = self.spots
        om_f = np.asarray(self.ds.omega[iscan], float)
        step = float(np.median(np.abs(np.diff(om_f))))
        fr = np.flatnonzero(np.abs(om_f - sp['omega'][ispot]) <= (sp['ho'][ispot] + extra + 0.5) * step)
        fr = fr[np.argsort(om_f[fr])]
        hs, hf = int(sp['hs'][ispot] + sp['bw'][ispot]), int(sp['hf'][ispot] + sp['bw'][ispot])
        s0, f0 = int(sp['s0'][ispot]), int(sp['f0'][ispot])
        out = np.zeros((len(fr), 2 * hs + 1, 2 * hf + 1), np.float32)
        rd = open_reader(self.scan_file(iscan), self.ds.limapath, self._good, self.reader)
        nfast = self._good.shape[1]
        try:
            for j, i in enumerate(fr):
                idx, vals = rd(int(i))
                r, c = idx // nfast, idx % nfast
                k = (np.abs(r - s0) <= hs) & (np.abs(c - f0) <= hf)
                out[j, r[k] - s0 + hs, c[k] - f0 + hf] = vals[k] * self.scale[iscan, i]
        finally:
            rd.close()
        return out, om_f[fr], (s0 - hs, f0 - hf)

    def show_spot(self, ispot, iscan, extra=2, ncols=8, vmax=None):
        """Plot the frames of one reflection in one scan, with the box (red) and annulus (blue) marked"""
        import matplotlib.pyplot as plt
        stack, oms, _ = self.spot_frames(ispot, iscan, extra)
        sp = self.spots
        n = len(stack)
        nr = int(np.ceil(n / ncols))
        fig, axs = plt.subplots(nr, ncols, figsize=(1.8 * ncols, 1.9 * nr + 0.6), squeeze=False)
        vmax = vmax or max(float(stack.max()), 1.0)
        bw, phs, phf = int(sp['bw'][ispot]), int(sp['hs'][ispot]), int(sp['hf'][ispot])
        hs, hf = phs + bw, phf + bw
        inbox = np.abs(oms - sp['omega'][ispot]) <= (sp['ho'][ispot] + 0.5) * np.median(np.abs(np.diff(oms)))
        for j, ax in enumerate(axs.ravel()):
            ax.set_axis_off()
            if j >= n:
                continue
            ax.imshow(stack[j], vmin=0, vmax=vmax, cmap='viridis', interpolation='nearest')
            if inbox[j]:
                ax.add_patch(plt.Rectangle((bw - 0.5, bw - 0.5), 2 * phf + 1, 2 * phs + 1,
                                           fill=False, ec='r', lw=1))
                ax.add_patch(plt.Rectangle((-0.5, -0.5), 2 * hf + 1, 2 * hs + 1, fill=False, ec='c', lw=1))
            ax.plot(sp['f'][ispot] - sp['f0'][ispot] + hf, sp['s'][ispot] - sp['s0'][ispot] + hs, 'w+', ms=6)
            ax.set_title('%.2f' % oms[j], fontsize=8)
        fig.suptitle('grain %d (%d %d %d) scan %d' % (sp['grain'][ispot], sp['h'][ispot], sp['k'][ispot], sp['l'][ispot], iscan))
        return fig
