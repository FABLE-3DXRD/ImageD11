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
import time
import multiprocessing

import numpy as np
import h5py
import numba

import ImageD11.transform
import ImageD11.unitcell
import ImageD11.columnfile

try:
    import hdf5plugin  # noqa: F401  (registers the bitshuffle / zstd filters for h5py)
except ImportError:
    hdf5plugin = None

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
        sigma_omega^2 = (frame * step)^2 + (phi_tot^2 + (2 tan(theta) eps_tot)^2) / sin(eta)^2
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
                 sin_eta_min=0.05, bg_estimator='clipped', p_clip=1e-3):
        self.psf, self.bandpass, self.phi0, self.phi0_t = psf, bandpass, phi0, phi0_t
        self.phi, self.eps, self.frame = phi, eps, frame
        self.nsigma, self.margin_px, self.margin_frames = nsigma, margin_px, margin_frames
        self.bg_factor, self.bg_min = bg_factor, int(bg_min)
        self.max_half_frames, self.sin_eta_min = int(max_half_frames), sin_eta_min
        if bg_estimator not in ('mean', 'clipped'):
            raise ValueError("bg_estimator must be 'mean' or 'clipped'")
        self.bg_estimator, self.p_clip = bg_estimator, float(p_clip)

    def __repr__(self):
        return ('ShoeboxModel(psf=%g px, bandpass=%g, phi0=%g deg, phi=%s deg, eps=%s, nsigma=%g, bg_factor=%g)' %
                (self.psf, self.bandpass, self.phi0, np.round(np.median(self.phi), 4), np.round(np.median(self.eps), 6),
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
        sig_w = np.sqrt(self.frame ** 2 + (np.degrees(rock) / omega_step) ** 2)
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
    plane) and the instrumental d-spacing spread, reported as `bandpass` = 2.355 x eps0) to the per-scan box moments of strong, clean reflections in a harvest
    file written by Harvester.run(outfile=...). The grain terms (model.phi, model.eps) are kept fixed.
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
            rows = np.sort(np.random.default_rng(seed).choice(rows, nmax, replace=False))
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

    def setp(p):
        trial.psf, trial.phi0, trial.phi0_t, trial.bandpass = p[0], p[1], p[2], p[3] * 2.355

    def resid(p):
        setp(p)
        a, b, c = trial.sigmas(sub, pars, step)
        return np.concatenate([(a - meas_s) / (meas_s + 0.3), (b - meas_f) / (meas_f + 0.3), (c - meas_w) / (meas_w + 0.3)])
    p0 = [model.psf, model.phi0, model.phi0 if model.phi0_t is None else model.phi0_t, model.bandpass / 2.355]
    r = scipy.optimize.least_squares(resid, p0, bounds=([0, 0, 0, 0], [10, 5, 5, 0.05]), loss='soft_l1', f_scale=0.3)
    setp(r.x)
    a, b, c = trial.sigmas(sub, pars, step)
    info = dict(nrows=len(rows), psf=r.x[0], phi0=r.x[1], phi0_t=r.x[2], eps0=r.x[3],
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
        codec = bslz4_to_sparse.detect_codec(self.d)
        self.dec = bslz4_to_sparse.chunk2sparse(np.ascontiguousarray(good, dtype=np.uint8), dtype=self.d.dtype,
                                                codec=codec)

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


def detector_position(tth, eta, origin_x, pars):
    """Corrected detector pixel (sc, fc) of rays (tth, eta) leaving the lab point (origin_x, 0, 0); origin_x
    may be an array (one per ray). Ray / detector-plane intersection with the plane from compute_xyz_lab."""
    p = dict(pars)
    rtth, reta = np.radians(tth), np.radians(eta)
    xyz = np.array([np.cos(rtth), -np.sin(rtth) * np.sin(reta), np.sin(rtth) * np.cos(reta)])
    pk = ImageD11.transform.compute_xyz_lab(np.array([(1, 0), (0, 1), (0, 0)], float).T, **p)
    d_s, d_f, d_o = pk[:, 0] - pk[:, 2], pk[:, 1] - pk[:, 2], pk[:, 2]
    origin = np.zeros((3, np.size(tth)))
    origin[0] = origin_x
    o = d_o[:, None] - origin                     # detector origin pixel, seen from each ray's origin
    norm = np.cross(d_s, d_f).dot(xyz)
    sc = np.einsum('ij,ij->j', np.cross(d_f[:, None], o, axis=0), xyz) / norm
    fc = np.einsum('ij,ij->j', np.cross(o, d_s[:, None], axis=0), xyz) / norm
    return sc, fc


def tth_eta_from_pixels(sc, fc, origin_x, pars):
    """tth, eta of corrected pixels (sc, fc) for diffraction origins (origin_x, 0, 0): compute_xyz_lab with
    a per-peak distance, as in point_by_point"""
    p = dict(pars)
    p['distance'] = pars['distance'] - np.asarray(origin_x, float)
    xyz = ImageD11.transform.compute_xyz_lab(np.array([sc, fc], float), **p)
    return ImageD11.transform.compute_tth_eta_from_xyz(xyz, None)


def invert_distortion(sc, fc, e2dx, e2dy, iterations=5):
    """Raw (s, f) pixel positions from distortion-corrected (sc, fc): corrected = raw + d(raw)"""
    if e2dx is None:
        return sc.copy(), fc.copy()
    s, f = sc.copy(), fc.copy()
    ns, nf = e2dx.shape
    for _ in range(iterations):
        i = np.clip(np.round(s).astype(int), 0, ns - 1)
        j = np.clip(np.round(f).astype(int), 0, nf - 1)
        s, f = sc - e2dy[i, j], fc - e2dx[i, j]
    return s, f


def correct_distortion(s, f, e2dx, e2dy):
    """Distortion-corrected (sc, fc) from raw (s, f), bilinear in the maps"""
    if e2dx is None:
        return s.copy(), f.copy()
    ns, nf = e2dx.shape
    s0 = np.clip(np.floor(s).astype(int), 0, ns - 2)
    f0 = np.clip(np.floor(f).astype(int), 0, nf - 2)
    ds_, df_ = np.clip(s - s0, 0, 1), np.clip(f - f0, 0, 1)

    def interp(m):
        return ((1 - ds_) * (1 - df_) * m[s0, f0] + ds_ * (1 - df_) * m[s0 + 1, f0] +
                (1 - ds_) * df_ * m[s0, f0 + 1] + ds_ * df_ * m[s0 + 1, f0 + 1])
    return s + interp(e2dy), f + interp(e2dx)


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
    s, f = invert_distortion(sc, fc, e2dx, e2dy)
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
    """Accumulate one frame into acc[spots] (see ACC). idx ascending flat pixel indices, vals counts.
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


# ------------------------------------------------------------------------------------- the harvester
_WORKER = {}


def _work_scan(iscan):
    try:
        return _WORKER['h']._harvest_scan(iscan)
    except Exception as e:          # e.g. an image file that is unreadable or still being restored
        return iscan, None, '%s: %r' % (_WORKER['h'].scan_file(iscan), e), 0, 0.0


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
    """

    def __init__(self, ds, grains, pars, ds_max, radius=None, mask=None, datapath=None, reader='auto',
                 monitor=True, box=None, e2dx=None, e2dy=None):
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
            return iscan, None, None, 0, 0.0
        frames, ptr, spot, dom = plan
        sp, b = self.spots, self.box
        ids = np.unique(spot)
        local = np.searchsorted(ids, spot).astype(np.int64)
        acc = np.zeros((len(ids), NACC))
        arr = dict((k, np.ascontiguousarray(sp[k][ids])) for k in ('s0', 'f0', 'hs', 'hf', 'bw', 'nbox', 'nann', 'S1s', 'S1f',
                                                                      'S2s', 'S2f'))
        t0 = time.time()
        rd = open_reader(self.scan_file(iscan), self.ds.limapath, self._good, self.reader)
        nfast = self._good.shape[1]
        try:
            for j, fr in enumerate(frames):
                idx, vals = rd(int(fr))
                a, z = ptr[j], ptr[j + 1]
                _integrate_frame(idx, vals, nfast, self._adj, local[a:z], arr['s0'], arr['f0'], arr['hs'], arr['hf'],
                                 arr['bw'], arr['nbox'], arr['nann'], arr['S1s'], arr['S1f'], arr['S2s'], arr['S2f'],
                                 dom[a:z], float(self.scale[iscan, fr]), b.bg_estimator == 'clipped', b.p_clip, acc)
        finally:
            rd.close()
        return iscan, ids, acc, len(frames), time.time() - t0

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
            for k in ACC:
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
            nproc = max(1, len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count())
        if nproc > 1:
            pool = multiprocessing.get_context('fork').Pool(nproc)
            it = pool.imap_unordered(_work_scan, scans)
        else:
            pool = None
            it = map(_work_scan, scans)
        try:
            for k, (iscan, ids, acc, nf, dt) in enumerate(it):
                if ids is None:
                    if acc is not None:
                        failed.append((iscan, acc))
                        if verbose:
                            print('skipped scan %d: %s' % (iscan, acc), flush=True)
                    continue
                nframes += nf
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
                if verbose and (k % 20 == 0 or k == len(scans) - 1):
                    el = time.time() - t0
                    print('%d/%d scans, %d frames, %.0f s, %.0f frames/s' % (k + 1, len(scans), nframes, el, nframes / max(el, 1e-9)),
                          flush=True)
        finally:
            if pool is not None:
                pool.close()
                pool.join()
        self.results = dict(total=tot, nscans=nscans)
        self.failed_scans = failed
        if failed and verbose:
            print('%d scans could not be read and were skipped: %s' % (len(failed), [f[0] for f in failed]))
        m = self.merged()
        if h5 is not None:
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
        sc, fc = correct_distortion(np.where(ok, out['s'], 0), np.where(ok, out['f'], 0), self.e2dx, self.e2dy)
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
