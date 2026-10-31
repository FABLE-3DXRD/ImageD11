"""
ImageD11.sinograms.profile_fit
==============================

Intensities and grain shapes from the per-(reflection, scan) rows of a harvest (ImageD11.sinograms.harvester).

A grain diffracts on the scans (dty) where it is in the beam. For a reflection at omega, the intensity versus dty is
the projection of the grain shape at that omega, times the reflection intensity. Here:

* project / shape_projector: that projection, for a density image in the reconstruction frame of
  ImageD11.sinograms.geometry, with a triangular beam footprint one ystep wide.
* row_estimates: per (reflection, scan) the box sum and the profile-weighted estimate of a harvest run with
  weights (wA / wB), each with a variance from the model (background + a model intensity), not from the observed
  counts: count-based weights pull the weighted means of weak reflections down.
* fit_dty: per reflection, the scale of the projected shape that fits the rows of the scans, with rows far above or
  below the profile rejected (another grain's peak in the box adds counts on a few scans only).
* intensities: the two together. Weak reflections (box-sum I / sigI < t_prof) take the profile-fitted intensity,
  the others the plain box sum, as in single-crystal practice: once counting statistics are good, an imperfect
  shape adds more error than it removes.
* mlem_shape: the grain shape from the box counts by Poisson MLEM with the reflection intensities in the forward
  model, box(y) = I_r [R x](y) + bg(y): reflections count in proportion to their photons. Use the strongest
  reflections only; points far above the model (other grains' peaks) are masked during the iterations.
"""

from __future__ import print_function, division

import numpy as np
import scipy.sparse

from ImageD11.sinograms import geometry


def _footprint(sx, sy, omegas, y0, ymin, ystep, ny):
    """(projection index, dty bin, voxel index, weight) of the triangular beam footprint of points (sx, sy)"""
    w = np.radians(np.asarray(omegas, float))[:, None]
    dty = y0 - sx[None, :] * np.sin(w) - sy[None, :] * np.cos(w)
    u = (dty - ymin) / ystep
    i0 = np.floor(u).astype(np.int64)
    f1 = u - i0
    p = np.repeat(np.arange(len(omegas)), len(sx)).reshape(len(omegas), len(sx))
    v = np.tile(np.arange(len(sx)), (len(omegas), 1))
    out = []
    for ii, ww in ((i0, 1 - f1), (i0 + 1, f1)):
        ok = (ii >= 0) & (ii < ny)
        out.append((p[ok], ii[ok], v[ok], ww[ok]))
    return [np.concatenate(x) for x in zip(*out)]


def project(rho, omegas, y0, ymin, ystep, ny, chunk=2e7):
    """(len(omegas), ny): projection of the density image rho (reconstruction frame, geometry.recon_to_sample) at
    each omega onto the dty bins ymin + i ystep (i = 0 .. ny - 1)"""
    ri, rj = np.nonzero(rho > 0)
    w = rho[ri, rj].astype(float)
    sx, sy = geometry.recon_to_sample(ri, rj, rho.shape, ystep)
    out = np.zeros((len(omegas), ny))
    if len(w) == 0:
        return out
    n = max(1, int(chunk // len(w)))
    for a in range(0, len(omegas), n):
        p, i, v, f = _footprint(sx, sy, omegas[a:a + n], y0, ymin, ystep, ny)
        o = np.zeros(out[a:a + n].size)
        np.add.at(o, p * ny + i, f * w[v])
        out[a:a + n] = o.reshape(-1, ny)
    return out


def shape_projector(support, omegas, y0, ymin, ystep, ny):
    """sparse (len(omegas) * ny, nvox) matrix projecting the voxels of the bool image support (row p * ny + i:
    projection p, dty bin i) and the voxel (row, column) indices"""
    ri, rj = np.nonzero(support)
    sx, sy = geometry.recon_to_sample(ri, rj, support.shape, ystep)
    p, i, v, f = _footprint(sx, sy, omegas, y0, ymin, ystep, ny)
    A = scipy.sparse.csr_matrix((f, (p * ny + i, v)), shape=(len(omegas) * ny, len(ri)))
    return A, (ri, rj)


def support_disc(translation, radius, shape, ystep):
    """bool image (reconstruction frame) of the disc of radius (sample units) around the grain position"""
    ri, rj = geometry.sample_to_recon(translation[0], translation[1], shape, ystep)
    ii, jj = np.mgrid[0:shape[0], 0:shape[1]]
    return (ii - ri) ** 2 + (jj - rj) ** 2 <= (radius / ystep) ** 2


def row_estimates(box, bg, var, nframes, nbox, w=None, I_model=None):
    """Per (scan, reflection) estimates and model variances. box, bg, var, nframes: the harvest per_scan columns
    as (nscan, nrefl) arrays; nbox: good pixels per box (nrefl,); w: dict of the profile-weighted columns (wA, wB,
    wM, wC, wVbg), or None; I_model: (nscan, nrefl) model intensity of each row (None: 0).
    Box sum: net = box - bg, variance bg + I_model + (var - box) (the background-estimate term), at least 1.
    Profile: wA / wB, variance (max(wM, wB / (nbox nframes)) + I_model wC + wVbg) / wB^2 (at least one background
    count expected over the box). Returns dict(net, v_sum, prof, v_prof) (prof, v_prof None without w)."""
    Im = 0.0 if I_model is None else np.maximum(I_model, 0)
    net = box - bg
    v_sum = np.maximum(np.maximum(bg, 0) + Im + np.maximum(var - box, 0), 1.0)
    out = dict(net=net, v_sum=v_sum, prof=None, v_prof=None)
    if w is not None:
        B = w['wB']
        ok = B > 0
        Bs = np.where(ok, B, 1.0)
        npf = np.maximum(nbox[None, :] * nframes, 1)
        wM = np.maximum(w['wM'], B / npf)
        out['prof'] = np.where(ok, w['wA'] / Bs, net)
        out['v_prof'] = np.where(ok, (wM + Im * w['wC'] + w['wVbg']) / Bs ** 2, v_sum)
    return out


def fit_dty(x, v, seen, shape, rej=4.0, niter=8):
    """Fit x (nscan, nrefl) with variances v by scale * shape (nscan, nrefl: the projected grain shape at each
    reflection's omega, on the same scans). Rows more than rej sigma ABOVE the fit are rejected, the worst one per
    reflection and iteration (another grain's peak only adds counts; rejecting all large residuals at once also
    throws out the good rows when one outlier has pulled the scale up). At least 4 rows are kept.
    Returns dict(I = scale x shape summed over the measured rows, sigI, scale, keep (rows used), chi2); I is NaN
    with fewer than 4 measured rows or no shape on them."""
    iv = 1.0 / np.maximum(v, 1e-12)
    K = seen.copy()
    cols = np.arange(x.shape[1])
    for it in range(niter):
        w2 = K * iv
        den = (shape * shape * w2).sum(0)
        sc = np.where(den > 0, (shape * x * w2).sum(0) / np.where(den > 0, den, 1), 0.0)
        res = (x - sc[None, :] * shape) * np.sqrt(iv)
        worst = np.argmax(np.where(K, res, -np.inf), axis=0)
        bad = (res[worst, cols] > rej) & (K.sum(0) > 4)
        if it == niter - 1 or not bad.any():
            break
        K[worst[bad], cols[bad]] = False
    pseen = (shape * seen).sum(0)
    ok = (seen.sum(0) >= 4) & (den > 0) & (pseen > 0)
    I = np.where(ok, sc * pseen, np.nan)
    sig = np.where(ok, pseen / np.sqrt(np.where(den > 0, den, 1)), np.nan)
    chi2 = (res ** 2 * K).sum(0) / np.maximum(K.sum(0) - 1, 1)
    return dict(I=I, sigI=sig, scale=sc, keep=K, chi2=chi2)


def intensities(box, bg, var, nframes, nbox, shape, w=None, t_prof=3.0, rej=4.0):
    """Per-reflection intensities of one grain from its harvest rows (see row_estimates) and the projected shape
    (nscan, nrefl). Box sums for the reflections with I / sigI >= t_prof, the profile-fitted intensity (dty-shape
    fit of the per-row best estimate, profile-weighted or box sum, whichever has the smaller variance) below it.
    Two passes: the row variances use the model intensities of the first. Returns dict(I, sigI, I_sum, sigI_sum,
    I_prof, sigI_prof, profile_used, model (fitted rows), keep, chi2)."""
    seen = nframes > 0
    Im = None
    for it in range(2):
        r = row_estimates(box, bg, var, nframes, nbox, w, Im)
        if r['prof'] is None:
            x, vx = r['net'], r['v_sum']
        else:
            use_w = r['v_prof'] < r['v_sum']
            x = np.where(use_w, r['prof'], r['net'])
            vx = np.where(use_w, r['v_prof'], r['v_sum'])
        f = fit_dty(x, vx, seen, shape, rej=rej)
        Im = np.where(seen, f['scale'][None, :] * shape, 0)
    I_sum = np.where(seen, r['net'], 0).sum(0)
    s_sum = np.sqrt(np.where(seen, r['v_sum'], 0).sum(0))
    weak = (I_sum < t_prof * s_sum) & np.isfinite(f['I'])
    return dict(I=np.where(weak, f['I'], I_sum), sigI=np.where(weak, f['sigI'], s_sum), I_sum=I_sum, sigI_sum=s_sum,
                I_prof=f['I'], sigI_prof=f['sigI'], profile_used=weak, model=Im, keep=f['keep'], chi2=f['chi2'])


def mlem_shape(box, bg, seen, I, omegas, support, y0, ymin, ystep, niter=40, outsig=5.0, prop=0.1,
               out_at=(5, 15)):
    """Grain shape (density on the bool image support, reconstruction frame) from the box counts of its reflections:
    box, bg, seen are (nrefl, ny) over all dty bins, I (nrefl,) the reflection intensities, omegas (nrefl,).
    Poisson MLEM for box = I [R x] + bg:  x <- x R^T[I box / model] / R^T[I], measured points only. At the
    iterations out_at, points more than outsig sigma above the model (sigma^2 = model + (prop x signal)^2) are
    masked. Returns (density image, number of masked points)."""
    ny = box.shape[1]
    A, (ri, rj) = shape_projector(support, omegas, y0, ymin, ystep, ny)
    Y, B, M = box.ravel().astype(float), bg.ravel().astype(float), seen.ravel().astype(float)
    Ir = np.repeat(np.asarray(I, float), ny)
    x = np.ones(len(ri))
    x *= (np.maximum(Y - B, 0) * M).sum() / max((Ir * A.dot(x) * M).sum(), 1e-30)
    den = A.T.dot(Ir * M)
    nout = 0
    for it in range(niter):
        sig_ = A.dot(x) * Ir
        model = sig_ + B
        if it in out_at:
            out = (M > 0) & ((Y - model) / np.sqrt(np.maximum(model, 1.0) + (prop * sig_) ** 2) > outsig)
            nout += int(out.sum())
            M[out] = 0
            den = A.T.dot(Ir * M)
        ratio = np.where(M > 0, Y / np.maximum(model, 1e-6), 0)
        x *= np.where(den > 0, A.T.dot(Ir * ratio) / np.where(den > 0, den, 1), 0)
    rho = np.zeros(support.shape)
    rho[ri, rj] = x
    return rho, nout
