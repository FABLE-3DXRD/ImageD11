"""
ImageD11.sinograms.absorption
=============================

Sample absorption for scanning 3DXRD. The sample cross-section comes from an MLEM reconstruction of the sinogram
of all the integrated intensity (every reflection's box sums, binned by omega and dty), so it does not depend on
the grain maps. The sample is taken to extend along z (the rotation axis) with that cross-section: the
diffracted ray's path is its chord in the cross-section divided by the horizontal fraction of its direction.
The 3D shape is not known; for small 2theta this is the beam chord through the sample at (omega, dty).

Lengths are in the units of dty (the reconstruction pixel is ystep); mu is per dty unit.
"""

from __future__ import print_function, division

import numpy as np
import numba
import scipy.ndimage

from ImageD11.sinograms import geometry, profile_fit


def sinogram(omega, spot, scan, net, omega_edges, ny):
    """(nomega, ny) sum of the rows net (one per (reflection, scan)) by the omega bin of their reflection
    (omega[spot]) and their scan (dty bin)"""
    io = np.clip(np.digitize(omega[spot], omega_edges) - 1, 0, len(omega_edges) - 2)
    return np.bincount(io * ny + scan, net, minlength=(len(omega_edges) - 1) * ny).reshape(-1, ny)


def sample_shape(sino, omegas, y0, ymin, ystep, shape, niter=60, threshold=0.1):
    """(density, mask) of the sample in the reconstruction frame (shape): Poisson MLEM of the sinogram (nomega, ny)
    over the disc covered by the dty range; mask = density > threshold x its 99th percentile, holes filled"""
    ny = sino.shape[1]
    support = profile_fit.support_disc(np.zeros(2), 0.5 * (ny - 1) * ystep, shape, ystep)
    A, (ri, rj) = profile_fit.shape_projector(support, omegas, y0, ymin, ystep, ny)
    Y = np.maximum(sino, 0).ravel()
    x = np.full(len(ri), Y.sum() / max(A.dot(np.ones(len(ri))).sum(), 1e-30))
    den = A.T.dot(np.ones(len(Y)))
    for it in range(niter):
        x *= np.where(den > 0, A.T.dot(Y / np.maximum(A.dot(x), 1e-12)) / np.where(den > 0, den, 1), 0)
    rho = np.zeros(shape)
    rho[ri, rj] = x
    mask = scipy.ndimage.binary_fill_holes(rho > threshold * np.percentile(x, 99))
    return rho, mask


@numba.njit(cache=True, parallel=True)
def _paths(x0, y0, uxi, uyi, uxo, uyo, hn, mask, Ai, Aj, step, lmax):
    n = x0.shape[0]
    out = np.zeros(n)
    ni, nj = mask.shape
    for r in numba.prange(n):
        Ls = 0.0
        for d in range(2):
            ux = uxi[r] if d == 0 else uxo[r]
            uy = uyi[r] if d == 0 else uyo[r]
            t = 0.0
            L = 0.0
            while t < lmax:
                i = int(np.rint(Ai[0] * (x0[r] + ux * t) + Ai[1] * (y0[r] + uy * t) + Ai[2]))
                j = int(np.rint(Aj[0] * (x0[r] + ux * t) + Aj[1] * (y0[r] + uy * t) + Aj[2]))
                if 0 <= i < ni and 0 <= j < nj and mask[i, j]:
                    L += step
                t += step
            Ls += L if d == 0 else L / hn[r]
        out[r] = Ls
    return out


def path_lengths(mask, ystep, x, y, omega, tth, eta, step=None):
    """path in the sample of each measurement: the incident beam (along +x lab) to the point (x, y) (sample frame)
    and the diffracted ray (tth, eta) from there, for the cross-section mask (reconstruction frame) extended along
    z. omega, tth, eta in degrees, one per measurement; x, y scalars or one per measurement."""
    omega = np.asarray(omega, float)
    x = np.broadcast_to(np.asarray(x, float), omega.shape).copy()
    y = np.broadcast_to(np.asarray(y, float), omega.shape).copy()
    ii, jj = np.mgrid[0:mask.shape[0], 0:mask.shape[1]]
    sx, sy = geometry.recon_to_sample(ii.ravel(), jj.ravel(), mask.shape, ystep)
    X = np.column_stack((sx, sy, np.ones(sx.size)))
    Ai = np.linalg.lstsq(X, ii.ravel(), rcond=None)[0]       # sample (x, y) -> reconstruction pixel
    Aj = np.linalg.lstsq(X, jj.ravel(), rcond=None)[0]
    w, t2, e = np.radians(omega), np.radians(tth), np.radians(eta)
    c, s = np.cos(w), np.sin(w)
    dx, dy = np.cos(t2), -np.sin(t2) * np.sin(e)              # horizontal part of the diffracted direction (lab)
    hn = np.hypot(dx, dy)
    # lab = R(omega) sample: x_l = x cos w - y sin w, y_l = x sin w + y cos w
    uxi, uyi = -c, s                                          # back along the incident beam, sample frame
    uxo, uyo = (dx * c + dy * s) / hn, (-dx * s + dy * c) / hn
    step = 0.5 * ystep if step is None else float(step)
    lmax = 2.0 * max(np.ptp(sx), np.ptp(sy))
    return _paths(x, y, uxi, uyi, uxo, uyo, hn, np.ascontiguousarray(mask, dtype=np.uint8), Ai, Aj, step, lmax)
