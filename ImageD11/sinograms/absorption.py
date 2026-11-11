"""
ImageD11.sinograms.absorption
=============================

Sample absorption for scanning 3DXRD. The sample cross-section comes from an MLEM reconstruction of the sinogram
of all the integrated intensity (every reflection's box sums, binned by omega and dty), so it does not depend on
the grain maps. The sample is taken to extend along z (the rotation axis) with that cross-section: the
diffracted ray's path is its chord in the cross-section divided by the horizontal fraction of its direction.
The 3D shape is not known; for small 2theta this is the beam chord through the sample at (omega, dty).

Lengths are in the units of dty (the reconstruction pixel is ystep); mu is per dty unit (see mu_per_dty_unit).
"""

from __future__ import print_function, division

import numpy as np
import numba
import scipy.ndimage

from ImageD11.sinograms import geometry, profile_fit


# length of one dty unit in cm
DTY_UNIT_CM = {'nm': 1e-7, 'um': 1e-4, 'mm': 0.1}


def mu_per_dty_unit(formula, density, wavelength, dty_unit):
    """linear absorption coefficient per dty unit of a compound (chemical formula, e.g. 'LaMn0.4Ga0.6O3';
    density in g/cm3) at the wavelength (Angstrom), from xraylib (total attenuation cross section).
    dty_unit: 'nm', 'um' or 'mm', the length of one unit of dty (and of the grain translations)."""
    import xraylib
    energy = 12.398419843320026 / wavelength          # keV
    return xraylib.CS_Total_CP(str(formula), energy) * float(density) * DTY_UNIT_CM[dty_unit]


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


@numba.njit(cache=True)
def _inside(A, x0, y0, ux, uy, n, tlo, thi):
    """(tlo, thi) narrowed to where the pixel coordinate A . (x0 + ux t, y0 + uy t, 1) is within a pixel of 0 .. n - 1"""
    a0 = A[0] * x0 + A[1] * y0 + A[2]
    a1 = A[0] * ux + A[1] * uy
    if a1 > 1e-12:
        return max(tlo, (-1.0 - a0) / a1), min(thi, (n - a0) / a1)
    if a1 < -1e-12:
        return max(tlo, (n - a0) / a1), min(thi, (-1.0 - a0) / a1)
    if a0 < -1.0 or a0 > n:
        return tlo, -1.0
    return tlo, thi


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
            # the samples t = k step with k kmin .. kmax: only those within a pixel of the grid can be inside it
            tlo, thi = _inside(Ai, x0[r], y0[r], ux, uy, ni, 0.0, lmax)
            tlo, thi = _inside(Aj, x0[r], y0[r], ux, uy, nj, tlo, thi)
            L = 0.0
            for k in range(int(np.ceil(tlo / step)), int(np.floor(thi / step)) + 1):
                t = k * step
                if t >= lmax:
                    break
                i = int(np.rint(Ai[0] * (x0[r] + ux * t) + Ai[1] * (y0[r] + uy * t) + Ai[2]))
                j = int(np.rint(Aj[0] * (x0[r] + ux * t) + Aj[1] * (y0[r] + uy * t) + Aj[2]))
                if 0 <= i < ni and 0 <= j < nj and mask[i, j]:
                    L += step
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
    ni, nj = mask.shape
    # geometry.recon_to_sample: sx = (i - ni // 2) ystep, sy = -(j - nj // 2) ystep
    Ai = np.array([1.0 / ystep, 0.0, ni // 2])                # sample (x, y) -> reconstruction pixel
    Aj = np.array([0.0, -1.0 / ystep, nj // 2])
    w, t2, e = np.radians(omega), np.radians(tth), np.radians(eta)
    c, s = np.cos(w), np.sin(w)
    dx, dy = np.cos(t2), -np.sin(t2) * np.sin(e)              # horizontal part of the diffracted direction (lab)
    hn = np.hypot(dx, dy)
    # lab = R(omega) sample: x_l = x cos w - y sin w, y_l = x sin w + y cos w
    uxi, uyi = -c, s                                          # back along the incident beam, sample frame
    uxo, uyo = (dx * c + dy * s) / hn, (-dx * s + dy * c) / hn
    step = 0.5 * ystep if step is None else float(step)
    lmax = 2.0 * max(ni - 1, nj - 1) * ystep                  # twice the larger extent of the grid
    return _paths(x, y, uxi, uyi, uxo, uyo, hn, np.ascontiguousarray(mask, dtype=np.uint8), Ai, Aj, step, lmax)
