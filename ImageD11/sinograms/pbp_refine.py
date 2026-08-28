"""Point-by-point refinement of orientation and strain.

Originally written by Axel Henningsson, with adaptations by James Ball.

https://www.nature.com/articles/s41598-024-71006-0

The method follows Henningsson et al., and the equation numbers below refer
to the supplementary material of that paper. This module refines what the
point-by-point indexer produced: for every voxel and every candidate
orientation it selects the peaks that saw that voxel, refines the unit cell
matrix, and then fits an elastic strain tensor.

Notation, and where each symbol lives in the code
-------------------------------------------------
    G           measured diffraction vector, in sample coordinates
                    G = U B G_hkl                                   eq (5)
                (eq (4) is the same in the lab frame, with the turntable
                rotation Omega applied; everything here works in sample
                coordinates, so Omega does not appear)
    G_hkl       integer Miller indices                       -> fit_hkl
    (UB)^-1     unit cell matrix inverse, "the UBI"          -> ubi_out
    U, B0       orientation, and the undeformed reference B  -> Umat, B0
    y + dy      distance from the beam to the voxel centroid -> g_yd

Peak selection
--------------
A peak is assigned to a voxel when the beam was illuminating that voxel as
the peak was recorded, i.e. when

    | y0 - sx sin(w) - sy cos(w) - dty |  <=  ystep

for the peak's own omega and dty. sx, sy is the voxel centroid in sample
coordinates and y0 the rotation axis offset. ystep is the scan step, which
is set to the measured beam FWHM.

Peaks recorded at different omega that share (h, k, l), the sign of eta and
the scan line are merged into one intensity-weighted peak before fitting.

Orientation refinement, eq (6) to (11)
--------------------------------------
Each merged peak carries a weight that falls off with how far the beam
passed from the voxel centre:

    w = r / (y + dy + r)                                     eq (6)

r is a regularisation length, so that w = 1 for a peak straight through the
centroid. The paper writes r as 1 um for a 3 um beam; it is beam_size / 3,
and lives in the code as `weight_reg`. This is the only absolute length in
the module -- everything else is a ratio of lengths and is therefore free of
any assumption about the dty motor units.

A peak counts towards the fit when its Miller indices come out close enough
to integers:

    || (UB)^-1 G - round((UB)^-1 G) ||_2  <  e_hkl            eq (8)

with e_hkl = hkl_tol_refine before merging and hkl_tol_refine_merged after.

The unit cell matrix is then the weighted least squares solution over all
such peaks, eq (11). The code solves the equivalent normal equations
directly, accumulating

    GtWWG = sum_n w^2 G G^T          GtWWH = sum_n w^2 G G_hkl^T

and getting (UB)^-1 from a 3x3 Cholesky rather than forming the full system
of eq (7) and (10).

Strain refinement, eq (12) to (20)
----------------------------------
Each merged peak is reduced to one scalar strain along its own scattering
direction, by comparing the measured G with the one the refined orientation
and an undeformed reference cell predict, G0 = U B0 G_hkl:

    eps = (G^T G0) / (G^T G) - 1                             eq (12)

The elastic strain tensor is the unknown that reproduces those numbers,

    eps = (G^T eps_tensor G) / (G^T G)                       eq (13)

which is linear in the six independent components, so many measurements form

    y = M s                                                  eq (14)
    s = [eps_xx, eps_yy, eps_zz, eps_xy, eps_xz, eps_yz]     eq (15)
    M_j = [k1^2, k2^2, k3^2, 2 k1 k2, 2 k1 k3, 2 k2 k3]      eq (16)
    kappa = G / ||G||_2                                      eq (17)

Weights combine the beam-proximity factor of eq (6) with the g-vector noise
sigma_g propagated through eq (12):

    w = ( r / (y + dy + r) ) * ( sigma_g^2 G0^T G0 / (G^T G)^2 )^(-1/2)
                                                             eq (18)

Measurements more than 3.5 standard deviations from the mean directional
strain are dropped, eq (19), and the tensor is the weighted least squares
solution, eq (20). As with the orientation fit the code accumulates the
normal equations (MtWWM, MtWWy) and solves a 6x6 Cholesky rather than
building M and W.

The whole procedure runs independently for every voxel and every candidate
orientation.

Where this code differs from the paper
--------------------------------------
Four differences, none of them large, all deliberate:

  * eq (10) and (11) fit the flattened UB by minimising the residual in
    G space. This code minimises in hkl space and obtains (UB)^-1 directly,
    which is what ImageD11's own indexing.refine does (Paciorek et al.,
    Acta A55 543, 1999). Least squares is not invariant under inverting the
    model, so the two differ -- by about 1e-8 in strain, four orders below
    the strain precision.

  * eq (18) as printed has (G^T G) where the code has (G^T G)^2. Propagating
    sigma_g through eq (12) gives sigma_eps = sigma_g / ||G||, which is what
    the code computes; the printed form reduces to a constant and would not
    weight anything. Confirmed typo with Axel.

  * eq (8) as printed omits the subtraction of the rounded indices. The code
    tests the residual, which is what the surrounding text describes.

  * the origin calculation in the original implementation compared the
    SQUARED hkl residual against an unsquared hkl_tol, so a nominal 0.05
    was really a tolerance of sqrt(0.05) = 0.224. Axel has confirmed this
    was a typo. Here hkl_tol_origins is squared before comparison, matching
    hkl_tol_refine, so 0.05 now means 0.05. Origins computed with earlier
    versions of this code, and with the original, used the looser value.

Two apparent typos in the supplement, for anyone reading along: eq (9) lists
UB23 twice and omits UB22, and eq (15) is labelled y where eq (14) and (20)
call it s.
"""

from __future__ import division, print_function

import os

from ImageD11.peakselect import mask_rings_by_ifrac

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import ctypes
import platform
import subprocess
import sys
import time

import h5py
import numba
import numpy as np
from skimage.filters import threshold_otsu
from skimage.morphology import convex_hull_image

import ImageD11.indexing
import ImageD11.sinograms.dataset
from ImageD11 import cImageD11, unitcell
from ImageD11.sinograms import geometry
from ImageD11.sinograms.point_by_point import PBPMap
from ImageD11.sinograms.roi_iradon import run_iradon
from ImageD11.sinograms.sinogram import save_array
from ImageD11.sinograms.tensor_map import unitcell_to_b
from ImageD11.sinograms.voxel_mask import (
    VoxelSinoMasker,
    choose_omega_bins,
    fill_voxel_idx,
)
from ImageD11.sinograms.voxel_mask import max_candidates as vm_max_candidates

# ignore Numba dot product performance warnings?
# import warnings
# import numba
# warnings.simplefilter('ignore', category=numba.core.errors.NumbaPerformanceWarning)


def _tune_malloc(mmap_threshold=1 << 30, trim_threshold=None):
    """
    Raise glibc's mmap threshold so numba's NRT allocations stay in the
    arena instead of being mmap'd and munmap'd per alloc/free.

    Only matters for allocations above the 128 kB default -- in refine_map
    that is compute_gve's (3, ng) temporaries, so ng > ~5500. The per-voxel
    scratch is allocated once per chunk and never hits this path.

    trim_threshold is left alone by default: raising it stops glibc
    returning freed memory to the OS, which is the wrong trade in a process
    handling hundreds of millions of peaks.

    Harmless if it fails (musl, macOS, Windows).
    """
    if platform.system() != "Linux":
        return False
    try:
        libc = ctypes.CDLL("libc.so.6", use_errno=True)
        M_TRIM_THRESHOLD = -1
        M_MMAP_THRESHOLD = -3
        ok = libc.mallopt(M_MMAP_THRESHOLD, ctypes.c_int(mmap_threshold))
        if trim_threshold is not None:
            ok &= libc.mallopt(M_TRIM_THRESHOLD, ctypes.c_int(trim_threshold))
        return bool(ok)
    except Exception:
        return False


_MALLOC_TUNED = _tune_malloc()


# ==========================================================================
#  Packed sort keys
#
#  key = ((((h+2048)*4096 + (k+2048))*4096 + (l+2048))*4 + s)*2**20 + dtyi'
#
#  h, k, l in (-2048, 2048); s in {0, 1, 2} for etasign {-1, 0, +1};
#  dtyi' = dtyi + 2**19.  Max magnitude 2**58, comfortably inside int64.
#  Monotone in (h, k, l, etasign, dtyi), so sorting packed keys reproduces
#  the lexicographic tuple sort used by count_unique_peaks
# ==========================================================================
_HOFF = 2048
_HMUL = 4096
_DOFF = 1 << 19
_DMUL = 1 << 20


@numba.njit(cache=True, inline="always")
def _etasign_code(e):
    if e < 0.0:
        return 0
    elif e == 0.0:
        return 1
    return 2

@numba.njit(cache=True, inline="always")
def _clamp_hkl(v):
    """|h| >= 2048 means a nonsense UBI. Clamping is safer than letting the
    packed int64 key wrap silently and mis-merge unrelated reflections."""
    if v > 2047:
        return 2047
    if v < -2047:
        return -2047
    return v

@numba.njit(cache=True, inline="always")
def _pack5(h, k, l, s, d):
    return ((((h + _HOFF) * _HMUL + (k + _HOFF)) * _HMUL + (l + _HOFF)) * 4 + s) * _DMUL + (d + _DOFF)


@numba.njit(cache=True, inline="always")
def _pack4(h, k, l, s):
    return (((h + _HOFF) * _HMUL + (k + _HOFF)) * _HMUL + (l + _HOFF)) * 4 + s


# ------------------------------- allocation-free argsort (heapsort) --------
@numba.njit(cache=True, inline="always")
def _siftdown(keys, idx, start, end):
    root = start
    while True:
        child = 2 * root + 1
        if child > end:
            return
        if child + 1 <= end and keys[idx[child]] < keys[idx[child + 1]]:
            child += 1
        if keys[idx[root]] < keys[idx[child]]:
            t = idx[root]
            idx[root] = idx[child]
            idx[child] = t
            root = child
        else:
            return


@numba.njit(cache=True)
def _argsort_into(keys, n, idx):
    """idx[0:n] <- argsort(keys[0:n]). No allocation."""
    for i in range(n):
        idx[i] = i
    for start in range(n // 2 - 1, -1, -1):
        _siftdown(keys, idx, start, n - 1)
    for end in range(n - 1, 0, -1):
        t = idx[0]
        idx[0] = idx[end]
        idx[end] = t
        _siftdown(keys, idx, 0, end - 1)


@numba.njit(cache=True)
def _label_groups(keys, n, srt, labels):
    """Assign a group label per entry, in ascending-key order. Returns ngroups."""
    if n == 0:
        return 0
    _argsort_into(keys, n, srt)
    lab = 0
    labels[srt[0]] = 0
    for t in range(1, n):
        if keys[srt[t]] != keys[srt[t - 1]]:
            lab += 1
        labels[srt[t]] = lab
    return lab + 1


# ==========================================================================
#  Small dense linear algebra (no LAPACK)
# ==========================================================================
@numba.njit(cache=True, inline="always")
def _det3(m):
    return (
        m[0, 0] * (m[1, 1] * m[2, 2] - m[1, 2] * m[2, 1])
        - m[0, 1] * (m[1, 0] * m[2, 2] - m[1, 2] * m[2, 0])
        + m[0, 2] * (m[1, 0] * m[2, 1] - m[1, 1] * m[2, 0])
    )


@numba.njit(cache=True, inline="always")
def _sv_extremes(m):
    """Largest and smallest singular value of a 3x3, with no allocation.
 
    They are the square roots of the extreme eigenvalues of m m^T, and a
    symmetric 3x3 eigenvalue problem has a closed form (Cardano, in the
    stable Deledalle/Smith arrangement). Squaring the matrix squares the
    condition number, so smin here is only trustworthy while cond(m) is below
    about 1e7 -- see _well_conditioned, which is what uses it.
    """
    a00 = m[0, 0] * m[0, 0] + m[0, 1] * m[0, 1] + m[0, 2] * m[0, 2]
    a11 = m[1, 0] * m[1, 0] + m[1, 1] * m[1, 1] + m[1, 2] * m[1, 2]
    a22 = m[2, 0] * m[2, 0] + m[2, 1] * m[2, 1] + m[2, 2] * m[2, 2]
    a01 = m[0, 0] * m[1, 0] + m[0, 1] * m[1, 1] + m[0, 2] * m[1, 2]
    a02 = m[0, 0] * m[2, 0] + m[0, 1] * m[2, 1] + m[0, 2] * m[2, 2]
    a12 = m[1, 0] * m[2, 0] + m[1, 1] * m[2, 1] + m[1, 2] * m[2, 2]
    q = (a00 + a11 + a22) / 3.0
    b00 = a00 - q
    b11 = a11 - q
    b22 = a22 - q
    p2 = (b00 * b00 + b11 * b11 + b22 * b22
          + 2.0 * (a01 * a01 + a02 * a02 + a12 * a12)) / 6.0
    if p2 <= 0.0:
        hi = q
        lo = q
    else:
        pp = np.sqrt(p2)
        d = (b00 * (b11 * b22 - a12 * a12)
             - a01 * (a01 * b22 - a12 * a02)
             + a02 * (a01 * a12 - b11 * a02)) / (2.0 * pp * pp * pp)
        if d < -1.0:
            d = -1.0
        if d > 1.0:
            d = 1.0
        phi = np.arccos(d) / 3.0
        hi = q + 2.0 * pp * np.cos(phi)
        lo = q + 2.0 * pp * np.cos(phi + 2.0943951023931953)   # + 2pi/3
        mid = 3.0 * q - hi - lo
        if lo > mid:
            lo = mid
        if hi < mid:
            hi = mid
    if lo < 0.0:
        lo = 0.0
    if hi < 0.0:
        hi = 0.0
    return np.sqrt(hi), np.sqrt(lo)


@numba.njit(cache=True)
def _chol_solve(A, B, n, nrhs, pivot_tol):
    """
    In-place Cholesky of the symmetric n x n block A[:n,:n] (lower triangle
    written), then solve A X = B for the nrhs columns of B[:n,:nrhs].
    Returns False if A is not positive definite to within pivot_tol.
    """
    for i in range(n):
        for j in range(i + 1):
            s = A[i, j]
            for k in range(j):
                s -= A[i, k] * A[j, k]
            if i == j:
                if s <= pivot_tol:
                    return False
                A[i, i] = np.sqrt(s)
            else:
                A[i, j] = s / A[j, j]
    for c in range(nrhs):
        for i in range(n):
            s = B[i, c]
            for k in range(i):
                s -= A[i, k] * B[k, c]
            B[i, c] = s / A[i, i]
        for i in range(n - 1, -1, -1):
            s = B[i, c]
            for k in range(i + 1, n):
                s -= A[k, i] * B[k, c]
            B[i, c] = s / A[i, i]
    return True

@numba.njit(cache=True)
def _well_conditioned(m, guard=1e6, cond_max=1e14):
    """Is m of rank 3 with cond(m) < cond_max? The worth_fitting test.
 
    np.linalg.svd on a 3x3 costs about 3 us and allocates LAPACK workspace,
    once per candidate UBI and inside a prange. _sv_extremes answers the same
    question for free, but only while cond(m) is well below 1e7. So: trust it
    when it says the matrix is comfortably conditioned (a decade of margin at
    guard = 1e6), and pay for LAPACK only in the ambiguous band, which real
    UBIs essentially never reach. Checked against pure LAPACK on 60000
    matrices from cond 1 to rank-deficient: no disagreements.
    """
    smax, smin = _sv_extremes(m)
    if smin > 0.0 and smax < smin * guard:
        return True
    sv = np.linalg.svd(m)[1]
    smax = sv[0]
    smin = sv[2]
    if smin <= 0.0:
        return False
    if not (smin > smax * 3.0 * 2.220446049250313e-16):
        return False
    return (smax / smin) < cond_max


# stuff we need to compute g-vectors
# from ImageD11.transform

def gvec_pars(pars):
    """Pack the sample<->lab part of the pars for the fast kernels.
 
    Returns (A, oc, invwvln, a_is_eye):
 
      A        (3, 3)  W(wedge) @ C(chi).  origin_lab = A @ R(omega) @ t and
                       g = R(omega)^T @ A^T @ k, matching ImageD11's
                       compute_grain_origins and compute_g_from_k.
      oc       (3, 3)  origin_lab[i] = cos(om) oc[i,0] + sin(om) oc[i,1]
                                       + oc[i,2]
      invwvln  1 / wavelength
      a_is_eye True when wedge == chi == 0, so A^T can be skipped
    """
    def g(k):
        return float(pars[k])
    wedge, chi = g("wedge"), g("chi")
    tx, ty, tz = g("t_x"), g("t_y"), g("t_z")
    cw, sw = np.cos(np.radians(wedge)), np.sin(np.radians(wedge))
    cc, sc_ = np.cos(np.radians(chi)), np.sin(np.radians(chi))
    W = np.array([[cw, 0.0, -sw], [0.0, 1.0, 0.0], [sw, 0.0, cw]], np.float64)
    C = np.array([[1.0, 0.0, 0.0], [0.0, cc, -sc_], [0.0, sc_, cc]], np.float64)
    A = np.ascontiguousarray(W @ C)
    oc = np.empty((3, 3), np.float64)
    oc[:, 0] = A[:, 0] * tx + A[:, 1] * ty
    oc[:, 1] = A[:, 1] * tx - A[:, 0] * ty
    oc[:, 2] = A[:, 2] * tz
    return A, oc, 1.0 / g("wavelength"), (wedge == 0.0 and chi == 0.0)



@numba.njit(cache=True)
def detector_rotation_matrix(tilt_x, tilt_y, tilt_z):
    r1 = np.array([[np.cos(tilt_z), -np.sin(tilt_z), 0.0],  # note this is r.h.
                   [np.sin(tilt_z), np.cos(tilt_z), 0.0],
                   [0.0, 0.0, 1.0]], np.float64)
    r2 = np.array([[np.cos(tilt_y), 0.0, np.sin(tilt_y)],
                   [0.0, 1.0, 0.0],
                   [-np.sin(tilt_y), 0.0, np.cos(tilt_y)]], np.float64)
    r3 = np.array([[1.0, 0.0, 0.0],
                   [0.0, np.cos(tilt_x), -np.sin(tilt_x)],
                   [0.0, np.sin(tilt_x), np.cos(tilt_x)]], np.float64)
    r2r1 = np.dot(np.dot(r3, r2), r1)
    return r2r1


@numba.njit(cache=True, inline="always")
def gvec_one(xl, yl, zl, xpos, so, co, A, oc, invw, a_is_eye):
    """One g-vector, no allocation and no trig.
 
    xl, yl, zl  peak position in the lab frame, at the nominal distance
    xpos        diffraction origin along lab x (the refined origin)
    so, co      sin, cos of that peak's omega
    """
    o0 = co * oc[0, 0] + so * oc[0, 1] + oc[0, 2] + xpos
    o1 = co * oc[1, 0] + so * oc[1, 1] + oc[1, 2]
    o2 = co * oc[2, 0] + so * oc[2, 1] + oc[2, 2]
    dx = xl - o0
    dy = yl - o1
    dz = zl - o2
    r = invw / np.sqrt(dx * dx + dy * dy + dz * dz)
    kx = dx * r - invw              # k_out - k_in, k_in = (1, 0, 0) / wvln
    ky = dy * r
    kz = dz * r
    if a_is_eye:
        ux = kx
        uy = ky
        uz = kz
    else:                           # u = A^T k = C(chi)^T W(wedge)^T k
        ux = A[0, 0] * kx + A[1, 0] * ky + A[2, 0] * kz
        uy = A[0, 1] * kx + A[1, 1] * ky + A[2, 1] * kz
        uz = A[0, 2] * kx + A[1, 2] * ky + A[2, 2] * kz
    return co * ux + so * uy, co * uy - so * ux, uz     # g = R(omega)^T u



@numba.njit(cache=True, parallel=True)
def gvectors_from_lab(xl, yl, zl, xpos, sinomega, cosomega,
                      A, oc, invw, a_is_eye, nblocks):
    """(N, 3) g-vectors, C-contiguous, in one parallel pass.
 
    Blocked rather than a bare prange over N, because refine_map sets
    numba.set_parallel_chunksize(1) and that would otherwise make every
    single peak its own task.
    """
    n = xl.shape[0]
    out = np.empty((n, 3), np.float64)
    if nblocks < 1:
        nblocks = 1
    bs = (n + nblocks - 1) // nblocks
    if bs < 1:
        bs = 1
    for b in numba.prange(nblocks):
        lo = b * bs
        hi = lo + bs
        if hi > n:
            hi = n
        for q in range(lo, hi):
            gx, gy, gz = gvec_one(xl[q], yl[q], zl[q], xpos[q],
                                  sinomega[q], cosomega[q],
                                  A, oc, invw, a_is_eye)
            out[q, 0] = gx
            out[q, 1] = gy
            out[q, 2] = gz
    return out


@numba.njit(cache=True)
def build_csr(cells, ncell):
    """Counting sort of item indices into `ncell` buckets.  O(N)."""
    counts = np.zeros(ncell + 1, dtype=np.int64)
    for p in range(cells.shape[0]):
        counts[cells[p] + 1] += 1
    start = np.cumsum(counts)
    order = np.empty(cells.shape[0], dtype=np.int64)
    fill = start[:-1].copy()
    for p in range(cells.shape[0]):
        c = cells[p]
        order[fill[c]] = p
        fill[c] += 1
    return start, order

# no longer needed as we go from xl, yl, zl
# @numba.njit(cache=True)
# def compute_gve(sc, fc, omega, xpos,
#                 distance, y_center, y_size, tilt_y, z_center, z_size, tilt_z, tilt_x,
#                 o11, o12, o21, o22,
#                 t_x, t_y, t_z, wedge, chi, wavelength):
#     this_distance = distance - xpos

#     tth, eta = compute_tth_eta(sc, fc, omega,
#                                y_center=y_center,
#                                y_size=y_size,
#                                tilt_y=tilt_y,
#                                z_center=z_center,
#                                z_size=z_size,
#                                tilt_z=tilt_z,
#                                tilt_x=tilt_x,
#                                distance=this_distance,
#                                o11=o11,
#                                o12=o12,
#                                o21=o21,
#                                o22=o22,
#                                t_x=t_x,
#                                t_y=t_y,
#                                t_z=t_z,
#                                wedge=wedge,
#                                chi=chi
#                                )

#     gve = compute_g_vectors(tth, eta,
#                             omega,
#                             wavelength,
#                             wedge=wedge,
#                             chi=chi)

#     return gve


@numba.njit(cache=True)
def _ubi_to_u_safe(ubi, U):
    """
    U from UBI via the metric tensor -- same arithmetic as
    ubi_to_unitcell and ubi_and_ucell_to_u, but it returns False instead of
    raising on a degenerate cell.

    We used to call np.linalg.inv on the metric tensor g = ubi.ubi^T. Since
    worth_fitting only requires cond(ubi) < 1e14, cond(g) can reach 1e28 and
    LAPACK's LU throws LinAlgError. Inside a prange that exception is not
    reliably propagated: the thread stops writing results and every remaining
    voxel it owned is silently left at NaN.
    """
    mt = np.dot(ubi, ubi.T)
    if mt[0, 0] <= 0.0 or mt[1, 1] <= 0.0 or mt[2, 2] <= 0.0:
        return False
    a = np.sqrt(mt[0, 0])
    b = np.sqrt(mt[1, 1])
    c = np.sqrt(mt[2, 2])
    if not (np.isfinite(a) and np.isfinite(b) and np.isfinite(c)):
        return False
    xa = mt[1, 2] / (b * c)
    xb = mt[0, 2] / (a * c)
    xg = mt[0, 1] / (a * b)
    if abs(xa) >= 1.0 or abs(xb) >= 1.0 or abs(xg) >= 1.0:
        return False           # arccos would give NaN
    ca = np.cos(np.radians(np.degrees(np.arccos(xa))))
    cb = np.cos(np.radians(np.degrees(np.arccos(xb))))
    cg = np.cos(np.radians(np.degrees(np.arccos(xg))))

    # det(g) = (abc)^2 * vf, so vf is a dimensionless cell-validity measure
    vf = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg
    if not (vf > 1e-10):
        return False

    g = np.empty((3, 3), dtype=np.float64)
    g[0, 0] = a * a
    g[0, 1] = a * b * cg
    g[0, 2] = a * c * cb
    g[1, 0] = a * b * cg
    g[1, 1] = b * b
    g[1, 2] = b * c * ca
    g[2, 0] = a * c * cb
    g[2, 1] = b * c * ca
    g[2, 2] = c * c

    det = _det3(g)
    gi = np.empty((3, 3), dtype=np.float64)
    gi[0, 0] = (g[1, 1] * g[2, 2] - g[1, 2] * g[2, 1]) / det
    gi[0, 1] = (g[0, 2] * g[2, 1] - g[0, 1] * g[2, 2]) / det
    gi[0, 2] = (g[0, 1] * g[1, 2] - g[0, 2] * g[1, 1]) / det
    gi[1, 0] = (g[1, 2] * g[2, 0] - g[1, 0] * g[2, 2]) / det
    gi[1, 1] = (g[0, 0] * g[2, 2] - g[0, 2] * g[2, 0]) / det
    gi[1, 2] = (g[0, 2] * g[1, 0] - g[0, 0] * g[1, 2]) / det
    gi[2, 0] = (g[1, 0] * g[2, 1] - g[1, 1] * g[2, 0]) / det
    gi[2, 1] = (g[0, 1] * g[2, 0] - g[0, 0] * g[2, 1]) / det
    gi[2, 2] = (g[0, 0] * g[1, 1] - g[0, 1] * g[1, 0]) / det

    if gi[0, 0] <= 0.0 or gi[1, 1] <= 0.0 or gi[2, 2] <= 0.0:
        return False
    astar = np.sqrt(gi[0, 0])
    bstar = np.sqrt(gi[1, 1])
    cstar = np.sqrt(gi[2, 2])
    yb = gi[0, 2] / (astar * cstar)
    yg = gi[0, 1] / (astar * bstar)
    if abs(yb) >= 1.0 or abs(yg) >= 1.0:
        return False
    betas = np.degrees(np.arccos(yb))
    gammas = np.degrees(np.arccos(yg))

    B = np.zeros((3, 3), dtype=np.float64)
    B[0, 0] = astar
    B[0, 1] = bstar * np.cos(np.radians(gammas))
    B[0, 2] = cstar * np.cos(np.radians(betas))
    B[1, 1] = bstar * np.sin(np.radians(gammas))
    B[1, 2] = -cstar * np.sin(np.radians(betas)) * ca
    B[2, 2] = 1.0 / c

    Ut = np.dot(B, ubi).T
    for i in range(3):
        for j in range(3):
            if not np.isfinite(Ut[i, j]):
                return False
            U[i, j] = Ut[i, j]
    return True


def build_partition(refine, omega_binsize="auto", omega_step=None,
                    verbose=True):
    """Partition the peaks into (omega bin, dty bin) cells.
 
    See voxel_mask.VoxelSinoMasker.partition for how the partitioning works.
 
    Path 1: already built this session, held as refine._partition
    Path 2: build it

    This is not saved to disk, because we don't need to memmap like for indexing.
    """
    icolf = refine.icolf
    n = icolf.nrows
 
    # id(icolf) in the key means setpeaks and loadpeaks invalidate this for
    # free: both build the new columnfile while self.icolf still holds the
    # old one, so the addresses cannot collide. addcolumn does NOT
    # invalidate, which is what we want -- get_origins adding xpos_refined
    # leaves omega, dty, dtyi and frm untouched.
    memkey = (omega_binsize, omega_step, id(icolf), icolf.nrows)
    held = getattr(refine, "_partition", None)
    if held is not None and held[0] == memkey:
        if verbose:
            print("reusing the partition built earlier this session")
        idx = dict(held[1])
        return idx
 
    omega = np.ascontiguousarray(icolf.omega, dtype=np.float64)
    dty   = np.ascontiguousarray(icolf.dty,   dtype=np.float64)
    rmax = float(np.hypot(refine.sx_grid, refine.sy_grid)[refine.mask].max())
 
    t0 = time.perf_counter()
    kw, iom = choose_omega_bins(refine, omega, omega_binsize, omega_step,
                                verbose)
    M = VoxelSinoMasker(omega, dty, refine.ystep, ymin=refine.ymin)
    M.partition(keep_columns=False, **kw)
    nom = int(M.sinomega_bins.size)
    ndty = int(M.dty_partitions.shape[1]) - 1
 
    # max |omega - bin centre|, which sets how far a ray can wander inside a
    # cell and hence the voxel-strip half width used by compute_origins.
    # From the unpermuted omega: keep_columns=False means the masker no
    # longer holds a permuted copy.
    ob = np.asarray(M.omega_bins, dtype=np.float64)
    if iom is None:                       # heuristic or explicit binsize
        iom = np.clip(np.round((omega - omega.min()) / M.omega_binsize
                               ).astype(np.int64), 0, nom - 1)
    dev = float(np.abs(omega - ob[iom]).max())
    ray_margin = 1.5 + rmax * np.radians(dev) / refine.ystep + 0.5
 
    if verbose:
        print("peak index: %d omega bins x %d dty bins = %d cells, %.1f "
              "peaks/cell (%.2f s)" % (nom, ndty, nom * ndty,
                                       n / float(nom * ndty),
                                       time.perf_counter() - t0))
        print("  binsize %s, max |omega - bin centre| %.5f deg, r_max %.1f "
              "-> ray half-width %.2f*ystep"
              % ("supplied frames" if M.omega_binsize is None
                 else "%.4f deg" % M.omega_binsize, dev, rmax, ray_margin))
        if n < nom * ndty:
            print("  WARNING: fewer peaks than cells (%.2f/cell); the index "
                  "is sparse and traversal will dominate."
                  % (n / float(nom * ndty)))
 
    idx = {"order": M.peak_ordering,
           "omega_partitions": M.omega_partitions,
           "dty_partitions": M.dty_partitions,
           "usin": M.sinomega_bins, "ucos": M.cosomega_bins,
           "nom": nom, "ndty": ndty, "ymin": float(M.ymin),
           "ray_margin": ray_margin, "dev": dev,
          }
 
    refine._partition = (memkey, dict(idx))
    return idx



def build_refine_inputs(refine, omega_binsize="auto", omega_step=None,
                        verbose=True):
    """Everything refine_map needs, on top of the partition.
 
    Returns {"partition": ..., "peaks": ..., "pstart": ..., "porder": ...,
             "pbpmap_ubis": ..., "nrj": ..., "gvec_pars": ...}
    """
    icolf = refine.icolf
    part = build_partition(refine, omega_binsize=omega_binsize,
                           omega_step=omega_step, verbose=verbose)
    order = part["order"]
 
    def col(name, dtype=np.float64):
        return np.ascontiguousarray(getattr(icolf, name)[order], dtype=dtype)
 
    # We only permute the columns we need for the g-vector computation (with
    # corrected xpos) and for the merge inside refine_map.
    #
    # xl, yl, zl replace sc, fc: they are the peak positions in the lab frame,
    # which is what the g-vector needs. compute_xyz_lab is affine in (sc, fc),
    # so the intensity-weighted mean of xl/yl/zl over a merge group equals
    # xl/yl/zl of the mean sc, fc -- refine_map can therefore average these
    # directly and never touch the detector transform.
    peaks_permuted = {
        "eta": col("eta"),
        "sum_intensity": col("sum_intensity"),
        "sinomega": col("sinomega"), "cosomega": col("cosomega"),
        "omega": col("omega"), "dty": col("dty"),
        # int64, not the float64 default: _pack5 builds keys up to 2**58 and
        # float64 has 53 bits of mantissa, so a float dtyi would collapse
        # every dty distinction in the merge
        "dtyi": col("dtyi", np.int64),
        "xpos": col("xpos_refined"),
        "xl": col("xl"), "yl": col("yl"), "zl": col("zl"),
    }
 
    # Remember these are permuted!
    pars = icolf.parameters.get_parameters()
    gpars = gvec_pars(pars)
 
    t0 = time.perf_counter()
    peaks_permuted["gve_all"] = gvectors_from_lab(      # (N, 3)
        peaks_permuted["xl"], peaks_permuted["yl"], peaks_permuted["zl"],
        peaks_permuted["xpos"],
        peaks_permuted["sinomega"], peaks_permuted["cosomega"],
        gpars[0], gpars[1], gpars[2], gpars[3],
        16 * max(numba.get_num_threads(), 1))
    if verbose:
        print("computed g-vectors in %.2f s" % (time.perf_counter() - t0))
 
    # ---- pbpmap CSR over (ri, rj), and (M, 3, 3) UBIs --------------------
    ri_col, rj_col = geometry.step_to_recon(refine.pbpmap.i, refine.pbpmap.j,
                                            refine.mask.shape)
    nri, nrj = refine.mask.shape
    pcells = (np.asarray(ri_col, dtype=np.int64) * nrj
              + np.asarray(rj_col, dtype=np.int64))
    pstart, porder = build_csr(pcells, nri * nrj)
 
    return {
        "partition": part,
        "peaks": peaks_permuted,
        "pstart": pstart,
        "porder": porder,
        "pbpmap_ubis": np.ascontiguousarray(
            refine.pbpmap.ubi.transpose(2, 0, 1)),
        "nrj": nrj,
        "gvec_pars": gpars,
    }




def grid_axes(sx_grid, sy_grid):
    """
    sx_grid/sy_grid come from np.meshgrid(sx, sy, indexing='ij'), so they are
    separable: sx_grid[i, j] == sx[i] and sy_grid[i, j] == sy[j]. The strip
    walk in compute_origins relies on that, so check it rather than assume.
    """
    sx_ax = np.ascontiguousarray(sx_grid[:, 0], dtype=np.float64)
    sy_ax = np.ascontiguousarray(sy_grid[0, :], dtype=np.float64)
    if not (np.allclose(sx_grid, sx_ax[:, None])
            and np.allclose(sy_grid, sy_ax[None, :])):
        raise ValueError("sx_grid/sy_grid are not a separable meshgrid")
    if len(sx_ax) < 2 or len(sy_ax) < 2:
        raise ValueError("grid needs at least 2 points per axis")
    dsx = sx_ax[1] - sx_ax[0]
    dsy = sy_ax[1] - sy_ax[0]
    if not (np.allclose(np.diff(sx_ax), dsx) and np.allclose(np.diff(sy_ax), dsy)):
        raise ValueError("sx_grid/sy_grid are not uniformly spaced")
    return sx_ax, sy_ax, float(dsx), float(dsy)


@numba.njit(cache=True, parallel=True)
def compute_origins(singlemap, sample_mask,
                    gve, sinomega, cosomega, dty,
                    sx_ax, sy_ax, dsx, dsy,
                    y0, ystep, hkl_tol, weight_reg,
                    omega_partitions, dty_partitions, nom, ndty,
                    usin, ucos, ymin, ray_margin,
                    max_cell, nchunks):

    """
    Compute the origin of diffraction for each peak.

    For each (omega frame, dtyi) cell we cast a ray through singlemap; where a
    voxel's orientation indexes the peak we accumulate its position, weighted
    by 1/(ydist+1). All arrays are in CSR (permuted) order, so cells are
    contiguous runs; the caller un-permutes the result.

    Three changes from the version that grouped with lexsort + np.split:

      - cells come from the CSR index instead of unique_with_counts and
        np.split, which built reflected Python lists per omega group

      - the ray is a line, so instead of building full-map ydist / ray_mask /
        final_mask / weights arrays per cell (~1.5 MB of allocation each, and
        above the mmap threshold) we solve the strip analytically: for fixed i
        the set of j with |dty - ygrid| <= ystep is a contiguous interval.
        Drive along whichever axis keeps the strip narrow.

      - the inner test is fused. orien.dot(g), np.round and the sum of squares
        allocated three small arrays for every (peak, voxel) pair.

    The CSR cell is only a pre-filter: the exact predicate
    |dty_p - (y0 - sx sin(w_p) - sy cos(w_p))| <= ystep is still applied per
    peak, with that peak's own sinomega/cosomega.
    """
    npk = sinomega.shape[0]
    NI = sx_ax.shape[0]
    NJ = sy_ax.shape[0]
    lx_modified = np.zeros(npk, dtype=np.float64)
    W = ystep * ray_margin

    # the very original function was checking tolsq. this was a bug.
    tolsq = hkl_tol * hkl_tol

    for chunk in numba.prange(nchunks):
        accx = np.zeros(max_cell, dtype=np.float64)
        accy = np.zeros(max_cell, dtype=np.float64)
        accw = np.zeros(max_cell, dtype=np.float64)

        for io in range(chunk, nom, nchunks):
            so = usin[io]
            co = ucos[io]
            drive_i = abs(co) >= abs(so)

            obase = omega_partitions[io]
            for k in range(ndty):
                plo = obase + dty_partitions[io, k]
                phi = obase + dty_partitions[io, k + 1]
                npc = phi - plo
                if npc == 0:
                    continue
                for t in range(npc):
                    accx[t] = 0.0
                    accy[t] = 0.0
                    accw[t] = 0.0

                dty_c = k * ystep + ymin      # nominal dty for this cell

                if drive_i:
                    for i in range(NI):
                        A = y0 - sx_ax[i] * so - dty_c
                        # |A + sy*co| <= W  ->  sy between the two roots
                        ja = (A - W) / co
                        jb = (A + W) / co
                        ta = (ja - sy_ax[0]) / dsy
                        tb = (jb - sy_ax[0]) / dsy
                        if ta > tb:
                            tmp = ta; ta = tb; tb = tmp
                        jlo = int(np.floor(ta))
                        jhi = int(np.ceil(tb))
                        if jlo < 0:
                            jlo = 0
                        if jhi > NJ - 1:
                            jhi = NJ - 1
                        for j in range(jlo, jhi + 1):
                            if not sample_mask[i, j]:
                                continue
                            sxv = sx_ax[i]
                            syv = sy_ax[j]
                            # (UB)^-1 for this voxel, from the single-valued map
                            u00 = singlemap[i, j, 0, 0]
                            if np.isnan(u00):
                                continue
                            u01 = singlemap[i, j, 0, 1]
                            u02 = singlemap[i, j, 0, 2]
                            u10 = singlemap[i, j, 1, 0]
                            u11 = singlemap[i, j, 1, 1]
                            u12 = singlemap[i, j, 1, 2]
                            u20 = singlemap[i, j, 2, 0]
                            u21 = singlemap[i, j, 2, 1]
                            u22 = singlemap[i, j, 2, 2]
                            for q in range(plo, phi):
                                # y + dy of eq (6): beam to voxel centroid
                                yd = abs(y0 - sxv * sinomega[q]
                                         - syv * cosomega[q] - dty[q])
                                if yd > ystep:
                                    continue
                                Gx = gve[q, 0]              # G, eq (5)
                                Gy = gve[q, 1]
                                Gz = gve[q, 2]
                                hf0 = u00 * Gx + u01 * Gy + u02 * Gz
                                hf1 = u10 * Gx + u11 * Gy + u12 * Gz
                                hf2 = u20 * Gx + u21 * Gy + u22 * Gz
                                dh0 = hf0 - np.round(hf0)   # eq (8) residual
                                dh1 = hf1 - np.round(hf1)
                                dh2 = hf2 - np.round(hf2)
                                if dh0 * dh0 + dh1 * dh1 + dh2 * dh2 < tolsq:
                                    # eq (6)
                                    w = weight_reg / (yd + weight_reg)
                                    accx[q - plo] += sxv * w
                                    accy[q - plo] += syv * w
                                    accw[q - plo] += w

                else:
                    for j in range(NJ):
                        B = y0 - sy_ax[j] * co - dty_c
                        ia = (B - W) / so
                        ib = (B + W) / so
                        ta = (ia - sx_ax[0]) / dsx
                        tb = (ib - sx_ax[0]) / dsx
                        if ta > tb:
                            tmp = ta; ta = tb; tb = tmp
                        ilo = int(np.floor(ta))
                        ihi = int(np.ceil(tb))
                        if ilo < 0:
                            ilo = 0
                        if ihi > NI - 1:
                            ihi = NI - 1
                        for i in range(ilo, ihi + 1):
                            if not sample_mask[i, j]:
                                continue
                            sxv = sx_ax[i]
                            syv = sy_ax[j]
                            # (UB)^-1 for this voxel, from the single-valued map
                            u00 = singlemap[i, j, 0, 0]
                            if np.isnan(u00):
                                continue
                            u01 = singlemap[i, j, 0, 1]
                            u02 = singlemap[i, j, 0, 2]
                            u10 = singlemap[i, j, 1, 0]
                            u11 = singlemap[i, j, 1, 1]
                            u12 = singlemap[i, j, 1, 2]
                            u20 = singlemap[i, j, 2, 0]
                            u21 = singlemap[i, j, 2, 1]
                            u22 = singlemap[i, j, 2, 2]
                            for q in range(plo, phi):
                                # y + dy of eq (6): beam to voxel centroid
                                yd = abs(y0 - sxv * sinomega[q]
                                         - syv * cosomega[q] - dty[q])
                                if yd > ystep:
                                    continue
                                Gx = gve[q, 0]              # G, eq (5)
                                Gy = gve[q, 1]
                                Gz = gve[q, 2]
                                hf0 = u00 * Gx + u01 * Gy + u02 * Gz
                                hf1 = u10 * Gx + u11 * Gy + u12 * Gz
                                hf2 = u20 * Gx + u21 * Gy + u22 * Gz
                                dh0 = hf0 - np.round(hf0)   # eq (8) residual
                                dh1 = hf1 - np.round(hf1)
                                dh2 = hf2 - np.round(hf2)
                                if dh0 * dh0 + dh1 * dh1 + dh2 * dh2 < tolsq:
                                    # eq (6)
                                    w = weight_reg / (yd + weight_reg)
                                    accx[q - plo] += sxv * w
                                    accy[q - plo] += syv * w
                                    accw[q - plo] += w

                for t in range(npc):
                    if accw[t] > 0.0:
                        q = plo + t
                        sxc = accx[t] / accw[t]
                        syc = accy[t] / accw[t]
                        lx_modified[q] = sxc * cosomega[q] - syc * sinomega[q]

    return lx_modified


@numba.njit(cache=True, parallel=True)
def count_ray_voxels(sample_mask, sx_ax, sy_ax, dsx, dsy, y0, ystep,
                     dty_partitions, nom, ndty,
                     usin, ucos, ymin, ray_margin):
    """Diagnostic: mean voxels visited per non-empty cell, and cell occupancy."""
    NI = sx_ax.shape[0]
    NJ = sy_ax.shape[0]
    W = ystep * ray_margin
    tot = np.zeros(nom, dtype=np.int64)
    ncell = np.zeros(nom, dtype=np.int64)
    for io in numba.prange(nom):
        so = usin[io]
        co = ucos[io]
        drive_i = abs(co) >= abs(so)
        for k in range(ndty):
            if dty_partitions[io, k + 1] == dty_partitions[io, k]:
                continue
            ncell[io] += 1
            dty_c = k * ystep + ymin
            n = 0
            if drive_i:
                for i in range(NI):
                    A = y0 - sx_ax[i] * so - dty_c
                    ta = ((A - W) / co - sy_ax[0]) / dsy
                    tb = ((A + W) / co - sy_ax[0]) / dsy
                    if ta > tb:
                        tmp = ta; ta = tb; tb = tmp
                    jlo = max(int(np.floor(ta)), 0)
                    jhi = min(int(np.ceil(tb)), NJ - 1)
                    for j in range(jlo, jhi + 1):
                        if sample_mask[i, j]:
                            n += 1
            else:
                for j in range(NJ):
                    B = y0 - sy_ax[j] * co - dty_c
                    ta = ((B - W) / so - sx_ax[0]) / dsx
                    tb = ((B + W) / so - sx_ax[0]) / dsx
                    if ta > tb:
                        tmp = ta; ta = tb; tb = tmp
                    ilo = max(int(np.floor(ta)), 0)
                    ihi = min(int(np.ceil(tb)), NI - 1)
                    for i in range(ilo, ihi + 1):
                        if sample_mask[i, j]:
                            n += 1
            tot[io] += n
    return tot.sum(), ncell.sum()


class PBPRefine:
    """Class to manage point-by-point refinement.
    The maps in this class are set up on a grid that matches the PBPMap reference frame
    This is in step space in geometry.py"""

    def __init__(self,
                 dset,
                 phase_name,
                 hkl_tol_origins=0.05,  # we forgot to square this before!
                 hkl_tol_refine=0.1,
                 hkl_tol_refine_merged=0.05,
                 ds_tol=0.005,
                 etacut=0.1,
                 ifrac=None,
                 forref=None,
                 y0=0.0,
                 min_grain_npks=6,
                 beam_size=None,
                 ):
        self.dset = dset
        self.phase_name = phase_name
        # peak selection parameters
        self.forref = forref
        self.ds_tol = ds_tol

        if ifrac is None:
            self.ifrac = 1.0 / len(self.dset.ybincens)
        else:
            self.ifrac = ifrac

        self.etacut = etacut
        # compute_origins parameters
        self.hkl_tol_origins = hkl_tol_origins
        # refinement parameters
        self.hkl_tol_refine = hkl_tol_refine
        self.hkl_tol_refine_merged = hkl_tol_refine_merged
        self.min_grain_npks = min_grain_npks
        # geometry stuff
        self.ystep = self.dset.ystep
        self.y0 = y0
        self.ybincens = self.dset.ybincens
        self.ymin = self.ybincens.min()

        # X-ray beam size, in the same units as the dty motor.
        #
        # The distance weighting in both the UBI and the strain fit uses
        # beam_size/3 as its regularisation length: w = r / (ydist + r).
        # Henningsson et al. eq (18) writes that as 1 um, for a 3 um beam.
        #
        # This is the only absolute length in the refinement. Everything else
        # (the dty window, the ray strip, the omega binsize heuristic) is a
        # ratio and cancels, and the g-vectors are in 1/wavelength regardless.
        # Left as a bare 1.0 it means "1 um" whatever the motor units, so on a
        # mm dataset the weighting goes flat and every peak on the ray counts
        # the same however far the beam passed from the voxel centre.
        #
        # Defaults to ystep, i.e. assumes the scan step matches the beam. Set
        # it explicitly if you know better -- it is a physical property of the
        # beamline, not of the scan.
        self.beam_size = self.ystep if beam_size is None else float(beam_size)

        # set default paths
        self.pbpmap_filename = self.dset.refmapfile  # pbp map input
        self.icolf_filename = self.dset.refpeaksfile  # icolf for refinement
        self.refinedmap_filename = self.dset.refoutfile  # refined pbp map output
        self.own_filename = self.dset.refmanfile  # myself as an H5

    def setmap(self, pbpmap):
        """Set an input PBPMap Python object to use"""
        self.pbpmap = pbpmap
        # get arrays of integer translations in step space
        # e.g -33, -32, -31, ... 32, 33
        si = np.arange(self.pbpmap.i.min(), self.pbpmap.i.max() + 1)
        sj = np.arange(self.pbpmap.j.min(), self.pbpmap.j.max() + 1)

        # convert these arrays to sample space
        sx, sy = geometry.step_to_sample(si, sj, self.ystep)

        # make a grid of these
        self.sx_grid, self.sy_grid = np.meshgrid(sx, sy, indexing='ij')

    def loadmap(self, filename=None, refined=False):
        """Load an existing input/refined PBPMap from h5. If you only have a .txt file, you want self.setmap()"""
        if filename is None:
            if refined:
                filename = self.refinedmap_filename
            else:
                filename = self.pbpmap_filename
        pmap = PBPMap(new=True)
        if not os.path.exists(filename):
            raise OSError("Can't find map on disk!")
        if refined:
            name = 'refinedmap'
        else:
            name = 'pbpmap'
        pmap = ImageD11.columnfile.colfile_from_hdf(filename, obj=pmap, name=name)
        if refined:
            self.refinedmap = pmap
            self.refinedmap_filename = filename
        else:
            self.setmap(pmap)
            self.pbpmap_filename = filename

    def savemap(self, filename=None, refined=False):
        """Save self.pbpmap/self.refinedmap to an H5 file"""
        if filename is None:
            if refined:
                filename = self.refinedmap_filename
            else:
                filename = self.pbpmap_filename
        if refined:
            name = 'refinedmap'
            obj = self.refinedmap
        else:
            name = 'pbpmap'
            obj = self.pbpmap
        ImageD11.columnfile.colfile_to_hdf(obj, filename, compression=None, name=name)
        if refined:
            self.refinedmap_filename = filename
        else:
            self.pbpmap_filename = filename

    def setpeaks(self, colf, keep_colf=True):
        """Set icolf to refine against.
        
        keep_colf: keep a reference to colf as self.colf (prevents `del colf` calls from notebook from working)
         """

        # Set the parameters for the peaks
        self.dset.update_colfile_pars(colf, phase_name=self.phase_name)

        uc = unitcell.unitcell_from_parameters(colf.parameters)
        uc.makerings(colf.ds.max(), self.ds_tol)
        if self.forref is None:
            self.forref = range(len(uc.ringds))

        sel, npks, hmax = mask_rings_by_ifrac(
            colf, self.ds_tol, colf.ds.max(), self.ifrac, uc,
            forref=self.forref, verbose=True, return_npks_hmax=True)
        
        # peaks we refine with
        isel = sel & (np.abs(np.sin(np.radians(colf.eta))) > self.etacut)

        dtyi = geometry.dty_to_dtyi(colf.dty, self.ystep, self.ybincens.min())
        colf.addcolumn(dtyi, "dtyi")

        # cache these to speed up selections later
        colf.addcolumn(np.sin(np.radians(colf.omega)), "sinomega")
        colf.addcolumn(np.cos(np.radians(colf.omega)), "cosomega")

        if keep_colf:
            self.colf = colf
        else:
            # iplot only needs these two, and only to draw. Keep a decimated
            # copy so the plot survives; the full colf is often tens of GB.
            step = max(1, colf.nrows // 2_000_000)
            self._plot_ds = np.array(colf.ds[::step])
            self._plot_I = np.array(colf.sum_intensity[::step])
            self.colf = None
        self.icolf = colf.copyrows(isel)
        
        self.npks = npks
        print(
            "Using for refinement:",
            self.icolf.nrows,
            "npks, forref",
            self.npks,
            self.forref,
        )
        self.uc = uc

        # unlike pbp index, partitioning happens in get_origins

    def savepeaks(self, icolf_filename=None):
        if icolf_filename is None:
            icolf_filename = self.icolf_filename
        else:
            self.icolf_filename = icolf_filename
        ImageD11.columnfile.colfile_to_hdf(self.icolf, icolf_filename, compression=None)

    def loadpeaks(self, icolf_filename=None):
        """Load icolf from disk"""
        if icolf_filename is None:
            icolf_filename = self.icolf_filename
        else:
            self.icolf_filename = icolf_filename
        self.icolf = ImageD11.columnfile.colfile_from_hdf(icolf_filename)
        # update pars
        self.dset.update_colfile_pars(self.icolf, phase_name=self.phase_name)

    def iplot(self, skip=1):
        """Same as PBP.iplot for now"""
        import pylab as pl

        if self.colf is not None:
            bg_ds, bg_I = self.colf.ds[::skip], self.colf.sum_intensity[::skip]
        else:
            bg_ds, bg_I = self._plot_ds, self._plot_I

        f, ax = pl.subplots(2, 1)
        ax[0].plot(bg_ds, bg_I, ",")
        ax[0].plot(self.icolf.ds[::skip], self.icolf.sum_intensity[::skip], ",")
        ax[0].set(yscale="log", ylabel="sum intensity", xlabel=r'$d^{*}~(\AA^{-1})$')
        histo = self.dset.sinohist(
            omega=self.icolf.omega, dty=self.icolf.dty, weights=self.icolf.sum_intensity
        )
        f.colorbar(
            ax[1].pcolormesh(
                self.dset.obinedges,
                self.dset.ybinedges,
                histo.T,
                norm=pl.matplotlib.colors.LogNorm(),
            ),
            ax=ax[1],
        )
        ax[1].set(ylabel="dty", xlabel=r'$\omega~(\degree)$')
        ax[0].vlines(self.uc.ringds, 1e4, 3e4, color='red')
        return f, ax

    def setmask(self, manual_threshold=None, doplot=False, use_icolf=True, use_singlemap=False):
        """Set a mask for choosing what to refine or not."""
        if use_singlemap:
            whole_sample_mask = ~np.isnan(self.singlemap[:, :, 0, 0])
            recon_man_mask = whole_sample_mask.astype(float)
            binary = recon_man_mask
            chull = recon_man_mask
        else:
            if use_icolf:
                dty = self.icolf.dty
                omega = self.icolf.omega
            else:
                dty = self.colf.dty
                omega = self.colf.omega
            whole_sample_sino, xedges, yedges = np.histogram2d(
                dty, omega, bins=[self.dset.ybinedges, self.dset.obinedges])
            shift, _ = geometry.sino_shift_and_pad(
                self.y0, len(self.ybincens), self.ybincens.min(), self.ystep)
            nthreads = cImageD11.cores_available()
            pad = self.sx_grid.shape[0] - whole_sample_sino.shape[0]
            whole_sample_recon = run_iradon(whole_sample_sino, self.dset.obincens,
                                            pad, shift, workers=nthreads)
            recon_man_mask = whole_sample_recon
            if manual_threshold is None:
                thresh = threshold_otsu(recon_man_mask)
            else:
                thresh = manual_threshold
            binary = recon_man_mask > thresh
            chull = convex_hull_image(binary)
            whole_sample_mask = chull

        if doplot:
            from matplotlib import pyplot as plt
            fig, axs = plt.subplots(1, 3, sharex=True, sharey=True,
                                    constrained_layout=True)
            axs[0].imshow(recon_man_mask, vmin=0, origin="lower")
            axs[1].imshow(binary, origin="lower")
            axs[2].imshow(chull, origin="lower")
            axs[0].set_title("Reconstruction")
            axs[1].set_title("Binarised threshold")
            axs[2].set_title("Convex hull")
            fig.supxlabel("<-- Y axis")
            fig.supylabel("Beam >")
            plt.show()

        self.mask = whole_sample_mask

    def setsingle(self, pbpmap, method='best', **method_args):
        """Function to go from a multi-valued PBPMap to a single-valued PBPMap"""
        if method == 'best':
            pbpmap.choose_best(**method_args)
            self.singlemap = pbpmap.best_ubi
        else:
            raise ValueError('Unsupported method!')

    def to_h5(self, filename=None, h5group='PBPRefine', save_peaks=True):
        # save to h5 file
        if filename is None:
            filename = self.own_filename
        else:
            self.own_filename = filename

        if save_peaks:
            print('Saving icolf to disk')
            self.savepeaks()
        print('Saving input map to disk')
        self.savemap()
        try:
            print('Saving output map to disk')
            self.savemap(refined=True)
        except AttributeError:
            print("Couldn't find a refined map to save, continuing")
            pass

        print('Saving myself to disk')
        with h5py.File(filename, "a") as hout:
            parent_group = hout.require_group(h5group)

            array_mappings = {
                'sx_grid': 'sx_grid',
                'sy_grid': 'sy_grid',
                'mask': 'mask',
                'singlemap': 'singlemap',
            }
            for array_name, attr_name in array_mappings.items():
                try:
                    save_array(parent_group, array_name, getattr(self, attr_name))
                except AttributeError:
                    continue

            parent_group.attrs['dsetfile'] = self.dset.dsfile

            filename_mappings = {
                'pbpmapfile': 'pbpmap_filename',
                'icolffile': 'icolf_filename',
                'refmapfile': 'refinedmap_filename',
            }
            for path_name, attr_name in filename_mappings.items():
                try:
                    parent_group.attrs[path_name] = getattr(self, attr_name)
                except AttributeError:
                    continue

            pars = ['phase_name', 'hkl_tol_origins', 'hkl_tol_refine',
                    'hkl_tol_refine_merged', 'ds_tol', 'etacut', 'ifrac',
                    'forref', 'y0', 'min_grain_npks','beam_size']
            for par in pars:
                try:
                    if getattr(self, par) is not None:
                        parent_group.attrs[par] = getattr(self, par)
                except AttributeError:
                    continue

    @classmethod
    def from_h5(cls, filename, h5group='PBPRefine',
                load_peaks=True, load_input=True, load_output=True):
        manager_filename = filename

        with h5py.File(filename, "r") as hin:
            parent_group = hin[h5group]

            arrays = {}
            for array_name in parent_group.keys():
                arrays[array_name] = parent_group[array_name][:]

            dsfile = parent_group.attrs['dsetfile']
            filename_mappings = {
                'pbpmapfile': 'pbpmap_filename',
                'icolffile': 'icolf_filename',
                'refmapfile': 'refinedmap_filename',
            }
            filenames = {}
            for attrkey, attrname in filename_mappings.items():
                try:
                    filenames[attrname] = parent_group.attrs[attrkey]
                except (AttributeError, KeyError):
                    continue

            pars = ['phase_name', 'hkl_tol_origins', 'hkl_tol_refine',
                    'hkl_tol_refine_merged', 'ds_tol', 'etacut', 'ifrac',
                    'forref', 'y0', 'min_grain_npks','beam_size']
            pars_dict = {}
            for par in pars:
                try:
                    pars_dict[par] = parent_group.attrs[par]
                except (AttributeError, KeyError):
                    continue

        dset = ImageD11.sinograms.dataset.load(dsfile)
        refine_obj = cls(dset=dset, **pars_dict)
        refine_obj.own_filename = manager_filename
        for filename_attr, filename in filenames.items():
            setattr(refine_obj, filename_attr, filename)
        for array_attr, array in arrays.items():
            setattr(refine_obj, array_attr, array)

        if load_peaks:
            print('Loading peaks')
            try:
                refine_obj.loadpeaks(refine_obj.icolf_filename)
            except AttributeError:
                pass

        if load_input:
            print('Loading input map')
            try:
                refine_obj.loadmap(refine_obj.pbpmap_filename)
            except AttributeError:
                pass

        if load_output:
            print('Loading output map')
            try:
                refine_obj.loadmap(refine_obj.refinedmap_filename, refined=True)
            except (AttributeError, OSError):
                pass
        return refine_obj

    def get_origins(self, nthreads=None, nchunks=None,
                    omega_binsize="auto", omega_step=None,
                    verbose=True):
        """Fill the xpos_refined column."""
        if 'xpos_refined' in self.icolf.titles:
            raise ValueError('We already have origins in self.icolf! Not recomputing')
 
        if nthreads is None:
            nthreads = max(cImageD11.cores_available() - 1, 1)
        numba.set_num_threads(nthreads)
        try:
            numba.set_parallel_chunksize(1)
        except AttributeError:
            pass
        
        # Build the partition. Doesn't permute the peaks!
        part = build_partition(self, omega_binsize=omega_binsize,
                               omega_step=omega_step, verbose=verbose)
        order = part["order"]
 
        # Permute un-corrected g-vectors into partition order
        # also sinomega, cosomega, dty
 
        gve = np.empty((self.icolf.nrows, 3), dtype=np.float64)
        for c, name in enumerate(("gx", "gy", "gz")):
            np.take(getattr(self.icolf, name), order, out=gve[:, c])
        sinomega = np.ascontiguousarray(self.icolf.sinomega[order])
        cosomega = np.ascontiguousarray(self.icolf.cosomega[order])
        dty = np.ascontiguousarray(self.icolf.dty[order])
        # everything we need for compute_origins is now in partition order
 
        sx_ax, sy_ax, dsx, dsy = grid_axes(self.sx_grid, self.sy_grid)
        singlemap = np.ascontiguousarray(self.singlemap, dtype=np.float64)
        mask = np.ascontiguousarray(self.mask).astype(np.bool_)
 
        max_cell = int(np.diff(part["dty_partitions"], axis=1).max())
        if nchunks is None:
            nchunks = nthreads
        nchunks = max(1, min(nchunks, part["nom"]))
 
        if verbose:
            nvox, ncells = count_ray_voxels(
                mask, sx_ax, sy_ax, dsx, dsy, self.y0, self.ystep,
                part["dty_partitions"],
                part["nom"], part["ndty"], part["usin"], part["ucos"],
                part["ymin"], part["ray_margin"])
            print("%d non-empty cells, max %d peaks in a cell" % (ncells, max_cell))
            print("mean %.0f voxels visited per ray (full map would be %d)"
                  % (nvox / max(ncells, 1), mask.size))
            print("running on %d threads" % nthreads)
 
        weight_reg = self.beam_size / 3.0
        if verbose:
            print("beam size %.4g, origin weight regularisation %.4g"
                  % (self.beam_size, weight_reg))
 
        t0 = time.perf_counter()
        lx_perm = compute_origins(
            singlemap, mask, gve, sinomega, cosomega, dty,
            sx_ax, sy_ax, dsx, dsy,
            self.y0, self.ystep, self.hkl_tol_origins, weight_reg,
            part["omega_partitions"], part["dty_partitions"],
            part["nom"], part["ndty"], part["usin"], part["ucos"],
            part["ymin"], part["ray_margin"],
            max_cell, nchunks)
        if verbose:
            print("origins took %.1f s" % (time.perf_counter() - t0))
 
        # un-permute back to icolf row order
        lx_modified = np.empty_like(lx_perm)
        lx_modified[order] = lx_perm
 
        self.icolf.addcolumn(lx_modified, 'xpos_refined')
        print('xpos_refined column added to self.icolf')
 
        return lx_modified

    def run_refine(self, points_step_space=None, npoints=None,
                   output_filename=None, use_cluster=False, pythonpath=None,
                   nthreads=None, nchunks=None, omega_binsize="auto",
                   omega_step=None, 
                   save=True, verbose=True,
                   cpus_per_task=16, time_h=1, partition="nice", mem_G=32):
        """
        Refine every UBI candidate in self.pbpmap.
 
        Computes the diffraction origins first if they are not already there,
        so on the cluster the worker does that too -- get_origins is the
        expensive part and there is no reason to run it on a submission node.
        Call get_origins() yourself first if you want to inspect them.
        """
        if use_cluster:
            if pythonpath is None:
                raise ValueError("Must supply pythonpath to run on cluster!")
            if points_step_space is not None or npoints is not None:
                raise ValueError("Choosing points to refine on the cluster is "
                                 "not implemented")
            # The worker does from_h5 -> loadpeaks, so the peaks and the maps
            # have to be on disk.
 
            print("Saving peaks and maps for the worker")
            self.to_h5()
            print("Making bash script")
            bash_script_path = prepare_refine_bash(
                self, pythonpath, output_filename=output_filename,
                cpus_per_task=cpus_per_task, time_h=time_h,
                partition=partition, mem_G=mem_G)
            print("Made bash script at", bash_script_path)
            print("Submitting to cluster")
            result = subprocess.run("sbatch {}".format(bash_script_path),
                                    capture_output=True, shell=True)
            print(result.stdout.decode("utf-8"))
            if result.returncode != 0:
                print(result.stderr.decode("utf-8"))
            return None
 
        # ---- origins, if we do not have them yet --------------------------
        if "xpos_refined" not in self.icolf.titles:
            if verbose:
                print("no xpos_refined column -- computing origins first")
            self.get_origins(nthreads=nthreads, nchunks=nchunks,
                             omega_binsize=omega_binsize,
                             omega_step=omega_step, verbose=verbose)
 
        if nthreads is None:
            nthreads = max(cImageD11.cores_available() - 1, 1)
 
        numba.set_num_threads(nthreads)
        try:
            numba.set_parallel_chunksize(1)
        except AttributeError:
            pass
        
        # If there's a good cache, build_refine_inputs will also get it from disk
        inputs = build_refine_inputs(self, omega_binsize=omega_binsize,
                                     omega_step=omega_step, verbose=verbose)
        part = inputs["partition"]
        pk = inputs["peaks"]            # permuted by the partition scheme
 
        # voxels to refine
        if points_step_space is None:
            pts = np.array(np.nonzero(self.mask)).T
            if npoints is not None:
                pts = pts[:npoints]
        else:
            pts = np.array([geometry.step_to_recon(si, sj, self.mask.shape)
                            for (si, sj) in points_step_space])
        points_ri = np.ascontiguousarray(pts[:, 0], dtype=np.int64)
        points_rj = np.ascontiguousarray(pts[:, 1], dtype=np.int64)
 
        # Points are divided into chunks based on how many threads we use
        # Each chunk is persistent so we don't malloc
        if nchunks is None:
            nchunks = nthreads
        nchunks = max(1, min(nchunks, len(points_ri)))
 
        # size the buffers
        maxlocal = int(vm_max_candidates(
            np.ascontiguousarray(self.sx_grid[self.mask], dtype=np.float64),
            np.ascontiguousarray(self.sy_grid[self.mask], dtype=np.float64),
            self.y0, self.ystep, part["ymin"],
            part["dty_partitions"],
            part["usin"], part["ucos"]))
        maxlocal = max(maxlocal, 1)
        if verbose:
            mb = maxlocal * 200 * nchunks / 1e6
            print("worst-case peaks on a single ray: %d" % maxlocal)
            print("scratch: ~%.0f MB across %d chunks (virtual; touched lazily)" % (mb, nchunks))
            print("launching numba parallel refinement on %d threads" % nthreads)
 
        gpk = inputs["gvec_pars"]        # (A, oc, 1/wavelength, a_is_eye)
        uc = unitcell.unitcell_from_parameters(self.icolf.parameters)
        B0 = unitcell_to_b(uc.lattice_parameters, np.eye(3))
 
        weight_reg = self.beam_size / 3.0
        if verbose:
            print("beam size %.4g, fit weight regularisation %.4g "
                  "(ystep %.4g)" % (self.beam_size, weight_reg, self.ystep))
 
        t0 = time.perf_counter()
        (ubis_m, eps_m, npks, nuniq,
         overflow, degenerate, gather_fail) = refine_map(
            points_ri, points_rj,
             # pbpmap CSR:
            inputs["pstart"], inputs["porder"], inputs["pbpmap_ubis"], inputs["nrj"],
            self.sx_grid, self.sy_grid, self.mask,
             # partition stuff:
            part["omega_partitions"], part["dty_partitions"],
            part["usin"], part["ucos"],
            part["ymin"],
            # permuted peaks:
            pk["xl"], pk["yl"], pk["zl"], pk["eta"], pk["sum_intensity"],
            pk["sinomega"], pk["cosomega"], pk["omega"], pk["dty"], pk["dtyi"],
            pk["xpos"], pk["gve_all"],
            self.ystep, self.y0, B0,
            gpk[0], gpk[1], gpk[2], gpk[3],
            float(self.hkl_tol_refine), float(self.hkl_tol_refine_merged),
            int(self.min_grain_npks), 1e-4, weight_reg,
            maxlocal, nchunks,
        )
        if verbose:
            print("refinement took %.1f s" % (time.perf_counter() - t0))
        if gather_fail.sum():
            print("ERROR: the peak gather overflowed its buffer at %d voxels. "
                  "maxlocal was wrong and those voxels were skipped."
                  % int(gather_fail.sum()))
        if degenerate.sum():
            print("note: %d UBIs gave a degenerate metric tensor; strain skipped "
                  "for those voxels (UBI/npks/nuniq still written)."
                  % int(degenerate.sum()))
        if overflow.sum():
            print("WARNING: %d peak assignments had |hkl| > 2047 and were clamped "
                  "in the merge key. That means some candidate UBIs are nonsense; "
                  "the affected voxels may merge peaks incorrectly."
                  % int(overflow.sum()))
 
        # back to the (3, 3, M) layout the rest of ImageD11 expects
        final_ubis = np.ascontiguousarray(ubis_m.transpose(1, 2, 0))
        final_eps = np.ascontiguousarray(eps_m.transpose(1, 2, 0))
 
        final_ubis[:, :, np.isnan(final_ubis[0, 0, :])] = np.eye(3)[..., np.newaxis]
        final_eps[:, :, np.isnan(final_eps[0, 0, :])] = 0
        npks[np.isnan(npks)] = 0
        nuniq[np.isnan(nuniq)] = 0
 
        output_map = PBPMap(new=True)
        output_map.nrows = final_ubis.shape[2]
        ub = final_ubis.reshape(9, final_ubis.shape[2])
        ep = final_eps.reshape(9, final_eps.shape[2])
        output_map.addcolumn(self.pbpmap.i, "i")
        output_map.addcolumn(self.pbpmap.j, "j")
        output_map.addcolumn(npks, "ntotal")
        output_map.addcolumn(nuniq, "nuniq")
        for n, name in enumerate(("ubi00", "ubi01", "ubi02", "ubi10", "ubi11",
                                  "ubi12", "ubi20", "ubi21", "ubi22")):
            output_map.addcolumn(ub[n], name)
        for n, name in enumerate(("eps00", "eps01", "eps02", "eps10", "eps11",
                                  "eps12", "eps20", "eps21", "eps22")):
            output_map.addcolumn(ep[n], name)
 
        self.refinedmap = output_map
        if save:
            if verbose:
                print("saving refined map to disk")
            self.savemap(filename=output_filename, refined=True)
        return output_map


def prepare_refine_bash(pbp_object, id11_code_path, output_filename=None,
                        cpus_per_task=16, time_h=1, partition="nice",
                        mem_G=32):
    ds = pbp_object.dset
    pbp_object_file = pbp_object.own_filename

    slurm_pbp_path = os.path.join(ds.analysispath, "slurm_pbp_refine")
    if not os.path.exists(slurm_pbp_path):
        os.mkdir(slurm_pbp_path)

    bash_script_path = os.path.join(
        slurm_pbp_path, ds.dsname + "_pbp_refine_slurm.sh")
    worker_path = os.path.abspath(__file__)
    if worker_path.endswith(".pyc"):
        worker_path = worker_path[:-1]
    outfile_path = os.path.join(
        slurm_pbp_path, ds.dsname + "_pbp_refine_slurm_%A.out")
    errfile_path = os.path.join(
        slurm_pbp_path, ds.dsname + "_pbp_refine_slurm_%A.err")
    log_path = os.path.join(
        slurm_pbp_path, ds.dsname + "_pbp_refine_slurm_$SLURM_JOB_ID.log")

    outname = "" if output_filename is None else " " + output_filename

    bash_script_string = """#!/bin/bash
#SBATCH --job-name=pbp_refine
#SBATCH --output={outfile_path}
#SBATCH --error={errfile_path}
#SBATCH --time={time_h}:00:00
#SBATCH --partition={partition}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem={mem_G}G

date
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH={id11_code_path}
python3 {worker_path} {pbp_object_file}{outname} > {log_path} 2>&1
date
""".format(outfile_path=outfile_path, errfile_path=errfile_path,
           time_h=time_h, partition=partition, cpus_per_task=cpus_per_task,
           mem_G=mem_G, id11_code_path=id11_code_path,
           worker_path=worker_path, outname=outname,
           pbp_object_file=os.path.abspath(pbp_object_file).replace("/mnt/storage", ""),
           log_path=log_path)

    with open(bash_script_path, "w") as f:
        f.write(bash_script_string)
    return bash_script_path


# ==========================================================================
#  The refinement kernel
# ==========================================================================
@numba.njit(cache=True, parallel=True)
def refine_map(
    # voxels to refine (reconstruction space)
    points_ri, points_rj,
    # pbpmap: CSR over (ri, rj) -> rows of pbpmap_ubis, which is (M, 3, 3)
    pstart, porder, pbpmap_ubis, nrj,
    # sample grids / mask
    sx_grid, sy_grid, mask,
    # peak data, PERMUTED into partition order, all C-contiguous
    omega_partitions, dty_partitions, usin, ucos, ymin,
    xl, yl, zl, eta, sum_intensity, sinomega, cosomega, omega, dty, dtyi, xpos,
    gve_all,                                   # (N, 3)
    # geometry: ystep/y0/B0, then the g-vector pack from gvec_pars()
    ystep, y0, B0,
    A, oc, invwvln, a_is_eye,
    # tolerances
    tol, merge_tol, min_grain_npks, sigma_g, weight_reg,
    # execution
    maxlocal, nchunks,
):
    M = pbpmap_ubis.shape[0]
    final_ubis = np.full((M, 3, 3), np.nan)
    final_eps = np.full((M, 3, 3), np.nan)
    final_npks = np.full(M, np.nan)
    final_nuniq = np.full(M, np.nan)
    overflow = np.zeros(nchunks, dtype=np.int64)
    degenerate = np.zeros(nchunks, dtype=np.int64)
    gather_fail = np.zeros(nchunks, dtype=np.int64)
 
    npoints = points_ri.shape[0]
    tolsq_assign = tol * tol
    tolsq_merge = merge_tol * merge_tol
 
    # paralellise over chunks
    # each chunk is persistent and does its own iteration over its voxels
    # this is MUCH better for memory allocation
    for chunk in numba.prange(nchunks):
        # ---- scratch, allocated once per chunk (not once per voxel) -------
        loc = np.empty(maxlocal, dtype=np.int64)          # peak indices on the ray
        loc_yd = np.empty(maxlocal, dtype=np.float64)     # and their y distances
        gvl = np.empty((maxlocal, 3), dtype=np.float64)   # their g-vectors
 
        sel = np.empty(maxlocal, dtype=np.int64)          # assigned subset
        hklsel = np.empty((maxlocal, 3), dtype=np.int64)
        keys = np.empty(maxlocal, dtype=np.int64)
        srt = np.empty(maxlocal, dtype=np.int64)
        labels = np.empty(maxlocal, dtype=np.int64)
 
        m_sI = np.empty(maxlocal, dtype=np.float64)
        # lab-frame peak positions, not sc/fc: compute_xyz_lab is affine in
        # (sc, fc), so the intensity-weighted mean of xl/yl/zl is exactly
        # xl/yl/zl of the intensity-weighted mean sc/fc. That takes the
        # detector transform out of the inner loop entirely.
        m_xl = np.empty(maxlocal, dtype=np.float64)
        m_yl = np.empty(maxlocal, dtype=np.float64)
        m_zl = np.empty(maxlocal, dtype=np.float64)
        m_om = np.empty(maxlocal, dtype=np.float64)
        m_dty = np.empty(maxlocal, dtype=np.float64)
        m_xp = np.empty(maxlocal, dtype=np.float64)
        m_eta = np.empty(maxlocal, dtype=np.float64)
 
        g_sel = np.empty(maxlocal, dtype=np.int64)        # merged peaks on the ray
        g_yd = np.empty(maxlocal, dtype=np.float64)
        gvm = np.empty((maxlocal, 3), dtype=np.float64)   # merged g-vectors
 
        fit_sel = np.empty(maxlocal, dtype=np.int64)
        fit_hkl = np.empty((maxlocal, 3), dtype=np.int64)
 
        # UBI fit
        GtWWG = np.empty((3, 3), dtype=np.float64)   # sum w^2 G G^T
        GtWWH = np.empty((3, 3), dtype=np.float64)   # sum w^2 G G_hkl^T
 
        # strain fit
        y_eps = np.empty(maxlocal, dtype=np.float64)   # eq (14) measurements
        Mj = np.empty(6, dtype=np.float64)             # eq (16) one row of M
        Wdiag = np.empty(maxlocal, dtype=np.float64)   # eq (18) weights
        MtWWM = np.empty((6, 6), dtype=np.float64)     # eq (20) normal matrix
        MtWWy = np.empty((6, 1), dtype=np.float64)     # eq (20) rhs, then s
        
        ubi_out = np.empty((3, 3), dtype=np.float64)
        Umat = np.empty((3, 3), dtype=np.float64)
        n_ovf = 0
        n_degen = 0
        n_gfail = 0
 
        # voxels within each chunk
        for pi in range(chunk, npoints, nchunks):
            ri = points_ri[pi]
            rj = points_rj[pi]
            if not mask[ri, rj]:
                continue
 
            cell = ri * nrj + rj
            plo = pstart[cell]
            phi = pstart[cell + 1]
            if phi == plo:
                continue
 
            any_ubi = False
            for q in range(plo, phi):
                if not np.isnan(pbpmap_ubis[porder[q], 0, 0]):
                    any_ubi = True
                    break
            if not any_ubi:
                continue
 
            xi0 = sx_grid[ri, rj]
            yi0 = sy_grid[ri, rj]
 
            # ---- gather the peaks whose ray passes through this voxel ------
            # voxel_mask.fill_voxel_idx is the authoritative selection: it
            # walks its own (omega bin, dty bin) partition with a per-bin dty
            # padding, then applies the exact predicate per peak. Returns -1
            # rather than raising if the buffer is short, because an exception
            # here would not escape the prange.
            # loc and loc_yd are idx_buffer and ydist_buffer
            nloc = fill_voxel_idx(
                xi0, yi0, y0, ystep, ymin,
                omega_partitions, dty_partitions, dty,
                sinomega, cosomega, usin, ucos,
                loc, loc_yd)
            if nloc < 0:
                n_gfail += 1          # maxlocal was wrong -- should be impossible
                continue
            if nloc == 0:
                continue
            for t in range(nloc):
                q = loc[t]
                gvl[t, 0] = gve_all[q, 0]
                gvl[t, 1] = gve_all[q, 1]
                gvl[t, 2] = gve_all[q, 2]
 
            # ---- loop over candidate UBIs at this voxel --------------------
            for q in range(plo, phi):
                row = porder[q]
                u00 = pbpmap_ubis[row, 0, 0]
                if np.isnan(u00):
                    continue
                u01 = pbpmap_ubis[row, 0, 1]
                u02 = pbpmap_ubis[row, 0, 2]
                u10 = pbpmap_ubis[row, 1, 0]
                u11 = pbpmap_ubis[row, 1, 1]
                u12 = pbpmap_ubis[row, 1, 2]
                u20 = pbpmap_ubis[row, 2, 0]
                u21 = pbpmap_ubis[row, 2, 1]
                u22 = pbpmap_ubis[row, 2, 2]
 
                # ---- assign with `tol` ------------------------------------
                nsel = 0
                for t in range(nloc):
                    Gx = gvl[t, 0]
                    Gy = gvl[t, 1]
                    Gz = gvl[t, 2]
                    hf0 = u00 * Gx + u01 * Gy + u02 * Gz    # (UB)^-1 G
                    hf1 = u10 * Gx + u11 * Gy + u12 * Gz
                    hf2 = u20 * Gx + u21 * Gy + u22 * Gz
                    hi0 = np.round(hf0)                     # nearest integers
                    hi1 = np.round(hf1)
                    hi2 = np.round(hf2)
                    dh0 = hf0 - hi0                         # eq (8) residual
                    dh1 = hf1 - hi1
                    dh2 = hf2 - hi2
                    if dh0 * dh0 + dh1 * dh1 + dh2 * dh2 < tolsq_assign:
                        sel[nsel] = t
                        hklsel[nsel, 0] = np.int64(hi0)
                        hklsel[nsel, 1] = np.int64(hi1)
                        hklsel[nsel, 2] = np.int64(hi2)
                        nsel += 1
 
                if nsel == 0:
                    continue
 
                # ---- merge on (h, k, l, etasign, dtyi) --------------------
                for t in range(nsel):
                    p = loc[sel[t]]
                    h0i = _clamp_hkl(hklsel[t, 0])
                    h1i = _clamp_hkl(hklsel[t, 1])
                    h2i = _clamp_hkl(hklsel[t, 2])
                    if (h0i != hklsel[t, 0] or h1i != hklsel[t, 1]
                            or h2i != hklsel[t, 2]):
                        n_ovf += 1
                    keys[t] = _pack5(h0i, h1i, h2i, _etasign_code(eta[p]), dtyi[p])
                nlab = _label_groups(keys, nsel, srt, labels)
 
                for t in range(nlab):
                    m_sI[t] = 0.0
                    m_xl[t] = 0.0
                    m_yl[t] = 0.0
                    m_zl[t] = 0.0
                    m_om[t] = 0.0
                    m_dty[t] = 0.0
                    m_xp[t] = 0.0
                    m_eta[t] = 0.0
                # accumulate in input order, to match np.bincount exactly
                for t in range(nsel):
                    lb = labels[t]
                    p = loc[sel[t]]
                    wI = sum_intensity[p]
                    m_sI[lb] += wI
                    m_xl[lb] += xl[p] * wI
                    m_yl[lb] += yl[p] * wI
                    m_zl[lb] += zl[p] * wI
                    m_om[lb] += omega[p] * wI
                    m_dty[lb] += dty[p] * wI
                    m_xp[lb] += xpos[p] * wI
                    m_eta[lb] += eta[p] * wI
                for t in range(nlab):
                    s = m_sI[t]
                    m_xl[t] /= s
                    m_yl[t] /= s
                    m_zl[t] /= s
                    m_om[t] /= s
                    m_dty[t] /= s
                    m_xp[t] /= s
                    m_eta[t] /= s
 
                # ---- re-apply the voxel mask, and build the merged
                #      g-vectors in the same pass. The mask test needs sin and
                #      cos of the merged omega, which is exactly what the
                #      g-vector needs, so they are computed once. Straight into
                #      the pre-allocated gvm: no per-candidate allocation.
                ng = 0
                for t in range(nlab):
                    orad = np.radians(m_om[t])
                    som = np.sin(orad)
                    com = np.cos(orad)
                    d = abs(y0 - xi0 * som - yi0 * com - m_dty[t])
                    if d <= ystep:
                        g_sel[ng] = t
                        g_yd[ng] = d
                        gx, gy, gz = gvec_one(m_xl[t], m_yl[t], m_zl[t],
                                              m_xp[t], som, com,
                                              A, oc, invwvln, a_is_eye)
                        gvm[ng, 0] = gx
                        gvm[ng, 1] = gy
                        gvm[ng, 2] = gz
                        ng += 1
                if ng == 0:
                    continue
 
                # ---- reassign the merged peaks, eq (8) with e_hkl ---------
                nfit = 0
                for t in range(ng):
                    Gx = gvm[t, 0]                          # G, eq (5)
                    Gy = gvm[t, 1]
                    Gz = gvm[t, 2]
                    hf0 = u00 * Gx + u01 * Gy + u02 * Gz    # (UB)^-1 G
                    hf1 = u10 * Gx + u11 * Gy + u12 * Gz
                    hf2 = u20 * Gx + u21 * Gy + u22 * Gz
                    hi0 = np.round(hf0)                     # nearest integers
                    hi1 = np.round(hf1)
                    hi2 = np.round(hf2)
                    dh0 = hf0 - hi0                         # eq (8) residual
                    dh1 = hf1 - hi1
                    dh2 = hf2 - hi2
                    if dh0 * dh0 + dh1 * dh1 + dh2 * dh2 < tolsq_merge:
                        fit_sel[nfit] = t
                        fit_hkl[nfit, 0] = np.int64(hi0)    # G_hkl
                        fit_hkl[nfit, 1] = np.int64(hi1)
                        fit_hkl[nfit, 2] = np.int64(hi2)
                        nfit += 1
 
                if nfit <= min_grain_npks:
                    continue
 
                # ---- weighted UBI fit via 3x3 normal equations ------------
                for a in range(3):
                    for b in range(3):
                        GtWWG[a, b] = 0.0
                        GtWWH[a, b] = 0.0
                for t in range(nfit):
                    ii = fit_sel[t]
                    # eq (6): 1 at the beam centre, beam_size/3 regularisation
                    w = weight_reg / (g_yd[ii] + weight_reg)
                    ww = w * w
                    Gx = gvm[ii, 0]                     # G, eq (5)
                    Gy = gvm[ii, 1]
                    Gz = gvm[ii, 2]
                    h0 = np.float64(fit_hkl[t, 0])      # G_hkl
                    h1 = np.float64(fit_hkl[t, 1])
                    h2 = np.float64(fit_hkl[t, 2])
                    GtWWG[0, 0] += ww * Gx * Gx
                    GtWWG[0, 1] += ww * Gx * Gy
                    GtWWG[0, 2] += ww * Gx * Gz
                    GtWWG[1, 1] += ww * Gy * Gy
                    GtWWG[1, 2] += ww * Gy * Gz
                    GtWWG[2, 2] += ww * Gz * Gz
                    GtWWH[0, 0] += ww * Gx * h0
                    GtWWH[0, 1] += ww * Gx * h1
                    GtWWH[0, 2] += ww * Gx * h2
                    GtWWH[1, 0] += ww * Gy * h0
                    GtWWH[1, 1] += ww * Gy * h1
                    GtWWH[1, 2] += ww * Gy * h2
                    GtWWH[2, 0] += ww * Gz * h0
                    GtWWH[2, 1] += ww * Gz * h1
                    GtWWH[2, 2] += ww * Gz * h2
                GtWWG[1, 0] = GtWWG[0, 1]
                GtWWG[2, 0] = GtWWG[0, 2]
                GtWWG[2, 1] = GtWWG[1, 2]
 
                scale = GtWWG[0, 0]
                if GtWWG[1, 1] > scale:
                    scale = GtWWG[1, 1]
                if GtWWG[2, 2] > scale:
                    scale = GtWWG[2, 2]
                ok = _chol_solve(GtWWG, GtWWH, 3, 3, 1e-14 * scale)
 
                if ok:
                    # GtWWH now holds X = (G^T W^T W G)^-1 G^T W^T W G_hkl,
                    # and (UB)^-1 = X^T, so that (UB)^-1 G ~= G_hkl.
                    for a in range(3):
                        for b in range(3):
                            ubi_out[a, b] = GtWWH[b, a]
 
                    finite = True
                    for a in range(3):
                        for b in range(3):
                            if not np.isfinite(ubi_out[a, b]):
                                finite = False
                    # rank 3, cond < 1e14, positive determinant -- same test
                    # as before, but without a LAPACK call per candidate
                    worth_fitting = (finite and _det3(ubi_out) > 0.0
                                     and _well_conditioned(ubi_out))
                else:
                    worth_fitting = False
 
                if not worth_fitting:
                    ubi_out[0, 0] = u00
                    ubi_out[0, 1] = u01
                    ubi_out[0, 2] = u02
                    ubi_out[1, 0] = u10
                    ubi_out[1, 1] = u11
                    ubi_out[1, 2] = u12
                    ubi_out[2, 0] = u20
                    ubi_out[2, 1] = u21
                    ubi_out[2, 2] = u22
 
                # ---- reassign with the (probably refined) UBI -------------
                v00 = ubi_out[0, 0]; v01 = ubi_out[0, 1]; v02 = ubi_out[0, 2]
                v10 = ubi_out[1, 0]; v11 = ubi_out[1, 1]; v12 = ubi_out[1, 2]
                v20 = ubi_out[2, 0]; v21 = ubi_out[2, 1]; v22 = ubi_out[2, 2]
 
                # how many peaks we fit strain with
                # ---- reassign with the refined (UB)^-1, eq (8) ------------
                nfit2 = 0
                for t in range(ng):
                    Gx = gvm[t, 0]
                    Gy = gvm[t, 1]
                    Gz = gvm[t, 2]
                    hf0 = v00 * Gx + v01 * Gy + v02 * Gz
                    hf1 = v10 * Gx + v11 * Gy + v12 * Gz
                    hf2 = v20 * Gx + v21 * Gy + v22 * Gz
                    hi0 = np.round(hf0)
                    hi1 = np.round(hf1)
                    hi2 = np.round(hf2)
                    dh0 = hf0 - hi0
                    dh1 = hf1 - hi1
                    dh2 = hf2 - hi2
                    if dh0 * dh0 + dh1 * dh1 + dh2 * dh2 < tolsq_merge:
                        fit_sel[nfit2] = t
                        fit_hkl[nfit2, 0] = np.int64(hi0)
                        fit_hkl[nfit2, 1] = np.int64(hi1)
                        fit_hkl[nfit2, 2] = np.int64(hi2)
                        nfit2 += 1
 
                # ---- unique peak count -----------------------------------
                for t in range(nfit2):
                    keys[t] = _pack4(_clamp_hkl(fit_hkl[t, 0]),
                                     _clamp_hkl(fit_hkl[t, 1]),
                                     _clamp_hkl(fit_hkl[t, 2]),
                                     _etasign_code(m_eta[g_sel[fit_sel[t]]]))
                nuniq = _label_groups(keys, nfit2, srt, labels)
 
                for a in range(3):
                    for b in range(3):
                        final_ubis[row, a, b] = ubi_out[a, b]
                final_npks[row] = nfit2
                final_nuniq[row] = nuniq
 
                if not worth_fitting or nfit2 == 0:
                    continue
 
                # ---- directional strain, Henningsson et al. eq (12)-(20) ---
                if not _ubi_to_u_safe(ubi_out, Umat):
                    n_degen += 1
                    continue          # cell is degenerate: no strain, no crash
                UB0 = Umat.dot(B0)    # U.B0, the undeformed reference cell
 
                mu_eps = 0.0
                for t in range(nfit2):
                    ii = fit_sel[t]
                    Gx = gvm[ii, 0]           # G, measured diffraction vector
                    Gy = gvm[ii, 1]
                    Gz = gvm[ii, 2]
                    h0 = np.float64(fit_hkl[t, 0])
                    h1 = np.float64(fit_hkl[t, 1])      # integer hkl
                    h2 = np.float64(fit_hkl[t, 2])
                    # G0 = U.B0.h, the model diffraction vector
                    G0x = UB0[0, 0] * h0 + UB0[0, 1] * h1 + UB0[0, 2] * h2
                    G0y = UB0[1, 0] * h0 + UB0[1, 1] * h1 + UB0[1, 2] * h2
                    G0z = UB0[2, 0] * h0 + UB0[2, 1] * h1 + UB0[2, 2] * h2
 
                    GtG0 = Gx * G0x + Gy * G0y + Gz * G0z
                    GtG = Gx * Gx + Gy * Gy + Gz * Gz
                    eps = (GtG0 / GtG) - 1.0            # eq (12)
                    y_eps[t] = eps
                    mu_eps += eps
 
                    # eq (18). The first factor down-weights peaks whose ray
                    # passes further from the voxel centroid; the second is
                    # sigma_g propagated onto eps through eq (12).
                    # there's a factor of GtG missing from the equation in the paper
                    # but this form is correct!
                    G0tG0 = G0x * G0x + G0y * G0y + G0z * G0z
                    if GtG != 0.0:
                        sigma_eps_i = np.sqrt(sigma_g * sigma_g * G0tG0
                                              / (GtG * GtG))
                    else:
                        sigma_eps_i = 1.0
                    Wdiag[t] = (weight_reg / (g_yd[ii] + weight_reg)) / sigma_eps_i
 
                mu_eps /= nfit2
                var_eps = 0.0
                for t in range(nfit2):
                    d = y_eps[t] - mu_eps
                    var_eps += d * d
                sigma_eps = np.sqrt(var_eps / nfit2)
                eps_hi = mu_eps + 3.5 * sigma_eps       # eq (19)
                eps_lo = mu_eps - 3.5 * sigma_eps
 
                wmax = 0.0
                for t in range(nfit2):
                    if y_eps[t] > eps_hi or y_eps[t] < eps_lo:
                        Wdiag[t] = 0.0                  # outlier, eq (19)
                    elif Wdiag[t] > wmax:
                        wmax = Wdiag[t]
                if wmax <= 0.0:
                    continue
                for t in range(nfit2):
                    Wdiag[t] /= wmax
 
                # eq (20): s = (M^T W^T W M)^-1 M^T W^T W y, accumulated
                # directly rather than materialising M and W.
                for a in range(6):
                    MtWWy[a, 0] = 0.0
                    for b in range(6):
                        MtWWM[a, b] = 0.0
                for t in range(nfit2):
                    w = Wdiag[t]
                    if w == 0.0:
                        continue
                    ii = fit_sel[t]
                    Gx = gvm[ii, 0]
                    Gy = gvm[ii, 1]
                    Gz = gvm[ii, 2]
                    Gnorm = np.sqrt(Gx * Gx + Gy * Gy + Gz * Gz)
                    k1 = Gx / Gnorm                     # kappa, eq (17)
                    k2 = Gy / Gnorm
                    k3 = Gz / Gnorm
                    Mj[0] = k1 * k1                     # eq (16)
                    Mj[1] = k2 * k2
                    Mj[2] = k3 * k3
                    Mj[3] = 2.0 * k1 * k2
                    Mj[4] = 2.0 * k1 * k3
                    Mj[5] = 2.0 * k2 * k3
                    ww = w * w
                    for a in range(6):
                        MtWWy[a, 0] += ww * Mj[a] * y_eps[t]
                        for b in range(a, 6):
                            MtWWM[a, b] += ww * Mj[a] * Mj[b]
                for a in range(6):
                    for b in range(a):
                        MtWWM[a, b] = MtWWM[b, a]
 
                mscale = MtWWM[0, 0]
                for a in range(1, 6):
                    if MtWWM[a, a] > mscale:
                        mscale = MtWWM[a, a]
                if not _chol_solve(MtWWM, MtWWy, 6, 1, 1e-14 * mscale):
                    continue
 
                # MtWWy now holds s = [eps_xx, eps_yy, eps_zz,
                #                      eps_xy, eps_xz, eps_yz], eq (15)
                eps_xx = MtWWy[0, 0]; eps_yy = MtWWy[1, 0]; eps_zz = MtWWy[2, 0]
                eps_xy = MtWWy[3, 0]; eps_xz = MtWWy[4, 0]; eps_yz = MtWWy[5, 0]
                final_eps[row, 0, 0] = eps_xx
                final_eps[row, 0, 1] = eps_xy
                final_eps[row, 0, 2] = eps_xz
                final_eps[row, 1, 0] = eps_xy
                final_eps[row, 1, 1] = eps_yy
                final_eps[row, 1, 2] = eps_yz
                final_eps[row, 2, 0] = eps_xz
                final_eps[row, 2, 1] = eps_yz
                final_eps[row, 2, 2] = eps_zz
 
        overflow[chunk] = n_ovf
        degenerate[chunk] = n_degen
        gather_fail[chunk] = n_gfail
 
    return (final_ubis, final_eps, final_npks, final_nuniq,
            overflow, degenerate, gather_fail)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print("usage: python3 pbp_refine.py <PBPRefine.h5> [output_map.h5]")
        return 1
    pbprefinefile = argv[0]
    output_filename = argv[1] if len(argv) > 1 else None

    nthreads = int(os.environ.get("SLURM_CPUS_PER_TASK",
                                  max(cImageD11.cores_available() - 1, 1)))
    print("Loading refinement object from", pbprefinefile)
    refine = PBPRefine.from_h5(pbprefinefile)
    print("Running refinement on", nthreads, "threads")
    refine.run_refine(output_filename=output_filename, nthreads=nthreads)
    print("Refinement complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())