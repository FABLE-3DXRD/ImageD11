"""Point-by-point refinement code. Originally written by Axel Henningsson with adaptations by James Ball"""
from __future__ import print_function, division

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import sys


import ctypes
import platform
import subprocess
import hashlib

import time
import numpy as np

import h5py
import numba

# ignore Numba dot product performance warnings?
import warnings
warnings.simplefilter('ignore', category=numba.core.errors.NumbaPerformanceWarning)

from skimage.filters import threshold_otsu
from skimage.morphology import convex_hull_image

from ImageD11 import cImageD11
from ImageD11 import unitcell
import ImageD11.sinograms.dataset
import ImageD11.indexing

from ImageD11.sinograms import geometry
from ImageD11.sinograms.roi_iradon import run_iradon
from ImageD11.sinograms.tensor_map import unitcell_to_b
from ImageD11.sinograms.sinogram import save_array
from ImageD11.sinograms.point_by_point import PBPMap
from ImageD11.sinograms.voxel_mask import (
    VoxelSinoMasker, fill_voxel_idx, max_candidates as vm_max_candidates)

# --------------------------------------------------------------------------
# glibc allocator tuning.  Numba's NRT calls malloc; blocks above the mmap
# threshold (128 kB by default) are mmap'd and munmap'd on every alloc/free,
# which means the kernel zero-fills every page again and every free triggers
# a TLB shootdown IPI to all cores.  Raising the threshold keeps them in the
# arena.  Harmless if it fails (musl, macOS, Windows).
# --------------------------------------------------------------------------
def _tune_malloc(threshold=1 << 30):
    if platform.system() != "Linux":
        return False
    try:
        libc = ctypes.CDLL("libc.so.6", use_errno=True)
        M_TRIM_THRESHOLD = -1
        M_MMAP_THRESHOLD = -3
        ok = libc.mallopt(M_TRIM_THRESHOLD, ctypes.c_int(threshold))
        ok &= libc.mallopt(M_MMAP_THRESHOLD, ctypes.c_int(threshold))
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
    """idx[0:n] <- argsort(keys[0:n]).  No allocation."""
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
    """Assign a group label per entry, in ascending-key order.  Returns ngroups."""
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


# stuff we need to compute g-vectors
# from ImageD11.transform
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


@numba.njit(cache=True)
def compute_grain_origins(omega, wedge, chi, t_x, t_y, t_z):
    w = np.radians(wedge)
    WI = np.array([[np.cos(w),         0, -np.sin(w)],
                   [0,           1.0,         0],
                   [np.sin(w),         0,  np.cos(w)]], np.float64)
    c = np.radians(chi)
    CI = np.array([[1,            0.0,         0],
                   [0,     np.cos(c), -np.sin(c)],
                   [0,     np.sin(c),  np.cos(c)]], np.float64)
    t = np.zeros((3, omega.shape[0]), np.float64)  # crystal translations
    # Rotations in reverse order compared to making g-vector
    # also reverse directions. this is trans at all zero to
    # current setting. gv is scattering vector to all zero
    om_r = np.radians(omega)
    # This is the real rotation (right handed, g back to k)
    t[0, :] = np.cos(om_r) * t_x - np.sin(om_r) * t_y
    t[1, :] = np.sin(om_r) * t_x + np.cos(om_r) * t_y
    t[2, :] = t_z
    if chi != 0.0:
        c = np.cos(np.radians(chi))
        s = np.sin(np.radians(chi))
        u = np.zeros(t.shape, np.float64)
        u[0, :] = t[0, :]
        u[1, :] = c * t[1, :] + -s * t[2, :]
        u[2, :] = s * t[1, :] + c * t[2, :]
        t = u
    if wedge != 0.0:
        c = np.cos(np.radians(wedge))
        s = np.sin(np.radians(wedge))
        u = np.zeros(t.shape, np.float64)
        u[0, :] = c * t[0, :] + -s * t[2, :]
        u[1, :] = t[1, :]
        u[2, :] = s * t[0, :] + c * t[2, :]
        t = u
    return t


@numba.njit(cache=True)
def compute_xyz_lab(sc, fc,
                    y_center=0., y_size=0., tilt_y=0.,
                    z_center=0., z_size=0., tilt_z=0.,
                    tilt_x=0.,
                    distance=0.,
                    o11=1.0, o12=0.0, o21=0.0, o22=-1.0):
    # Matrix for the tilt rotations
    r2r1 = detector_rotation_matrix(tilt_x, tilt_y, tilt_z)
    # Peak positions in 3D space
    #  - apply detector orientation
    peaks_on_detector = np.stack((sc, fc))
    peaks_on_detector[0, :] = (peaks_on_detector[0, :] - z_center) * z_size
    peaks_on_detector[1, :] = (peaks_on_detector[1, :] - y_center) * y_size
    #
    detector_orientation = [[o11, o12], [o21, o22]]
    flipped = np.dot(np.array(detector_orientation, np.float64),
                     peaks_on_detector)

    vec = np.stack((np.zeros(flipped.shape[1]), flipped[1, :], flipped[0, :]))

    # Position of diffraction spots in 3d space after detector tilts about
    # the beam centre on the detector
    rotvec = np.dot(r2r1, vec)
    # Now add the distance (along x)
    rotvec[0, :] = rotvec[0, :] + distance
    return rotvec


@numba.njit(cache=True)
def compute_tth_eta_from_xyz(peaks_xyz, omega,
                             t_x=0.0, t_y=0.0, t_z=0.0,
                             wedge=0.0,  # Wedge == theta on 4circ
                             chi=0.0):  # last line is for laziness -

    s1 = peaks_xyz - compute_grain_origins(omega, wedge, chi, t_x, t_y, t_z)

    # CHANGED to HFP convention 4-9-2007
    eta = np.degrees(np.arctan2(-s1[1, :], s1[2, :]))
    s1_perp_x = np.sqrt(s1[1, :] * s1[1, :] + s1[2, :] * s1[2, :])
    tth = np.degrees(np.arctan2(s1_perp_x, s1[0, :]))
    return tth, eta


@numba.njit(cache=True)
def compute_tth_eta(sc, fc, omega,
                    y_center=0., y_size=0., tilt_y=0.,
                    z_center=0., z_size=0., tilt_z=0.,
                    tilt_x=0.,
                    distance=0.,
                    o11=1.0, o12=0.0, o21=0.0, o22=-1.0,
                    t_x=0.0, t_y=0.0, t_z=0.0,
                    wedge=0.0,
                    chi=0.0):
    peaks_xyz = compute_xyz_lab(
        sc, fc,
        y_center=y_center, y_size=y_size, tilt_y=tilt_y,
        z_center=z_center, z_size=z_size, tilt_z=tilt_z,
        tilt_x=tilt_x,
        distance=distance,
        o11=o11, o12=o12, o21=o21, o22=o22)

    tth, eta = compute_tth_eta_from_xyz(
        peaks_xyz, omega,
        t_x=t_x, t_y=t_y, t_z=t_z,
        wedge=wedge,
        chi=chi)

    return tth, eta


@numba.njit(cache=True)
def compute_k_vectors(tth, eta, wvln):
    """
    generate k vectors - scattering vectors in laboratory frame
    """
    tth = np.radians(tth)
    eta = np.radians(eta)
    c = np.cos(tth / 2)  # cos theta
    s = np.sin(tth / 2)  # sin theta
    ds = 2 * s / wvln
    k = np.zeros((3, tth.shape[0]), np.float64)
    # x - along incident beam
    k[0, :] = -ds * s  # this is negative x
    # y - towards door
    k[1, :] = -ds * c * np.sin(eta)  # CHANGED eta to HFP convention 4-9-2007
    # z - towards roof
    k[2, :] = ds * c * np.cos(eta)
    return k


@numba.njit(cache=True)
def compute_g_from_k(k, omega, wedge=0, chi=0):
    """
    Compute g-vectors with cached k-vectors
    """
    om = np.radians(omega)
    # G-vectors - rotate k onto the crystal axes
    g = np.zeros((3, k.shape[1]), np.float64)
    t = np.zeros((3, k.shape[1]), np.float64)

    if wedge != 0.0:
        c = np.cos(np.radians(wedge))
        s = np.sin(np.radians(wedge))
        t[0, :] = c * k[0, :] + s * k[2, :]
        t[1, :] = k[1, :]
        t[2, :] = -s * k[0, :] + c * k[2, :]
        k = t.copy()
    if chi != 0.0:
        c = np.cos(np.radians(chi))
        s = np.sin(np.radians(chi))
        t[0, :] = k[0, :]
        t[1, :] = c * k[1, :] + s * k[2, :]
        t[2, :] = -s * k[1, :] + c * k[2, :]
        k = t.copy()
    # This is the reverse rotation (left handed, k back to g)
    g[0, :] = np.cos(om) * k[0, :] + np.sin(om) * k[1, :]
    g[1, :] = -np.sin(om) * k[0, :] + np.cos(om) * k[1, :]
    g[2, :] = k[2, :]
    return g


@numba.njit(cache=True)
def compute_g_vectors(tth,
                      eta,
                      omega,
                      wvln,
                      wedge=0.0,
                      chi=0.0):
    """
    Generates spot positions in reciprocal space from
      twotheta, wavelength, omega and eta
    Assumes single axis vertical
    ... unless a wedge angle is specified
    """
    k = compute_k_vectors(tth, eta, wvln)
    return compute_g_from_k(k, omega, wedge, chi)


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


@numba.njit(cache=True)
def compute_gve(sc, fc, omega, xpos,
                distance, y_center, y_size, tilt_y, z_center, z_size, tilt_z, tilt_x,
                o11, o12, o21, o22,
                t_x, t_y, t_z, wedge, chi, wavelength):
    this_distance = distance - xpos

    tth, eta = compute_tth_eta(sc, fc, omega,
                               y_center=y_center,
                               y_size=y_size,
                               tilt_y=tilt_y,
                               z_center=z_center,
                               z_size=z_size,
                               tilt_z=tilt_z,
                               tilt_x=tilt_x,
                               distance=this_distance,
                               o11=o11,
                               o12=o12,
                               o21=o21,
                               o22=o22,
                               t_x=t_x,
                               t_y=t_y,
                               t_z=t_z,
                               wedge=wedge,
                               chi=chi
                               )

    gve = compute_g_vectors(tth, eta,
                            omega,
                            wavelength,
                            wedge=wedge,
                            chi=chi)

    return gve


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


# ==========================================================================
#  INDEX CACHE
#
#  Split in two, because get_origins needs the bucket index but runs *before*
#  xpos_refined exists -- it is what produces it.
#
#    PBPPeakIndex : the (omega frame, dtyi) CSR and the permutation.
#                   Depends on omega, dty, dtyi, ystep, y0, mask shape.
#                   Used by both get_origins and run_refine.
#    PBPGveCache  : gve_all in permuted order.
#                   Additionally depends on sc, fc, xpos_refined and the
#                   geometry parameters. Only run_refine needs it.
#
#  Not cached: the permuted columns and the pbpmap CSR, which are a
#  fancy-index and a counting sort away and cost less than writing them out.
# ==========================================================================
_PEAK_H5GROUP = "PBPPeakIndex"
_GVE_H5GROUP = "PBPGveCache"
_INDEX_VERSION = 4


def _hash_cols(h, icolf, names):
    for name in names:
        if name not in icolf.titles:
            h.update(("%s:absent;" % name).encode())
            continue
        arr = np.ascontiguousarray(getattr(icolf, name))
        h.update(name.encode())
        h.update(arr.dtype.str.encode())
        h.update(np.int64(arr.shape[0]).tobytes())
        h.update(arr.tobytes())


def _peak_index_fingerprint(refine, omega_binsize="auto", omega_step=None):
    """Inputs build_peak_index depends on. Deliberately excludes xpos_refined."""
    h = hashlib.blake2b(digest_size=16)
    h.update(("v%d peak" % _INDEX_VERSION).encode())
    _hash_cols(h, refine.icolf, ("omega", "dty", "dtyi", "frm"))
    h.update(repr((float(refine.ystep), float(refine.y0),
                   None if omega_step is None else float(omega_step),
                   omega_binsize,
                   tuple(refine.mask.shape))).encode())
    return h.hexdigest()


def _gve_fingerprint(refine, peak_fingerprint):
    """Everything the g-vectors depend on, on top of the peak index."""
    h = hashlib.blake2b(digest_size=16)
    h.update(("v%d gve " % _INDEX_VERSION).encode())
    h.update(peak_fingerprint.encode())
    _hash_cols(h, refine.icolf, ("sc", "fc", "xpos_refined"))
    pars = refine.icolf.parameters.get_parameters()
    for k in sorted(pars):
        h.update(("%s=%r;" % (k, pars[k])).encode())
    return h.hexdigest()


def save_peak_index(idx, filename, fingerprint):
    with h5py.File(filename, "a") as hout:
        if _PEAK_H5GROUP in hout:
            del hout[_PEAK_H5GROUP]
        g = hout.create_group(_PEAK_H5GROUP)
        g.attrs["fingerprint"] = fingerprint
        g.attrs["version"] = _INDEX_VERSION
        for k in ("nom", "ndty"):
            g.attrs[k] = int(idx[k])
        for k in ("ymin", "ray_margin", "dev"):
            g.attrs[k] = float(idx[k])
        for k in ("order", "omega_partitions", "dty_partitions", "usin", "ucos"):
            g.create_dataset(k, data=idx[k])


def load_peak_index(filename, fingerprint):
    if filename is None or not os.path.exists(filename):
        return None
    try:
        with h5py.File(filename, "r") as hin:
            if _PEAK_H5GROUP not in hin:
                return None
            g = hin[_PEAK_H5GROUP]
            if g.attrs.get("version", -1) != _INDEX_VERSION:
                return None
            if g.attrs.get("fingerprint", "") != fingerprint:
                return None
            return dict(
                nom=int(g.attrs["nom"]), ndty=int(g.attrs["ndty"]),
                ymin=float(g.attrs["ymin"]),
                ray_margin=float(g.attrs["ray_margin"]),
                dev=float(g.attrs["dev"]),
                order=g["order"][:],
                omega_partitions=g["omega_partitions"][:],
                dty_partitions=g["dty_partitions"][:],
                usin=g["usin"][:], ucos=g["ucos"][:],
            )
    except (OSError, KeyError):
        return None


def save_gve_cache(gve_all, filename, fingerprint):
    with h5py.File(filename, "a") as hout:
        if _GVE_H5GROUP in hout:
            del hout[_GVE_H5GROUP]
        g = hout.create_group(_GVE_H5GROUP)
        g.attrs["fingerprint"] = fingerprint
        g.attrs["version"] = _INDEX_VERSION
        g.create_dataset("gve_all", data=gve_all)


def load_gve_cache(filename, fingerprint):
    if filename is None or not os.path.exists(filename):
        return None
    try:
        with h5py.File(filename, "r") as hin:
            if _GVE_H5GROUP not in hin:
                return None
            g = hin[_GVE_H5GROUP]
            if g.attrs.get("version", -1) != _INDEX_VERSION:
                return None
            if g.attrs.get("fingerprint", "") != fingerprint:
                return None
            return g["gve_all"][:]
    except (OSError, KeyError):
        return None


def default_index_filename(refine):
    base = getattr(refine, "icolf_filename", None)
    if base is None:
        raise ValueError("no icolf_filename to derive an index cache path from")
    return os.path.splitext(base)[0] + "_pbpindex.h5"


# ==========================================================================
#  Python driver
# ==========================================================================
def _omega_frame_index(refine, omega, omega_step=None, verbose=True):
    """
    Map each peak onto an omega *frame*, not onto its exact encoder readback.

    omega is read back from the rotation encoder, so nominally-identical frames
    differ in the last few digits and np.unique returns roughly one value per
    frame *per dty scan*. That makes the (omega, dtyi) index degenerate: far
    too many cells, almost all empty, and the offset array stops fitting in
    cache. Returns (iomega, omega_centres).
    """
    dset = getattr(refine, "dset", None)
    if omega_step is None and dset is not None and getattr(dset, "obinedges", None) is not None:
        edges = np.asarray(dset.obinedges, dtype=np.float64)
        cens = np.asarray(dset.obincens, dtype=np.float64)
        om = omega % 360.0 if edges[0] >= -1e-9 and omega.min() < -1e-9 else omega
        iomega = np.digitize(om, edges) - 1
        np.clip(iomega, 0, len(cens) - 1, out=iomega)
        if verbose:
            print("omega binned onto %d frames from dset.obinedges" % len(cens))
        return iomega.astype(np.int64), cens

    if omega_step is not None:
        o0 = float(omega.min())
        iomega = np.round((omega - o0) / float(omega_step)).astype(np.int64)
        iomega -= iomega.min()
        cens = o0 + np.arange(iomega.max() + 1) * float(omega_step)
        if verbose:
            print("omega binned onto %d frames with omega_step=%g"
                  % (len(cens), omega_step))
        return iomega, cens

    cens, iomega = np.unique(omega, return_inverse=True)
    return iomega.astype(np.int64), cens


def heuristic_omega_binsize(ystep, rmax):
    """
    Omega bin width such that a voxel at r_max moves ~2*ystep per bin step.

    From Axel Henningsson's voxel_mask.VoxelSinoMasker._heuristic_omega_binsize.
    Coarser bins mean fewer outer iterations per voxel but a wider dty window;
    this factor of 2 is his empirical compromise. His version bounds the voxel
    radius as sqrt(2)*max(|dty|); we have the real r_max from the mask, so we
    use that.
    """
    if rmax <= 0.0:
        return 180.0
    return 2.0 * np.degrees(ystep / rmax)


def _frame_index_from_frm(refine, omega, verbose=True):
    """
    Omega bin index via the per-peak frame number.

    frm says which frame each peak came from. We do NOT assume the
    frame-in-scan index tracks omega: a zigzag scan, or a rotation that runs
    continuously across dty rows, breaks that and there is nothing in the file
    that promises otherwise. Instead we bin each *frame's* omega once, straight
    out of dset.omega, and index that table with frm. Exact for any scan
    ordering, and it digitizes nscans*nframes values rather than one per peak.

    Returns (iomega, obincens, dev) or None if unusable.
    """
    icolf = refine.icolf
    if "frm" not in icolf.titles:
        return None
    dset = getattr(refine, "dset", None)
    om_img = getattr(dset, "omega", None)
    edges = getattr(dset, "obinedges", None)
    cens = getattr(dset, "obincens", None)
    if om_img is None or edges is None or cens is None:
        if verbose:
            print("frm present but dset.omega/obinedges missing; "
                  "falling back to omega binning")
        return None
    om_img = np.asarray(om_img, dtype=np.float64)
    edges = np.asarray(edges, dtype=np.float64)
    cens = np.asarray(cens, dtype=np.float64)
    frm = np.round(np.asarray(icolf.frm, dtype=np.float64)).astype(np.int64)
    if frm.size == 0 or frm.min() < 0 or frm.max() >= om_img.size:
        if verbose:
            print("frm values fall outside dset.omega (%d frames); "
                  "falling back to omega binning" % om_img.size)
        return None

    flat = om_img.ravel()
    if edges[0] >= -1e-9 and flat.min() < -1e-9:
        flat = flat % 360.0
    iom_of_frame = np.digitize(flat, edges) - 1
    np.clip(iom_of_frame, 0, len(cens) - 1, out=iom_of_frame)
    iomega = iom_of_frame[frm]
    dev = float(np.abs(omega - cens[iomega]).max())

    if verbose:
        live = int((np.bincount(iomega, minlength=len(cens)) > 0).sum())
        note = ""
        if om_img.ndim == 2 and "dtyi" in icolf.titles:
            # does the scan row implied by frm agree with the dty bin?
            iscan = frm // om_img.shape[1]
            db = np.asarray(icolf.dtyi, dtype=np.int64)
            note = (", scan row agrees with dty bin for %.1f%% of peaks"
                    % (100.0 * float((iscan == db - db.min()).mean())))
        print("omega bins via frm: %d bins (%d with peaks), %.1f peaks/bin, "
              "|omega - bin centre| max %.5f deg%s"
              % (len(cens), live, frm.size / max(live, 1), dev, note))
    return iomega, cens, dev


def build_peak_index(refine, omega_binsize="auto", omega_step=None,
                     cache_filename=None, use_cache=True, verbose=True):
    """
    Partition the peaks into (omega bin, dty bin) cells.

    The partition itself is voxel_mask.VoxelSinoMasker -- that is the
    authoritative implementation. This wrapper adds what is specific to us:
    the frm/obinedges frame decode, and the on-disk cache.

    Shared by get_origins and run_refine. Needs only omega, dty and dtyi, so
    it can be built before xpos_refined exists.
    """
    icolf = refine.icolf
    n = icolf.nrows

    if use_cache and cache_filename is None:
        try:
            cache_filename = default_index_filename(refine)
        except ValueError:
            use_cache = False

    fingerprint = (_peak_index_fingerprint(refine, omega_binsize, omega_step)
                   if use_cache else None)
    cached = load_peak_index(cache_filename, fingerprint) if use_cache else None

    dtyi = np.ascontiguousarray(icolf.dtyi, dtype=np.int64)
    dbin = dtyi - int(dtyi.min())

    if cached is not None:
        if verbose:
            print("loaded peak index from %s (%d omega bins, %d dty bins)"
                  % (cache_filename, cached["nom"], cached["ndty"]))
        idx = dict(cached)
        idx["dbin"] = dbin
        idx["cache_filename"] = cache_filename
        return idx

    omega = np.ascontiguousarray(icolf.omega, dtype=np.float64)
    dty = np.ascontiguousarray(icolf.dty, dtype=np.float64)
    rmax = float(np.hypot(refine.sx_grid, refine.sy_grid)[refine.mask].max())

    # our frame decode, handed to Axel's partitioner as explicit bins
    kw = {}
    got = None
    if omega_step is not None:
        kw = dict(omega_binsize=float(omega_step))
    elif omega_binsize not in (None, "auto"):
        kw = dict(omega_binsize=float(omega_binsize))   # explicit wins
    else:
        got = _frame_index_from_frm(refine, omega, verbose)   # frm is the default
        if got is not None:
            kw = dict(omega_bin_indices=got[0], omega_bins=got[1])
        elif omega_binsize == "auto":
            kw = dict(omega_binsize=None)               # his heuristic

    t0 = time.perf_counter()
    M = VoxelSinoMasker(omega, dty, refine.ystep)
    M.partition(**kw)                           # inplace=False: icolf untouched
    nom = int(M.sinomega_bins.size)
    ndty = int(M.dty_partitions.shape[1]) - 1

    # max |omega - bin centre|, which sets how far a ray can wander inside a
    # cell and hence the voxel-strip half width used by compute_origins
    ob = np.asarray(M.omega_bins, dtype=np.float64)
    dev = 0.0
    for _i in range(nom):
        _lo = int(M.omega_partitions[_i])
        _hi = int(M.omega_partitions[_i + 1])
        if _hi > _lo:
            _d = float(np.abs(M.omega[_lo:_hi] - ob[_i]).max())
            if _d > dev:
                dev = _d
    ray_margin = 1.5 + rmax * np.radians(dev) / refine.ystep + 0.5

    
    if verbose:
        print("peak index: %d omega bins x %d dty bins = %d cells, %.1f "
              "peaks/cell (%.2f s)" % (nom, ndty, nom * ndty,
                                       n / float(nom * ndty),
                                       time.perf_counter() - t0))
        print("  binsize %s, max |omega - bin centre| %.5f deg, r_max %.1f "
              "-> ray half-width %.2f*ystep"
              % ("supplied frames" if got is not None
                 else "%.4f deg" % M.omega_binsize, dev, rmax, ray_margin))
        if n < nom * ndty:
            print("  WARNING: fewer peaks than cells (%.2f/cell); the index is "
                  "sparse and traversal will dominate." % (n / float(nom * ndty)))

    idx = dict(order=M.peak_ordering,
               omega_partitions=M.omega_partitions,
               dty_partitions=M.dty_partitions,
               usin=M.sinomega_bins, ucos=M.cosomega_bins,
               nom=nom, ndty=ndty, ymin=float(M.ymin),
               ray_margin=ray_margin, dev=dev,
               dbin=dbin, cache_filename=cache_filename)

    if use_cache:
        refine.index_filename = cache_filename
        try:
            save_peak_index(idx, cache_filename, fingerprint)
            if verbose:
                print("saved peak index to %s" % cache_filename)
        except OSError as e:
            print("could not write index cache (%s), continuing" % e)
    return idx


def build_indexes(refine, omega_binsize="auto", omega_step=None,
                  cache_filename=None, use_cache=True, verbose=True):
    """
    Everything refine_map needs: the shared peak index, the permuted columns,
    the g-vectors, and the pbpmap CSR. Returns a dict you can hand back via
    run_refine(index_cache=...).
    """
    icolf = refine.icolf
    idx = build_peak_index(refine, omega_binsize=omega_binsize, omega_step=omega_step,
                           cache_filename=cache_filename,
                           use_cache=use_cache, verbose=verbose)
    cache_filename = idx.pop("cache_filename")
    dbin = idx.pop("dbin")
    order = idx["order"]

    # ---- permuted columns: cheap, always rebuilt --------------------------
    def col(name, dtype=np.float64):
        return np.ascontiguousarray(getattr(icolf, name)[order], dtype=dtype)

    peaks = dict(
        sc=col("sc"), fc=col("fc"), eta=col("eta"),
        sum_intensity=col("sum_intensity"),
        sinomega=col("sinomega"), cosomega=col("cosomega"),
        omega=col("omega"), dty=col("dty"),
        dtyi=np.ascontiguousarray(dbin[order]),
        xpos=col("xpos_refined"),
    )

    gve_all = None
    if use_cache and cache_filename is not None:
        gfp = _gve_fingerprint(refine,
            _peak_index_fingerprint(refine, omega_binsize, omega_step))
        gve_all = load_gve_cache(cache_filename, gfp)
        if gve_all is not None and verbose:
            print("loaded g-vectors from %s" % cache_filename)
    else:
        gfp = None

    if gve_all is None:
        pars = icolf.parameters.get_parameters()
        t0 = time.perf_counter()
        gve = compute_gve(
            peaks["sc"], peaks["fc"], peaks["omega"], peaks["xpos"],
            distance=float(pars["distance"]), y_center=float(pars["y_center"]),
            y_size=float(pars["y_size"]), tilt_y=float(pars["tilt_y"]),
            z_center=float(pars["z_center"]), z_size=float(pars["z_size"]),
            tilt_z=float(pars["tilt_z"]), tilt_x=float(pars["tilt_x"]),
            o11=float(pars["o11"]), o12=float(pars["o12"]),
            o21=float(pars["o21"]), o22=float(pars["o22"]),
            t_x=float(pars["t_x"]), t_y=float(pars["t_y"]), t_z=float(pars["t_z"]),
            wedge=float(pars["wedge"]), chi=float(pars["chi"]),
            wavelength=float(pars["wavelength"]),
        )
        gve_all = np.ascontiguousarray(gve.T)          # (N, 3)
        if verbose:
            print("computed g-vectors in %.2f s" % (time.perf_counter() - t0))
        if gfp is not None:
            try:
                save_gve_cache(gve_all, cache_filename, gfp)
            except OSError as e:
                print("could not write gve cache (%s), continuing" % e)
    peaks["gve_all"] = gve_all

    # ---- pbpmap CSR over (ri, rj), and (M, 3, 3) UBIs --------------------
    ri_col, rj_col = geometry.step_to_recon(refine.pbpmap.i, refine.pbpmap.j,
                                            refine.mask.shape)
    nri, nrj = refine.mask.shape
    pcells = (np.asarray(ri_col, dtype=np.int64) * nrj
              + np.asarray(rj_col, dtype=np.int64))
    pstart, porder = build_csr(pcells, nri * nrj)

    idx["peaks"] = peaks
    idx["pstart"] = pstart
    idx["porder"] = porder
    idx["pbpmap_ubis"] = np.ascontiguousarray(refine.pbpmap.ubi.transpose(2, 0, 1))
    idx["nrj"] = nrj
    return idx


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


# ==========================================================================
#  Origins
# ==========================================================================
@numba.njit(cache=True, parallel=True)
def compute_origins(singlemap, sample_mask,
                    gve, sinomega, cosomega, dty,
                    sx_ax, sy_ax, dsx, dsy,
                    y0, ystep, hkl_tol,
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
                            o00 = singlemap[i, j, 0, 0]
                            if np.isnan(o00):
                                continue
                            o01 = singlemap[i, j, 0, 1]
                            o02 = singlemap[i, j, 0, 2]
                            o10 = singlemap[i, j, 1, 0]
                            o11 = singlemap[i, j, 1, 1]
                            o12 = singlemap[i, j, 1, 2]
                            o20 = singlemap[i, j, 2, 0]
                            o21 = singlemap[i, j, 2, 1]
                            o22 = singlemap[i, j, 2, 2]
                            for q in range(plo, phi):
                                yd = abs(y0 - sxv * sinomega[q]
                                         - syv * cosomega[q] - dty[q])
                                if yd > ystep:
                                    continue
                                g0 = gve[q, 0]
                                g1 = gve[q, 1]
                                g2 = gve[q, 2]
                                hf0 = o00 * g0 + o01 * g1 + o02 * g2
                                hf1 = o10 * g0 + o11 * g1 + o12 * g2
                                hf2 = o20 * g0 + o21 * g1 + o22 * g2
                                e0 = hf0 - np.round(hf0)
                                e1 = hf1 - np.round(hf1)
                                e2 = hf2 - np.round(hf2)
                                if e0 * e0 + e1 * e1 + e2 * e2 < hkl_tol:
                                    w = 1.0 / (yd + 1.0)  # TODO: this is 1.0 um hard-coded!
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
                            o00 = singlemap[i, j, 0, 0]
                            if np.isnan(o00):
                                continue
                            o01 = singlemap[i, j, 0, 1]
                            o02 = singlemap[i, j, 0, 2]
                            o10 = singlemap[i, j, 1, 0]
                            o11 = singlemap[i, j, 1, 1]
                            o12 = singlemap[i, j, 1, 2]
                            o20 = singlemap[i, j, 2, 0]
                            o21 = singlemap[i, j, 2, 1]
                            o22 = singlemap[i, j, 2, 2]
                            for q in range(plo, phi):
                                yd = abs(y0 - sxv * sinomega[q]
                                         - syv * cosomega[q] - dty[q])
                                if yd > ystep:
                                    continue
                                g0 = gve[q, 0]
                                g1 = gve[q, 1]
                                g2 = gve[q, 2]
                                hf0 = o00 * g0 + o01 * g1 + o02 * g2
                                hf1 = o10 * g0 + o11 * g1 + o12 * g2
                                hf2 = o20 * g0 + o21 * g1 + o22 * g2
                                e0 = hf0 - np.round(hf0)
                                e1 = hf1 - np.round(hf1)
                                e2 = hf2 - np.round(hf2)
                                if e0 * e0 + e1 * e1 + e2 * e2 < hkl_tol:
                                    w = 1.0 / (yd + 1.0)  # TODO: this is 1.0 um hard-coded!
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
                     omega_partitions, dty_partitions, nom, ndty,
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
                 hkl_tol_origins=0.05,
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
        self.index_filename = default_index_filename(self)  # CSR / gve cache

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

    def setpeaks(self, colf, icolf_filename=None, del_existing=True, prompt_del=True):
        """Similar to PBP.setpeaks for now"""
        if del_existing and prompt_del:
            print('I will delete an existing refined peaks H5 file if I find it on disk!')
            print('Waiting 10 seconds for you to interrupt this if you are unhappy')
            print('To disable this prompt, set prompt_del=False when calling setpeaks()')
            time.sleep(10)
            print('Continuing')

        # Set the parameters for the peaks
        self.dset.update_colfile_pars(colf, phase_name=self.phase_name)

        uc = unitcell.unitcell_from_parameters(colf.parameters)
        uc.makerings(colf.ds.max(), self.ds_tol)
        # peaks that are on rings
        sel = np.zeros(colf.nrows, bool)
        # rings to use for indexing
        isel = np.zeros(colf.nrows, bool)
        npks = 0
        if self.forref is None:
            self.forref = range(len(uc.ringds))
        hmax = 0
        for i in range(len(uc.ringds)):
            ds = uc.ringds[i]
            hkls = uc.ringhkls[ds]
            if i in self.forref:
                rm = abs(colf.ds - ds) < self.ds_tol
                if rm.sum() == 0:
                    continue
                icut = np.max(colf.sum_intensity[rm]) * self.ifrac
                rm = rm & (colf.sum_intensity > icut)
                isel |= rm
                npks += len(hkls)
                print(
                    i,
                    "%.4f" % (ds),
                    hkls[-1],
                    len(hkls),
                    rm.sum(),
                    "used, sum_intensity>",
                    icut,
                )
                sel |= rm
                for hkl in hkls:
                    hmax = max(np.abs(hkl).max(), hmax)
            else:
                print(i, "%.4f" % (ds), hkls[-1], len(hkls), "skipped")

        isel = isel & ((np.abs(np.sin(np.radians(colf.eta)))) > self.etacut)

        colf.addcolumn(isel, "isel")  # peaks selected for refinement
        dtyi = geometry.dty_to_dtyi(colf.dty, self.ystep, self.ybincens.min())
        colf.addcolumn(dtyi, "dtyi")

        # cache these to speed up selections later
        colf.addcolumn(np.sin(np.radians(colf.omega)), "sinomega")
        colf.addcolumn(np.cos(np.radians(colf.omega)), "cosomega")

        # peaks that are on any rings
        self.colf = colf
        self.icolf = colf.copy()
        self.icolf.filter(sel)
        self.npks = npks
        print(
            "Using for refinement:",
            self.icolf.nrows,
            "npks, forref",
            self.npks,
            self.forref,
        )
        self.uc = uc
        self.savepeaks(icolf_filename=icolf_filename, del_existing=del_existing)

    def savepeaks(self, icolf_filename=None, del_existing=False):
        if icolf_filename is None:
            icolf_filename = self.icolf_filename
        else:
            self.icolf_filename = icolf_filename
        if os.path.exists(icolf_filename) and del_existing:
            os.remove(icolf_filename)
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

        f, ax = pl.subplots(2, 1)
        ax[0].plot(self.colf.ds[::skip], self.colf.sum_intensity[::skip], ",")
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

    def to_h5(self, filename=None, h5group='PBPRefine'):
        # save to h5 file
        if filename is None:
            filename = self.own_filename
        else:
            self.own_filename = filename

        print('Saving icolf to disk')
        self.savepeaks(del_existing=False)
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
                'indexfile': 'index_filename',
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
                'indexfile': 'index_filename',
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

    # ----------------------------------------------------------------------
    def get_origins(self, save_peaks_after=True, nthreads=None, nchunks=None,
                    omega_binsize="auto", omega_step=None, index_cache=None,
                    cache_filename=None, use_cache=True, verbose=True,
                    guess_speed=None, guess_npks=None):
        """
        Fill the xpos_refined column.

        guess_speed / guess_npks are accepted and ignored -- the timing probe
        they drove was there because this used to take hours.
        """
        if 'xpos_refined' in self.icolf.titles:
            raise ValueError('We already have origins in self.icolf! Not recomputing')
        if guess_speed is not None or guess_npks is not None:
            print("note: guess_speed/guess_npks are ignored")

        if nthreads is None:
            nthreads = max(cImageD11.cores_available() - 1, 1)
        numba.set_num_threads(nthreads)
        try:
            numba.set_parallel_chunksize(1)
        except AttributeError:
            pass

        idx = index_cache if index_cache is not None else build_peak_index(
            self, omega_binsize=omega_binsize, omega_step=omega_step,
            cache_filename=cache_filename, use_cache=use_cache, verbose=verbose)
        order = idx["order"]

        # g-vectors as the columnfile already has them (no xpos correction --
        # that is what we are about to compute), permuted into CSR order
        gve = np.ascontiguousarray(
            np.column_stack((self.icolf.gx, self.icolf.gy, self.icolf.gz))[order])
        sinomega = np.ascontiguousarray(self.icolf.sinomega[order])
        cosomega = np.ascontiguousarray(self.icolf.cosomega[order])
        dty = np.ascontiguousarray(self.icolf.dty[order])

        sx_ax, sy_ax, dsx, dsy = grid_axes(self.sx_grid, self.sy_grid)
        singlemap = np.ascontiguousarray(self.singlemap, dtype=np.float64)
        mask = np.ascontiguousarray(self.mask).astype(np.bool_)

        max_cell = int(np.diff(idx["dty_partitions"], axis=1).max())
        if nchunks is None:
            nchunks = nthreads
        nchunks = max(1, min(nchunks, idx["nom"]))

        if verbose:
            nvox, ncells = count_ray_voxels(
                mask, sx_ax, sy_ax, dsx, dsy, self.y0, self.ystep,
                idx["omega_partitions"], idx["dty_partitions"],
                idx["nom"], idx["ndty"], idx["usin"], idx["ucos"],
                idx["ymin"], idx["ray_margin"])
            print("%d non-empty cells, max %d peaks in a cell" % (ncells, max_cell))
            print("mean %.0f voxels visited per ray (full map would be %d)"
                  % (nvox / max(ncells, 1), mask.size))
            print("running on %d threads" % nthreads)

        t0 = time.perf_counter()
        lx_perm = compute_origins(
            singlemap, mask, gve, sinomega, cosomega, dty,
            sx_ax, sy_ax, dsx, dsy,
            self.y0, self.ystep, self.hkl_tol_origins,
            idx["omega_partitions"], idx["dty_partitions"],
            idx["nom"], idx["ndty"], idx["usin"], idx["ucos"],
            idx["ymin"], idx["ray_margin"],
            max_cell, nchunks)
        if verbose:
            print("origins took %.1f s" % (time.perf_counter() - t0))

        # un-permute back to icolf row order
        lx_modified = np.empty_like(lx_perm)
        lx_modified[order] = lx_perm

        self.icolf.addcolumn(lx_modified, 'xpos_refined')
        print('xpos_refined column added to self.icolf')
        if save_peaks_after:
            print('Saving self.icolf to disk with new column')
            self.savepeaks()
        return lx_modified

    # ----------------------------------------------------------------------
    def run_refine(self, points_step_space=None, npoints=None,
                   output_filename=None, use_cluster=False, pythonpath=None,
                   nthreads=None, nchunks=None, omega_binsize="auto", omega_step=None,
                   index_cache=None, cache_filename=None, use_cache=True,
                   save=True, verbose=True,
                   cpus_per_task=16, time_h=1, partition="nice", mem_G=32):
        """
        Refine every UBI candidate in self.pbpmap.

        use_cluster=True builds the index locally, writes it to the cache
        alongside self, and submits a slurm job that reloads both. The worker
        then does almost nothing but the kernel, so ask for a short walltime.
        """
        if "xpos_refined" not in self.icolf.titles:
            raise ValueError("Run get_origins() first -- no xpos_refined column")

        if use_cluster:
            if pythonpath is None:
                raise ValueError("Must supply pythonpath to run on cluster!")
            if points_step_space is not None or npoints is not None:
                raise ValueError("Choosing points to refine on the cluster is "
                                 "not implemented")
            # build the index here so the worker loads it instead of redoing it
            if verbose:
                print("Building index before submission")
            build_indexes(self, omega_binsize=omega_binsize, omega_step=omega_step,
                          cache_filename=cache_filename, use_cache=True,
                          verbose=verbose)
            print("Saving everything to disk")
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

        if nthreads is None:
            nthreads = max(cImageD11.cores_available() - 1, 1)
        numba.set_num_threads(nthreads)
        try:
            numba.set_parallel_chunksize(1)
        except AttributeError:
            pass

        idx = index_cache if index_cache is not None else build_indexes(
            self, omega_binsize=omega_binsize, omega_step=omega_step,
            cache_filename=cache_filename, use_cache=use_cache, verbose=verbose)
        pk = idx["peaks"]

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

        if nchunks is None:
            nchunks = nthreads
        nchunks = max(1, min(nchunks, len(points_ri)))

        maxlocal = int(vm_max_candidates(
            np.ascontiguousarray(self.sx_grid[self.mask], dtype=np.float64),
            np.ascontiguousarray(self.sy_grid[self.mask], dtype=np.float64),
            self.y0, self.ystep, idx["ymin"],
            idx["omega_partitions"], idx["dty_partitions"],
            idx["usin"], idx["ucos"]))
        maxlocal = max(maxlocal, 1)
        if verbose:
            mb = maxlocal * 200 * nchunks / 1e6
            print("worst-case peaks on a single ray: %d" % maxlocal)
            print("scratch: ~%.0f MB across %d chunks (virtual; touched lazily)" % (mb, nchunks))
            print("launching numba parallel refinement on %d threads" % nthreads)

        pars = self.icolf.parameters.get_parameters()
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
            idx["pstart"], idx["porder"], idx["pbpmap_ubis"], idx["nrj"],
            self.sx_grid, self.sy_grid, self.mask,
            idx["omega_partitions"], idx["dty_partitions"],
            idx["nom"], idx["ndty"], idx["usin"], idx["ucos"],
            idx["ymin"],
            pk["sc"], pk["fc"], pk["eta"], pk["sum_intensity"],
            pk["sinomega"], pk["cosomega"], pk["omega"], pk["dty"], pk["dtyi"],
            pk["xpos"], pk["gve_all"],
            self.ystep, self.y0, B0,
            float(pars["distance"]), float(pars["y_center"]), float(pars["y_size"]),
            float(pars["tilt_y"]), float(pars["z_center"]), float(pars["z_size"]),
            float(pars["tilt_z"]), float(pars["tilt_x"]),
            float(pars["o11"]), float(pars["o12"]), float(pars["o21"]), float(pars["o22"]),
            float(pars["t_x"]), float(pars["t_y"]), float(pars["t_z"]),
            float(pars["wedge"]), float(pars["chi"]), float(pars["wavelength"]),
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
    omega_partitions, dty_partitions, nom, ndty, usin, ucos, ymin,
    sc, fc, eta, sum_intensity, sinomega, cosomega, omega, dty, dtyi, xpos,
    gve_all,                                   # (N, 3)
    # geometry
    ystep, y0, B0,
    distance, y_center, y_size, tilt_y, z_center, z_size, tilt_z, tilt_x,
    o11, o12, o21, o22, t_x, t_y, t_z, wedge, chi, wavelength,
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
        m_sc = np.empty(maxlocal, dtype=np.float64)
        m_fc = np.empty(maxlocal, dtype=np.float64)
        m_om = np.empty(maxlocal, dtype=np.float64)
        m_dty = np.empty(maxlocal, dtype=np.float64)
        m_xp = np.empty(maxlocal, dtype=np.float64)
        m_eta = np.empty(maxlocal, dtype=np.float64)

        g_sel = np.empty(maxlocal, dtype=np.int64)        # merged peaks on the ray
        g_yd = np.empty(maxlocal, dtype=np.float64)
        c_sc = np.empty(maxlocal, dtype=np.float64)
        c_fc = np.empty(maxlocal, dtype=np.float64)
        c_om = np.empty(maxlocal, dtype=np.float64)
        c_xp = np.empty(maxlocal, dtype=np.float64)

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

        # ---- cyclic distribution over voxels ------------------------------
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
                    m_sc[t] = 0.0
                    m_fc[t] = 0.0
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
                    m_sc[lb] += sc[p] * wI
                    m_fc[lb] += fc[p] * wI
                    m_om[lb] += omega[p] * wI
                    m_dty[lb] += dty[p] * wI
                    m_xp[lb] += xpos[p] * wI
                    m_eta[lb] += eta[p] * wI
                for t in range(nlab):
                    s = m_sI[t]
                    m_sc[t] /= s
                    m_fc[t] /= s
                    m_om[t] /= s
                    m_dty[t] /= s
                    m_xp[t] /= s
                    m_eta[t] /= s

                # ---- re-apply the voxel mask to the merged peaks ----------
                ng = 0
                for t in range(nlab):
                    orad = np.radians(m_om[t])
                    d = abs(y0 - xi0 * np.sin(orad) - yi0 * np.cos(orad) - m_dty[t])
                    if d <= ystep:
                        g_sel[ng] = t
                        g_yd[ng] = d
                        c_sc[ng] = m_sc[t]
                        c_fc[ng] = m_fc[t]
                        c_om[ng] = m_om[t]
                        c_xp[ng] = m_xp[t]
                        ng += 1
                if ng == 0:
                    continue

                # merged g-vectors
                gvm = compute_gve(
                    c_sc[:ng], c_fc[:ng], c_om[:ng], c_xp[:ng],
                    distance=distance, y_center=y_center, y_size=y_size,
                    tilt_y=tilt_y, z_center=z_center, z_size=z_size,
                    tilt_z=tilt_z, tilt_x=tilt_x,
                    o11=o11, o12=o12, o21=o21, o22=o22,
                    t_x=t_x, t_y=t_y, t_z=t_z,
                    wedge=wedge, chi=chi, wavelength=wavelength,
                )  # (3, ng)

                # ---- reassign the merged peaks, eq (8) with e_hkl ---------
                nfit = 0
                for t in range(ng):
                    Gx = gvm[0, t]                          # G, eq (5)
                    Gy = gvm[1, t]
                    Gz = gvm[2, t]
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
                    Gx = gvm[0, ii]                     # G, eq (5)
                    Gy = gvm[1, ii]
                    Gz = gvm[2, ii]
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
                    dt = _det3(ubi_out) if finite else 0.0
                    sv = np.linalg.svd(ubi_out)[1] if finite else np.zeros(3)
                    smax = sv[0]
                    smin = sv[2]
                    rank3 = smin > smax * 3.0 * 2.220446049250313e-16
                    cnd = 1e300 if smin <= 0.0 else smax / smin
                    worth_fitting = finite and rank3 and (cnd < 1e14) and (dt > 0.0)
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
                    Gx = gvm[0, t]
                    Gy = gvm[1, t]
                    Gz = gvm[2, t]
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
                    Gx = gvm[0, ii]           # G, measured diffraction vector
                    Gy = gvm[1, ii]
                    Gz = gvm[2, ii]
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
                    Gx = gvm[0, ii]
                    Gy = gvm[1, ii]
                    Gz = gvm[2, ii]
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
    refine.run_refine(output_filename=output_filename, nthreads=nthreads,
                      cache_filename=getattr(refine, "index_filename", None))
    print("Refinement complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())