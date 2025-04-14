#!/usr/bin/env python
# coding: utf-8

# Try to build the point-by-point mapping code ...
from __future__ import print_function, division

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import sys
from ImageD11 import cImageD11

cImageD11.check_multiprocessing(patch=True)  # forkserver
try:
    import threadpoolctl
except ImportError:
    threadpoolctl = None

import multiprocessing
import time, random
import numpy as np
import subprocess

import h5py
import numba

# ignore Numba dot product performance warnings?
import warnings

warnings.simplefilter('ignore', category=numba.core.errors.NumbaPerformanceWarning)

import pprint

from skimage.filters import threshold_otsu
from skimage.morphology import convex_hull_image

from ImageD11 import sym_u, unitcell, parameters
import ImageD11.sinograms.dataset
import ImageD11.indexing
from ImageD11.columnfile import columnfile
from ImageD11.sinograms import geometry
from ImageD11.sinograms.roi_iradon import run_iradon
from ImageD11.sinograms.tensor_map import unitcell_to_b
from ImageD11.sinograms.sinogram import save_array


# GOTO - find somewhere!
# Everything written in Numba could move to C?
@numba.njit(boundscheck=True, cache = True)
def hkluniq(ubi, gx, gy, gz, eta, m, tol, hmax):
    """count uniq peaks - move to cImageD11 in the future"""
    # index from -hmax to hmax
    hcount = np.zeros((2 * hmax + 1, 2 * hmax + 1, 2 * hmax + 1, 2), dtype=np.int32)
    tcount = 0
    for i in range(len(gx)):
        if ~m[i]:
            continue
        H = ubi[0][0] * gx[i] + ubi[0][1] * gy[i] + ubi[0][2] * gz[i]
        K = ubi[1][0] * gx[i] + ubi[1][1] * gy[i] + ubi[1][2] * gz[i]
        L = ubi[2][0] * gx[i] + ubi[2][1] * gy[i] + ubi[2][2] * gz[i]
        ih = np.round(H)
        ik = np.round(K)
        il = np.round(L)
        dh2 = (H - ih) ** 2 + (K - ik) ** 2 + (L - il) ** 2
        if dh2 < tol * tol:
            # fixme sign(yl)
            if abs(ih) < hmax and abs(ik) < hmax and abs(il) < hmax:
                if eta[i] > 0:
                    s = 1
                else:
                    s = 0
                hcount[int(ih) + hmax, int(ik) + hmax, int(il) + hmax, s] = 1
            tcount += 1
    ucount = 0
    for i in range(hcount.size):
        if hcount.flat[i] > 0:
            ucount += 1
    return tcount, ucount

# stuff we need to compute g-vectors
# from ImageD11.transform
@numba.njit(cache = True)
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


@numba.njit(cache = True)
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


@numba.njit(cache = True)
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
    # logging.debug("detector_orientation = "+str(detector_orientation))
    flipped = np.dot(np.array(detector_orientation, np.float64),
                     peaks_on_detector)
    #
    # vec = np.array([np.zeros(flipped.shape[1]),  # place detector at zero,
    #                 # sample at -dist
    #                 flipped[1, :],             # x in search, frelon +z
    #                 flipped[0, :]], np.float64)     # y in search, frelon -y

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
    # TOOD: don't need this, use xl, yl, zl
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
def count_unique_peaks(hkl, etasign, dtyi):
    N = hkl.shape[1]  # Number of entries
    indices = np.zeros(N, dtype=np.int64)

    # Step 1: Create a list of tuples (h, k, l, etasign, dtyi)
    combined = []
    for i in range(N):
        key = (hkl[0, i], hkl[1, i], hkl[2, i], etasign[i], dtyi[i])
        combined.append((key, i))

    # Step 2: Sort lexicographically based on the tuple
    combined.sort()

    # Step 3: Create the unique indices array
    unique_map = {}
    rank = 0

    for sorted_key, original_index in combined:
        if sorted_key not in unique_map:
            unique_map[sorted_key] = rank
            rank += 1
        indices[original_index] = unique_map[sorted_key]

    return indices


@numba.njit(cache=True)
def count_unique_peaks_no_dtyi(hkl, etasign):
    N = hkl.shape[1]  # Number of entries
    indices = np.zeros(N, dtype=np.int64)

    # Step 1: Create a list of tuples (h, k, l, etasign, dtyi)
    combined = []
    for i in range(N):
        key = (hkl[0, i], hkl[1, i], hkl[2, i], etasign[i])
        combined.append((key, i))

    # Step 2: Sort lexicographically based on the tuple
    combined.sort()

    # Step 3: Create the unique indices array
    unique_map = {}
    rank = 0

    for sorted_key, original_index in combined:
        if sorted_key not in unique_map:
            unique_map[sorted_key] = rank
            rank += 1
        indices[original_index] = unique_map[sorted_key]

    return indices


@numba.njit(cache=True)
def merge(hkl, etasign, dtyi, sum_intensity, sc, fc, omega, dty, xpos_refined, eta):
    """
    merge peaks with the same (h, k, l, etasign, dtyi)
    weighted by sum_intensity
    """

    labels = count_unique_peaks(hkl, etasign, dtyi)
    wI = sum_intensity
    sI = np.bincount(labels, weights=wI)
    merged_sum_intensity = sI
    merged_sc = np.bincount(labels, weights=sc * wI) / sI
    merged_fc = np.bincount(labels, weights=fc * wI) / sI
    merged_omega = np.bincount(labels, weights=omega * wI) / sI
    merged_dty = np.bincount(labels, weights=dty * wI) / sI
    merged_xpos_refined = np.bincount(labels, weights=xpos_refined * wI) / sI
    merged_eta = np.bincount(labels, weights=eta * wI) / sI

    return merged_sum_intensity, merged_sc, merged_fc, merged_omega, merged_dty, merged_xpos_refined, merged_eta


@numba.njit(cache=True)
def get_voxel_idx(y0, xi0, yi0, sinomega, cosomega, dty, ystep):
    """
    get peaks at xi0, yi0
    basically just geometry.dty_values_grain_in_beam_sincos
    """

    ydist = np.abs(y0 - xi0 * sinomega - yi0 * cosomega - dty)
    idx = np.where(ydist <= ystep)[0]

    return idx, ydist


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
def weighted_lstsq_ubi_fit(ydist, gve, hkl):
    # run the weighted fit
    # a.T @ gve = h =>  gve.T @ a = h.T => a = np.linalg.pinv(gve.T) @ h.T, same for b and c
    w = (1. / (ydist + 1)).reshape(gve.shape[1], 1)
    a = w * gve.T
    b = w * hkl.T
    m, n = a.shape[-2:]
    rcond = np.finfo(b.dtype).eps * max(n, m)
    ubifitT, residuals, rank, sing_vals = np.linalg.lstsq(a, b, rcond=rcond)
    ubifit = ubifitT.T

    return w, ubifit, residuals, rank, sing_vals


@numba.njit(cache=True)
def gve_norm(gve):
    norms = np.zeros(gve.shape[1])
    for i in range(gve.shape[1]):
        gv = gve[:, i]
        norms[i] = np.sqrt(np.sum(np.power(gv, 2)))

    return norms


@numba.njit(cache=True)
def divide_where(arr1, arr2, out, wherearr):
    """
    Do arr1/arr2.
    In locations where wherearr == 0, return out instead
    """
    div = np.divide(arr1, arr2)
    return np.where(wherearr != 0, div, out)


@numba.njit(cache=True)
def ubi_to_unitcell(ubi):
    # fast numba version, can't use guvec version from tensor_map.py here unfortunately
    mt = np.dot(ubi, ubi.T)
    G = mt
    a, b, c = np.sqrt(np.diag(G))
    al = np.degrees(np.arccos(G[1, 2] / b / c))
    be = np.degrees(np.arccos(G[0, 2] / a / c))
    ga = np.degrees(np.arccos(G[0, 1] / a / b))
    return np.array([a, b, c, al, be, ga])


@numba.njit(cache=True)
def ubi_and_ucell_to_u(ubi, ucell):
    # compute B
    a, b, c = ucell[:3]
    ralpha, rbeta, rgamma = np.radians(ucell[3:])  # radians
    ca = np.cos(ralpha)
    cb = np.cos(rbeta)
    cg = np.cos(rgamma)
    g = np.full((3, 3), np.nan, np.float64)
    g[0, 0] = a * a
    g[0, 1] = a * b * cg
    g[0, 2] = a * c * cb
    g[1, 0] = a * b * cg
    g[1, 1] = b * b
    g[1, 2] = b * c * ca
    g[2, 0] = a * c * cb
    g[2, 1] = b * c * ca
    g[2, 2] = c * c
    gi = np.linalg.inv(g)
    astar, bstar, cstar = np.sqrt(np.diag(gi))
    betas = np.degrees(np.arccos(gi[0, 2] / astar / cstar))
    gammas = np.degrees(np.arccos(gi[0, 1] / astar / bstar))

    B = np.zeros((3, 3))

    B[0, 0] = astar
    B[0, 1] = bstar * np.cos(np.radians(gammas))
    B[0, 2] = cstar * np.cos(np.radians(betas))
    B[1, 0] = 0.0
    B[1, 1] = bstar * np.sin(np.radians(gammas))
    B[1, 2] = -cstar * np.sin(np.radians(betas)) * ca
    B[2, 0] = 0.0
    B[2, 1] = 0.0
    B[2, 2] = 1.0 / c

    u = np.dot(B, ubi).T
    return u


@numba.njit(cache=True)
def nb_choose_best(i, j, u, n, NY, ubiar,
                   minpeaks=6):
    # map of the unique scores
    uniq = np.ones((NY, NY), dtype='q')
    uniq.fill(minpeaks)  # peak cutorr
    npk = np.zeros((NY, NY), dtype='q')
    ubi = np.zeros((NY, NY, 3, 3), dtype='d')
    ubi.fill(np.nan)
    # store the index of the best row for each voxel
    best_rows = np.full(uniq.shape, -1, dtype=np.int64)
    # iterate through rows
    for k in range(i.size):
        # get the i and j values in the row
        ip = i[k]
        jp = j[k]
        # u[k] is the ubi in the row
        # if u in this row is better than the uniq at this voxel
        # (if no uniq at this voxel, it'll be minpeaks because we filled)
        if u[k] > uniq[ip, jp]:
            # uniq at this point becomes the uniq at this row
            uniq[ip, jp] = u[k]
            # npk at this point becomes the npk at this row
            npk[ip, jp] = n[k]
            # best_rows at this point becomes the index of this row
            best_rows[ip, jp] = k
            # ubi at this point becomes the ubi at this row
            for ii in range(3):
                for jj in range(3):
                    ubi[ip, jp, ii, jj] = ubiar[ii, jj, k]
    return uniq, npk, ubi, best_rows


@numba.njit(cache=True)
def nb_inv(mats, imats):
    for i in range(mats.shape[0]):
        for j in range(mats.shape[1]):
            if np.isnan(mats[i, j, 0, 0]):
                imats[i, j] = np.nan
            else:
                try:
                    imats[i, j] = np.linalg.inv(mats[i, j])
                except:
                    print(i, j, mats[i, j])
                    break
                    imats[i, j] = 42.


@numba.njit(cache=True)
def nb_inv_3d(mats, imats):
    for i in range(mats.shape[0]):
        for j in range(mats.shape[1]):
            for k in range(mats.shape[2]):
                if np.isnan(mats[i, j, k, 0, 0]):
                    imats[i, j, k] = np.nan
                else:
                    try:
                        imats[i, j, k] = np.linalg.inv(mats[i, j, k])
                    except:
                        imats[i, j, k] = np.nan


class PBPMap(columnfile):
    """Columnfile for point-by-point indexing and refinement results.
    The results of point-by-point indexing and refinement can be multi-valued
    e.g there can be more than one UBI found at a given map pixel
    Therefore this class stores pixel results in the following way:
    
    #  i  j  ntotal  nuniq  ubi00  ubi01  ubi02  ubi10  ubi11  ubi12  ubi20  ubi21  ubi22
    
    We have an extra method to get all the results at one pixel: self.get_pixel_mask()
    There's another concept: selecting a single best UBI at each pixel
    This will fill in self.best_UBI, self.best_npks, self.best_nuniq
    
    """

    # inherit init
    def __init__(self, *args, **kwargs):
        super(PBPMap, self).__init__(*args, **kwargs)

    def get_pixel_mask(self, i, j):
        """Return mask to columnfile with those specified i and j values"""
        return (self.i == i) & (self.j == j)

    @property
    def i_shift(self):
        return (self.i - self.i.min()).astype(int)

    @property
    def j_shift(self):
        return (self.j - self.j.min()).astype(int)

    @property
    def NI(self):
        return int(self.i_shift.max() - self.i_shift.min()) + 1

    @property
    def NJ(self):
        return int(self.j_shift.max() - self.j_shift.min()) + 1

    @property
    def NY(self):
        return max(self.NI, self.NJ)

    @property
    def ubi(self):
        """ubi computed on-demand from reshaping ubi columns.
        shape is (3, 3, nrows) to match previous PBPMap definition"""
        return np.vstack((
            self.ubi00, self.ubi01, self.ubi02, self.ubi10, self.ubi11, self.ubi12, self.ubi20, self.ubi21,
            self.ubi22)).reshape(3, 3, self.nrows)

    @property
    def eps(self):
        """eps computed on-demand from reshaping eps columns
        shape is (3, 3, nrows) to match previous PBPMap definition"""
        return np.vstack((
            self.eps00, self.eps01, self.eps02, self.eps10, self.eps11, self.eps12, self.eps20, self.eps21,
            self.eps22)).reshape(3, 3, self.nrows)

    @property
    def UBI(self):
        return self.ubi

    def plot_nuniq_hist(self):
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(figsize=(7,7), constrained_layout=True)
        ax.hist(self.nuniq, bins=np.arange(0.5, np.max(self.nuniq) + 0.51, 1))
        ax.set_xlabel('Unique spots per pixel')
        ax.set_ylabel('Count')
        plt.show()

    def choose_best(self, minpeaks=6):
        self.best_nuniq, self.best_npks, self.best_ubi, self.best_rows = nb_choose_best(
            self.i_shift, self.j_shift, self.nuniq.astype(int), self.ntotal.astype(int), self.NY, self.ubi, minpeaks)
        # see if we can get a self.best_eps too
        if 'eps00' in self.titles:
            best_eps = np.full_like(self.best_ubi, np.nan)
            all_eps = self.eps
            for i in range(best_eps.shape[0]):
                for j in range(best_eps.shape[1]):
                    if self.best_rows[i, j] > 0:
                        best_eps[i, j, :, :] = all_eps[:, :, self.best_rows[i, j]]
            self.best_eps = best_eps

    def plot_best(self, minpeaks=6):
        from matplotlib import pyplot as plt
        fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10, 5), constrained_layout=True)
        r = np.where(self.best_nuniq > minpeaks, self.best_nuniq, 0)
        axs[0].imshow(self.best_nuniq, origin="lower")
        axs[0].set_title('All pixels')
        axs[1].imshow(r, origin="lower")
        axs[1].set_title('Best pixels')
        fig.supxlabel('< Lab Y axis')
        fig.supylabel('Lab X axis')
        plt.show()


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
                 min_grain_npks=6
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

        # set default paths
        self.pbpmap_filename = self.dset.refmapfile  # pbp map input
        self.icolf_filename = self.dset.refpeaksfile  # icolf for refinement
        self.refinedmap_filename = self.dset.refoutfile  # refined pbp map output
        self.own_filename = self.dset.refmanfile  # myself as an H5

    def setmap(self, pbpmap):
        """Set an input PBPMap Python object to use"""
        # set up the grid to refine on
        # set up grids of X and Y positions in the sample reference frame
        # get an array of the integer I and J step values from the pmap
        # this is in step space
        # e.g -33, -32, -31, ... 32, 33
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
            # use the default from the init
            if refined:
                filename = self.refinedmap_filename
            else:
                filename = self.pbpmap_filename
        # make an empty pmap object
        pmap = PBPMap(new=True)
        # fill it using HDF load method
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
            # use the default from the init
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
            time.sleep(10)
            print('Continuing')
            print('To disable this prompt, set prompt_del=False when calling setpeaks()')

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
            # use the default one from the init
            icolf_filename = self.icolf_filename
        else:
            # update the default with user-supplied value
            self.icolf_filename = icolf_filename
        if os.path.exists(icolf_filename) and del_existing:
            os.remove(icolf_filename)
        ImageD11.columnfile.colfile_to_hdf(self.icolf, icolf_filename, compression=None)

    def loadpeaks(self, icolf_filename=None):
        """Load icolf from disk"""
        if icolf_filename is None:
            # use the default from the init
            icolf_filename = self.icolf_filename
        else:
            # update self with user-supplied value
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
        """Set a mask for choosing what to refine or not.
        You can choose whether to use self.colf (all peaks) or self.icolf (selected peaks)
        At the moment it does an iradon on the sinogram of all the 2D peaks in self.colf"""
        if use_singlemap:
            # take all non-nan singlemap values
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
            whole_sample_sino, xedges, yedges = np.histogram2d(dty, omega,
                                                               bins=[self.dset.ybinedges, self.dset.obinedges])
            shift, _ = geometry.sino_shift_and_pad(self.y0, len(self.ybincens), self.ybincens.min(), self.ystep)
            nthreads = cImageD11.cores_available()
            # make sure the shape is the same as sx_grid
            pad = self.sx_grid.shape[0] - whole_sample_sino.shape[0]
            whole_sample_recon = run_iradon(whole_sample_sino, self.dset.obincens, pad, shift, workers=nthreads)

            # we should be able to easily segment this using scikit-image
            recon_man_mask = whole_sample_recon

            # we can also override the threshold if we don't like it:
            # manual_threshold = 0.025

            if manual_threshold is None:
                thresh = threshold_otsu(recon_man_mask)
            else:
                thresh = manual_threshold

            binary = recon_man_mask > thresh

            chull = convex_hull_image(binary)

            whole_sample_mask = chull

        if doplot:
            from matplotlib import pyplot as plt
            fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, constrained_layout=True)
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
        # invoke the PBPMap method
        # then the single-valued map can be accessed from 
        if method == 'best':
            pbpmap.choose_best(**method_args)
            self.singlemap = pbpmap.best_ubi
        else:
            raise ValueError('Unsupported method!')

    def to_h5(self, filename=None, h5group='PBPRefine'):
        # save to h5 file
        if filename is None:
            # use init default
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

            # array name in h5: attribute name on object
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

            # path name in h5: attribute name on object
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

            # other pars we need for refinement
            pars = ['phase_name', 'hkl_tol_origins', 'hkl_tol_refine', 'hkl_tol_refine_merged', 'ds_tol',
                    'etacut', 'ifrac', 'forref', 'y0', 'min_grain_npks']

            for par in pars:
                try:
                    if getattr(self, par) is not None:
                        parent_group.attrs[par] = getattr(self, par)
                except AttributeError:
                    continue

    @classmethod
    def from_h5(cls, filename, h5group='PBPRefine'):
        # load the stuff in
        # then make an object
        manager_filename = filename

        with h5py.File(filename, "r") as hin:
            parent_group = hin[h5group]

            # dict to hold the arrays we loaded
            arrays = {}
            for array_name in parent_group.keys():
                arrays[array_name] = parent_group[array_name][:]  # load the array from disk

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

            pars = ['phase_name', 'hkl_tol_origins', 'hkl_tol_refine', 'hkl_tol_refine_merged', 'ds_tol',
                    'etacut', 'ifrac', 'forref', 'y0', 'min_grain_npks']
            pars_dict = {}
            for par in pars:
                try:
                    pars_dict[par] = parent_group.attrs[par]
                except (AttributeError, KeyError):
                    continue
        # load the dataset
        dset = ImageD11.sinograms.dataset.load(dsfile)
        refine_obj = cls(dset=dset, **pars_dict)
        refine_obj.own_filename = manager_filename
        for filename_attr, filename in filenames.items():
            setattr(refine_obj, filename_attr, filename)
        for array_attr, array in arrays.items():
            setattr(refine_obj, array_attr, array)

        # load the stuff we found
        print('Loading peaks')
        try:
            refine_obj.loadpeaks(refine_obj.icolf_filename)
        except AttributeError:
            pass

        print('Loading input map')
        try:
            refine_obj.loadmap(refine_obj.pbpmap_filename)
        except AttributeError:
            pass

        print('Loading output map')
        try:
            refine_obj.loadmap(refine_obj.refinedmap_filename, refined=True)
        except (AttributeError, OSError):
            pass
        return refine_obj

    def get_origins(self, guess_speed=True, guess_npks=10000, save_peaks_after=True):
        if 'xpos_refined' in self.icolf.titles:
            raise ValueError('We already have origins in self.icolf! Not recomputing')

        print('Getting gvecs...')
        # compute gves
        gve = np.column_stack((self.icolf.gx, self.icolf.gy, self.icolf.gz))

        nthreads = max( cImageD11.cores_available() - 1, 1 )
        numba.set_num_threads(nthreads)

        if guess_speed:
            guess_npks = min(guess_npks, self.icolf.nrows)
            print('Running test function twice for the first', guess_npks, 'peaks to guess speed...')
            print('This is because we get speed advantages with peaks on consecutive frames')

            idx = np.arange(guess_npks)
            sorter_reduced = np.lexsort((self.icolf.dty[idx], self.icolf.omega[idx]))

            print('First time to trigger NJIT')
            _ = compute_origins(self.singlemap, self.mask,
                                gve[idx], self.icolf.sinomega[idx], self.icolf.cosomega[idx], self.icolf.omega[idx],
                                self.icolf.dty[idx],
                                self.sx_grid, self.sy_grid,
                                self.y0, self.ystep, self.hkl_tol_origins,
                                sorter_reduced)

            print('Now we time')
            start = time.perf_counter()
            _ = compute_origins(self.singlemap, self.mask,
                                gve[idx], self.icolf.sinomega[idx], self.icolf.cosomega[idx], self.icolf.omega[idx],
                                self.icolf.dty[idx],
                                self.sx_grid, self.sy_grid,
                                self.y0, self.ystep, self.hkl_tol_origins,
                                sorter_reduced)
            end = time.perf_counter()
            time_taken = end - start
            time_for_one_peak = time_taken / guess_npks
            scaling_factor = 0.093  # experimental, from observation of guess vs reality
            time_for_all_peaks = time_for_one_peak * self.icolf.nrows * scaling_factor
            print('I estimate roughly', time_for_all_peaks, 'seconds for all peaks in self.icolf')
            print("That's", time_for_all_peaks / 60 / 60, 'hours')

        print('Lexsort...')
        # get the indices that sort all the arrays by omega then dty
        sorter = np.lexsort((self.icolf.dty, self.icolf.omega))

        print('Running numba computation on', nthreads, 'threads, may take a while!')

        lx_modified = compute_origins(self.singlemap, self.mask,
                                      gve, self.icolf.sinomega, self.icolf.cosomega, self.icolf.omega, self.icolf.dty,
                                      self.sx_grid, self.sy_grid,
                                      self.y0, self.ystep, self.hkl_tol_origins,
                                      sorter)

        self.icolf.addcolumn(lx_modified, 'xpos_refined')
        print('xpos_refined column added to self.icolf')
        if save_peaks_after:
            print('Saving self.icolf to disk with new column')
            self.savepeaks()

    def run_refine(self, points_step_space=None, npoints=None, output_filename=None, use_cluster=False,
                   pythonpath=None):
        """
        Call ImageD11.sinograms.point_by_point.refine_map to refine
        If you say use_cluster, will save the refinement object to disk
        Then submit a job to slurm to refine
        """
        if use_cluster:
            if pythonpath is None:
                raise ValueError('Must supply pythonpath to run refinement on cluster!')
            if points_step_space is not None or npoints is not None:
                raise ValueError("Choosing points to refine on the cluster is not yet implemented :(")
            print('Saving everything to disk')
            self.to_h5()
            print('Making bash script')
            bash_script_path = prepare_refine_bash(self, pythonpath, output_filename=output_filename)
            print('Made bash script at', bash_script_path)
            print('Submitting to cluster')
            submit_command = "sbatch {}".format(bash_script_path)
            sbatch_submit_result = subprocess.run(submit_command, capture_output=True, shell=True).stdout.decode(
                "utf-8")
            print(sbatch_submit_result)
        else:
            # get columns for refine.pbpmap
            # which contain the ri, rj points
            # right now it's indexed in step space
            ri_col, rj_col = geometry.step_to_recon(self.pbpmap.i, self.pbpmap.j, self.mask.shape)

            if points_step_space is None:
                # get the list of peak indices masks in reconstruction space
                points_recon_space = np.array(np.nonzero(self.mask)).T
                if npoints is not None:
                    points_recon_space = points_recon_space[:npoints]
                    # print(points_recon_space)
            else:
                # convert input points to reconstruction space, then pass to the function
                points_recon_space = [geometry.step_to_recon(si, sj, self.mask.shape) for (si, sj) in points_step_space]

            # columnfile by [3, 3, (ri, rj)]
            all_pbpmap_ubis = self.pbpmap.ubi

            pars = self.icolf.parameters.get_parameters()

            dummy_var = np.eye(3)
            uc = unitcell.unitcell_from_parameters(self.icolf.parameters)
            B0 = unitcell_to_b(uc.lattice_parameters, dummy_var)

            nthreads = max( cImageD11.cores_available() - 1, 1 )
            numba.set_num_threads(nthreads)
            print('Launching Numba parallel refinement on', nthreads, 'threads')

            final_ubis, final_eps, final_npks, final_nuniq = refine_map(points_recon_space, all_pbpmap_ubis, ri_col,
                                                                        rj_col,
                                                                        self.sx_grid, self.sy_grid, self.mask,
                                                                        self.icolf.sc, self.icolf.fc, self.icolf.eta,
                                                                        self.icolf.sum_intensity, self.icolf.sinomega,
                                                                        self.icolf.cosomega, self.icolf.omega,
                                                                        self.icolf.dty,
                                                                        self.icolf.dtyi, self.icolf.xpos_refined,
                                                                        self.ystep, self.y0,
                                                                        B0,
                                                                        pars['distance'], pars['y_center'],
                                                                        pars['y_size'],
                                                                        pars['tilt_y'],
                                                                        pars['z_center'], pars['z_size'],
                                                                        pars['tilt_z'],
                                                                        pars['tilt_x'],
                                                                        pars['o11'], pars['o12'], pars['o21'],
                                                                        pars['o22'],
                                                                        pars['t_x'], pars['t_y'], pars['t_z'],
                                                                        pars['wedge'],
                                                                        pars['chi'],
                                                                        pars['wavelength'],
                                                                        tol=self.hkl_tol_refine,
                                                                        merge_tol=self.hkl_tol_refine_merged,
                                                                        min_grain_npks=int(self.min_grain_npks)
                                                                        )
            # now we match how the indexer returns dodgy values
            # mask nan 3x3 entires to identity matrix
            final_ubis[:, :, np.isnan(final_ubis[0, 0, :])] = np.eye(3)[..., np.newaxis]
            # replace nans with 0 in the eps
            final_eps[:, :, np.isnan(final_eps[0, 0, :])] = 0
            # replace nans with 0 in the npks
            final_npks[np.isnan(final_npks)] = 0
            # replace nans with 0 in the nuniq
            final_nuniq[np.isnan(final_nuniq)] = 0

            # make a PBPMap to contain the output, then save to disk
            output_map = PBPMap(new=True)
            output_map.nrows = final_ubis.shape[2]
            final_ubis_cols = final_ubis.reshape(9, final_ubis.shape[2])
            final_eps_cols = final_eps.reshape(9, final_eps.shape[2])
            ubi00, ubi01, ubi02, ubi10, ubi11, ubi12, ubi20, ubi21, ubi22 = final_ubis_cols
            eps00, eps01, eps02, eps10, eps11, eps12, eps20, eps21, eps22 = final_eps_cols

            # add the columns
            output_map.addcolumn(self.pbpmap.i, 'i')
            output_map.addcolumn(self.pbpmap.j, 'j')
            output_map.addcolumn(final_npks, 'ntotal')
            output_map.addcolumn(final_nuniq, 'nuniq')
            output_map.addcolumn(ubi00, 'ubi00')
            output_map.addcolumn(ubi01, 'ubi01')
            output_map.addcolumn(ubi02, 'ubi02')
            output_map.addcolumn(ubi10, 'ubi10')
            output_map.addcolumn(ubi11, 'ubi11')
            output_map.addcolumn(ubi12, 'ubi12')
            output_map.addcolumn(ubi20, 'ubi20')
            output_map.addcolumn(ubi21, 'ubi21')
            output_map.addcolumn(ubi22, 'ubi22')
            output_map.addcolumn(eps00, 'eps00')
            output_map.addcolumn(eps01, 'eps01')
            output_map.addcolumn(eps02, 'eps02')
            output_map.addcolumn(eps10, 'eps10')
            output_map.addcolumn(eps11, 'eps11')
            output_map.addcolumn(eps12, 'eps12')
            output_map.addcolumn(eps20, 'eps20')
            output_map.addcolumn(eps21, 'eps21')
            output_map.addcolumn(eps22, 'eps22')

            self.refinedmap = output_map

            print('Saving refined map to disk')
            self.savemap(filename=output_filename, refined=True)


def prepare_refine_bash(pbp_object, id11_code_path, output_filename):
    ds = pbp_object.dset
    pbp_object_file = pbp_object.own_filename

    slurm_pbp_path = os.path.join(ds.analysispath, "slurm_pbp_refine")

    if not os.path.exists(slurm_pbp_path):
        os.mkdir(slurm_pbp_path)

    bash_script_path = os.path.join(slurm_pbp_path, ds.dsname + '_pbp_refine_slurm.sh')
    python_script_path = os.path.join(id11_code_path, "ImageD11/nbGui/S3DXRD/run_pbp_refine.py")
    outfile_path = os.path.join(slurm_pbp_path, ds.dsname + '_pbp_refine_slurm_%A.out')
    errfile_path = os.path.join(slurm_pbp_path, ds.dsname + '_pbp_refine_slurm_%A.err')
    log_path = os.path.join(slurm_pbp_path,
                            ds.dsname + '_pbp_refine_slurm_$SLURM_JOB_ID.log')

    # python 2 version
    bash_script_string = """#!/bin/bash
#SBATCH --job-name=pbp_refine
#SBATCH --output={outfile_path}
#SBATCH --error={errfile_path}
#SBATCH --time=48:00:00
#SBATCH --partition=nice-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64

#
date
source /cvmfs/hpc.esrf.fr/software/packages/linux/x86_64/jupyter-slurm/latest/envs/jupyter-slurm/bin/activate
echo OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH={id11_code_path} python3 {python_script_path} {pbp_object_file} > {log_path} 2>&1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH={id11_code_path} python3 {python_script_path} {pbp_object_file} > {log_path} 2>&1
date
        """.format(outfile_path=outfile_path,
                   errfile_path=errfile_path,
                   python_script_path=python_script_path,
                   id11_code_path=id11_code_path,
                   pbp_object_file=pbp_object_file,
                   log_path=log_path)

    with open(bash_script_path, "w") as bashscriptfile:
        bashscriptfile.writelines(bash_script_string)

    return bash_script_path


@numba.njit(cache=True)
def unique_with_counts(a):
    # https://github.com/numba/numba/pull/2959/commits/6657709adcf128d97eaaf8371d9106d48e2360ba
    b = np.sort(a.ravel())
    unique = list(b[:1])
    counts = [1 for _ in unique]
    for x in b[1:]:
        if x != unique[-1]:
            unique.append(x)
            counts.append(1)
        else:
            counts[-1] += 1
    return np.array(unique), np.array(counts)


@numba.njit(parallel=True,cache=True)
def compute_origins(singlemap, sample_mask,
                    gve, sinomega, cosomega, omega, dty,
                    sx_grid, sy_grid,
                    y0, ystep, hkl_tol,
                    sorter):
    """
    Compute origins of diffraction for each peak
    For each dty, omega value, we cast a ray through singlemap
    Along the ray, we get the UBIs in singlemap
    If the UBI indexes the peak, we consider it
    For each peak, we then work out the weighted origin position in x
    Adapted from Axel codebase
    """
    # set up empty array to hold corrected lx positions
    lx_modified = np.zeros_like(sinomega)

    # the game is to loop over omega as the outside loop
    # then loop over dty as the inside loop
    # this way we can compute ygrid from sinomega and cosomega as minimally as possible
    # to do this, we have the variable sorter
    # which is an index to all the cf columns that sorts in the order we want

    # get the sorted omegas
    sorted_omegas = omega[sorter]
    # get the unique omegas and the counts
    unique_omegas, omega_counts = unique_with_counts(sorted_omegas)
    # get the grouped indices to the peaks by incrementing omega value
    omega_idx_groups = np.split(sorter, omega_counts.cumsum())
    # omega_idx_groups is a list of peak indices, grouped and sorted by unique omega value
    # e.g omega_idx_groups[0] will contain all peak indices where omega = -180.00
    # and omega_idx_groups[1] will contain all peak indices where omega = -179.00

    # iterate over omega index groups
    # we do this in parallel for speed using numba.prange
    for omega_idx_group_idx in numba.prange(len(omega_idx_groups)):
        omega_idx_group = omega_idx_groups[omega_idx_group_idx]
        if len(omega_idx_group) > 0:
            # get an index that can get us the sinomega and cosomega value for this group
            omega_idx = omega_idx_group[0]

            so = sinomega[omega_idx]
            co = cosomega[omega_idx]

            # geometry.dty_values_grain_in_beam_sincos
            # for each grid point (sx, sy) get a dty value that brings the grid point
            # into the beam at this omega angle

            ygrid = y0 - sx_grid * so - sy_grid * co

            # now within this group we need to iterate over dty
            # get the unique dtys and the counts in this group
            dty_in_group = dty[omega_idx_group]

            unique_dtys, dty_counts = unique_with_counts(dty_in_group)
            # make index groups for dty
            dty_idx_groups = np.split(omega_idx_group, dty_counts.cumsum())

            for dty_idx_group in dty_idx_groups:
                if len(dty_idx_group) > 0:
                    this_dty = dty[dty_idx_group[0]]

                    # work out ydist etc
                    ydist = np.abs(this_dty - ygrid)
                    ray_mask = ydist <= ystep
                    final_mask = ray_mask * sample_mask
                    n_ray_px = final_mask.sum()

                    if n_ray_px > 0:
                        # get final_mask non-zero indices
                        ii, jj = np.nonzero(final_mask)
                        # meaningful weights
                        weights = (1 / (ydist + 1))

                        # now iterate over the peaks
                        for peak_idx in dty_idx_group:
                            # get the g-vector
                            g_vector = gve[peak_idx]

                            sx_weight_prod_cumsum = 0
                            sy_weight_prod_cumsum = 0
                            weight_cumsum = 0

                            for i, j in zip(ii, jj):
                                # active pixel
                                # get the orientation at this pixel
                                orien = singlemap[i, j, :, :]
                                # does orientation index the peak?
                                hklf = orien.dot(g_vector)
                                hklr = np.round(hklf)
                                err = ((hklf - hklr) ** 2).sum()

                                if err < hkl_tol:
                                    sx_weight_prod_cumsum += sx_grid[i, j] * weights[i, j]
                                    sy_weight_prod_cumsum += sy_grid[i, j] * weights[i, j]
                                    weight_cumsum += weights[i, j]

                            if weight_cumsum > 0:
                                sx_grid_corrected = sx_weight_prod_cumsum / weight_cumsum
                                sy_grid_corrected = sy_weight_prod_cumsum / weight_cumsum

                                lx_corrected = sx_grid_corrected * co - sy_grid_corrected * so
                                lx_modified[peak_idx] = lx_corrected

    return lx_modified


@numba.njit(parallel=True, cache=True)
def refine_map(refine_points, all_pbpmap_ubis, ri_col, rj_col, sx_grid, sy_grid, mask,  # refinement stuff
               sc, fc, eta, sum_intensity, sinomega, cosomega, omega, dty, dtyi, xpos,  # icolf columns
               ystep, y0,
               B0,
               distance, y_center, y_size, tilt_y, z_center, z_size, tilt_z, tilt_x,
               o11, o12, o21, o22,
               t_x, t_y, t_z, wedge, chi, wavelength,
               tol=0.1, merge_tol=0.05, min_grain_npks=6, gnoise_std=1e-4
               ):
    """
    Refine UBIs in all_pbpmap_ubis
    Adapted from Axel codebase

    refine_points: list or Nx2 array of points in reconstruction space (ri, rj) to refine at
    all_pbpmap_ubis: 3x3xM array of many-valued UBIs from point-by-point map
    ri_col, rj_col: length-M array of ri and rj coordinates (reconstruction space) for each UBI in all_pbpmap_ubis
    sx_grid, sy_grid: same as before (sample space sx, sy positions)
    a bunch of icolf columns
    ystep, y0: same as before
    B0: B matrix of dzero_cell
    then a bunch of pars from the parameter file
    then assignment tolerances:
    tol is the assignment hkl tolerance for non-merged peaks
    merge_tol is the assignment hkl tolerance for merged peaks
    """
    # make a non-awkward container for results
    final_ubis = np.full_like(all_pbpmap_ubis, np.nan)
    final_eps = np.full_like(all_pbpmap_ubis, np.nan)
    final_npks = np.full(all_pbpmap_ubis.shape[2], np.nan)
    final_nuniq = np.full(all_pbpmap_ubis.shape[2], np.nan)

    # iterate through point indices
    # they are in reconstruction space

    npoints = len(refine_points)

    all_idx = np.arange(len(ri_col))

    for refine_idx in numba.prange(npoints):
        ri, rj = refine_points[refine_idx]
        if mask[ri, rj]:
            # mask is valid at this pixel
            # print('at ri rj', ri, rj)

            # mask all_ubis by the pbpmap points
            pbpmap_mask = (ri_col == ri) & (rj_col == rj)
            pbpmap_idx = all_idx[pbpmap_mask]

            # get ubis at this point
            ubis_here = all_pbpmap_ubis[:, :, pbpmap_mask]
            # check that not all ubis are nan
            if np.all(np.isnan(ubis_here[0, 0, :])):
                continue

            # get xi0, xi0 at this point
            xi0 = sx_grid[ri, rj]
            yi0 = sy_grid[ri, rj]

            # get a mask to the peaks at this point
            # this is basically geometry.dtyimask_from_sincos
            # but we already have x, y

            # boolean masking is very slow for large arrays, even in Numba
            # instead, we get the indices we want
            # see more: https://stackoverflow.com/questions/46041811/performance-of-various-numpy-fancy-indexing-methods-also-with-numba
            local_idx, _ = get_voxel_idx(y0, xi0, yi0, sinomega, cosomega, dty, ystep)

            dtyi_local = dtyi[local_idx]
            sum_intensity_local = sum_intensity[local_idx]
            sc_local = sc[local_idx]
            fc_local = fc[local_idx]
            omega_local = omega[local_idx]
            dty_local = dty[local_idx]
            xpos_local = xpos[local_idx]
            eta_local = eta[local_idx]

            # get g-vectors for this voxel
            gve_voxel = compute_gve(sc_local, fc_local, omega_local, xpos_local,
                                    distance=distance, y_center=y_center, y_size=y_size, tilt_y=tilt_y,
                                    z_center=z_center, z_size=z_size, tilt_z=tilt_z, tilt_x=tilt_x,
                                    o11=o11, o12=o12, o21=o21, o22=o22,
                                    t_x=t_x, t_y=t_y, t_z=t_z, wedge=wedge, chi=chi, wavelength=wavelength)

            # iterate through the ubis at this voxel
            for ubi_idx in np.arange(ubis_here.shape[2]):
                ubi = ubis_here[:, :, ubi_idx]
                if np.isnan(ubi[0, 0]):
                    continue

                # we're scoring and assigning one UBI to a bunch of gves
                # all we need is to generate a mask

                tolsq = tol ** 2
                hklf = ubi.dot(gve_voxel)  # hkl float
                hkli = np.round(hklf).astype(np.int32)  # hkl int
                err = hklf - hkli
                sumsq = (err ** 2).sum(axis=0)
                grain_peak_mask = sumsq < tolsq

                grain_npks = np.sum(grain_peak_mask)
                hkl = hkli[:, grain_peak_mask]
                etasign = np.sign(eta_local[grain_peak_mask])

                # merge the assigned peaks in omega
                merged_sum_intensity, merged_sc, merged_fc, merged_omega, merged_dty, merged_xpos_refined, merged_eta = merge(
                    hkl,
                    etasign,
                    dtyi_local[
                        grain_peak_mask],
                    sum_intensity_local[
                        grain_peak_mask],
                    sc_local[
                        grain_peak_mask],
                    fc_local[
                        grain_peak_mask],
                    omega_local[
                        grain_peak_mask],
                    dty_local[
                        grain_peak_mask],
                    xpos_local[
                        grain_peak_mask],
                    eta_local[grain_peak_mask])

                merged_sinomega = np.sin(np.radians(merged_omega))
                merged_cosomega = np.cos(np.radians(merged_omega))

                # re-compute voxel peak mask with merged peaks
                local_idx_of_grain, ydist = get_voxel_idx(y0, xi0, yi0, merged_sinomega, merged_cosomega, merged_dty,
                                                          ystep)

                merged_sc_grain = merged_sc[local_idx_of_grain]
                merged_fc_grain = merged_fc[local_idx_of_grain]
                merged_omega_grain = merged_omega[local_idx_of_grain]
                merged_xpos_refined_grain = merged_xpos_refined[local_idx_of_grain]
                ydist_grain = ydist[local_idx_of_grain]
                merged_eta_grain = merged_eta[local_idx_of_grain]

                # re-compute merged g-vectors
                gve_voxel_merged = compute_gve(merged_sc_grain, merged_fc_grain, merged_omega_grain,
                                               merged_xpos_refined_grain,
                                               distance=distance, y_center=y_center, y_size=y_size, tilt_y=tilt_y,
                                               z_center=z_center, z_size=z_size, tilt_z=tilt_z, tilt_x=tilt_x,
                                               o11=o11, o12=o12, o21=o21, o22=o22,
                                               t_x=t_x, t_y=t_y, t_z=t_z, wedge=wedge, chi=chi, wavelength=wavelength)

                # reassign merged g-vectors with smaller merge_tol
                tolsq = merge_tol ** 2
                hklf = ubi.dot(gve_voxel_merged)
                hkli = np.round(hklf).astype(np.int32)
                err = hklf - hkli
                sumsq = (err ** 2).sum(axis=0)
                grain_peak_mask = sumsq < tolsq
                grain_npks = np.sum(grain_peak_mask)

                gve_grain = gve_voxel_merged[:, grain_peak_mask]
                hkl = hkli[:, grain_peak_mask]
                grain_ydist = ydist_grain[grain_peak_mask]

                # now we're ready to refine!

                if grain_npks > min_grain_npks:
                    w, ubifit, residuals, rank, sing_vals = weighted_lstsq_ubi_fit(grain_ydist, gve_grain, hkl)
                    
                    # check the quality of the fit
                    worth_fitting = (ubifit is not None) and (rank == 3) and (np.linalg.cond(ubifit) < 1e14) and (
                            np.linalg.det(ubifit) > 0) and (np.linalg.matrix_rank(ubifit) == 3)

                    # do we like the quality?
                    if worth_fitting:
                        ubi_out = ubifit
                    else:
                        ubi_out = ubi.copy()

                    # now reassign the (probably refined) UBI

                    tolsq = merge_tol ** 2
                    hklf = ubi_out.dot(gve_voxel_merged)
                    hkli = np.round(hklf).astype(np.int32)
                    err = hklf - hkli
                    sumsq = (err ** 2).sum(axis=0)
                    grain_peak_mask = sumsq < tolsq
                    grain_npks = np.sum(grain_peak_mask)

                    # compute the number of unique peaks
                    hkl_uniqcount = hkli[:, grain_peak_mask]
                    etasign_uniqcount = np.sign(merged_eta_grain[grain_peak_mask])
                    labels = count_unique_peaks_no_dtyi(hkl_uniqcount, etasign_uniqcount)
                    grain_nuniq = np.unique(labels).shape[0]

                    final_ubis[:, :, pbpmap_idx[ubi_idx]] = ubi_out
                    final_npks[pbpmap_idx[ubi_idx]] = grain_npks
                    final_nuniq[pbpmap_idx[ubi_idx]] = grain_nuniq

                    # if we like the fit quality, redetermine the g-vectors from the new peak mask
                    # then fit the strain tensor

                    if worth_fitting:
                        gve_grain_strainfit = gve_voxel_merged[:, grain_peak_mask]
                        hkl = hkli[:, grain_peak_mask]
                        ydist = ydist_grain[grain_peak_mask]
                        
                        ucell = ubi_to_unitcell(ubifit)
                        U = ubi_and_ucell_to_u(ubifit, ucell)
                        
                        gve0 = U.dot(B0).dot(hkl.astype(np.float64))
                        gTg0 = np.sum(gve_grain_strainfit * gve0, axis=0)
                        gTg = np.sum(gve_grain_strainfit * gve_grain_strainfit, axis=0)
                        directional_strain = (gTg0 / gTg) - 1

                        kappa = gve_grain_strainfit / gve_norm(gve_grain_strainfit)
                        kx, ky, kz = kappa

                        M = np.column_stack((kx * kx, ky * ky, kz * kz, 2 * kx * ky, 2 * kx * kz, 2 * ky * kz))

                        w = (1. / (ydist + 1)).reshape(gve_grain_strainfit.shape[1], 1)
                        # The noise in the directional strain now propagates according to the linear transform
                        a = np.sum(gve0 * (gnoise_std ** 2) * gve0, axis=0)
                        strain_noise_std = np.sqrt(divide_where(a, gTg ** 2, out=np.ones_like(gTg), wherearr=gTg))
                        w = w * (1. / strain_noise_std.reshape(w.shape))

                        w[directional_strain > np.mean(directional_strain) + np.std(
                            directional_strain) * 3.5] = 0  # outliers
                        w[directional_strain < np.mean(directional_strain) - np.std(
                            directional_strain) * 3.5] = 0  # outliers

                        w = w / np.max(w)
                        a = w * M
                        b = w.flatten() * directional_strain
                        m, n = a.shape[-2:]
                        rcond = np.finfo(b.dtype).eps * max(n, m)
                        eps_vec = np.linalg.lstsq(a, b, rcond=rcond)[0].T
                        sxx, syy, szz, sxy, sxz, syz = eps_vec
                        eps_tensor = np.array([[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]])

                        final_eps[:, :, pbpmap_idx[ubi_idx]] = eps_tensor

    return final_ubis, final_eps, final_npks, final_nuniq


def get_local_gv(si, sj, ystep, omega, sinomega, cosomega, xl, yl, zl):
    """
    Given a point (si, sj) in step space, compute the x values in the lab frame for each omega angle
    Then subtract them from xl to account for the change in the diffraction origin point.
    Then re-compute the g-vectors.
    """
    # compute the sample coordinates for this pixel
    sx, sy = geometry.step_to_sample(si, sj, ystep)

    # compute the position in the lab frame
    # geometry.sample_to_lab_sincos
    # lx = sxr
    # sxr = sx * cosomega - sy * sinomega
    # lx = sx * cosomega - sy * sinomega
    # we only care about the x axis
    x_offset = sx * cosomega - sy * sinomega
    # compute local offset
    xl_local = xl - x_offset
    # stack xl, yl, zl
    xyz_local = np.column_stack((xl_local, yl, zl))
    # hold computed g-vectors
    gv = np.empty((len(cosomega), 3))
    # compute g-vectors from xyz
    ImageD11.cImageD11.compute_gv(xyz_local,
                                  omega,
                                  parglobal.get('omegasign'),
                                  parglobal.get('wavelength'),
                                  parglobal.get('wedge'),
                                  parglobal.get('chi'),
                                  np.array((0., 0., 0.)),
                                  gv)

    gx = gv[:, 0]
    gy = gv[:, 1]
    gz = gv[:, 2]
    return gv, gx, gy, gz


def idxpoint(
        si,
        sj,
        isel,
        omega,
        sinomega,
        cosomega,
        dtyi,
        xl,
        yl,
        zl,
        eta,
        ystep=2.00,
        y0=0.0,
        ymin=-2.00,
        minpks=1000,
        hkl_tol=0.1,
        ds_tol=0.005,
        cosine_tol=np.cos(np.radians(90 - 0.1)),
        forgen=None,
        hmax=-1,  # auto
        uniqcut=0.75,
):
    """
    Indexing function called at one point in space.
    Selects peaks from the sinogram and attempts to index.
    More of this could be indexing.py I guess.

    si, sj = position in step space (origin is rotation axis)

    Effectively a columnfile:
        isel = column of rings to use
        sinomega = column of sinomega
        cosomega = column of cosomega
        dtyi = column of dty on integer bins
        gx/gy/gz/eta = peak co-ordinates

    ystep, y0 = compute the dtyi bins to select peaks (dtymask function above)

    forgen = the rings to use for INDEXING / orientation searching (a speedup effect)

    minpks = number of peaks to find with rings
    hkl_tol = ImageD11 tolerance of indexed == | h-int(h) | < hkl_tol
    cosine_tol = ImageD11 tolerance for finding pairs of peaks

    uniqcut = unique peaks cutoff (uniq> ucut * max are returned)

    returns a list of :
        (npks, nuniq, ubi)
    """
    cImageD11.cimaged11_omp_set_num_threads(1)
    # args    isel, omega, dtyi, gx, gy, gz
    # ring mask
    mr = isel
    # mask all the peaks by dtyi and omega (for scoring)

    m = geometry.dtyimask_from_step_sincos(si, sj, sinomega, cosomega, dtyi, y0, ystep, ymin)
    # combine masks together (for indexing)
    local_mask = m & mr

    # we need to index using peaks[m & mr]
    # then score with peaks[m]
    # so we need to compute g-vectors for peaks[m]
    # then mask those for indexing via [mr]

    # compute g-vectors for peaks[m]
    gv, gx, gy, gz = get_local_gv(si, sj, ystep, omega[m], sinomega[m], cosomega[m], xl[m], yl[m], zl[m])
    eta_local = eta[m]

    # mask g-vectors to [m & mr]
    gv_local_mask = gv[local_mask[m], :]

    # index with masked g-vectors
    ind = ImageD11.indexing.indexer(
        unitcell=ucglobal,
        wavelength=parglobal.get("wavelength"),
        gv=gv_local_mask,
    )
    ind.minpks = minpks
    ind.hkl_tol = hkl_tol
    ind.cosine_tol = cosine_tol  # degrees
    ind.ds_tol = ds_tol
    try:
        ind.assigntorings()
    except ValueError:
        return [
            (0, 0, np.eye(3)),
        ]
    for ind.ring_1 in forgen:
        for ind.ring_2 in forgen:
            ind.find()
            try:
                ind.scorethem()
            except (TypeError, ValueError):
                continue
    if ImageD11.indexing.loglevel <= 3:
        sys.stdout.flush()
    #flake8: global symglobal
    ind.ubis = [sym_u.find_uniq_u(ubi, symglobal) for ubi in ind.ubis]
    # count the unique peaks per grain
    uniqs = np.empty((len(ind.ubis), 2), int)
    for i, ubi in enumerate(ind.ubis):
        # we are scoring with peaks[m]
        # here, gx, gy, gz are already computed via peaks[m]
        # so no mask required (passing ones)
        uniqs[i] = hkluniq(ubi, gx, gy, gz, eta_local, np.ones_like(eta_local, dtype=bool), hkl_tol, hmax)
    if len(ind.ubis) == 0:
        return [
            (0, 0, np.eye(3)),
        ]
    if len(ind.ubis) == 1:
        return [
            (uniqs[0][0], uniqs[0][1], ind.ubis[0]),
        ]
    order = np.lexsort(uniqs.T)[::-1]  # decreasing ...
    best = uniqs[order[0]][1]  # number of unique spots
    last = int(best * uniqcut)
    grains = []
    for i in order:
        ntotal, nuniq = uniqs[i]
        if nuniq < last:
            break
        grains.append((ntotal, nuniq, ind.ubis[i]))
    return grains


def proxy(args):
    """Wrapper function for multiprocessing"""
    i, j, idxopts = args
    sm = colglobal
    g = idxpoint(
        i,
        j,
        sm["isel"],
        sm["omega"],
        sm["sinomega"],
        sm["cosomega"],
        sm["dtyi"],
        sm["xl"],
        sm["yl"],
        sm["zl"],
        sm["eta"],
        **idxopts
    )
    return i, j, g


# one of these per process:
ucglobal = None
symglobal = None
parglobal = None
colglobal = None


def initializer(parfile, phase_name, symmetry, colfile, loglevel=3):
    global ucglobal, symglobal, parglobal, colglobal
    if threadpoolctl is not None:
        threadpoolctl.threadpool_limits(limits=1)
    parglobal = parameters.read_par_file(parfile, phase_name=phase_name)
    ucglobal = unitcell.unitcell_from_parameters(parglobal)
    symglobal = sym_u.getgroup(symmetry)()
    colglobal = ImageD11.columnfile.mmap_h5colf(colfile)
    colglobal.isel = colglobal.isel.astype(bool)
    ImageD11.indexing.loglevel = loglevel


class PBP:
    def __init__(
            self,
            parfile,
            dset,
            hkl_tol=0.01,
            fpks=0.7,
            ds_tol=0.005,
            etacut=0.1,
            ifrac=None,
            gmax=5,
            cosine_tol=np.cos(np.radians(90 - 0.1)),
            y0=0,
            symmetry="cubic",
            loglevel=3,
            foridx=None,
            forgen=None,
            uniqcut=0.75,
            phase_name=None
    ):
        """
        parfile = ImageD11 parameter file (for the unit cell + geometry)
        dsname = name of dset file, or object with "ybincens" array
        """
        self.parfile = parfile
        self.dset = dset
        self.hkl_tol = hkl_tol
        self.fpks = fpks
        self.ds_tol = ds_tol
        self.etacut = etacut
        self.symmetry = symmetry
        self.foridx = foridx
        self.forgen = forgen
        self.gmax = gmax
        self.ybincens = np.array(dset.ybincens)
        if ifrac is None:
            self.ifrac = 1.0 / len(self.ybincens)
        else:
            self.ifrac = ifrac
        self.ystep = dset.ystep
        self.y0 = y0
        self.ymin = self.ybincens.min()
        self.uniqcut = uniqcut
        self.cosine_tol = cosine_tol
        self.loglevel = loglevel
        self.phase_name = phase_name

    def setpeaks(self, colf, icolf_filename=None):
        """
        This is called on the main parent process

        The peaks will be written to file
        """
        if icolf_filename is None:
            icolf_filename = self.dset.icolfile
        # Load the peaks
        colf.parameters.loadparameters(self.parfile, phase_name=self.phase_name)
        if "ds" not in colf.titles:
            colf.updateGeometry()  # for ds
        #
        uc = unitcell.unitcell_from_parameters(colf.parameters)
        uc.makerings(colf.ds.max(), self.ds_tol)
        # peaks that are on rings
        sel = np.zeros(colf.nrows, bool)
        # rings to use for indexing
        isel = np.zeros(colf.nrows, bool)
        npks = 0
        if self.foridx is None:
            self.foridx = range(len(uc.ringds))
        hmax = 0
        for i in range(len(uc.ringds)):
            ds = uc.ringds[i]
            hkls = uc.ringhkls[ds]
            if i in self.foridx:
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
        if self.forgen is None:
            self.forgen = self.foridx

        isel = isel & ((np.abs(np.sin(np.radians(colf.eta)))) > self.etacut)

        colf.addcolumn(isel, "isel")  # peaks selected for indexing

        dtyi = geometry.dty_to_dtyi(colf.dty, self.ystep, self.ymin)
        colf.addcolumn(dtyi, "dtyi")

        # cache these to speed up selections later
        colf.addcolumn(np.sin(np.radians(colf.omega)), "sinomega")
        colf.addcolumn(np.cos(np.radians(colf.omega)), "cosomega")

        # peaks that are on any rings
        self.colf = colf
        self.icolf = colf.copy()
        self.icolf.filter(sel)
        self.npks = npks
        if self.fpks < 1:
            self.minpks = int(npks * self.fpks)
        else:
            self.minpks = self.fpks
        print(
            "Using for indexing:",
            self.icolf.nrows,
            "npks, minpks, forgen",
            self.npks,
            self.minpks,
            self.forgen,
        )
        self.hmax = hmax
        if os.path.exists(icolf_filename):
            os.remove(icolf_filename)
        ImageD11.columnfile.colfile_to_hdf(self.icolf, icolf_filename, compression=None)
        self.icolf_filename = icolf_filename
        self.uc = uc

    def iplot(self, skip=1):
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
        ax[0].vlines(self.uc.ringds, 1e4, 3e4, color='red', label=self.phase_name)
        return f, ax

    def point_by_point(
            self,
            grains_filename,
            icolf_filename=None,  # use self
            nprocs=None,
            gridstep=1,
            debugpoints=None,
            loglevel=3,
    ):
        """
        grains_filename = output file
        icolf_filename = hdf5file to write. Allows mmap to be used.
        """
        self.loglevel = loglevel
        start = time.time()
        # FIXME - look at dataset ybincens and roi_iradon output_size ?
        if nprocs is None:
            nprocs = cImageD11.cores_available()

        idxopt = {
            "ystep": self.ystep,
            "y0": self.y0,
            "ymin": self.ymin,
            "minpks": self.minpks,
            "hkl_tol": self.hkl_tol,
            "ds_tol": self.ds_tol,
            "cosine_tol": self.cosine_tol,
            "forgen": self.forgen,
            "hmax": self.hmax,
            "uniqcut": self.uniqcut,
        }
        pprint.pprint(idxopt)

        if debugpoints is None:
            args = [(i, j, idxopt) for (i, j) in geometry.step_grid_from_ybincens(self.ybincens, self.ystep, gridstep, self.y0)]
            random.shuffle(args)
        else:
            args = [(i, j, idxopt) for i, j in debugpoints]
        gmap = {}
        t0 = time.time()
        # main process is writing
        with open(grains_filename, "w") as gout:
            gout.write(
                "#  i  j  ntotal  nuniq  ubi00  ubi01  ubi02  ubi10  ubi11  ubi12  ubi20  ubi21  ubi22\n"
            )
            done = 0
            ng = 0
            with multiprocessing.Pool(
                    nprocs,
                    initializer=initializer,
                    initargs=(
                            self.parfile,
                            self.phase_name,
                            self.symmetry,
                            self.icolf_filename,
                            self.loglevel,
                    ),
            ) as p:
                for i, j, g in p.imap_unordered(proxy, args):
                    gmap[i, j] = g
                    done += 1
                    ng += len(g)
                    if done % nprocs == 0:
                        sys.stdout.flush()  # before!
                        dt = time.time() - t0
                        print(
                            "Done %7.3f %%, average grains/point %6.2f, %.3f /s/point, total %.1f /s"
                            % (
                                100.0 * done / len(args),
                                float(ng) / done,
                                dt / done,
                                dt,
                            ),
                            end="\r",
                        )
                    sys.stdout.flush()
                    for k in range(min(self.gmax, len(g))):
                        # only output gmax grains per point max
                        n, u, ubi = g[k]
                        # print(g[k])
                        gout.write("%d  %d  %d  %d  " % (i, j, n, u))
                        gout.write(("%f " * 9) % tuple(ubi.ravel()) + "\n")
                        gout.flush()

        end = time.time()
        print(end - start, "seconds", (end - start) / len(args), "s per point")
        print(end - t0, "seconds", (end - t0) / len(args), "s per point without setup")


if __name__ == "__main__":
    print("Enter main thread", time.ctime())
    dset = ImageD11.sinograms.dataset.load(sys.argv[1])
    ImageD11.cImageD11.cimaged11_omp_set_num_threads(
        ImageD11.cImageD11.cores_available()
    )
    cf2d = dset.get_cf_2d()
    ImageD11.cImageD11.cimaged11_omp_set_num_threads(1)

    print("Test run with defaults that will not work for you!!!")

    cf2d.filter(cf2d.Number_of_pixels > 15)

    fac = 1.5

    pbpmanager = PBP(
        dset.parfile,
        dset,
        hkl_tol=np.sqrt(4 * 4 + 4 * 4) * np.sin(np.radians(dset.ostep * fac)),
        fpks=0.7,
        ds_tol=0.005,
        etacut=0.1,
        ifrac=1e-3,
        gmax=20,
        cosine_tol=np.cos(np.radians(90 - dset.ostep * fac)),
        y0=0,  # check?
        symmetry="cubic",
        foridx=[0, 1, 2, 4, 5, 10],
        forgen=[5, 1],
        uniqcut=0.6,
    )

    pbpmanager.setpeaks(cf2d, sys.argv[2])
    pbpmanager.point_by_point(sys.argv[3], loglevel=3, gridstep=1)
    print("Exit main thread", time.ctime())
