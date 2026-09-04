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

import numba

# ignore Numba dot product performance warnings?
import warnings

warnings.simplefilter('ignore', category=numba.core.errors.NumbaPerformanceWarning)

import pprint

from ImageD11 import sym_u, unitcell, parameters
import ImageD11.sinograms.dataset
import ImageD11.indexing
from ImageD11.columnfile import columnfile
from ImageD11.peakselect import mask_rings_by_ifrac
from ImageD11.sinograms import geometry
from ImageD11.sinograms.voxel_mask import (
    VoxelSinoMasker, fill_voxel_idx, max_candidates as vm_max_candidates,
    default_index_filename,
    save_peak_index_cache, load_peak_index_cache,
    choose_omega_bins,
    _CACHE_VERSION)

# GOTO - find somewhere!
# Everything written in Numba could move to C?
@numba.njit(boundscheck=True, cache=True)
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


@numba.njit(cache=True)
def select_point_peaks(si, sj, ystep, y0, ymin,
                       omega_partitions, dty_partitions,
                       sinomega, cosomega, dty, dtyi,
                       sinomega_bins, cosomega_bins,
                       idx_buf, ydist_buf, out, relax):
    """Peaks that put step-space point (si, sj) in the beam.

    Wraps fill_voxel_idx() but with an optional stricter predicate (relax=False)
    The original geometry.dtyimask_from_step_sincos() predicate is exact equality of the dtyi bin
    The voxel_mask predicate is |.| <= ystep, which is the two nearest rows -- exactly 2x the peaks
     
    relax=False reproduces geometry.dtyimask_from_step_sincos exactly:
        dty_to_dtyi(y0 - sx sin(w) - sy cos(w)) == dtyi

    Fills idx_buf, ydist_buf and out (sel_buf) inplace


    Parameters
    ----------
    si, sj
        Sample indices in step space
    ystep, y0, ymin
        Passed to fill_voxel_idx()
    omega_partitions, dty_partitions
        See VoxelSinoMasker documentation
    sinomega, cosomega, dty, dtyi
        Columnfile columns, needed for the stricter predicate when relax=False
    sinomega_bins, cosomega_bins
        See VoxelSinoMasker documentation
    idx_buf, ydist_buf
        Buffers for fill_voxel_idx() to write into
    out
        sel_buf: stricter version of idx_buf when relax=False
    relax
        Use voxel_mask predicate (|.| <= ystep) instead of exact equality of the dtyi bin

    Returns
    -------
    m
        Number of peaks written into out, or -1 if nothing was written because the buffers were too small
    """

    sx = si * ystep
    sy = -sj * ystep          # geometry.step_to_sample: note the sign on sy
 
    n = fill_voxel_idx(sx, sy, y0, ystep, ymin,
                       omega_partitions, dty_partitions, dty,
                       sinomega, cosomega, sinomega_bins, cosomega_bins,
                       idx_buf, ydist_buf)
    if n < 0:  # buffer too small, catch and return -1 to the caller
        return -1
 
    # TODO: relax=True does not seem to be used or tested yet?
    if relax:  # match voxel_mask predicate, so just return the peaks that fill_voxel_idx() found
        for t in range(n):
            out[t] = idx_buf[t]
        return n
 
    m = 0
    inv_ystep = 1.0 / ystep
    for t in range(n):
        q = idx_buf[t]
        dty_calc = y0 - sx * sinomega[q] - sy * cosomega[q]
        # round() here matches numpy's round-half-to-even, as dty_to_dtyi uses
        if int(round((dty_calc - ymin) * inv_ystep)) == dtyi[q]:
            out[m] = q
            m += 1
    return m


def idxpoint(si, sj,
    omega, sinomega, cosomega, dty, dtyi, xl, yl, zl, eta,
    omega_partitions, dty_partitions,
    sinomega_bins, cosomega_bins,
    idx_buf, ydist_buf, sel_buf,
    ystep=2.00, y0=0.0, ymin=-2.00, minpks=1000, hkl_tol=0.1,
    ds_tol=0.005, cosine_tol=np.cos(np.radians(90 - 0.1)), forgen=None, hmax=-1,
    uniqcut=0.75, relax_mask=False):
    """Indexing function called at one point in space.
    
    Selects peaks from the sinogram and attempts to index.

    Parameters
    ----------
    si, sj
        Voxel indices in step space
    omega, sinomega, cosomega, dty, dtyi, xl, yl, zl, eta
        Columnfile columns, already permuted to partition order
    omega_partitions
        omega_partitions[i] is the index of the first (sorted) peak in omega bin i
    dty_partitions
        dty_partitions[i,j] is the index of the first (sorted) peak in omega bin i and dty bin j
    sinomega_bins, cosomega_bins
        For speed
    idx_buf
        Buffer - will be filled with the indices of the peaks that put (si, sj) in the beam
    ydist_buf
        Buffer - will be filled with the ydist of the peaks that put (si, sj) in the beam
    sel_buf
        Buffer - used if relax_mask=False to select the peaks that put (si, sj) in the beam and satisfy the stricter predicate
    ystep, y0, ymin
        Passed to fill_voxel_idx()
    minpks, optional
        number of peaks to find with rings, by default 1000
    hkl_tol, optional
        ImageD11 tolerance of indexed == | h-int(h) | < hkl_tol, by default 0.1
    ds_tol, optional
        Ring assignment tolerance, by default 0.005
    cosine_tol, optional
        Cosine tolerance for finding pairs of peaks, by default np.cos(np.radians(90 - 0.1))
    forgen, optional
        Rings to generate orientations, by default None
    hmax, optional
        Maximum hkl value, by default -1
    uniqcut, optional
        Unique peaks cutoff (uniq> ucut * max are returned), by default 0.75
    relax_mask, optional
        Use previously-used stricter predicate for ydist, by default False

    Returns
    -------
        list of (npks, nuniq, ubi)

    Raises
    ------
    ValueError
        If peak buffer too small for this point. This shouldn't happen if max_canidates was computed correctly
    """
    cImageD11.cimaged11_omp_set_num_threads(1)  # ty: ignore[unresolved-attribute]

    # fill idx_buf, ydist_buf, sel_buf
    # idx_buf and sel_buf are peak indices into the partitioned columns
    npk = select_point_peaks(
        si, sj, ystep, y0, ymin,
        omega_partitions, dty_partitions,
        sinomega, cosomega, dty, dtyi,          # the mmap'd columns, already permuted
        sinomega_bins, cosomega_bins,
        idx_buf, ydist_buf, sel_buf, relax_mask)
    if npk < 0:
        raise ValueError(
            "peak buffer too small at (%d, %d) -- maxlocal was computed for "
            "the standard grid; this point is outside it" % (si, sj))
    if npk == 0:
        return [(0, 0, np.eye(3))]
    # sel_buf is sized for max_candidates, but only the first npk entries are valid
    sel = sel_buf[:npk] # indices to icolf

    gv, gx, gy, gz = get_local_gv(
        si, sj, ystep,
        omega[sel], sinomega[sel], cosomega[sel],
        xl[sel], yl[sel], zl[sel]
    )
    eta_local = eta[sel]

    # index with masked g-vectors
    ind = ImageD11.indexing.indexer(
        unitcell=ucglobal,
        wavelength=parglobal.get("wavelength"),
        gv=gv,
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
        i, j,
        sm["omega"], sm["sinomega"], sm["cosomega"], sm["dty"], sm["dtyi"],
        sm["xl"], sm["yl"], sm["zl"], sm["eta"],
        omega_partitions=partglobal["omega_partitions"],
        dty_partitions=partglobal["dty_partitions"],
        sinomega_bins=partglobal["usin"],
        cosomega_bins=partglobal["ucos"],
        idx_buf=bufglobal["idx"],
        ydist_buf=bufglobal["ydist"],
        sel_buf=bufglobal["sel"],
        **idxopts
    )
    return i, j, g

# one of these per process:
ucglobal = None
symglobal = None
parglobal = None
colglobal = None
partglobal = None
bufglobal = None

def initializer(parfile, phase_name, symmetry, colfile, index_filename, loglevel=3):
    global ucglobal, symglobal, parglobal, colglobal, partglobal, bufglobal
    try:
        if threadpoolctl is not None:
            threadpoolctl.threadpool_limits(limits=1)
        parglobal = parameters.read_par_file(parfile, phase_name=phase_name)
        ucglobal = unitcell.unitcell_from_parameters(parglobal)
        symglobal = sym_u.getgroup(symmetry)()
        colglobal = ImageD11.columnfile.mmap_h5colf(colfile)
        ImageD11.indexing.loglevel = loglevel
 
        # mmap, not read: identical and read-only in every worker, so one
        # mapping through the page cache replaces one copy per process
        partglobal = load_peak_index_cache(index_filename, mmap=True)
        if partglobal is None:
            raise RuntimeError(
                "peak index %s is missing or stale -- rerun setpeaks()"
                % index_filename)
 
        n = partglobal["maxlocal"]
        if n is None:
            raise RuntimeError(
                "peak index %s has no maxlocal -- it predates _CACHE_VERSION "
                "%d. Rerun setpeaks()." % (index_filename, _CACHE_VERSION))

        bufglobal = {
            "idx": np.empty(n, np.int64),
            "ydist": np.empty(n, np.float64),
            "sel": np.empty(n, np.int64)}
    except Exception:
        # Pool discards initializer exceptions and respawns the worker, so a
        # failure here is an invisible infinite loop that looks like a
        # slowdown. Make it loud.
        import traceback
        print("INITIALIZER FAILED in pid %d" % os.getpid(), file=sys.stderr)
        traceback.print_exc()
        sys.stderr.flush()
        raise


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

        self.icolf_filename = None      # set by setpeaks
        self.index_filename = None      # set by setpeaks

    def _build_index(self, verbose=True):
        """Take the partition setpeaks built, size the buffers with
        max_candidates, and write it to disk ready for mmap in the workers.
        Returns the partition dictionary.
        """
        self.index_filename = default_index_filename(self)
 
        # The partition is already built in self._masker; we only need to
        # size the buffers and write it out for the workers.
        M = self._masker
        nom = int(M.sinomega_bins.size)          # number of omega bins
        ndty = int(M.dty_partitions.shape[1]) - 1  # bin edges - 1
 
        # Worst case peaks any point can pull out, so the workers size their
        # buffers exactly. Uses the refiner's |.| <= ystep predicate, which
        # is a superset of bin equality.
        pts = np.asarray(geometry.step_grid_from_ybincens(
            self.ybincens, self.ystep, 1, self.y0))
        sx, sy = geometry.step_to_sample(pts[:, 0], pts[:, 1], self.ystep)
        maxlocal = int(vm_max_candidates(
            np.ascontiguousarray(sx, dtype=np.float64),
            np.ascontiguousarray(sy, dtype=np.float64),
            self.y0, self.ystep, M.ymin,
            M.dty_partitions,
            M.sinomega_bins, M.cosomega_bins))
        maxlocal = max(maxlocal, 1)
 
        # This is what gets mem-mapped in the workers, so they can size their
        # buffers and select peaks. Also what makes the cluster path work.
        idx = {
            "order": M.peak_ordering,
            "omega_partitions": M.omega_partitions,
            "dty_partitions": M.dty_partitions,
            "usin": M.sinomega_bins, "ucos": M.cosomega_bins,
            "nom": nom, "ndty": ndty, "ymin": float(M.ymin),
            "ray_margin": 0.0, "dev": 0.0,     # refiner-only, unused here
            "maxlocal": maxlocal,
        }
        if verbose:
            print("indexing partition: %d omega bins x %d dty bins = %d "
                  "cells, %.1f peaks/cell, worst-case %d peaks/point"
                  % (nom, ndty, nom * ndty,
                     self.icolf.nrows / float(nom * ndty), maxlocal))
            if self.icolf.nrows < nom * ndty:
                print("  WARNING: fewer peaks than cells -- the index is "
                      "sparse and traversal will dominate. Check whether the "
                      "frm decode ran; the heuristic over-refines without it.")
        try:
            save_peak_index_cache(idx, self.index_filename)
            if verbose:
                print("saved peak index to %s" % self.index_filename)
        except OSError as e:
            print("could not write index cache (%s), continuing" % e)
        return idx

    def setpeaks(self, colf, icolf_filename=None, keep_colf=True):
        """
        Given a cf_2d in RAM (colf), make an icolf in RAM which contains the selection of peaks we want to index
        """
        # Load parameters if not already loaded
        colf.parameters.loadparameters(self.parfile, phase_name=self.phase_name)
        if "ds" not in colf.titles:
            colf.updateGeometry()  # for ds

        uc = unitcell.unitcell_from_parameters(colf.parameters)
        uc.makerings(colf.ds.max(), self.ds_tol)
        self.uc = uc

        if self.foridx is None:
            self.foridx = range(len(uc.ringds))

        sel, npks, hmax = mask_rings_by_ifrac(
            colf, self.ds_tol, colf.ds.max(), self.ifrac, uc,
            forref=self.foridx, verbose=True, return_npks_hmax=True)
        isel = sel & (np.abs(np.sin(np.radians(colf.eta))) > self.etacut)

        dtyi = geometry.dty_to_dtyi(colf.dty, self.ystep, self.ymin)
        colf.addcolumn(dtyi, "dtyi")

        # cache these to speed up selections later
        colf.addcolumn(np.sin(np.radians(colf.omega)), "sinomega")
        colf.addcolumn(np.cos(np.radians(colf.omega)), "cosomega")

        if keep_colf:
            self.colf = colf
        else:
            # iplot only needs these two, and only to draw. Keep a decimated
            # copy so the plot survives; the full colf is often tens of GB.
            step = max(1, colf.nrows // 2000000)
            self._plot_ds = np.array(colf.ds[::step])
            self._plot_I = np.array(colf.sum_intensity[::step])
            self.colf = None
        self.icolf = colf.copyrows(isel)

        omega = np.ascontiguousarray(self.icolf.omega, dtype=np.float64)
        dty = np.ascontiguousarray(self.icolf.dty, dtype=np.float64)

        # Choose omega bins for the partitioning of the icolf
        kw, _ = choose_omega_bins(self, omega, verbose=True)

        # Now build the partition, for fast selection of peaks at each point in the scan.
        masker = VoxelSinoMasker(omega, dty, self.ystep, ymin=self.ymin)
        masker.partition(keep_columns=False, **kw)
        self._masker = masker          # _build_index reuses it

        # Reorder the icolf to match the partitioning, so that the partition indices are valid for the icolf
        self.icolf.reorder(masker.peak_ordering)

        # now we can compute npks too
        if self.fpks < 1:
            self.minpks = int(npks * self.fpks)
        else:
            self.minpks = self.fpks
        print(
            "Using for indexing:",
            self.icolf.nrows,
            "npks, minpks, forgen, foridx",
            npks,
            self.minpks,
            self.forgen,
            self.foridx
        )
        self.hmax = hmax

        # now save the peaks to disk
        # Delete existing columnfile if present
        if icolf_filename is None:
            icolf_filename = self.dset.icolfile
        if os.path.exists(icolf_filename):
            os.remove(icolf_filename)
        ImageD11.columnfile.colfile_to_hdf(self.icolf, icolf_filename, compression=None)

        self.icolf_filename = icolf_filename

        # Delete existing partition if present
        index_filename = default_index_filename(self)
        if os.path.exists(index_filename):
            os.remove(index_filename)

        # Prepare the partition index for memory mapping in the workers
        # It is actually mem-mapped in the initializer
        self._build_index()


    def iplot(self, skip=1):
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
        ax[0].vlines(self.uc.ringds, 1e4, 3e4, color='red', label=self.phase_name)
        return f, ax

    def write_config(self, config_path, cpus_per_chunk):
        with open(config_path, 'w') as f:
            f.write("dset_path={}\n".format(self.dset.dsfile))
            f.write("parfile={}\n".format(self.parfile))
            f.write("phase_name={}\n".format(self.phase_name))
            f.write("symmetry={}\n".format(self.symmetry))
            f.write("icolf_filename={}\n".format(self.icolf_filename))
            f.write("index_filename={}\n".format(self.index_filename))
            f.write("y0={}\n".format(self.y0))
            f.write("hkl_tol={}\n".format(self.hkl_tol))
            f.write("ds_tol={}\n".format(self.ds_tol))
            f.write("cosine_tol={}\n".format(self.cosine_tol))
            f.write("minpks={}\n".format(self.minpks))
            f.write("hmax={}\n".format(self.hmax))
            f.write("forgen={}\n".format(','.join(map(str, self.forgen))))
            f.write("uniqcut={}\n".format(self.uniqcut))
            f.write("nprocs={}\n".format(cpus_per_chunk))

    def submit_slurm_chunks(self, grains_prefix, id11_code_path, gridstep=1, n_chunks=4, cpus_per_chunk=64, time_h=48, partition="nice-long", mem_G=32, debugpoints=None):
        ds = self.dset
        slurm_pbp_path = os.path.join(ds.analysispath, "slurm_pbp")

        if not os.path.exists(slurm_pbp_path):
            os.mkdir(slurm_pbp_path)

        bash_script_path = os.path.join(slurm_pbp_path, ds.dsname + '_pbp_recon_slurm_' + self.phase_name + '.sh')
        python_script_path = os.path.join(id11_code_path, "ImageD11/nbGui/S3DXRD/run_pbp_recon_chunk.py")
        outfile_path = os.path.join(slurm_pbp_path, ds.dsname + '_pbp_recon_slurm_' + self.phase_name + '_%A_%a.out')
        errfile_path = os.path.join(slurm_pbp_path, ds.dsname + '_pbp_recon_slurm_' + self.phase_name + '_%A_%a.err')
        config_path = os.path.join(slurm_pbp_path, ds.dsname + '_pbp_recon_slurm_config_' + self.phase_name + ".txt")
        
        self.write_config(config_path, cpus_per_chunk=cpus_per_chunk)
        
        # Automatically generate all points
        if debugpoints is None:
            all_points = geometry.step_grid_from_ybincens(self.ybincens, self.ystep, gridstep, self.y0)
        else:
            all_points = np.array(debugpoints)
        
        # Split into chunks, one per array task
        chunks = np.array_split(all_points, n_chunks)
        
        # Save each chunk to text file using np.savetxt
        for idx, chunk in enumerate(chunks):
            np.random.shuffle(chunk)  # shuffle points within each chunk
            chunk_prefix = os.path.join(slurm_pbp_path, ds.dsname + "_chunk_" + self.phase_name + "_")
            chunk_suffix = ".txt"
            chunk_file = "{}{}{}".format(chunk_prefix, idx, chunk_suffix)
            np.savetxt(chunk_file, chunk, fmt='%d')
        
        # Create single sbatch script using job array
        sbatch_content = """#!/bin/bash

#SBATCH --job-name=pbp_scanning
#SBATCH --output={outfile_path}
#SBATCH --error={errfile_path}
#SBATCH --time={time_h}:00:00
#SBATCH --partition={partition}
#SBATCH --cpus-per-task={cpus_per_chunk}
#SBATCH --array=0-{array_end}
#SBATCH --mem={mem_G}G
export NUMBA_CACHE_DIR=/tmp/numba_${{SLURM_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}
mkdir -p $NUMBA_CACHE_DIR
CHUNK_FILE={chunk_prefix}${{SLURM_ARRAY_TASK_ID}}{chunk_suffix}
OMP_NUM_THREADS=1 PYTHONPATH={id11_code_path} python {python_script_path} \
{config_path} $CHUNK_FILE {grains_prefix}${{SLURM_ARRAY_TASK_ID}}.txt
""".format(
            outfile_path=outfile_path,
            errfile_path=errfile_path,
            time_h=time_h,
            partition=partition,
            cpus_per_chunk=cpus_per_chunk,
            array_end=n_chunks - 1,
            mem_G=mem_G,
            chunk_prefix=chunk_prefix,
            chunk_suffix=chunk_suffix,
            id11_code_path=id11_code_path,
            python_script_path=python_script_path,
            config_path=config_path,
            grains_prefix=grains_prefix
        )
        
        with open(bash_script_path, 'w') as f:
            f.write(sbatch_content)

        # output file paths
        grains_files = ['{}{}.txt'.format(grains_prefix, chunk) for chunk in range(n_chunks)]
        
        return bash_script_path, grains_files
    
    def point_by_point(
            self,
            grains_filename,
            icolf_filename=None,  # use self
            index_filename=None,  # use self
            nprocs=None,
            gridstep=1,
            debugpoints=None,
            loglevel=3,
    ):
        """
        grains_filename = output file
        icolf_filename = hdf5file to write. Allows mmap to be used.
        index_filename = hdf5file for voxel_mask index to mmap
        """
        if icolf_filename is not None:
            self.icolf_filename = icolf_filename

        if index_filename is not None:
            self.index_filename = index_filename
        
        self.loglevel = loglevel
        
        start = time.time()
        # FIXME - look at dataset ybincens and roi_iradon output_size ?
        if nprocs is None:
            nprocs = cImageD11.cores_available()

        print("nprocs=%s cores_available=%s SLURM_CPUS_PER_TASK=%s cache_dir=%s"
          % (nprocs, cImageD11.cores_available(),
             os.environ.get("SLURM_CPUS_PER_TASK"),
             os.environ.get("NUMBA_CACHE_DIR")), flush=True)

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
        # gmap = {}

        # Warm the numba on-disk cache in the parent. With forkserver the
        # compiled code is not inherited, but a warm cache turns a compile
        # into a load in each worker. The argument TYPES must match what the
        # worker passes: a read-only memmap is a distinct numba type from a
        # writable array, so this loads the index and columns the same way
        # initializer does rather than using dummy arrays.
        _part = load_peak_index_cache(self.index_filename, mmap=True)
        if _part is None:
            raise RuntimeError(
                "peak index %s is missing or the wrong version -- rerun "
                "setpeaks()" % self.index_filename)

        _sm = ImageD11.columnfile.mmap_h5colf(self.icolf_filename)
        _n = _part["maxlocal"]
        _t0 = time.time()
        select_point_peaks(
            0, 0,
            float(self.ystep), float(self.y0), float(self.ymin),
            _part["omega_partitions"], _part["dty_partitions"],
            _sm["sinomega"], _sm["cosomega"], _sm["dty"], _sm["dtyi"],
            _part["usin"], _part["ucos"],
            np.empty(_n, np.int64), np.empty(_n, np.float64),
            np.empty(_n, np.int64),
            False)
        print("numba warmup %.1f s, cache hits %d misses %d"
              % (time.time() - _t0,
                 sum(select_point_peaks.stats.cache_hits.values()),
                 sum(select_point_peaks.stats.cache_misses.values())),
                 flush=True)
        for _f in (select_point_peaks, fill_voxel_idx):
            print("  %s: hits %d misses %d"
                  % (_f.__name__,
                     sum(_f.stats.cache_hits.values()),
                     sum(_f.stats.cache_misses.values())), flush=True)
        del _part, _sm
        
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
                            self.index_filename,
                            self.loglevel,
                    ),
            ) as p:
                for i, j, g in p.imap_unordered(proxy, args):
                    # gmap[i, j] = g  # unused anyway
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
