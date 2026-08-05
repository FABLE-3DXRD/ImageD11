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

from ImageD11 import sym_u, unitcell, parameters
import ImageD11.sinograms.dataset
import ImageD11.indexing
from ImageD11.columnfile import columnfile
from ImageD11.sinograms import geometry
from ImageD11.sinograms.voxel_mask import VoxelSinoMasker, fill_voxel_idx

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
    """
    Peaks that put step-space point (si, sj) in the beam.
 
    relax=False reproduces geometry.dtyimask_from_step_sincos exactly:
        dty_to_dtyi(y0 - sx sin(w) - sy cos(w)) == dtyi
    i.e. only the single nearest scan row at each omega, |ddty| <= ystep/2.
 
    relax=True uses the refiner's predicate instead:
        |y0 - sx sin(w) - sy cos(w) - dty| <= ystep
    which is the two nearest scan rows, |ddty| <= ystep. Measured on a 275-row
    scan this is exactly 2.00x the peaks, all of the extra ones sitting 0.5 to
    1.0 ystep out, i.e. one scan row away.
 
    These are different questions, not a tolerance to tune:
      - bin equality asks "was this peak measured at the scan position that
        puts this voxel in the beam", which is what you want when *deciding*
        what sits at (si, sj)
      - the relaxed form asks "within one step of it", which is what you want
        when *fitting* an orientation you have already assigned, and the
        refiner pairs it with a 1/(ydist+1) weight so the off-row peaks count
        for less. The indexer has no such weighting -- every peak counts the
        same in ind.scorethem() and in hkluniq.
 
    All peak arrays must be in partition order. Returns the number written
    into `out`, or -1 if the buffers were too small (never raises -- this is
    called per point in worker processes and a raise here is a silent skip).
    """
    sx = si * ystep
    sy = -sj * ystep          # geometry.step_to_sample: note the sign on sy
 
    n = fill_voxel_idx(sx, sy, y0, ystep, ymin,
                       omega_partitions, dty_partitions, dty,
                       sinomega, cosomega, sinomega_bins, cosomega_bins,
                       idx_buf, ydist_buf)
    if n < 0:
        return -1
 
    if relax:
        # fill_voxel_idx already applied |dty_calc - dty| <= ystep
        for t in range(n):
            out[t] = idx_buf[t]
        return n
 
    # tighten to dty_to_dtyi(dty_calc) == dtyi. round() here matches numpy's
    # round-half-to-even, which is what geometry.dty_to_dtyi uses.
    m = 0
    inv_ystep = 1.0 / ystep
    for t in range(n):
        q = idx_buf[t]
        dty_calc = y0 - sx * sinomega[q] - sy * cosomega[q]
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
    """
    Indexing function called at one point in space.
    Selects peaks from the sinogram and attempts to index.
    More of this could be indexing.py I guess.

    si, sj = position in step space (origin is rotation axis)

    Effectively a columnfile:
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
    npk = select_point_peaks(
        si, sj, ystep, y0, ymin,
        omega_partitions, dty_partitions,
        sinomega, cosomega, dty, dtyi,
        sinomega_bins, cosomega_bins,
        idx_buf, ydist_buf, sel_buf, relax_mask)
    if npk <= 0:
        # npk < 0 means the buffer was short, which max_candidates should
        # have prevented; either way there is nothing to index here
        return [(0, 0, np.eye(3))]
    sel = sel_buf[:npk]
 
    # gathers are now O(npk), not O(N_peaks)
    gv, gx, gy, gz = get_local_gv(si, sj, ystep,
                                  omega[sel], sinomega[sel], cosomega[sel],
                                  xl[sel], yl[sel], zl[sel])
    eta_local = eta[sel]
    gv_local_mask = gv            # isel is all-True within the icolf


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
        i, j,
        sm["omega"], sm["sinomega"], sm["cosomega"], sm["dty"], sm["dtyi"],
        sm["xl"], sm["yl"], sm["zl"], sm["eta"],
        partglobal["omega_partitions"], partglobal["dty_partitions"],
        partglobal["sinomega_bins"], partglobal["cosomega_bins"],
        bufglobal["idx"], bufglobal["ydist"], bufglobal["sel"],
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
 
 
def initializer(parfile, phase_name, symmetry, colfile, loglevel=3,
                buffer_size=None):
    global ucglobal, symglobal, parglobal, colglobal, partglobal, bufglobal
    if threadpoolctl is not None:
        threadpoolctl.threadpool_limits(limits=1)
    parglobal = parameters.read_par_file(parfile, phase_name=phase_name)
    ucglobal = unitcell.unitcell_from_parameters(parglobal)
    symglobal = sym_u.getgroup(symmetry)()
    colglobal = ImageD11.columnfile.mmap_h5colf(colfile)
    ImageD11.indexing.loglevel = loglevel
 
    with h5py.File(colfile, "r") as hin:
        g = hin["PBPPartition"]
        partglobal = {k: g[k][:] for k in
                      ("omega_partitions", "dty_partitions",
                       "sinomega_bins", "cosomega_bins")}
        partglobal["ymin"] = float(g.attrs["ymin"])
 
    # scratch, allocated once per worker rather than once per point
    if buffer_size is None:
        buffer_size = max(1024, colglobal.nrows // 10)
    bufglobal = dict(idx=np.empty(buffer_size, np.int64),
                     ydist=np.empty(buffer_size, np.float64),
                     sel=np.empty(buffer_size, np.int64))



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
        Given a cf_2d in RAM (colf), make an icolf in RAM which contains the selection of peaks we want to index
        """
        # Load parameters if not already loaded
        colf.parameters.loadparameters(self.parfile, phase_name=self.phase_name)
        if "ds" not in colf.titles:
            colf.updateGeometry()  # for ds
        #
        uc = unitcell.unitcell_from_parameters(colf.parameters)
        uc.makerings(colf.ds.max(), self.ds_tol)
        self.uc = uc
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
        self.icolf = colf.copyrows(isel)

        masker = VoxelSinoMasker(np.ascontiguousarray(self.icolf.omega),
                                 np.ascontiguousarray(self.icolf.dty),
                                 self.ystep)
        masker.partition()                       # heuristic picks the omega bins
        order = masker.peak_ordering
     
        # put the peaks in partition order so the workers do no reordering
        self.icolf = self.icolf.copyrows(order)
        self.partition = dict(
            omega_partitions=masker.omega_partitions,
            dty_partitions=masker.dty_partitions,
            sinomega_bins=masker.sinomega_bins,
            cosomega_bins=masker.cosomega_bins,
            ymin=float(masker.ymin),
        )
        print("indexing partition: %d omega bins (%.4f deg) x %d dty bins"
              % (masker.sinomega_bins.size, masker.omega_binsize,
                 masker.dty_partitions.shape[1] - 1))

        # how many peaks did we see?
        self.npks = npks

        # now we can compute npks too
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

        # now save the peaks to disk
        if icolf_filename is None:
            icolf_filename = self.dset.icolfile
        if os.path.exists(icolf_filename):
            os.remove(icolf_filename)
        ImageD11.columnfile.colfile_to_hdf(self.icolf, icolf_filename, compression=None)

        with h5py.File(icolf_filename, "a") as hout:
            g = hout.require_group("PBPPartition")
            for k, v in self.partition.items():
                if np.isscalar(v):
                    g.attrs[k] = v
                else:
                    if k in g:
                        del g[k]
                    g.create_dataset(k, data=v)
        
        self.icolf_filename = icolf_filename

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

    def write_config(self, config_path, cpus_per_chunk):
        with open(config_path, 'w') as f:
            f.write("dset_path={}\n".format(self.dset.dsfile))
            f.write("parfile={}\n".format(self.parfile))
            f.write("phase_name={}\n".format(self.phase_name))
            f.write("symmetry={}\n".format(self.symmetry))
            f.write("icolf_filename={}\n".format(self.icolf_filename))
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
            nprocs=None,
            gridstep=1,
            debugpoints=None,
            loglevel=3,
    ):
        """
        grains_filename = output file
        icolf_filename = hdf5file to write. Allows mmap to be used.
        """
        if icolf_filename is not None:
            self.icolf_filename = icolf_filename
        
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
        # gmap = {}
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
