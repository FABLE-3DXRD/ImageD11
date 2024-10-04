#!/usr/bin/env python
# coding: utf-8

# Try to build the point-by-point mapping code ...
from __future__ import print_function, division

import os
os.environ["OMP_NUM_THREADS"] = "1"
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
import pprint
from matplotlib import pyplot as plt
from skimage.filters import threshold_otsu
from skimage.morphology import convex_hull_image

from ImageD11 import sym_u, unitcell, parameters
import ImageD11.sinograms.dataset
import ImageD11.indexing
from ImageD11.columnfile import columnfile
from ImageD11.sinograms.geometry import dty_to_dtyi, dtyimask_from_sincos, step_to_sample, sample_to_lab_sincos, \
    dtyi_to_dty
from ImageD11.sinograms.roi_iradon import run_iradon

# GOTO - find somewhere!
# Everything written in Numba could move to C?
@numba.njit(boundscheck=True)
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


@numba.njit
def nb_choose_best(i, j, u, n, NY, ubiar,
                   minpeaks=6):
    # map of the unique scores
    uniq = np.ones((NY, NY), dtype='q')
    uniq.fill(minpeaks)  # peak cutorr
    npk = np.zeros((NY, NY), dtype='q')
    ubi = np.zeros((NY, NY, 3, 3), dtype='d')
    ubi.fill(np.nan)
    for k in range(i.size):
        ip = i[k]
        jp = j[k]
        #        if ip == 96 and jp == 510:
        #            print(ip,jp,k, ubiar[:,:,k])
        if u[k] > uniq[ip, jp]:
            uniq[ip, jp] = u[k]
            npk[ip, jp] = n[k]
            for ii in range(3):
                for jj in range(3):
                    ubi[ip, jp, ii, jj] = ubiar[ii, jj, k]
    return uniq, npk, ubi


@numba.njit
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


@numba.njit
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
    def __init__(self, *kwargs):
        super(PBPMap, self).__init__(*kwargs)

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
    def UBI(self):
        return self.ubi

    def plot_nuniq_hist(self):
        fig, ax = plt.subplots()
        ax.hist(self.nuniq, bins=np.arange(0.5, np.max(self.nuniq) + 0.51, 1))
        ax.set_xlabel('Unique spots per pixel')
        ax.set_ylabel('Count')
        plt.show()

    def choose_best(self, minpeaks=6):
        self.best_nuniq, self.best_npks, self.best_ubi = nb_choose_best(
            self.i_shift, self.j_shift, self.nuniq.astype(int), self.ntotal.astype(int), self.NY, self.ubi, minpeaks)

    def plot_best(self, minpeaks=6):
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
                 hkl_tol=0.01,
                 fpks=0.7,
                 ds_tol=0.005,
                 etacut=0.1,
                 ifrac=None,
                 forref=None,
                 y0=0.0,
                 phase_name=None,

                 ):
        self.dset = dset
        self.phase_name = phase_name
        # peak selection parameters
        self.fpks = fpks
        self.forref = forref
        self.ds_tol = ds_tol

        if ifrac is None:
            self.ifrac = 1.0 / len(self.dset.ybincens)
        else:
            self.ifrac = ifrac

        self.etacut = etacut
        # refinement parameters
        self.hkl_tol = hkl_tol
        # geometry stuff
        self.ystep = self.dset.ystep
        self.y0 = y0

    def setmap(self, pbpmap, pbpmap_filename=None):
        """Set a memory-mapped input PBPMap to use"""
        if pbpmap_filename is None:
            # get from the dset
            pbpmap_filename = self.dset.refmapfile
        self.pbpmap = pbpmap
        ImageD11.columnfile.colfile_to_hdf(self.pbpmap, pbpmap_filename, compression=None)
        self.pbpmap_filename = pbpmap_filename

        # set up the grid to refine on
        # set up grids of X and Y positions in the sample reference frame
        # get an array of the integer I and J step values from the pmap
        # this is in step space
        # e.g -33, -32, -31, ... 32, 33

        # get arrays of integer translations in step space
        # e.g -33, -32, -31, ... 32, 33
        si = np.arange(self.pbpmap.i.min(), self.pbpmap.i.max() + 1)
        sj = np.arange(self.pbpmap.j.min(), self.pbpmap.j.max() + 1)

        # convert these arrays to sample space
        sx, sy = step_to_sample(si, sj, self.ystep)

        # make a grid of these
        self.sx_grid, self.sy_grid = np.meshgrid(sx, sy, indexing='ij')

    def setpeaks(self, colf, icolf_filename=None):
        """Similar to PBP.setpeaks for now"""
        if icolf_filename is None:
            icolf_filename = self.dset.refpeaksfile
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
        dtyi = dty_to_dtyi(colf.dty, self.ystep)
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
            "Using for refinement:",
            self.icolf.nrows,
            "npks, forref",
            self.npks,
            self.forref,
        )

        self.savepeaks(filename=icolf_filename, del_existing=True)

    def savepeaks(self, filename=None, del_existing=False):
        # if no name supplied
        if filename is None:
            # try to use one in self
            filename_to_use = self.icolf_filename
        else:
            # use supplied one and change the parameter in self
            filename_to_use = filename
            self.icolf_filename = filename
        if os.path.exists(filename_to_use) and del_existing:
            os.remove(filename_to_use)
        ImageD11.columnfile.colfile_to_hdf(self.icolf, filename_to_use, compression=None)

    def loadpeaks(self, filename=None):
        """Load icolf from disk"""
        if filename is None:
            # try to use one in self
            filename_to_use = self.icolf_filename
        else:
            # use supplied one and change the parameter in self
            filename_to_use = filename
            self.icolf_filename = filename
        self.icolf = ImageD11.columnfile.colfile_from_hdf(filename_to_use)

    def iplot(self, skip=1):
        """Same as PBP.iplot for now"""
        import pylab as pl

        f, ax = pl.subplots(2, 1)
        ax[0].plot(self.colf.ds[::skip], self.colf.sum_intensity[::skip], ",")
        ax[0].plot(self.icolf.ds[::skip], self.icolf.sum_intensity[::skip], ",")
        ax[0].set(yscale="log", ylabel="sum intensity", xlabel="d-star")
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
        ax[1].set(ylabel="dty", xlabel="omega")
        return f, ax

    def setmask(self, manual_threshold=None, doplot=False):
        """Set a mask for choosing what to refine or not.
        At the moment it does an iradon on the sinogram of all the 2D peaks in self.colf"""
        whole_sample_sino, xedges, yedges = np.histogram2d(self.colf.dty, self.colf.omega,
                                                           bins=[self.dset.ybinedges, self.dset.obinedges])
        shift = -self.y0 / self.ystep
        nthreads = len(os.sched_getaffinity(os.getpid()))
        whole_sample_recon = run_iradon(whole_sample_sino, self.dset.obincens, 0, shift, workers=nthreads)

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

    def get_origins(self, guess_speed=True, guess_npks=10000):
        print('Getting gvecs...')
        # compute gves
        gve = np.column_stack((self.icolf.gx, self.icolf.gy, self.icolf.gz))

        nthreads = len(os.sched_getaffinity(os.getpid())) - 1
        numba.set_num_threads(nthreads)

        if guess_speed:
            print('Running test function twice for the first 1000 peaks to guess speed...')
            print('This is because we get speed advantages with peaks on consecutive frames')

            idx = np.arange(guess_npks)
            sorter_reduced = np.lexsort((self.icolf.dty[idx], self.icolf.omega[idx]))

            print('First time to trigger NJIT')
            _ = compute_origins(self.singlemap, self.mask,
                                gve[idx], self.icolf.sinomega[idx], self.icolf.cosomega[idx], self.icolf.omega[idx], self.icolf.dty[idx],
                                self.sx_grid, self.sy_grid,
                                self.y0, self.ystep, self.hkl_tol,
                                sorter_reduced)

            print('Now we time')
            start = time.perf_counter()
            _ = compute_origins(self.singlemap, self.mask,
                                gve[idx], self.icolf.sinomega[idx], self.icolf.cosomega[idx], self.icolf.omega[idx], self.icolf.dty[idx],
                                self.sx_grid, self.sy_grid,
                                self.y0, self.ystep, self.hkl_tol,
                                sorter_reduced)
            end = time.perf_counter()
            time_taken = end - start
            time_for_one_peak = time_taken/guess_npks
            scaling_factor = 1/6  # experimental, from observation of guess vs reality
            time_for_all_peaks = time_for_one_peak * self.icolf.nrows * scaling_factor
            print('I estimate roughly', time_for_all_peaks, 'seconds for all peaks in self.icolf')
            print("That's", time_for_all_peaks/60/60, 'hours')

        print('Lexsort...')
        # get the indices that sort all the arrays by omega then dty
        sorter = np.lexsort((self.icolf.dty, self.icolf.omega))

        print('Running numba computation on', nthreads, 'threads, may take a while!')
        
        lx_modified = compute_origins(self.singlemap, self.mask,
                                      gve, self.icolf.sinomega, self.icolf.cosomega, self.icolf.omega, self.icolf.dty,
                                      self.sx_grid, self.sy_grid,
                                      self.y0, self.ystep, self.hkl_tol,
                                      sorter)

        self.icolf.addcolumn(lx_modified, 'xpos_refined')


@numba.njit
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


@numba.njit(parallel=True)
def compute_origins(singlemap, sample_mask,
                    gve, sinomega, cosomega, omega, dty,
                    sx_grid, sy_grid,
                    y0, ystep, hkl_tol,
                    sorter):
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
            # get an index that can get us the sinomega and cosomega value for thsi group
            omega_idx = omega_idx_group[0]

            so = sinomega[omega_idx]
            co = cosomega[omega_idx]

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


def idxpoint(
        i,
        j,
        isel,
        omega,
        sinomega,
        cosomega,
        dtyi,
        sc,
        fc,
        ystep=2.00,
        y0=0.0,
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

    i,j = integer step from the rotation axis center

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
    mr = isel
    # select the peaks using the geometry module
    m = dtyimask_from_sincos(i, j, sinomega, cosomega, dtyi, y0, ystep)
    # now we correct g-vectors for origin point

    # compute the sample coordinates for this pixel
    sx, sy = step_to_sample(i, j, ystep)
    dty = dtyi_to_dty(dtyi, ystep)

    # compute x distance offset for this pixel
    lx, _ = sample_to_lab_sincos(sx, sy, y0, dty, sinomega, cosomega)
    parlocal = parglobal.get_parameters().copy()
    parlocal['distance'] = parlocal['distance'] - lx

    # compute tth eta for all peaks with corrected distance at this pixel
    tth, eta = ImageD11.transform.compute_tth_eta((sc, fc), **parlocal)

    # compute gvectors for all peaks with corrected distance at this pixel
    gx, gy, gz = ImageD11.transform.compute_g_vectors(tth, eta,
                                                      omega,
                                                      parlocal['wavelength'],
                                                      wedge=parlocal['wedge'],
                                                      chi=parlocal['chi'])

    # mask g-vectors to this pixel for indexing
    inds = np.mgrid[0: len(gx)][m & mr]
    gv = np.empty((len(inds), 3), "d")
    gv[:, 0] = gx[inds]
    gv[:, 1] = gy[inds]
    gv[:, 2] = gz[inds]

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
            except TypeError:
                continue
    if ImageD11.indexing.loglevel <= 3:
        sys.stdout.flush()
    global symglobal
    ind.ubis = [sym_u.find_uniq_u(ubi, symglobal) for ubi in ind.ubis]
    # count the unique peaks per grain
    uniqs = np.empty((len(ind.ubis), 2), int)
    for i, ubi in enumerate(ind.ubis):
        # pass in the dty mask for peak counting with all peaks
        uniqs[i] = hkluniq(ubi, gx, gy, gz, eta, m, hkl_tol, hmax)
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
        sm["sc"],
        sm["fc"],
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

        # colf.addcolumn(np.round((colf.dty - self.y0) / self.ystep).astype(int), "dtyi")
        # dtyi = dty_to_dtyi(colf.dty - self.y0, self.ystep)
        dtyi = dty_to_dtyi(colf.dty, self.ystep)
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

    def iplot(self, skip=1):
        import pylab as pl

        f, ax = pl.subplots(2, 1)
        ax[0].plot(self.colf.ds[::skip], self.colf.sum_intensity[::skip], ",")
        ax[0].plot(self.icolf.ds[::skip], self.icolf.sum_intensity[::skip], ",")
        ax[0].set(yscale="log", ylabel="sum intensity", xlabel="d-star")
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
        ax[1].set(ylabel="dty", xlabel="omega")
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
            nprocs = len(os.sched_getaffinity(os.getpid()))

        idxopt = {
            "ystep": self.ystep,
            "y0": self.y0,
            "minpks": self.minpks,
            "hkl_tol": self.hkl_tol,
            "ds_tol": self.ds_tol,
            "cosine_tol": self.cosine_tol,
            "forgen": self.forgen,
            "hmax": self.hmax,
            "uniqcut": self.uniqcut,
        }
        pprint.pprint(idxopt)

        # nlo = np.floor((self.ybincens.min() - self.y0) / self.ystep).astype(int)
        # nhi = np.ceil((self.ybincens.max() - self.y0) / self.ystep).astype(int)

        nlo = np.floor((self.ybincens.min()) / self.ystep).astype(int)
        nhi = np.ceil((self.ybincens.max()) / self.ystep).astype(int)

        rng = range(nlo, nhi + 1, gridstep)

        if debugpoints is None:
            args = [(i, j, idxopt) for i in rng for j in rng]
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
