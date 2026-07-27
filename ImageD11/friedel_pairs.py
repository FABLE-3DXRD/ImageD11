# friedel_pairs.py
# Match Friedel reflection pairs in a columnfile
#
#    Friedel relationships
#    ----------------------
#    omega pairs : (h,k,l), -(h,k,l) reflexion 180° apart in ω:              g -> -g, eta -> 180 - eta
#    eta pairs   : "entry" and "exit" (h,k,l),-(h,k,l) reflexions at ω ± θ:  g -> -g, eta -> 180 + eta
#
# Works both for box-beam and scanning-3DXRD
# Two modes of pairing: 
#       global    -> split columnfile in half and match the two mirror subsets
#       by chunks -> split columnfile in chunks along dty (scans) or omega+dty (frames) and match mirror chunks
#                    (for scanning-3DXRD data only)
#
# Friedel pair matching handled by FriedelPairIndexer obj. 
# For the chunk method, columnfile splitting handled by PeakSubsets obj. 

# written by Jean-Baptiste Jacob, jean-baptiste.jacob@esrf.fr
# Last updated: 26 May 2026


import os, sys, time
import logging
import numpy as np
import multiprocessing as mp
from collections import namedtuple
from contextlib import contextmanager
from tqdm import tqdm

import scipy.spatial
from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt
from ImageD11 import columnfile, transform

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"



#  ──────────────────────────────────────────────────────────────────────────────────────────
# PeakSubset: peak sorting and subset selection for pairing
# ──────────────────────────────────────────────────────────────────────────────────────────
Pair_dty       = namedtuple('Pair', ['yi_hi', 'yi_lo', 'dty_hi', 'dty_lo'])
Pair_omega_dty = namedtuple('Pair', ['yi_A', 'yi_B', 'omi_A', 'omi_B', 'dty_A', 'dty_B', 'omega_A', 'omega_B'])
Pair_eta       = namedtuple('Pair', ['ei_hi', 'ei_lo', 'eta_hi', 'eta_lo'])

def _make_dty_pair(yi_hi, yi_lo, yc):
    return Pair_dty(
        yi_hi=yi_hi,
        yi_lo=yi_lo,
        dty_hi=float(yc[yi_hi]),
        dty_lo=float(yc[yi_lo]))

def _make_omega_dty_pair(oi_A, yi_A, oi_B, yi_B, yc, oc):
    return Pair_omega_dty(
        yi_A=yi_A,
        yi_B=yi_B,
        omi_A=oi_A,
        omi_B=oi_B,
        dty_A=float(yc[yi_A]),
        dty_B=float(yc[yi_B]),
        omega_A=float(oc[oi_A]),
        omega_B=float(oc[oi_B]))

def _make_eta_pair(ei_hi, ei_lo, ec):
    return Pair_eta(
                ei_hi=ei_hi,
                ei_lo=ei_lo,
                eta_hi=float(ec[ei_hi]),
                eta_lo=float(ec[ei_lo]))

    
class PeakSubsets:
    """
    Make list of paired subsets (subset_1, subset_2) from whole dataset for matching Friedel pairs. 
    Subsets are either:
        - pairs of mirror dty scans
        - pairs of mirror omega-dty frames (omega-pairs)
        - pairs of mirror eta bins (eta, eta+180°) (eta-pairs)
    Columnfile is sorted, either in sinogram order (dty, omega) for omega-pairs or in eta order for eta-pairs 
    and subset indices are stored in lookup tables (LUT) for fast selection. 
    """
    def __init__(self, cf, ds, n_eta_bins=360, y0=None):
        self.columnfile = cf
        self.dataset = ds
        self.is_half_acquisition = cf.omega.max() - cf.omega.min() < 181
        self.scans_LUT = None
        self.frames_LUT = None
        self.eta_bins_LUT = None
        self.scans_subsets = None
        self.frames_subsets = None
        self.eta_bins_subsets = None
        self.valid_scans_subsets = None
        self.valid_frames_subsets = None
        self.valid_eta_bins_subsets = None
        
        self._get_y0(y0)
        self._set_aliases()
        self._compute_subsets(n_eta_bins)
        
        if self.is_half_acquisition:
            logger.warning('half-rotations acquisition: cannot pair dty scans and frames!')

    def _set_aliases(self):
        self.obincens = self.dataset.obincens
        self.obinedges = self.dataset.obinedges
        self.ybincens = self.dataset.ybincens
        self.ybinedges = self.dataset.ybinedges

    def _compute_subsets(self, n_eta_bins):
        if not self.is_half_acquisition:
            self.get_scan_subsets()
            self.get_frames_subsets()
        self.get_eta_bins_subsets(n_eta_bins)

    def _reset_LUTs(self):
        for attr in list(self.__dict__.keys()):
            if attr.endswith('_LUT'):
                setattr(self, attr, None)

    def _get_y0(self, y0):
        if y0 is None:
            y0 = 0.5 * (self.dataset.ymax + self.dataset.ymin)
        self.dataset.y0 = y0
        self.dataset.save()
        self.dataset.correct_bins_for_half_scan(self.dataset.y0)
                
    # --- Internal methods ---

    def get_scan_subsets(self):
        """
        Find mirror pairs of dty scans on each side of the central y0:
            scan hi :  dty = y0 + Δy
            scan lo :  dty = y0 - Δy

        Sets ds.scans_subsets: list of Pair namedtuples ('Pair', ['yi_hi', 'yi_lo', 'dty_hi', 'dty_lo'])
        where yi_* are indices into ds.ybincens.
        """
        if self.is_half_acquisition:
            raise ValueError(
                "Half-rotation acquisition: cannot pair frames in omega.")

        yc       = self.ybincens.copy()
        n_y      = len(yc)
        y0 = self.dataset.y0
        
        pairs = []
        if n_y % 2 == 1:
            # Odd number of dty steps → central bin exists
            yi0 = n_y // 2
            for yi in range(0, yi0 + 1):
                pairs.append(_make_dty_pair(yi0 + yi, yi0 - yi, yc))
        else:
            # Even number of dty steps → no central bin
            half = n_y // 2
            for yi_hi in range(half, n_y):
                yi_lo = n_y - 1 - yi_hi
                pairs.append(_make_dty_pair(yi_hi, yi_lo, yc))

        self.scans_subsets = pairs
        logger.info("[get_scans_subsets] %d pairs of dty scans (y0=%s)", len(pairs), y0)

    def get_frames_subsets(self):
        """
        Find mirror pairs of frames (A, B) in omega and dty:
            frame A :  dty = y0 + Δy ,  omega = ω
            frame B :  dty = y0 - Δy ,  omega = (ω + 180) mod 360

        Sets ds.frames_subsets: list of Pair namedtuples
        ('Pair', ['yi_A','yi_B','omi_A','omi_B', 'dty_A','dty_B','omega_A','omega_B'])
        with indices yi*, omi* into ds.ybincens and ds.obincens
        """
        yc   = self.ybincens.copy()
        y0   = self.dataset.y0
        oc   = self.obincens.copy()
        n_y  = len(yc)
        n_o  = len(oc)
        k    = n_o // 2          # index shift for a 180-deg omega step

        # Sanity checks
        if self.is_half_acquisition:
            raise ValueError(
                "Half-rotation acquisition: cannot pair frames in omega.")
        if n_o % 2 != 0:
            raise ValueError(
                "omega_bin_centers must have an even number of bins, "
                "got {n}. Switch to dty-scan-wise pairing.".format(n=n_o))
        if not np.allclose((oc[k:] - oc[:k]) % 360, 180.0, atol=self.dataset.ostep / 2):
            raise ValueError(
                "ds.obincens[i + k] is not ds.obincens[i] + 180 deg for all i. "
                "Check that omega covers exactly 360 deg with uniform spacing.")

        pairs = []

        if n_y % 2 == 1:
            # Odd number of dty steps → central bin exists
            yi0 = n_y // 2
            # Off-centre bins: pair mirror scans (yi0+yi, yi0-yi)
            for yi in range(1, yi0 + 1):
                yi_hi = yi0 + yi
                yi_lo = yi0 - yi
                for oi in range(n_o):
                    oi_lo = (oi + k) % n_o
                    pairs.append(_make_omega_dty_pair(oi, yi_hi, oi_lo, yi_lo, yc, oc))

            # Central bin: pair lower-omega half with upper half (same yi)
            for oi in range(k):
                pairs.append(_make_omega_dty_pair(oi, yi0, oi + k, yi0, yc, oc))

        else:
            # Even number of dty steps → no central bin
            half = n_y // 2
            for yi_hi in range(half, n_y):
                yi_lo = n_y - 1 - yi_hi
                for oi in range(n_o):
                    oi_lo = (oi + k) % n_o
                    pairs.append(_make_omega_dty_pair(oi, yi_hi, oi_lo, yi_lo, yc, oc))

        self.frames_subsets = pairs

        n_y0  = sum(1 for p in pairs if p.yi_A == p.yi_B)
        n_off = len(pairs) - n_y0
        logger.info(
            "[get_frames_subsets] %d pairs of 2D frames (y0=%s): "
            "%d at y0 (omega-only), %d off-y0.",
            len(pairs), y0, n_y0, n_off)

    def _get_eta_bins(self, n_bins):
        """ compute eta bins """
        emin, emax = -180, 180
        eta_step = (emax - emin) / n_bins
        self.ebinedges = np.linspace(emin, emax, n_bins+1)
        self.ebincens  = np.linspace(emin + eta_step / 2, emax - eta_step / 2, len(self.ebinedges) - 1)
    
    def get_eta_bins_subsets(self, n_bins):
        """
        Find mirror pairs of eta bins:
            bin lo :  eta_lo
            bin hi :  eta_lo + 180

        Sets ds.eta_bins_subsets: list of Pair namedtuples ('Pair', ['ei_hi', 'ei_lo', 'eta_hi', 'eta_lo'])
        where ei_* are indices into self.ebincens.
        """
        self._get_eta_bins(n_bins)
        ec = self.ebincens
        eb = self.ebinedges
        n_eta = len(ec)
        eta_half = n_eta // 2          # index shift for a 180-deg omega step
        pairs = []
        for ei_lo in range(0, eta_half):
            ei_hi = (ei_lo + eta_half ) % n_eta
            pairs.append(_make_eta_pair(ei_hi, ei_lo, ec))
        self.eta_bins_subsets = pairs

        logger.info(
            "[get_eta_bins_subsets] %d pairs of eta bins ",
            len(pairs))
    

    def sort_by_sinogram(self):
        """
        Sort columnfile in sinogram space (dty primary, omega secondary)
        and build frame-level and scan-level lookup tables for quick peak selection
        by omega and dty bins.

        Sets
        self.scans_LUT : dict { yi, slice(lo, hi), n_peaks, total_intensity) }
        self.frames_LUT : dict { (omi, yi), slice(lo, hi), n_peaks, total_intensity) }
        
        y_i / om_y_i indices match those in self.scans_subsets and self.frames_subsets
        Frames / scans with zero peaks are absent from the dict.
        """
        oc  = self.obincens
        yc  = self.ybincens
        obe = self.obinedges
        ybe = self.ybinedges
        self._reset_LUTs()
        BinEntry = namedtuple('BinEntry', ['slice', 'npeaks', 'sumI'])

        # 1. Bin assignment
        i_omega = np.clip(
            np.searchsorted(obe, self.columnfile.omega, side="right") - 1,
            0, len(oc) - 1).astype(np.int32)          # int32 saves memory vs int64
        j_dty = np.clip(
            np.searchsorted(ybe, self.columnfile.dty,  side="right") - 1,
            0, len(yc) - 1).astype(np.int32)

        # 2. Sort (skip if already done)
        if self.columnfile.sortedby == "dty_omega":
            logger.info("[sort_by_sinogram] already sorted, %d peaks", self.columnfile.nrows)
            i_omega_s, j_dty_s = i_omega, j_dty
        else:
            logger.info("[sort_by_sinogram] sorting %d peaks…", self.columnfile.nrows)

            # Encode two int32 indices into one int64 key → single argsort pass,
            n_omega   = len(oc)
            keys      = j_dty.astype(np.int64) * n_omega + i_omega  
            order     = np.argsort(keys, kind="stable")               

            self.columnfile.reorder(order)
            self.columnfile.sortedby = "dty_omega"
            i_omega_s = i_omega[order]
            j_dty_s   = j_dty[order]
    
        # 3. Build frames_LUT
        first_idx_f, counts_f = _rle2(i_omega_s, j_dty_s)
        # Vectorised group sums — single reduceat call, no Python loop over groups
        group_sums_f = np.add.reduceat(self.columnfile.sum_intensity, first_idx_f)

        self.frames_LUT = {
            (int(i_omega_s[lo]), int(j_dty_s[lo])): BinEntry(
                slice=slice(int(lo), int(lo + cnt)),
                npeaks=int(cnt),
                sumI=float(gs))
            for lo, cnt, gs in zip(first_idx_f, counts_f, group_sums_f)}

        # 4. Build scans_LUT
        first_idx_s, counts_s = _rle1(j_dty_s)
        group_sums_s = np.add.reduceat(self.columnfile.sum_intensity, first_idx_s)

        self.scans_LUT = {
            int(j_dty_s[lo]): BinEntry(
                slice=slice(int(lo), int(lo + cnt)),
                npeaks=int(cnt),
                sumI=float(gs))
            for lo, cnt, gs in zip(first_idx_s, counts_s, group_sums_s) }

        # 5. Statistics
        n_frames_total = len(oc) * len(yc)
        n_scans_total  = len(yc)
        logger.info("%d non-empty frames (%d empty).",
                    len(self.frames_LUT), n_frames_total - len(self.frames_LUT))
        logger.info("%d non-empty scans  (%d empty).",
                    len(self.scans_LUT),  n_scans_total  - len(self.scans_LUT))


    def sort_by_eta(self):
        """
        Sort columnfile in eta bins order and build lookup table for quick peak selection

        Sets
        self.eta_bins_LUT : dict { ei, slice(lo, hi), n_peaks, total_intensity) }
        ei indices match those in self.eta_bins_subsets
        """
        # initialize LUTs
        ec = self.ebincens
        eb = self.ebinedges
        self._reset_LUTs()
        BinEntry = namedtuple('BinEntry', ['slice', 'npeaks', 'sumI'])  

        # skip sorting if already done
        if self.columnfile.sortedby == "eta":
            logger.info("[sort_by_eta] peakfile already sorted – %d peaks",
                        self.columnfile.nrows)
            # 1. bin by eta
            i_eta = np.clip(
                    np.searchsorted(eb, self.columnfile.eta, side="right") - 1, 0, len(ec) - 1)
            i_eta_s = i_eta

        # 1. bin by eta and sort
        else:
            logger.info("[sort_by_eta] sorting columnfile in eta order  –  %d peaks to sort, may take some time... ",
                        self.columnfile.nrows)
        
            i_eta = np.clip(
                    np.searchsorted(eb, self.columnfile.eta, side="right") - 1, 0, len(ec) - 1)
            order = np.argsort(i_eta)
            self.columnfile.reorder(order)
            self.columnfile.sortedby = "eta"
            i_eta_s = i_eta[order]  # sorted array

        # 2. Build eta_bins_LUT: Mapping (ei) -> BinEntry
        first_idx, counts = _rle1(i_eta_s)
        group_sums = np.add.reduceat(self.columnfile.sum_intensity, first_idx)

        self.eta_bins_LUT = {
            int(i_eta_s[lo]): BinEntry(
                slice=slice(int(lo), int(lo + cnt)),
                npeaks=int(cnt),
                sumI=float(gs))
            for lo, cnt, gs in zip(first_idx, counts, group_sums)}
        
        # 3. Statistics
        n_bins_total = len(ec)
        n_bins_empty = n_bins_total - len(self.eta_bins_LUT)
        logger.info(
            "%d non-empty bins indexed (%d bins have no peaks).",
             len(self.eta_bins_LUT), n_bins_empty)
        
    
    def find_valid_subsets(self):
        """
        compute boolean arrays identifying paired subsets where both members contain
        at least one peak in corresponding lookup table (frames_LUT or scans_LUT).
        Sets new attributes self.valid_scans_subsets / self.valid_frames_subsets
        """

        # which types of subsets: eta_bins / scans + frames. depends on sorting order
        cf = self.columnfile
        subset_type_list = ['eta_bins'] if cf.sortedby == 'eta' else ['scans', 'frames']
        
        for subset_type in subset_type_list:
            # Check if LUT is defined in PeakSubset
            LUT_name = str(subset_type + '_LUT')
            LUT = getattr(self, LUT_name, None)
            if LUT is None:
                raise RuntimeError(
                    "Lookup table {LUT_name} not found in PeakSubset. "
                    "Run sort_by_sinogram() or sort_by_eta() first.".format(LUT_name=LUT_name))
             # Get corresponding list of subsets
            pairs = getattr(self, subset_type+'_subsets', None)
            if pairs is None:
                raise RuntimeError(
                    "paired subsets list not found")
            
            valid_keys = set(LUT.keys())
            valid = np.zeros(len(pairs), dtype=bool)

            if subset_type == 'frames' :
                for i, p in enumerate(pairs):
                    if (p.omi_A, p.yi_A) in valid_keys and (p.omi_B, p.yi_B) in valid_keys:
                        valid[i] = True
            elif  subset_type == 'scans':
                for i, p in enumerate(pairs):
                    if p.yi_hi in valid_keys and p.yi_lo in valid_keys:
                        valid[i] = True
            else:
                for i, p in enumerate(pairs):
                    if p.ei_hi in valid_keys and p.ei_lo in valid_keys:
                        valid[i] = True

            nv = valid.sum()
            logger.info(
                "[find_valid_subsets]  %d pairs checked in %s\n"
                "  found %d subset pairs where both members have peaks "
                "(%.2f%% valid).",
                len(valid), LUT_name, nv, nv / len(valid) * 100)
            setattr(self, 'valid_'+subset_type+'_subsets', valid)

            
    def select_pair(self, pair_id, subset_type = 'scans', return_as = 'idx'):
        """
        select a pair of subsets from one of the pairs list"""
        
        LUT        = getattr(self, subset_type+'_LUT', None)
        pairs_list = getattr(self, subset_type+'_subsets', None)
        cf = self.columnfile

        if LUT is None:
            raise RuntimeError(
            "Lookup table not found. run sort_by_sinogram() / sort_by_eta_bins(), or check input name.")

        # select pair in pairs list
        pair = pairs_list[pair_id]
        # select peaks indices in lookup table
        if subset_type == 'frames':
            idx1 = (pair.omi_A, pair.yi_A)   
            idx2 = (pair.omi_B, pair.yi_B)   
        elif subset_type == 'scans':
            idx1 = pair.yi_hi   
            idx2 = pair.yi_lo
        else:
            idx1 = pair.ei_hi   
            idx2 = pair.ei_lo

        slc_1 = LUT[idx1].slice
        slc_2 = LUT[idx2].slice

        # check if any slice is missing data
        missing = []
        if slc_1 is None: missing.append("idx1:{}".format(idx1))
        if slc_2 is None: missing.append("idx2:{}".format(idx2))

        if len(missing) > 0:
            return None, None

        # retour output
        if return_as == 'idx':
            idx1 = np.asarray(slc_1) if not isinstance(slc_1, slice) else np.arange(*slc_1.indices(cf.nrows))
            idx2 = np.asarray(slc_2) if not isinstance(slc_2, slice) else np.arange(*slc_2.indices(cf.nrows))
            return idx1, idx2
        elif return_as == 'mask':
            m_1 = np.full(cf.nrows, False, dtype=bool)
            m_2 = m_1.copy()
            m_1[slc_1] = True
            m_2[slc_2] = True
            return m_1, m_2 
        else:
            def _make_cf(slc):
                return columnfile.colfile_from_dict(
                {t: cf.getcolumn(t)[slc] for t in cf.titles})
            return _make_cf(slc_1), _make_cf(slc_2)

    def _extended_pair_selec(self, pair_id, subset_type='eta_bins'):
        """ window selection of peaks around central bin k: [k-1:k+1]."""
        # central bin
        pid = pair_id
        pairs_list = getattr(self, subset_type+'_subsets', None)
                             
        b1, b2 = self.select_pair(pid, subset_type, return_as = 'idx')
        # neighbour bins. specific cases at boundaries
        if pid == 0:
            b1_sup, b2_sup = self.select_pair(pid+1, subset_type, return_as = 'idx')
            b1 = np.concatenate((b1, b1_sup))
            b2 = np.concatenate((b2, b2_sup))
        elif pid == len(pairs_list) -1:
            b1_inf, b2_inf = self.select_pair(pid-1, subset_type, return_as = 'idx')
            b1 = np.concatenate((b1_inf, b1))
            b2 = np.concatenate((b2_inf, b2))
        else:
            b1_inf, b2_inf = self.select_pair(pid-1, subset_type, return_as = 'idx')
            b1_sup, b2_sup = self.select_pair(pid+1, subset_type, return_as = 'idx')
            b1 = np.concatenate((b1_inf, b1, b1_sup))
            b2 = np.concatenate((b2_inf, b2, b2_sup))
        # return sorted indices
        return np.sort(b1), np.sort(b2)

    def check_symmetry(self):
        """
        Plot N_peaks and total_intensity vs. position for mirror dty scans
        (y0+Δy, y0-Δy) or mirror eta bins (eta, eta+180°) and check their correlation.

        If alignment is correct, values should match between mirror pairs.
        Significant mismatch suggests an incorrect y0 (for dty), sample movement,
        or a beam issue during scanning.
        """
        # ── auto-detect mode ─────────────────────────────────────────────────
        mode = None
        for attr in ['scans_LUT', 'eta_bins_LUT']:
            if getattr(self, attr, None) is not None:
                mode = attr.replace('_LUT', '')
            else:
                continue
        if mode is None:
            raise RuntimeError(
                    'No LUT found. Run sort_by_sinogram() or sort_by_eta() first.')
        
        self.find_valid_subsets()
        
        # ── mode-specific setup ───────────────────────────────────────────────
        if mode == 'scans':
            valid    = self.valid_scans_subsets
            id_hi    = np.array([p.yi_hi for p in self.scans_subsets])
            id_lo    = np.array([p.yi_lo for p in self.scans_subsets])
            lut      = self.scans_LUT
            y0       = self.dataset.y0
            x_hi     = np.abs(self.ybincens[id_hi] - y0)
            x_lo     = np.abs(self.ybincens[id_lo] - y0)
            xlabel   = '|dty - y0|'
            lo_label = 'lo-side (dty - y0 < 0)'
            hi_label = 'hi-side (dty - y0 > 0)'
            title    = 'dty symmetry - ' + str(self.dataset.dsname)
            #logger.info("lo:{}, hi:{}",x_lo, x_hi)

        else:
            valid    = self.valid_eta_bins_subsets
            id_hi    = np.array([p.ei_hi for p in self.eta_bins_subsets])
            id_lo    = np.array([p.ei_lo for p in self.eta_bins_subsets])
            lut      = self.eta_bins_LUT
            x_hi     = self.ebincens[id_hi]
            x_lo     = self.ebincens[id_lo] + 180
            xlabel   = 'eta (deg) [+180 for lo-side]'
            lo_label = 'lo-side (eta < 0)'
            hi_label = 'hi-side (eta > 0)'
            title    = 'eta symmetry - ' + str(self.dataset.dsname)

        # ── data extraction ───────────────────────────────────────────────────
        npks_hi    = np.array([lut[i].npeaks if v else np.nan for i, v in zip(id_hi, valid)])
        npks_lo    = np.array([lut[i].npeaks if v else np.nan for i, v in zip(id_lo, valid)])
        sumI_hi    = np.array([lut[i].sumI   if v else np.nan for i, v in zip(id_hi, valid)])
        sumI_lo    = np.array([lut[i].sumI   if v else np.nan for i, v in zip(id_lo, valid)])

        log_sumI_hi = np.log10(sumI_hi)
        log_sumI_lo = np.log10(sumI_lo)
        ratio_npks  = npks_hi / npks_lo
        ratio_sumI  = log_sumI_hi / log_sumI_lo

        corr_npks = np.corrcoef(npks_hi[valid],    npks_lo[valid])[0, 1]
        corr_sumI = np.corrcoef(log_sumI_hi[valid], log_sumI_lo[valid])[0, 1]

        # ── plot ──────────────────────────────────────────────────────────────
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
        axes = axes.flatten()

        axes[0].plot(x_hi[1:], npks_hi[1:],    '.-', label=hi_label)
        axes[0].plot(x_lo[1:], npks_lo[1:],    '.-', label=lo_label)
        axes[0].set_ylabel('npks per scan')
        axes[0].set_title('Number of peaks | Corr = {:.3f}'.format(corr_npks))
        axes[0].legend()

        axes[1].plot(x_hi[1:], log_sumI_hi[1:], '.-')
        axes[1].plot(x_lo[1:], log_sumI_lo[1:], '.-')
        axes[1].set_title('Total intensity | Corr = {:.3f}'.format(corr_sumI))
        axes[1].set_ylabel('log10(sumI) per scan')

        axes[2].plot(x_hi[1:], ratio_npks[1:], 'k.-')
        axes[2].axhline(1, color='gray', linestyle='--')
        axes[2].set_ylabel('Ratio npks (hi/lo)')
        axes[2].set_xlabel(xlabel)

        axes[3].plot(x_hi[1:], ratio_sumI[1:], 'k.-')
        axes[3].axhline(1, color='gray', linestyle='--')
        axes[3].set_ylabel('Ratio log10(sumI) (hi/lo)')
        axes[3].set_xlabel(xlabel)

        fig.suptitle(title, fontsize=14)
        plt.tight_layout()

        return fig, axes


def _rle1(a):
    """Return (first_indices, counts) for a single sorted index array."""
    n = len(a)
    if n == 0:
        return np.empty(0, np.intp), np.empty(0, np.intp)
    boundary        = np.empty(n + 1, dtype=np.bool_)
    boundary[0]     = True
    boundary[n]     = True
    boundary[1:n]   = a[1:] != a[:-1]
    starts          = np.flatnonzero(boundary[:-1])
    counts          = np.diff(np.flatnonzero(boundary))
    return starts, counts

def _rle2(a, b):
    """Return (first_indices, counts) for sorted parallel index arrays."""
    n = len(a)
    if n == 0:
        return np.empty(0, np.intp), np.empty(0, np.intp)
    boundary        = np.empty(n + 1, dtype=np.bool_)
    boundary[0]     = True
    boundary[n]     = True
    boundary[1:n]   = (a[1:] != a[:-1]) | (b[1:] != b[:-1])
    starts          = np.flatnonzero(boundary[:-1])   # first_idx
    counts          = np.diff(np.flatnonzero(boundary))
    return starts, counts



# ──────────────────────────────────────────────────────────────────────────────────────────
#  FriedelPairIndexer
# ──────────────────────────────────────────────────────────────────────────────────────────

_shared_cf = None

class FriedelPairIndexer:
    """
    Match Friedel reflection pairs in a columnfile
    
    Friedel relationships
    ----------------------
    omega pairs : (h,k,l), -(h,k,l) reflexion 180° apart in ω:              g -> -g, eta -> 180 - eta
    eta pairs   : "entry" and "exit" (h,k,l),-(h,k,l) reflexions at ω ± θ:  g -> -g, eta -> 180 + eta

    Parameters 
    ----------
    cf          : ImageD11 columnfile 
    ds          : ImageD11.sinograms.dataset object
    tol_gv      : float — g-vector tolerance
    tol_eta     : float — eta angle tolerance (degrees)
    tol_logI    : float — log10(intensity) tolerance (np.inf to ignore)
    weights     : dict | None — User-defined weights for scaling search-space dimensions (Optional)
    n_steps     : int — number of distance increments in iterative KDTree search
    PeakSubsets : PeakSubsets object | None — paired subsets for friedel pair search by chunks
    """

    def __init__(self, cf, ds,
                 tol_gv=0.02, tol_eta=0.5, tol_logI=np.inf,
                 weights=None, n_steps=5):

        self.cf      = cf
        self.ds      = ds
        self.tol_gv  = tol_gv
        self.tol_eta = tol_eta
        self.tol_logI = tol_logI
        self.weights  = weights
        self.n_steps  = n_steps
        self.PeakSubsets = None   # filled by set_peak_subsets()
        self.outputs     = []     # filled by run_pairing() / match_*_pairs()
        self.logger = logging.getLogger(__name__ + '.FriedelPairIndexer')
        self.logger.setLevel(logging.INFO)

    def __repr__(self):
        return (
            'FriedelPairIndexer('
            'cf: {cf_rows} peaks, '
            'tol_gv={tol_gv}, tol_eta={tol_eta}, tol_logI={tol_logI}, '
            'n_steps={n_steps}, '
            'outputs={n_out} chunk(s))'.format(
                cf_rows=self.cf.nrows if self.cf is not None else None,
                tol_gv=self.tol_gv,
                tol_eta=self.tol_eta,
                tol_logI=self.tol_logI,
                n_steps=self.n_steps,
                n_out=len(self.outputs)))

    def set_peak_subsets(self, n_eta_bins=360, y0=None):
        self.logger.info('--- SET PEAKSUBSETS ---')
        self.PeakSubsets = PeakSubsets(self.cf, self.ds, n_eta_bins, y0)

    def sort_peak_subsets(self, pair_type='omega'):
        if pair_type == 'omega':
            if self.PeakSubsets.is_half_acquisition:
                raise ValueError("Half acquisition. Cannot find omega pairs.")
            self.PeakSubsets.sort_by_sinogram()
        else:
            self.PeakSubsets.sort_by_eta() 
        self.PeakSubsets.find_valid_subsets()

    def reset_outputs(self, pair_id_col=None):
        self.outputs = []
        # reset the eta/omega_pair_id column in cf
        if pair_id_col is not None:
            if pair_id_col in self.cf.titles:
                initvals =  np.full(self.cf.nrows, -1, dtype=int)
                self.cf.addcolumn(initvals, pair_id_col)               

    # ------------------------------------------------------------------
    #  High-level API — Friedel pairing functions
    # ------------------------------------------------------------------
    def match_friedel_pairs(self, pair_type = 'omega',
                            drop_unpaired=False,
                            filter_mode = 'relaxed',
                            doplot=True,
                            timeout = 120):
        """
        Find Friedel pairs using a KDTree approach. The columnfile is split
        into two symmetric parts, and the symmetric partner is flipped to match
        with the reference subset. Two KDtrees are build with these subset and
        an iterative nearerst-neighbor search is performed to find the pairs
        within a certain distance defined by the tolerance thresholds. 

        Works both for omega-pairs and eta-pairs.
        WARNING: do not use on large (>10M peaks) sanning-3DXRD peakfiles: very slow!
                 Use the chunk approach instead (self.match_omega_pairs_by_chunks)
        
        Parameters
        ----------
        pair_type     : str; 'omega' or 'eta'
        drop_unpaired : bool; If True non-paired peaks are filtered out from self.cf 
        timeout       : float; max execution time (s) before timeout for _run_pairing.
        filter_mode   : str; 'strict' or 'relaxed'
                        controls how threshold filtering is performed during the search
                        strict is slower but allows higher completeness for pair matching
                        see also: self.run_pairing
        Returns
        -------
        adds new colum '<pair_type>_pair_id' to self.cf
        Returns self.cf with the new pair_id column
        -------
        """
        self.logger.info('--- MATCH FRIEDEL PAIRS [%s pairs] ---', pair_type)
        
        # -- 1. Initialization: sort and split
        self.reset_outputs()
        self.cf.sortby(pair_type)
        
        if pair_type == 'eta':
            idx1 = np.nonzero(self.cf.eta < 0)[0]
            idx2 = np.nonzero(self.cf.eta > 0)[0]
        else:
            if self.cf.omega.max() - self.cf.omega.min() > 181:
                idx1 = np.nonzero(self.cf.omega%360 > 180)[0]
                idx2 = np.nonzero(self.cf.omega%360 < 180)[0]
            else:
                raise ValueError("Half acquisition. Cannot find omega pairs.")

        self.cf.addcolumn(np.full(self.cf.nrows, -1, dtype=int), pair_type+'_pair_id')
        
        # -- 2. Search-space calibration -----------------------------
        self.logger.info('--- SEARCH SPACE CALIBRATION [pilot pairing] ---')

        # initial dist_max from tol params
        dist_max = _physical_to_normalised_cutoff(
            self.tol_gv, self.tol_eta, self.tol_logI, self.weights)
        
        # optimal weights for rescaling
        _, rescaling_wts = self.estimate_search_scales(
            idx1, idx2, 0.2 * dist_max, pair_type=pair_type)
        
        _w = {'gx': 1., 'gy': 1., 'gz': 1., 'eta': 1., 'I': 1.}
        if self.weights is not None:
            _w.update(self.weights)
            
        weights_eff = {k: _w[k] * rescaling_wts[k] for k in _w}
        # updated dist_max
        dist_max = _physical_to_normalised_cutoff(
                    self.tol_gv, self.tol_eta, self.tol_logI, weights_eff)

        self.logger.info('  Effective weights  : %s',
            {k: '{:.3f}'.format(v) for k, v in weights_eff.items()})
        self.logger.info('  Rescaled dist_max : %.4f', dist_max)

        # -- 3. pair matching ----------------------------------------
        verbose = self.logger.getEffectiveLevel() < 30
        self.run_pairing(idx1, idx2,
                         pair_type   = pair_type,
                         filter_mode = filter_mode,
                         weights     = weights_eff,
                         t_max = timeout)
        
        self.merge_outputs(pair_type     = pair_type,
                           drop_broken   = True,
                           drop_unpaired = drop_unpaired)
       
        if doplot:
            _ = plot_pair_distances(self.cf, pair_type = pair_type, bins=50, log_scale=False)
        return self.cf
        

    def match_friedel_pairs_by_chunks(self,
                                     chunk_type ='scans', 
                                     filter_mode='relaxed',
                                     n_eta_bins = 360,
                                     drop_unpaired = False,
                                     extended_bin_search=False,
                                     n_workers=-1,
                                     timeout = 120, 
                                     reset_psub = False,
                                     doplot = True):
        """
        Find Friedel pairs in a large columnfile using a chunked KDTree approach.
        Instead of building KDTrees over the whole columnfile, pairs are searched
        independently within mirror subsets (mirror dty scans or omega-dty frames),
        then merged into the master columnfile.
        
        Works only for scanning-3DXRD datasets with a dty column. Only search for
        omega-pairs. eta-pair search not implemented yet. 

        Parameters
        ----------
        chunk_type    : str 'frames', 'scans', 'eta_bins' ( default: 'scans')
                        For omega pairs, use 'frames' or 'scans' for friedel pair search at 
                        different levels of chunking granularity.
                        For eta pairs, use 'eta_bins' 
        n_eta_bins    : int, even nb. number of bins for eta chunks.  
        drop_unpaired : bool (default False)
        reset_psub    : bool (default False). reset self.Peaksubset if True. 
        n_workers     : int number of parallel worker processes.
                        -1 = use all available cores
        filter_mode   : str; 'strict' or 'relaxed' (default strict)
                        see also : self.match_friedel_pairs, self.run_pairing
                        
        Returns
        -------
        adds new colum '<pair_type>_pair_id' to self.cf
        Returns self.cf with the new pair_id column
        """
        global _shared_cf

        pair_type = 'eta' if chunk_type == 'eta_bins' else 'omega'
        self.logger.info('--- MATCH FRIEDEL PAIRS [%s pairs] ---', pair_type)

        # ── 1. Checks + init ───────────────────────────────────────────────────
        self.logger.info('--- INITIALIZATION ---')
        self.reset_outputs(pair_type+'_pair_id')

        pairs_list, valid_pairs, LUT = self.initialize_peaksubsets(chunk_type, n_eta_bins, reset_psub)
        Psub  = self.PeakSubsets
        n_valid = valid_pairs.sum()

        dist_max = _physical_to_normalised_cutoff(
            self.tol_gv, self.tol_eta, self.tol_logI, self.weights)
        
        self.logger.info('  %d pairs of subsets, %d valid', len(pairs_list), n_valid)

        # ── 2. Search-space calibration ──────────────────────────────────────
        self.logger.info('--- SEARCH SPACE CALIBRATION [pilot pairing] ---')

       # pilot_chunk = 'eta_bins' if chunk_type == 'eta_bins' else 'scans'
        
        # pick a pair near the middle 
        s0 = 0
        while not valid_pairs[s0]:
            s0 += 1
        idx1_pilot, idx2_pilot = Psub.select_pair(s0, chunk_type, return_as='idx')

        # find weight factors for optimal rescaling
        try:
            _, rescaling_wts = self.estimate_search_scales(
                idx1_pilot, idx2_pilot, 0.5 * dist_max, pair_type=pair_type)
        except RuntimeError:
            rescaling_wts = {'gx': 1., 'gy': 1., 'gz': 1., 'eta': 1., 'I': 1.}
            
        # combine user weights with auto-rescaling
        _w = {'gx': 1., 'gy': 1., 'gz': 1., 'eta': 1., 'I': 1.}
        if self.weights is not None:
            _w.update(self.weights)
        weights_eff = {k: _w[k] * rescaling_wts[k] for k in _w}
        dist_max    = _physical_to_normalised_cutoff(
            self.tol_gv, self.tol_eta, self.tol_logI, weights_eff)

        self.logger.info('  Effective weights : %s',
            {k: '{:.3f}'.format(v) for k, v in weights_eff.items()})
        self.logger.info('  Rescaled dist_max : %.4f', dist_max)

        #─ ─ 3. Friedel pair matching (parallelize across chunks) ────────────────
        self.logger.info('--- FRIEDEL PAIR MATCHING ---')

        # Resolve number of workers
        n_cores = len(os.sched_getaffinity(os.getpid())) - 1
        if n_workers == -1:
            n_workers = n_cores
        else:
            n_workers = max(1, min(n_workers, n_cores))
        verbose = (n_workers == 1)
        
        # Build worker argument list 
        worker_args = []
        
        for pid in range(len(pairs_list)):
            if not valid_pairs[pid]:
                continue
            
            # extended bin search: 3 contiguous bins grouped in selection
            if extended_bin_search:
                # select only even indices subsets -> 1-bin overlap between each subset 
                if pid % 2 != 0:
                    continue
                # skip last bin: already in the last extended selection
                if (pid >= len(pairs_list)-1):
                    continue

            if extended_bin_search and chunk_type!='frames':
                try:
                    idx1, idx2 = self.PeakSubsets._extended_pair_selec(pid, chunk_type)
                except KeyError as e: # with extended selection one bounding pair may end up being non-valid pair. in this case, use regular selection instead
                    idx1, idx2 = self.PeakSubsets.select_pair(pid, chunk_type, return_as='idx')
            else:
                idx1, idx2 = self.PeakSubsets.select_pair(pid, chunk_type, return_as='idx')

            worker_args.append((
                pid, idx1, idx2,
                pair_type,
                filter_mode,
                weights_eff,
                self.tol_gv, self.tol_eta, self.tol_logI,
                self.n_steps, timeout, verbose,
            ))

        self.logger.info('  %d chunks  |  %d worker(s)  |  %d cores available',
            len(worker_args), n_workers, n_cores)

        _shared_cf = self.cf

        try:
            if n_workers == 1:
                # serial path — useful for debugging, avoids fork overhead
                raw_results = [_worker_run_pairing(a) for a in tqdm(worker_args)]
            else:
                ctx  = mp.get_context('fork')
                with ctx.Pool(processes=n_workers,
                             initializer=_worker_initializer) as pool:
                    raw_results = list(
                    tqdm(
                            pool.imap_unordered(
                                _worker_run_pairing, worker_args),
                            total=len(worker_args),
                        )
                    )
        finally:
            # always release the global reference, even if an exception occurs
            _shared_cf = None

        # Restore chunk order and store outputs
        raw_results.sort(key=lambda x: x[0])
        self.outputs = [out for _, out in raw_results]

        # ── 4. Merge ─────────────────────────────────────────────────────────
        self.logger.info('--- MERGING OUTPUTS ---')
        deduplicate = extended_bin_search == True
        # unique pair identifier in omega_pair_id. -1 for unpaired peaks
        self.merge_outputs(pair_type=pair_type,
                           drop_broken=True,
                           drop_unpaired=drop_unpaired,
                           deduplicate = deduplicate)

        if doplot:
            _ = plot_pair_distances(self.cf, pair_type = pair_type, bins=50, log_scale=False)
        return self.cf


    
    def match_quadruplets(self):
        """
        Find quadruplets of peaks satisfying the full Friedel relation in both
        the omega and eta pairing directions.

        A valid quadruplet is a 2x2 block in the (omega_pair_id, eta_pair_id)
        grid:

                  eta_pair A    eta_pair B
        omega_pair 0:  peak p00  ,  peak p01
        omega_pair 1:  peak p10  ,  peak p11

        Both match_omega_pairs() and match_eta_pairs() must have been called
        first so that 'omega_pair_id' and 'eta_pair_id' columns exist in self.cf.

        Writes a new column 'quadruplet_id' to self.cf:
            -1  : peak not part of any quadruplet
            k   : peak belongs to quadruplet k  (each value appears exactly 4 times)

        Returns
        -------
        quadruplets : (M, 4) int array
                      Each row is [p_om0_eta0, p_om0_eta1, p_om1_eta0, p_om1_eta1]
                      with cf-level indices, arranged as the 2x2 block read row-major.
        """
        cf = self.cf

        # ── 0. Sanity checks ──────────────────────────────────────────────────
        for col in ('omega_pair_id', 'eta_pair_id'):
            if col not in cf.titles:
                raise RuntimeError(
                    "'{}' column not found. "
                    "Run match_{}_pairs() first.".format(col, col.split('_')[0]))

        om_col  = cf.getcolumn('omega_pair_id')   # (nrows,)
        eta_col = cf.getcolumn('eta_pair_id')     # (nrows,)

        # ── 1. Select peaks that belong to both an omega- and an eta-pair ─────
        self.logger.info('--- QUADRUPLET SEARCH ---')

        both_mask  = (om_col > -1) & (eta_col > -1)
        n_both     = both_mask.sum()
        self.logger.info('  Peaks with both omega- and eta-pair id : %d/%d',
            n_both, cf.nrows)

        if n_both == 0:
            self.logger.warning('  no doubly-paired peaks found — returning empty.')
            cf.addcolumn(np.full(cf.nrows, -1, dtype=int), 'quadruplet_id')
            return np.empty((0, 4), dtype=int)

        cf_indices = np.where(both_mask)[0]           # global cf positions
        om_ids     = om_col[cf_indices]               # omega_pair_id for each
        eta_ids    = eta_col[cf_indices]              # eta_pair_id   for each

        # ── 2. Build signature: omega_pair_id -> frozenset({eta_a, eta_b}) ────
        # sort by om_id so peers are adjacent
        sort_om      = np.argsort(om_ids, kind='stable')
        om_sorted    = om_ids[sort_om]
        eta_sorted   = eta_ids[sort_om]
        cf_sorted    = cf_indices[sort_om]
        
        unique_om, counts_om = np.unique(om_sorted, return_counts=True)

        # keep only omega-pairs where exactly 2 doubly-paired peaks exist
        valid_om_mask = counts_om == 2
        if not valid_om_mask.all():
            self.logger.info('  NOTE: %d omega-pair(s) skipped '
                   '(do not have exactly 2 doubly-paired peaks)',
                   (~valid_om_mask).sum())

        unique_om = unique_om[valid_om_mask]

        # positions of the two peers in the sorted arrays
        starts = np.searchsorted(om_sorted, unique_om, side='left')   # (K,)

        eta_a   = eta_sorted[starts]        # eta_pair_id of first  peer
        eta_b   = eta_sorted[starts + 1]   # eta_pair_id of second peer
        cf_a    = cf_sorted[starts]         # cf index    of first  peer
        cf_b    = cf_sorted[starts + 1]    # cf index    of second peer

        # discard degenerate cases (both peers share the same eta-pair)
        degenerate   = eta_a == eta_b
        if degenerate.any():
            self.logger.info('  NOTE: %d omega-pair(s) skipped '
                   '(both peaks share the same eta_pair_id)',
                   degenerate.sum())

        keep         = ~degenerate
        unique_om    = unique_om[keep]
        eta_a        = eta_a[keep];   eta_b  = eta_b[keep]
        cf_a         = cf_a[keep];    cf_b   = cf_b[keep]

        # canonical signature: always store the smaller eta_id first
        swap         = eta_a > eta_b
        eta_a[swap], eta_b[swap] = eta_b[swap].copy(), eta_a[swap].copy()
        cf_a[swap],  cf_b[swap]  = cf_b[swap].copy(),  cf_a[swap].copy()

        # ── 3. Group omega-pairs by their (eta_a, eta_b) signature ───────────
        #
        #  Encode the pair of eta ids as a single 64-bit integer for fast
        #  sorting/grouping (safe as long as eta_pair_id < 2^31 ~ 2 billion).
        # ─────────────────────────────────────────────────────────────────────
        max_eta_id  = int(eta_col.max()) + 1
        signature   = eta_a.astype(np.int64) * max_eta_id + eta_b.astype(np.int64)
        sort_sig    = np.argsort(signature, kind='stable')
        sig_sorted  = signature[sort_sig]
        om_s        = unique_om[sort_sig]
        cf_a_s      = cf_a[sort_sig]
        cf_b_s      = cf_b[sort_sig]

        # find runs of identical signatures -> same (eta_a,eta_b)
        boundaries  = np.flatnonzero(np.diff(sig_sorted)) + 1
        run_starts  = np.concatenate([[0], boundaries])
        run_ends    = np.concatenate([boundaries, [len(sig_sorted)]])
        run_lengths = run_ends - run_starts

        n_groups    = (run_lengths >= 2).sum()
        self.logger.info('  Omega-pair groups sharing an eta signature : %d', n_groups)

        # ── 4. Emit quadruplets ───────────────────────────────────────────────
        quadruplet_list = []

        for rs, rl in zip(run_starts, run_lengths):
            if rl < 2:
                continue                        # singleton — no match

            # take pairs within the group; warn if group > 2 (ambiguous)
            if rl > 2:
                self.logger.debug('  NOTE: signature group of size %d — '
                       'taking first pair only.', rl)

            # take the first two omega-pairs in the group
            i0, i1 = rs, rs + 1

            # recover the 4 cf indices arranged as the 2x2 block:
            # row 0: omega-pair om_s[i0]  -> (cf_a_s[i0], cf_b_s[i0])
            # row 1: omega-pair om_s[i1]  -> (cf_a_s[i1], cf_b_s[i1])
            # columns correspond to eta_a and eta_b (same for both rows)
            quadruplet_list.append([
                cf_a_s[i0], cf_b_s[i0],
                cf_a_s[i1], cf_b_s[i1]])

        # ── 5. Write quadruplet_id column ─────────────────────────────────────
        quad_arr = np.array(quadruplet_list, dtype=int)   # (M, 4)
        n_quads  = len(quad_arr)

        quad_id_col = np.full(cf.nrows, -1, dtype=int)

        if n_quads > 0:
            flat_indices = quad_arr.ravel()                      # (4*M,)
            flat_ids     = np.repeat(np.arange(n_quads), 4)     # [0,0,0,0,1,1,1,1,...]
            quad_id_col[flat_indices] = flat_ids

        if 'quadruplet_id' in cf.titles:
            cf.getcolumn('quadruplet_id')[:] = quad_id_col
        else:
            cf.addcolumn(quad_id_col, 'quadruplet_id')

        # ── 6. Validation ─────────────────────────────────────────────────────
        self.logger.info('--- VALIDATION ---')

        assigned = quad_id_col[quad_id_col > -1]
        if len(assigned) > 0:
            counts     = np.bincount(assigned)
            bad        = (counts != 4).sum()
            if bad:
                self.logger.warning('  %d quadruplet_id(s) do not appear '
                       'exactly 4 times.', bad)
            else:
                self.logger.info('  OK: all quadruplet_ids appear exactly 4 times.')

        n_peaks_in  = len(flat_indices) if n_quads > 0 else 0
        I_total     = cf.sum_intensity.sum()
        I_quad      = cf.sum_intensity[quad_id_col > -1].sum() if n_quads > 0 else 0.0

        self.logger.info('  Quadruplets found : %d', n_quads)
        self.logger.info('  Peaks assigned    : %d / %d (%.1f%%)',
            n_peaks_in, cf.nrows, 100.0 * n_peaks_in / cf.nrows)
        self.logger.info('  Intensity covered : %.1f%%', 100.0 * I_quad / I_total)

        return quad_arr

    
    # ------------------------------------------------------------------
    #  High-level helpers
    # ------------------------------------------------------------------
    def initialize_peaksubsets(self, chunk_type='scans', n_eta_bins=None, reset=False):
        """ Initialization and checks for match_friedel_pairs_by_chunks(): """

        # identify chunking scheme and pair_type
        assert chunk_type in ['frames','scans', 'eta_bins'], "chunk_type must be 'frames', 'scans' or 'eta_bins' " 
        self.logger.info('Chunk type for pairing: %s', chunk_type)
        pair_type = 'eta' if chunk_type=='eta_bins' else 'omega'

        # initialize peaksubset and sort accordingly with pair_type
        if self.PeakSubsets is None or reset:
            logger.info('reset Peaksubsets with %d eta bins', n_eta_bins)
            self.PeakSubsets = None
            self.set_peak_subsets(n_eta_bins)
            
        self.sort_peak_subsets(pair_type)
        Psub = self.PeakSubsets
        
        LUT         = getattr(Psub, chunk_type+'_LUT', None)
        pairs_list  = getattr(Psub, chunk_type+'_subsets', None)
        valid_pairs = getattr(Psub, 'valid_'+chunk_type+'_subsets', None)

        if (pair_type) == 'omega' and Psub.is_half_acquisition:
            raise ValueError(
                    "Half acquisition. Cannot find omega pairs.")
        if LUT is None:
            raise RuntimeError(
                'Lookup table not found in PeakSubsets')
        if pairs_list is None:
            raise RuntimeError(
                'Pairs list not found in PeakSubsets')
        if valid_pairs is None:
            raise RuntimeError(
                'No list of valid subset in PeakSubsets.'
                'Run self.PeakSubsets.find_valid_subsets()') 
            
        return pairs_list, valid_pairs, LUT


        
    def estimate_search_scales(self, idx1, idx2, dist_cutoff, pair_type='omega'):
        """
        Pilot pairing on a small subset to estimate typical pair distances
        along each search-space dimension and derive normalised weights.

        Parameters
        ----------
        idx1, idx2   : index arrays into self.cf (reference / partner subset)
        dist_cutoff  : float — initial normalised cutoff for the pilot KDTree
        verbose      : bool

        Returns
        -------
        scales       : (5,) array of median distances per dimension
        weights_dict : dict of normalised weights keyed by dimension name
        """
        self.logger.debug('pilot search: searching pairs within %.4f'
                  ' dist tolerance radius', dist_cutoff)
        
        cf = self.cf
        s1, _ = _search_space(cf, idx1, weights=None)
        s2, _ = _search_space(cf, idx2, weights=None, flip=pair_type)

        # iterative KDTree search — widen dist_cutoff in steps until enough pairs are found
        n_total   = min(len(idx1), len(idx2))
        min_pairs = int(0.5 * n_total) if n_total < 2000 else 1000
        step_size = dist_cutoff / 10

        dists = rows = cols = np.array([])
        for i in range(1, 11):
            cutoff_i = step_size * i
            dij = _compute_dist_matrix(s1, s2, cutoff_i)
            dists, rows, cols = _clean_dist_matrix(dij, verbose=False)
            self.logger.info('  pilot search: cutoff=%.4f  pairs found: %d',
                cutoff_i, len(rows))
            if len(rows) >= min_pairs:
                break

        if len(rows) < min_pairs:
            raise RuntimeError(
                'Only {found} pilot pairs found after 10 iterations '
                '(cutoff={cutoff:.4f}, threshold={thr}). '
                'Check peak overlap between the two frames, or increase '
                'dist_cutoff.'.format(
                    found=len(rows), cutoff=dist_cutoff, thr=min_pairs))
        
        # per-dimension residuals
        d_gv, d_eta, d_logI = _physical_pair_distance( cf, idx1[rows], idx2[cols], pair_type)

        # scale factors: median absolute deviation — robust to outliers
        d_gv_med  = np.nanmedian(d_gv)
        d_eta_med = np.nanmedian(d_eta)
        d_logI_med = np.nanmedian(d_logI)

        scales = np.array([d_gv_med, d_gv_med, d_gv_med, d_eta_med, d_logI_med])
        scales = np.where(scales < 1e-3 * np.median(scales), np.inf, scales) # near zero -> inf 
        d_eucl = np.linalg.norm(scales[np.isfinite(scales)])

        weights     = 1.0 / scales
        weights    /= weights.max()
        keys        = 'gx,gy,gz,eta,I'.split(',')
        weights_dict = {k: w for k, w in zip(keys, weights)}

        self.logger.info('Typical (median) pair distances in non-weighted search space:')
        self.logger.info('  gvecs: %.4f', d_gv_med)
        self.logger.info('  eta:   %.4f deg', d_eta_med)
        self.logger.info('  logI:  %.4f', d_logI_med)
        self.logger.info('  Euclidean dist: %.4f', d_eucl)
        return scales, weights_dict
        

    def run_pairing(self, idx1, idx2,
                    pair_type='omega',
                    filter_mode = 'strict',
                    weights=None,
                    t_max = 120):
        """ class-level launcher for _run_pairing"""

        verbose = self.logger.getEffectiveLevel() < 30

        results =  _run_pairing(self.cf, idx1, idx2,
                                pair_type   = pair_type,
                                filter_mode = filter_mode,
                                weights     = weights,
                                tol_gv      = self.tol_gv,
                                tol_eta     = self.tol_eta,
                                tol_logI    = self.tol_logI,
                                n_steps     = self.n_steps,
                                verbose     = verbose,
                                t_max       = t_max)
        
        self.outputs.append(results)

        
    def merge_outputs(self, cf_target=None, pair_type='omega',
                      drop_broken=False, drop_unpaired = False, deduplicate=False):
        """
        Merge all dicts of pair_ids in self.outputs into self.cf.
        For chunk-based matching, pair-id labels are made globally unique
        across chunks.

        Parameters
        ----------
        cf_target     : columnfile to write into (defaults to self.cf)
        pair_type     : 'omega' or 'eta'
        drop_broken   : bool — if True, remove singleton pair_ids from cf
        drop_unpaired : bool — if True, remove non-paired peaks from cf
        """
        if cf_target is None:
            cf_target = self.cf

        # add pair_id column (or reset it if already there)
        label_col_name = '{}_pair_id'.format(pair_type)
        cf_target.addcolumn(np.full(cf_target.nrows, -1, dtype=int), label_col_name)
        label_col = cf_target.getcolumn(label_col_name)

        self.logger.info('--- MERGING OUTPUTS  [%s-pair] ---', pair_type)
        self.logger.info('  %d chunk(s) to merge  (%d peaks total)',
            len(self.outputs), cf_target.nrows)

        valid_outputs = [o for o in self.outputs if len(o['idx']) > 0]

        if valid_outputs:
            # ── collect all candidates ────────────────────────────────────────
            # idx_ref and idx_partner are the two interleaved entries per pair
            cand_dist  = np.concatenate([o['d_gv'][0::2]    for o in valid_outputs])
            cand_i1    = np.concatenate([o['idx'][0::2]     for o in valid_outputs])
            cand_i2    = np.concatenate([o['idx'][1::2]     for o in valid_outputs])

            n_cand = len(cand_dist)
            self.logger.info('  %d candidate pairs from %d chunks',
                             n_cand, len(valid_outputs))

            # ── deduplicate exact duplicates first (same i1, i2 pair) ────────
            # sort by (i1, i2) then keep the one with smallest gvec_dist
            if deduplicate:
                lex_order  = np.lexsort((cand_dist, cand_i2, cand_i1))
                cand_dist  = cand_dist[lex_order]
                cand_i1    = cand_i1[lex_order]
                cand_i2    = cand_i2[lex_order]
                # first occurrence after lex sort has the smallest dist for each (i1,i2)
                pair_keys        = np.stack([cand_i1, cand_i2], axis=1)
                _, unique_pairs  = np.unique(pair_keys, axis=0, return_index=True)
                cand_dist  = cand_dist[unique_pairs]
                cand_i1    = cand_i1[unique_pairs]
                cand_i2    = cand_i2[unique_pairs]

                n_unique = len(cand_dist)
                self.logger.info('  %d unique candidate pairs after deduplication',
                                 n_unique)
            
                # ── greedy matching: sort by quality, accept if both peaks free ───
                quality_order = np.argsort(cand_dist)
                cand_dist     = cand_dist[quality_order]
                cand_i1       = cand_i1[quality_order]
                cand_i2       = cand_i2[quality_order]

                free = np.ones(cf_target.nrows, dtype=bool)
                both_free  = free[cand_i1] & free[cand_i2]
            
                # iterate in quality order and fix the sequential dependency with a small loop
                # (sequential because accepting pair k changes freedom for pair k+1)
                accepted   = np.zeros(n_unique, dtype=bool)

                for k in range(n_unique):
                    i1k, i2k = cand_i1[k], cand_i2[k]
                    if free[i1k] and free[i2k]:
                        accepted[k]  = True
                        free[i1k]    = False
                        free[i2k]    = False

                n_accepted = accepted.sum()
                self.logger.info('  %d pairs accepted by greedy matching', n_accepted)

                # ── assign globally unique pair_ids and scatter into label_col ───
                pair_ids         = np.arange(n_accepted)
                accepted_i1      = cand_i1[accepted]
                accepted_i2      = cand_i2[accepted]
            else:
                n_accepted   = len(cand_i1)
                pair_ids     = np.arange(n_accepted)
                accepted_i1  = cand_i1
                accepted_i2  = cand_i2
        
            label_col[accepted_i1] = pair_ids
            label_col[accepted_i2] = pair_ids
            next_id = n_accepted

        # ── validation ────────────────────────────────────────────────────
        paired  = label_col[label_col != -1]
        n_total = cf_target.nrows
        I_total = cf_target.sum_intensity.sum()

        if len(paired) > 0:
            counts     = np.bincount(paired)
            singletons = (counts == 1).sum()

            if singletons > 0:
                self.logger.warning('  %d pair_id(s) appear only once — broken pairs',
                    singletons)
                if drop_broken:
                    # identify broken pair_ids and reset to -1
                    broken_ids  = np.where(counts == 1)[0]
                    broken_mask = np.isin(label_col, broken_ids)
                    label_col[broken_mask] = -1
                    self.logger.info('  %d peaks reassigned to -1', broken_mask.sum())
            else:
                self.logger.info('  OK: all pair_ids appear exactly twice')

        paired_mask = label_col != -1
        n_paired    = paired_mask.sum()
        I_paired    = cf_target.sum_intensity[paired_mask].sum()
        self.logger.info('  %d pairs, %d peaks paired (%.1f%%)',
            int(n_paired//2), n_paired, 100.0 * n_paired / n_total)
        self.logger.info('  %.1f%% of total intensity matched',
            100.0 * I_paired / I_total)

        if drop_unpaired:
            self.logger.info('FILTERING: %d peaks not matching in any pair -> remove from cf',
                (~paired_mask).sum())
            cf_target.filter(paired_mask)
            self.logger.info('Done')
            

# ─────────────────────────────────────────────
#  Friedel Pair Indexing Helpers
# ─────────────────────────────────────────────
def _wrap_eta(eta):
    return (eta + 180.0) % 360.0 - 180.0    # maps to (-180, 180]

def _search_space(cf, idx, weights=None, flip=None):
    """
    Builds the 5D search space (gx, gy, gz, eta, logI) for friedel pairs.
    User weights allow fine-tuning contributions per dimension. 

    Parameters
    ----------
    idx     : index array into cf
    weights : dict, optional keys 'gx','gy','gz','eta','I'  (default all 1.0)
              user-defined scaling factors for search dimensions
              > 1 increase importance, < 1 decrease it, ~0 mutes it

    flip    : None | 'omega' | 'eta'.
              Flip partner space for omega or eta-pairs search
    Returns
    -------
    space   : (N, 5) search space rescaled according to input weights
    weights : rescaling weights dict
    """
    _w = {'gx': 1., 'gy': 1., 'gz': 1., 'eta': 1., 'I': 1.}
    if weights is not None:
        _w.update(weights)

    gx  = cf.gx[idx]
    gy  = cf.gy[idx]
    gz  = cf.gz[idx]
    eta = cf.eta[idx]
    logI = np.log10(np.maximum(cf.sum_intensity[idx], 1e-10))

    if flip == 'omega':
        gx, gy, gz = -gx, -gy, -gz
        eta = _wrap_eta(180.0 - eta)
    elif flip == 'eta':
        gx, gy, gz = -gx, -gy, -gz
        eta = _wrap_eta(180.0 + eta)

    space = np.column_stack([
        gx   * _w['gx'],
        gy   * _w['gy'],
        gz   * _w['gz'],
        eta  * _w['eta'],
        logI * _w['I']])

    return space, _w

# slow whith large dist_cutoff
#def _compute_dist_matrix__(space_ref, space_partner, dist_cutoff):
#    """
#    Build sparse distance matrix between reference and partner search spaces.
#    """
#    tree_ref     = scipy.spatial.cKDTree(space_ref)
#    tree_partner = scipy.spatial.cKDTree(space_partner)
#    dij = csr_matrix(tree_ref.sparse_distance_matrix(tree_partner, dist_cutoff))
#    return dij

def _compute_dist_matrix(space_ref, space_partner, dist_cutoff):
    """
    Find nearest neighbour in space_partner for each point in space_ref,
    within dist_cutoff. Returns as a sparse matrix. 
    """
    tree_partner = scipy.spatial.cKDTree(space_partner)
    # k=1: only the single nearest neighbour
    dists, col_ind = tree_partner.query(
        space_ref,
        k=1,
        distance_upper_bound=dist_cutoff,
        workers=1)     

    # filter out points with no neighbour within cutoff
    valid    = np.isfinite(dists)
    row_ind  = np.where(valid)[0]
    col_ind  = col_ind[valid]
    dists    = dists[valid]

    dij = csr_matrix(
        (dists, (row_ind, col_ind)),
        shape=(len(space_ref), len(space_partner)))
    return dij

def _clean_dist_matrix(dij, verbose=True):
    """
    Clean sparse distance matrix to ensure 1-to-1 peak matching.
    Returns row,col indices in dij and corresponding pair distances (dist)
    """
    if dij.nnz == 0:
        return np.array([]), np.array([], dtype=int), np.array([], dtype=int)

    nnz_col = np.count_nonzero(dij.count_nonzero(axis=0))
    nnz_row = np.count_nonzero(dij.count_nonzero(axis=1))
    N_candidates = min(nnz_col, nnz_row)

    dij_inv      = dij.copy()
    dij_inv.data = 1.0 / dij_inv.data

    # pass 1: for each column keep only the closest row
    row_ind  = np.array(dij_inv.argmax(axis=0)).ravel()
    col_ind  = np.arange(dij_inv.shape[1])
    max_vals = np.array(dij_inv[row_ind, col_ind]).ravel()
    valid    = max_vals > 0
    dij_inv  = csr_matrix(
        (max_vals[valid], (row_ind[valid], col_ind[valid])),
        shape=dij_inv.shape)

    # pass 2: for each row keep only the closest column
    col_ind  = np.array(dij_inv.argmax(axis=1)).ravel()
    row_ind  = np.arange(dij_inv.shape[0])
    max_vals = np.array(dij_inv[row_ind, col_ind]).ravel()
    valid    = max_vals > 0
    dij_inv  = csr_matrix(
        (max_vals[valid], (row_ind[valid], col_ind[valid])),
        shape=dij_inv.shape)

    dij_inv.eliminate_zeros()
    rows, cols = dij_inv.nonzero()
    dists      = 1.0 / dij_inv.data

    if verbose:
        print('{} pairs kept out of {} possible matches'.format(dij_inv.nnz, N_candidates))
    return dists, rows, cols

def _physical_pair_distance(cf, idx1, idx2, pair_type='omega'):
    """
    Physical distance in g-vector norm, eta and log(intensity) between matched pairs.

    Parameters
    ----------
    idx1 / idx2 : indices of paired subsets in cf
    pair_type   : 'eta' or 'omega'

    Returns
    -------
    d_gvec, d_eta, d_logI
    """
    # g-vector distance :(hkl),-(hkl) pairs: g + -g ~ 0 
    dgx = cf.gx[idx1] + cf.gx[idx2]
    dgy = cf.gy[idx1] + cf.gy[idx2]
    dgz = cf.gz[idx1] + cf.gz[idx2]
    d_gvec = np.sqrt(dgx**2 + dgy**2 + dgz**2)
    # eta distance
    if pair_type == 'omega':
        diff = cf.eta[idx1] - (180-cf.eta[idx2])   # absolute difference without modulo 
        d_eta = np.abs( _wrap_eta(diff) )          # wrap to [-180,180)
    else:
        diff = cf.eta[idx1] - (180+cf.eta[idx2])
        d_eta = np.abs( _wrap_eta(diff) )
    # logI distance (intensity mismatch)
    logI_1 = np.log10(cf.sum_intensity[idx1])
    logI_2 = np.log10(cf.sum_intensity[idx2])
    d_logI = np.abs(logI_1 - logI_2)
    
    return d_gvec, d_eta, d_logI

def _filter_pairs_physical(cf, idx1, idx2, pair_type, tol_gvec=0.05, tol_eta=1.0, tol_logI=1.0):
    """
    Filter matched pairs using physical tolerances.
    Return bool mask on idx1,idx2 for pairs that pass all tolerances
    """
    d_gvec, d_eta, d_logI = _physical_pair_distance(cf, idx1, idx2, pair_type)
    mask = (d_gvec < tol_gvec) & (d_eta < tol_eta) & (d_logI < tol_logI)
    return mask

def _wrong_omegas(cf, idx1, idx2, tol=1e-4):
    """ remove dubious pairs for eta pairing"""
    diff = cf.omega[idx1] - cf.omega[idx2]
    wrong =  np.abs(_wrap_eta(diff)) < tol
    return wrong

def _physical_to_normalised_cutoff(tol_gvec, tol_eta, tol_logI, weights=None):
    """
    Convert physical tolerances to a normalised dist_cutoff for KDTree
    search in the rescaled search space.

    Returns
    -------
    dist_max : float — normalised cutoff for KDTree
    """
    _w = {'gx': 1., 'gy': 1., 'gz': 1., 'eta': 1., 'I': 1.}
    if weights is not None:
        _w.update(weights)
    # normalised tolerance per dimension
    tol_gvec_n = tol_gvec * _w['gx']
    tol_eta_n  = tol_eta  * _w['eta']
    tol_logI_n = tol_logI * _w['I']

    # L2 norm of per-dimension tolerances, ignoring inf (muted dimensions)
    components = np.array([tol_gvec_n, tol_gvec_n, tol_gvec_n, tol_eta_n, tol_logI_n])
    components = components[np.isfinite(components)]
    dist_max   = np.sqrt(np.sum(components**2))
    return dist_max

def _run_pairing(cf, idx1, idx2,
                 pair_type='omega',
                 filter_mode='relaxed',
                 weights=None,
                 tol_gv=0.05,
                 tol_eta=0.5,
                 tol_logI=np.inf,
                 n_steps=5,
                 t_max = 120,
                 verbose=True):
    """
    Find Friedel pairs between two subsets of cf defined by idx1, idx2.
    Does NOT write to cf: returns a result dict

    Parameters
    ----------
    cf          : ImageD11 columnfile (read-only)
    idx1, idx2  : int arrays — indices into cf defining the two subsets
    pair_type   : 'omega' or 'eta'
    weights     : dict | None — search-space dimension weights,
                  keys: 'gx', 'gy', 'gz', 'eta', 'I'
    tol_gv      : float — g-vector tolerance (Å⁻¹)
    tol_eta     : float — eta angle tolerance (degrees)
    tol_logI    : float — log10(intensity) tolerance (np.inf to ignore)
    n_steps     : int   — number of distance increments in iterative search
    t_max       : max execution time (seconds) — stop early if exceeded
    verbose     : bool
    filter_mode: 'strict' or 'relaxed' (default: strict) 
            "strict"  -> tolerance filter applied at each KDtree search iteration.
                         slower but less likely to miss pairs 
            "relaxed" -> pairs are accumulated without filtering at each KDTree
                        search iterations. tolerance filter is applied once after
                        search loop is completed.
                        Faster but may miss some pairs. 
   
    result dict: keys 'idx', 'pair_id'
                 idx     : (2*N,), interleaved pairs indices [r0,p0,r1,p1,...]
                 pair_id : (2*N,), unique pair identifier
    """
    _w     = weights if weights is not None else weights

    def _print(msg):
        if verbose:
            print(msg)

    def _section(title):
        if verbose:
            print('\n{}'.format('-' * 50))
            print('  {}'.format(title))
            print('{}'.format('-' * 50))

    def _progress_bar(n_paired, n_total, width=30):
        frac   = n_paired / n_total if n_total > 0 else 0
        filled = int(width * frac)
        bar    = '#' * filled + '.' * (width - filled)
        return '  [{bar}] {pct:.1f}%  ({n}/{tot})'.format(
            bar=bar, pct=100 * frac, n=n_paired, tot=n_total)

    # ── 1. Initialisation ─────────────────────────────────────────────
    _section('FRIEDEL PAIR SEARCH  [{}]'.format(pair_type))

    free1    = np.ones(len(idx1), dtype=bool)
    free2    = np.ones(len(idx2), dtype=bool)
    partner1 = np.full(len(idx1), -1, dtype=int)
    dist_max = _physical_to_normalised_cutoff(tol_gv, tol_eta, tol_logI, _w)
    dist_steps = np.linspace(dist_max / n_steps, dist_max, n_steps)
    n_total    = min(len(idx1), len(idx2))
        
    _print('  Subsets    : {} + {} peaks'.format(len(idx1), len(idx2)))
    _print('  Tolerances : gv={}  eta={} deg  logI={}'.format(
            tol_gv, tol_eta, tol_logI))
    _print('  dist_max   : {:.4f}  ({} steps)'.format(dist_max, n_steps))

    # ── 2. Iterative matching ─────────────────────────────────────────
    _section('ITERATIVE PAIR MATCHING [{} steps]'.format(len(dist_steps)) )

    t_start  = time.perf_counter()
    timed_out = False

    for i, dist_cutoff in enumerate(dist_steps):

         # ── timeout check ───────────────
        elapsed = time.perf_counter() - t_start
        if elapsed > t_max:
            timed_out = True
            logger.warning(
                '_run_pairing timed out after %.1f s at step %d/%d '
                '(cutoff=%.4f). Returning %d pairs matched so far. '
                'Consider tightening tolerances or adjusting weights.',
                elapsed, i + 1, n_steps, dist_cutoff, (~free1).sum())
            break
        
        if free1.sum() < 1 or free2.sum() < 1:
            _print('  No free peaks remaining — stopping early.')
            break
            
        # ── indices into global cf ─────────
        cur_idx1 = idx1[free1]
        cur_idx2 = idx2[free2]

         # ── search space construction ─────
        space_1, _ = _search_space(cf, cur_idx1, _w)
        space_2, _ = _search_space(cf, cur_idx2, _w, flip=pair_type)

        # ── KDTree distance matrix
        dij = _compute_dist_matrix(space_1, space_2, dist_cutoff)
        if dij.nnz == 0:
            _print('  step {i}/{n}  cutoff={c:.4f}  -> no candidates'.format(
                i=i+1, n=n_steps, c=dist_cutoff))
            continue

        _, local_rows, local_cols = _clean_dist_matrix(dij, verbose=False)
        if len(local_rows) == 0:
            _print('  step {i}/{n}  cutoff={c:.4f}  -> no unique matches'.format(
                i=i+1, n=n_steps, c=dist_cutoff))
            continue
                
        # ── tolerance filtering (strict filter_mode)
        if filter_mode == 'strict':
            tol_mask = _filter_pairs_physical(
                            cf, cur_idx1[local_rows], cur_idx2[local_cols],
                            pair_type, tol_gv, tol_eta, tol_logI)
            
            wrong_om = _wrong_omegas(cf, idx1[local_rows], idx2[local_cols],
                                 tol=1e-4)
            tol_mask &= ~wrong_om

            n_new = tol_mask.sum()
            if n_new > 0:
                local_rows = local_rows[tol_mask]
                local_cols = local_cols[tol_mask]
        else:
            n_new = len(local_rows)
        if n_new == 0:
            _print('  step {i}/{n}  cutoff={c:.4f}  -> 0 pairs after filtering'.format(
                i=i+1, n=n_steps, c=dist_cutoff))
            continue

        # ── map local indices into global cf indices
        cf_idx1 = cur_idx1[local_rows]
        cf_idx2 = cur_idx2[local_cols]
        pos1    = np.searchsorted(idx1, cf_idx1)
        pos2    = np.searchsorted(idx2, cf_idx2)

        partner1[pos1] = cf_idx2
        free1[pos1]    = False
        free2[pos2]    = False

        n_paired_so_far = (~free1).sum()
        _print('  step {i}/{n}  cutoff={c:.4f}  -> +{new} pairs'.format(
            i=i+1, n=n_steps, c=dist_cutoff, new=n_new))
        _print(_progress_bar(n_paired_so_far, n_total))

    # ── 3. Assemble output ────────────────────────────────────────────
    _section('OUTPUT')

    paired_mask = ~free1
    n_pairs = paired_mask.sum()             
    _print('Total pairs found : {} / {}'.format(n_pairs, n_total))
    _print(_progress_bar(n_pairs, n_total))


    # relaxed filter_mode: apply tolerance filters once after loop
    tol_mask = _filter_pairs_physical(
            cf, idx1[paired_mask], partner1[paired_mask],
            pair_type, tol_gv, tol_eta, tol_logI)
        
    wrong_om = _wrong_omegas(cf, idx1[paired_mask], partner1[paired_mask],
                                 tol=1e-4)
    tol_mask &= ~wrong_om
        
    n_pairs = tol_mask.sum()
    idx1_paired = idx1[paired_mask][tol_mask]
    idx2_paired = partner1[paired_mask][tol_mask]
            
    _print('Applying final tolerance filters')
    _print('Total pairs after filtering : {} / {}'.format(n_pairs, n_total))
    _print(_progress_bar(n_pairs, n_total))
        
    if n_pairs < 1:
        _print('  WARNING: No pairs found — returning empty output.')
        result = {'idx'      : np.array([], dtype=int),
                  'pair_id'  : np.array([], dtype=int),
                  'd_gv'     : np.array([], dtype=float)}
        return result

    local_ids   = np.arange(n_pairs)
    # interleave [ref0, partner0, ref1, partner1, ...]
    idx_out    = np.empty(2 * n_pairs, dtype=int)
    pid_out    = np.empty(2 * n_pairs, dtype=int)
    dgv_out    = np.empty(2 * n_pairs, dtype=float)
    d_gv, _, _ = _physical_pair_distance(cf,idx1_paired,idx2_paired,pair_type)
    
    idx_out[0::2]    = idx1_paired;  idx_out[1::2]    = idx2_paired
    pid_out[0::2]    = local_ids;    pid_out[1::2]    = local_ids
    dgv_out[0::2]    = d_gv     ;    dgv_out[1::2]    = d_gv
        
    _print('\n{}\n'.format('-' * 50))

    result = {'idx'    : idx_out,
              'pair_id': pid_out,
              'd_gv'   : dgv_out}
    
    return result


def _worker_run_pairing(args):
    """ Module-level worker function for parallelised chunk pairing """
    pid, idx1, idx2, pair_type, filter_mode, weights, \
        tol_gv, tol_eta, tol_logI, n_steps, t_max, verbose = args
    
    out = _run_pairing(
        _shared_cf,
        idx1, idx2,
        pair_type=pair_type,
        filter_mode=filter_mode,
        weights=weights,
        tol_gv=tol_gv,
        tol_eta=tol_eta,
        tol_logI=tol_logI,
        n_steps=n_steps,
        t_max=t_max,
        verbose=verbose)
    return pid, out

def _worker_initializer():
    """Re-trigger numba JIT warm-up in the child process."""
    import numba
    numba.get_num_threads()  # forces TBB reinit in the child
    
# ─────────────────────────────────────────────
#  Other tools to work with Friedel pairs
# ─────────────────────────────────────────────
def get_pairs(cf, pair_type='omega'):
    """ returns (idx1, idx2) : index positions of pairs in cf """
    pair_id_name = pair_type + '_pair_id'
    pair_id_col = cf.getcolumn(pair_id_name)        
    paired = np.flatnonzero(pair_id_col > -1)
    if cf.sortedby == pair_id_name:
        idx1  = paired[::2]
        idx2  = paired[1::2]
    else:
    # sort so that pair_ids run [0,0,1,1,2,2,...] 
        order = np.argsort(pair_id_col[paired], kind='stable')
        idx1  = paired[order][::2]
        idx2  = paired[order][1::2]
    return idx1, idx2
    

def get_pair_distances(cf, pair_type='omega'):
    """
    Physical distance in g-vector norm, eta and log(intensity) between matched pairs.

    Returns 
    ---------
    distances    : dict with keys 'd_gvec', 'd_eta', 'd_logI'
    (idx1, idx2) : index positions of pairs in cf
    """
    idx1, idx2 = get_pairs(cf, pair_type)
    d_gv, d_eta, d_logI = _physical_pair_distance(cf, idx1, idx2, pair_type=pair_type)
    return {'d_gv': d_gv, 'd_eta':d_eta, 'd_logI':d_logI}, (idx1, idx2)

    
def plot_pair_distances(cf, pair_type='omega', bins=50, log_scale=False, **kwargs):
    """
    Plot distributions of pair distances along each physical dimension.

    Parameters
    ----------
    cf        : ImageD11 columnfile. Must contain a '{pair_type}_pair_id' column.
    pair_type : 'omega' or 'eta', used for column name lookup and title
    bins      : number of histogram bins
    log_scale : bool, log y-axis
    """
    dists, _ = get_pair_distances(cf, pair_type)
        
    # ── plot ─────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    fig.suptitle(
        '{} pair distance distributions  (N={} pairs)'.format(
            pair_type, len(dists['d_gv']) // 2),
        fontsize=12, fontweight='bold')

    datasets = [
        (dists['d_gv']  ,   r'gvec distance',             'steelblue'),
        (dists['d_eta'] ,  r'$\Delta\eta$ ($^\circ$)',    'darkorange'),
        (dists['d_logI'], r'$\Delta \log_{10}(\rm{Ints})$', 'seagreen')]

    for ax, (data, xlabel, color) in zip(axes, datasets):
        finite = data[np.isfinite(data)]
        ax.hist(finite, bins=bins, color=color, alpha=0.8,
                edgecolor='white', linewidth=0.4)
        ax.axvline(np.median(finite), color='k', lw=1.5, ls='--',
                   label='median={:.3f}'.format(np.median(finite)))
        ax.axvline(np.mean(finite),   color='k', lw=1.0, ls=':',
                   label='mean={:.3f}'.format(np.mean(finite)))
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel('counts', fontsize=10)
        ax.legend(fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)
        if log_scale:
            ax.set_yscale('log')

    plt.tight_layout()
    return fig, axes


def locate_eta_pairs(cf, pairs, ds=None, y0=0.):
    """
    Fit the centre of mass position of eta-pairs and write results
    back into cf as new columns 'sx' and 'sy'.
    Works only for scanning-3DXRD. For Box-beam, the relocation problems 
    has two more degrees of freedom and cannot be solved with only one pair. 

    Parameters
    ----------
    cf    : ImageD11 columnfile
    pairs : (idx1, idx2) arrays of cf indices for each pair
    y0    : float, beam centre offset
    """
    i1, i2 = pairs
    r  = np.radians(cf.omega)
    so = np.sin(r)
    co = np.cos(r)

    y = np.transpose((cf.dty[i1] - y0,
                      cf.dty[i2] - y0))

    R = np.transpose(((so[i1], co[i1]),
                      (so[i2], co[i2])), axes=(2, 0, 1))

    sx_pairs, sy_pairs = np.linalg.solve(R, y).T

    # ── write into cf. Overwrite pre-existing sx, sy
    sx = np.full(cf.nrows, np.nan)
    sy = np.full(cf.nrows, np.nan)
    sx[i1] = sx_pairs;  sx[i2] = sx_pairs
    sy[i1] = sy_pairs;  sy[i2] = sy_pairs

    for name, col in (('sx', sx), ('sy', sy)):
        cf.addcolumn(col, name)

    return sx, sy


def locate_omega_pairs(cf, pairs, ds=None, y0=None):
    """
    Fit the centre of mass position of omega-pairs and write results
    back into cf as new columns 'sx' and 'sy'.
    Works only for scanning-3DXRD. For Box-beam, the relocation problems 
    has two more degrees of freedom and cannot be solved with only one pair. 
    
    For omega-pairs the linear system used in 'locate_eta_pairs' is near-singular
    (paired peaks are ~180 degrees apart in omega), so the position is recovered from
    two-theta and dty:
    dx  : x-shift from rotation centre, derived from tth asymmetry
    ylab: y-shift in lab frame, taken as dty - y0
    Then rotate (dx, ylab) back into sample coordinates using omega.
    xlab, ylab and omega are averaged to make sure both peaks in each pair have the same (sx, sy) coordinates
    """
    i1, i2 = pairs
    L = cf.parameters.get('distance')
    if ds is not None and ds.dtymotor == 'diffty':    # diffty is in mm on TDXRD sation
        L /= 1000
    if hasattr(ds, 'y0'):
        y0 = ds.y0
    elif y0 is None:
        y0 = 0.5 * (ds.ymax + ds.ymin)
  
    # ── dx from two-theta asymmetry — 
    tantth1 = np.tan(np.radians(cf.tth[i1]))
    tantth2 = np.tan(np.radians(cf.tth[i2]))
    dx      = L * (tantth1 - tantth2) / (tantth1 + tantth2)

    # ── ylab from averaged dty — 
    dy   = (cf.dty[i1] - cf.dty[i2]) / 2                      
    # i1 sees +dx, +dy;  i2 sees -dx, -dy
    xlab_i1 =  dx;   xlab_i2 = -dx
    ylab_i1 = y0 + dy;  ylab_i2 = y0 - dy

    # ── averaged rotation matrix 
    r  = np.radians(cf.omega)
    co, so = np.cos(r), np.sin(r)
    co_av  = (co[i1] - co[i2]) / 2                            
    so_av  = (so[i1] - so[i2]) / 2                           

    # ── lab -> sample rotation ───
    #   sx =  co_av * xlab + so_av * ylab
    #   sy = -so_av * xlab + co_av * ylab
    sx_i1 =  co_av * xlab_i1 + so_av * ylab_i1    
    sy_i1 = -so_av * xlab_i1 + co_av * ylab_i1    

    sx_i2 = -co_av * xlab_i2 - so_av * ylab_i2    # (-R applied)
    sy_i2 =  so_av * xlab_i2 - co_av * ylab_i2    

    # sx_i1 == sx_i2 and sy_i1 == sy_i2 by construction — average for numerical safety
    sx_pairs = 0.5 * (sx_i1 + sx_i2)
    sy_pairs = 0.5 * (sy_i1 + sy_i2)

    # ── scatter into cf-length arrays ───
    sx = np.full(cf.nrows, np.nan)
    sy = np.full(cf.nrows, np.nan)
    sx[i1] = sx_pairs;  sx[i2] = sx_pairs
    sy[i1] = sy_pairs;  sy[i2] = sy_pairs

    # ── write to cf ──
    for name, col in (('sx', sx), ('sy', sy)):
        if name in cf.titles:
            cf.getcolumn(name)[:] = col
        else:
            cf.addcolumn(col, name)

    return sx, sy


def update_geometry_fpairs(cf, ds=None, add_xyz_lab=False, relocate_pairs=True):
    """
    Update diffraction geometry for (scanning)-3DXRD data using Friedel-pairs (omega-pair) symmetry correction:
    similar to cf.updateGeometry() but uses Friedel-pairs to correct for the offset from rotation-center.
    Details in Jacob et al. 2024, https://doi.org/10.1107/S1600576724009634

    Update the following columns in the ColumnFile:
    'tth', 'eta', 'ds', 'gx', 'gy', 'gz'
    If ``add_xyz_lab`` True:
    'xl', 'yl', 'zl' also updated
    If ``relocate_pairs`` True: 
    add 'sx' and 'sy' coordinates computed from omega pairs

    Parameters
    ----------
   -  cf : ImageD11 ColumnFile 
   -  ds : ImageD11.sinogram dataset (optional).
   -  add_xyz_lab : bool, optional. Add corrected lab-frame coordinates to cf
   - relocate_fpairs
    """
    if cf.parameters is None:
        raise ValueError(
            'No parameters in cf'
        )
    pars = cf.parameters
    assert "omega_pair_id" in cf.titles, 'no omega pairs in cf. Compute them first'
    if "sc" in cf.titles and "fc" in cf.titles:
            sc, fc = cf.sc, cf.fc
    elif "xc" in cf.titles and "yc" in cf.titles:
            sc, fc = cf.xc, cf.yc
    else:
        raise Exception("columnfile file misses xc/yc or sc/fc")
    
    # xyz coordinates in non-corrected geometry.
    comp = transform.Ctransform(pars.parameters)
    xyz_non_corr = comp.sf2xyz(sc, fc)
    i1, i2 = get_pairs(cf, pair_type = 'omega')
    
    # SOURCE RELOCATION: use raw (non-corrected) tth
    if relocate_pairs:
        tth_raw = comp.xyz2geometry( xyz_non_corr, cf.omega, tx=0, ty=0, tz=0)[:,0]
        cf.addcolumn(tth_raw, 'tth')
        _ = locate_omega_pairs(cf, (i1,i2), ds)
        
    # FRIEDEL PAIR CORRECTION
    # Back-rotate half of the peaks to give a more intuitive geometry where paired peaks are aligned with the diffraction source origin
    R = np.diag([-1,-1,1])     
    xyz_corr = 1/2 * (xyz_non_corr[i1] - xyz_non_corr[i2].dot(R))
   
    # merge arrays to go back to original colf size
    xyz_corr_full = np.full( (cf.nrows, 3), np.nan, dtype=xyz_non_corr.dtype)
    xyz_corr_full[i1] = xyz_corr
    xyz_corr_full[i2] = xyz_corr.dot(-R)
    del xyz_corr

    # compute geometry with corrected peaks
    out = comp.xyz2geometry( xyz_corr_full, cf.omega, tx=0, ty=0, tz=0 )
    for i,name in enumerate(("tth", "eta", "ds", "gx", "gy", "gz")):
        cf.addcolumn( out[:,i], name )
    if add_xyz_lab:
        for i,name in enumerate(['xl','yl','zl']):
            cf.addcolumn( xyz_corr_full[:,i], name )
        

        
        




