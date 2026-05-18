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

# TO DO:
#  - add updage_geometry function to compute corrected coordinates (gvecs, tth, ds, etc.).

# written by Jean-Baptiste Jacob, jean-baptiste.jacob@esrf.fr
# Last updated: 18 May 2026


import os, sys
import logging
from contextlib import contextmanager
import numpy as np
from collections import namedtuple
from tqdm import tqdm
import multiprocessing as mp

import scipy.spatial
from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt
from ImageD11 import columnfile, transform

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)



#  ──────────────────────────────────────────────────────────────────────────────────────────
# PeakSubset: peak sorting and subset selection for pairing
# ──────────────────────────────────────────────────────────────────────────────────────────
Pair_dty = namedtuple('Pair', ['yi_hi', 'yi_lo', 'dty_hi', 'dty_lo'])
Pair_omega_dty = namedtuple('Pair', ['yi_A', 'yi_B', 'omi_A', 'omi_B', 'dty_A', 'dty_B', 'omega_A', 'omega_B'])

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

class PeakSubsets:
    """
    Make list of paired subsets (subset_1, subset_2) from whole dataset for matching Friedel pairs. 
    Subsets are either pairs of mirror dty scans or mirror omega-dty frames (omega-pairs)
    columnfile is sorted by sinogram order (dty, omega) and subset indices are stored in
    lookup tables (LUT) for fast selection. 
    """
    def __init__(self, cf, ds):
        self.columnfile = cf
        self.dataset = ds
        self.frames_LUT = None
        self.scans_LUT = None
        self.scans_subsets = None
        self.frames_subsets = None
        self.valid_scans_subsets = None
        self.valid_frames_subsets = None
        self._set_aliases()
        self._compute_subsets()

    def _set_aliases(self):
        self.obincens = self.dataset.obincens
        self.obinedges = self.dataset.obinedges
        self.ybincens = self.dataset.ybincens
        self.ybinedges = self.dataset.ybinedges

    def _compute_subsets(self):
        self.get_scan_subsets()
        self.get_frames_subsets()
                
    # --- Internal methods ---

    def get_scan_subsets(self, y0=None):
        """
        Find mirror pairs of dty scans on each side of the central y0:
            scan hi :  dty = y0 + Δy
            scan lo :  dty = y0 - Δy

        Sets ds.dty_pairs: list of Pair namedtuples ('Pair', ['yi_hi', 'yi_lo', 'dty_hi', 'dty_lo'])
        where yi_* are indices into ds.ybincens.
        """
        if y0 is None:
            y0 = self.dataset.y0
        if y0 != 0:
            self.dataset.correct_bins_for_half_scan(y0)
            y0 = self.dataset.y0

        yc       = self.ybincens.copy()
        n_y      = len(yc)
        yi0      = n_y // 2

        pairs = []
        for yi in range(0, yi0 + 1):
            pairs.append(_make_dty_pair(yi0 + yi, yi0 - yi, yc))

        self.scans_subsets = pairs
        logger.info("[get_scans_subsets] %d pairs of dty scans (y0=%s)", len(pairs), y0)

    def get_frames_subsets(self, y0=None):
        """
        Find mirror pairs of frames (A, B) in omega and dty:
            frame A :  dty = y0 + Δy ,  omega = ω
            frame B :  dty = y0 - Δy ,  omega = (ω + 180) mod 360

        Sets ds.omega_dty_pairs: list of Pair namedtuples
        ('Pair', ['yi_A','yi_B','omi_A','omi_B', 'dty_A','dty_B','omega_A','omega_B'])
        with indices yi*, omi* into ds.ybincens and ds.obincens
        """
        if y0 is None:
            y0 = self.dataset.y0
        if y0 != 0 and self.dataset.ybin_real_mask is None:
            self.dataset.correct_bins_for_half_scan(y0)

        yc       = self.ybincens.copy()
        oc       = self.obincens.copy()
        n_y      = len(yc)
        n_o      = len(oc)
        k        = n_o // 2          # index shift for a 180-deg omega step

        # Sanity checks
        if self.dataset.omax - self.dataset.omin < 181:
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

    
    def sort_by_sinogram(self):
        """
        Sort columnfile in sinogram space (dty primary, omega secondary)
        and build frame-level and scan-level lookup tables for quick peak selection
        by omega and dty bins. 

        Sets
        self.scans_LUT : dict { yi, slice(lo, hi), n_peaks, total_intensity) }
        self.frames_index : dict { (omi, yi), slice(lo, hi), n_peaks, total_intensity) }
        
        y_i / om_y_i indices match those in self.scans_subsets and self.frames_subsets
        Frames / scans with zero peaks are absent from the dict.
        """
        oc  = self.obincens
        yc  = self.ybincens
        obe = self.obinedges
        ybe = self.ybinedges

        if (self.columnfile.sortedby == "dty_omega") and self.frames_LUT:
            n_frames_total = len(oc) * len(yc)
            n_frames_empty = n_frames_total - len(self.frames_LUT)
            n_scans_total  = len(yc)
            n_scans_empty  = n_scans_total - len(self.scans_LUT)
            
            logger.info(
                "[sort_by_sinogram] peakfile already sorted\n%d peaks sorted; "
                "%d non-empty frames indexed (%d frames have no peaks).",
                self.columnfile.nrows, len(self.frames_LUT), n_frames_empty)
            return
            

        # Initialize scans / frames LUT
        setattr(self, 'frames_LUT', {})
        setattr(self, 'scans_LUT' , {})
        BinEntry = namedtuple('BinEntry', ['slice', 'npeaks', 'sumI'])  # 
        
        # 1. Binning and Sorting by dty and omega
        i_omega = np.clip(
            np.searchsorted(obe, self.columnfile.omega, side="right") - 1, 0, len(oc) - 1)
        j_dty   = np.clip(
            np.searchsorted(ybe, self.columnfile.dty,   side="right") - 1, 0, len(yc) - 1)

        # Sort in-place: dty primary, omega secondary
        order = np.lexsort((i_omega, j_dty))
        self.columnfile.reorder(order)
        self.columnfile.sortedby = "dty_omega"

        i_omega_s = i_omega[order]
        j_dty_s   = j_dty[order]

        # 2. Build frames_LUT: Mapping (i, j) -> BinEntry
        bin_pairs = np.column_stack([i_omega_s, j_dty_s])
        _, first_idx, counts = np.unique(bin_pairs, axis=0, return_index=True, return_counts=True)

        self.frames_LUT = {}
        for lo, cnt in zip(first_idx, counts):
            om_y_i = (int(i_omega_s[lo]), int(j_dty_s[lo]))
            slc    = slice(int(lo), int(lo + cnt))
            Itot   = self.columnfile.sum_intensity[slc].sum()
    
            # Direct dictionary assignment using the tuple key
            self.frames_LUT[om_y_i] = BinEntry(slice=slc, npeaks=int(cnt), sumI=Itot)


        # 3. Build scans_LUT: Mapping j -> BinEntry
        _, first_idx_dty, counts_dty = np.unique(j_dty_s, return_index=True, return_counts=True)

        self.scans_LUT = {}
        for lo, cnt in zip(first_idx_dty, counts_dty):
            y_i  = int(j_dty_s[lo])
            slc  = slice(int(lo), int(lo + cnt))
            Itot = self.columnfile.sum_intensity[slc].sum()
    
            # Direct dictionary assignment using the integer key
            self.scans_LUT[y_i] = BinEntry(slice=slc, npeaks=int(cnt), sumI=Itot)

        # 4. Statistics
        n_frames_total = len(oc) * len(yc)
        n_frames_empty = n_frames_total - len(self.frames_LUT)
        n_scans_total  = len(yc)
        n_scans_empty  = n_scans_total - len(self.scans_LUT)
        
        logger.info(
            "[sort_by_sinogram]\n%d peaks sorted; "
            "%d non-empty frames indexed (%d frames have no peaks).",
            self.columnfile.nrows, len(self.frames_LUT), n_frames_empty)
        logger.info(
            "%d non-empty scans indexed (%d scans have no peaks).",
            len(self.scans_LUT), n_scans_empty)


    def find_valid_subsets(self):
        """
        compute boolean arrays marking paired subsets where both members contain
        at least one peak in corresponding lookup table (frames_LUT or scans_LUT).
        Sets new attributes self.valid_scans_subsets / self.valid_frames_subsets
        """
        for subset_type in ['scans', 'frames']:
            # Check if LUT is defined in PeakSubset
            LUT_name = str(subset_type + '_LUT')
            LUT = getattr(self, LUT_name, None)
            if LUT is None:
                raise RuntimeError(
                    "Lookup table {LUT_name} not found in PeakSubset. "
                    "Run sort_by_sinogram() first.".format(LUT_name=LUT_name))
             # Get corresponding list of subsets
            if LUT_name == 'frames_LUT':
                pairs = self.frames_subsets
            elif LUT_name == 'scans_LUT':
                pairs = self.scans_subsets
            else:
                raise ValueError("Unknown LUT: {LUT}. Expected 'frames_LUT' or 'scans_LUT'.".format(LUT=LUT_name))
            
            valid_keys = set(LUT.keys())
            valid = np.zeros(len(pairs), dtype=bool)

            if 'frames' in LUT_name:
                for i, p in enumerate(pairs):
                    if (p.omi_A, p.yi_A) in valid_keys and (p.omi_B, p.yi_B) in valid_keys:
                        valid[i] = True
            else:
                for i, p in enumerate(pairs):
                    if p.yi_hi in valid_keys and p.yi_lo in valid_keys:
                        valid[i] = True

            nv = valid.sum()
            logger.info(
                "[find_valid_subsets]  %d pairs checked in %s\n"
                "  found %d subset pairs where both members have peaks "
                "(%.2f%% valid).",
                len(valid), LUT_name, nv, nv / len(valid) * 100)

            setattr(self, 'valid_'+subset_type+'_subsets', valid)

    def select_pair(self, pair_id: int, subset_type = 'scans', return_as = 'idx'):
        """
        select a pair of subsets from one of the pairs list"""
        
        LUT        = getattr(self, subset_type+'_LUT', None)
        pairs_list = getattr(self, subset_type+'_subsets', None)
        cf = self.columnfile

        if LUT is None:
            raise RuntimeError(
            "Lookup table not found. run sort_by_sinogram(), or check input name.")

        # select pair in pairs list
        pair = pairs_list[pair_id]
        # select peaks indices in lookup table
        if subset_type == 'frames':
            idx1 = (pair.omi_A, pair.yi_A)   
            idx2 = (pair.omi_B, pair.yi_B)   
        else :
            idx1 = pair.yi_hi   
            idx2 = pair.yi_lo   
        
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

            
    def check_dty_symmetry(self, saveplot=False):
        """
        Plot N_peaks and total_intensity vs. abs(dty) for mirror dty_scans
        (y0+Δy, y0-Δy) and check their correlation.

        If alignment is correct, values should match between mirror scans.
        Significant mismatch suggests an incorrect y0, sample movement, or
        a beam issue during scanning.
        """

        if not hasattr(self, 'scans_subsets'):
            self.get_scan_subsets()
        if not hasattr(self, 'scans_LUT'):
            self.sort_by_sinogram()
        if 'valid_pairs_subsets' not in self.__dict__.keys():
            valid = self.valid_scans_subsets
        
        # paired dty scans on the hi and lo side
        scan_id_hi = np.array( [p.yi_hi for p in self.scans_subsets] )
        scan_id_lo = np.array( [p.yi_lo for p in self.scans_subsets] )

        # x-axis: dty position (flipped for lo side)
        y0 = self.dataset.y0
        x_hi = self.ybincens[scan_id_hi] - y0
        x_lo = y0 - self.ybincens[scan_id_lo]

        # y-axis: n_peaks and sum_intensity, NaN for invalid pairs
        npks_hi = np.array([self.scans_LUT[yi].npeaks if v else np.nan for yi,v in zip(scan_id_hi, valid)])
        npks_lo = np.array([self.scans_LUT[yi].npeaks if v else np.nan for yi,v in zip(scan_id_lo, valid)])
        sumI_hi = np.array([self.scans_LUT[yi].sumI if v else np.nan for yi,v in zip(scan_id_hi, valid)])
        sumI_lo = np.array([self.scans_LUT[yi].sumI if v else np.nan for yi,v in zip(scan_id_lo, valid)])

        log_sumI_hi = np.log10(sumI_hi)
        log_sumI_lo = np.log10(sumI_lo)
        ratio_npks  = np.array(npks_hi) / np.array(npks_lo)
        ratio_sumI  = log_sumI_hi / log_sumI_lo
        # correlation coeff between hi and lo sides
        corr_npks = np.corrcoef(np.array(npks_hi)[valid],
                                np.array(npks_lo)[valid])[0, 1]
        corr_sumI = np.corrcoef(log_sumI_hi[valid],
                                log_sumI_lo[valid])[0, 1]

        # plot figure
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
        axes = axes.flatten()

        # Panel 1: n_peaks per scan
        axes[0].plot(x_hi[1:], npks_hi[1:], '.-', label='hi-side (dty - y0 > 0)')
        axes[0].plot(x_lo[1:], npks_lo[1:], '.-', label='lo-side (dty - y0 < 0)')
        axes[0].set_ylabel('npks per scan')
        axes[0].set_title('Number of peaks per scan | Corr = {:.3f}'.format(corr_npks))
        axes[0].legend()
        # Panel 2: log(total intensity) per scan
        axes[1].plot(x_hi[1:], log_sumI_hi[1:], '.-')
        axes[1].plot(x_lo,     log_sumI_lo,      '.-')
        axes[1].set_title('Total intensity per scan | Corr = {:.3f}'.format(corr_sumI))
        axes[1].set_ylabel('log10(sumI) per scan')
        # Panel 3: ratio of n_peaks (hi/lo)
        axes[2].plot(x_hi[1:], ratio_npks[1:], 'k.-')
        axes[2].axhline(1, color='gray', linestyle='--')
        axes[2].set_ylabel('Ratio npks (hi/lo)')
        axes[2].set_xlabel('|dty - y0| (flipped for lo side)')
        # Panel 4: ratio of log(sumI) (hi/lo)
        axes[3].plot(x_hi[1:], ratio_sumI[1:], 'k.-')
        axes[3].axhline(1, color='gray', linestyle='--')
        axes[3].set_ylabel('Ratio log10(sumI) (hi/lo)')
        axes[3].set_xlabel('|dty - y0| (flipped for lo side)')

        fig.suptitle('dty alignment - ' + str(self.dataset.dsname), fontsize=14)
        plt.tight_layout()

        if saveplot:
            fname = os.path.join(
                self.dataset.analysispath, self.dataset.dsname + '_dty_alignment.svg')
            fig.savefig(fname, format='svg')
        return fig



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

    def set_peak_subsets(self):
        self.logger.info('--- SET PEAKSUBSETS ---')
        self.PeakSubsets = PeakSubsets(self.cf, self.ds)
        self.logger.info('Sorting columnfile in sinogram order. May take a bit of time...')
        self.PeakSubsets.sort_by_sinogram()
        self.PeakSubsets.find_valid_subsets()

    def reset_outputs(self, pair_id=None):
        self.outputs = []
        # reset the eta/omega_pair_id column in cf
        if pair_id is not None:
            if pair_id in self.cf.titles:
                initvals =  np.full(self.cf.nrows, -1, dtype=int)
                self.cf.addcolumn(initvals, pair_id)               

    # ------------------------------------------------------------------
    #  High-level API — Friedel pairing functions
    # ------------------------------------------------------------------
    def match_friedel_pairs(self, pair_type = 'omega', drop_unpaired=False, filter_mode = 'strict', doplot=True):
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
            idx1 = np.nonzero(self.cf.omega%360 > 180)[0]
            idx2 = np.nonzero(self.cf.omega%360 < 180)[0]

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
                         weights     = weights_eff)
        
        self.merge_outputs(pair_type     = pair_type,
                           drop_broken   = True,
                           drop_unpaired = drop_unpaired)
       
        if doplot:
            _ = plot_pair_distances(self.cf, pair_type = pair_type, bins=50, log_scale=False)
        return self.cf
        

    def match_omega_pairs_by_chunks(self, chunk_type='scans', 
                                    filter_mode='strict',
                                    drop_unpaired = False,
                                    n_workers=1,
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
        chunk_type    : str 'frames', 'scans' (default: 'scans')
        drop_unpaired : bool (default False)
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

        # ── 1. Checks + init ───────────────────────────────────────────────────
        self.logger.info('--- INITIALIZATION ---')
        self.reset_outputs()

        pairs_list, valid_pairs, LUT = self.checks_before_pairing(chunk_type)
        Psub  = self.PeakSubsets
        n_valid = valid_pairs.sum()

        dist_max = _physical_to_normalised_cutoff(
            self.tol_gv, self.tol_eta, self.tol_logI, self.weights)
        
        self.logger.info('  %d pairs of subsets, %d valid', len(pairs_list), n_valid)

        # ── 2. Search-space calibration ──────────────────────────────────────
        self.logger.info('--- SEARCH SPACE CALIBRATION [pilot pairing] ---')

        # pick a pair near the middle — avoids edge scans with fewer peaks
        s0 = len(pairs_list) // 2
        while not valid_pairs[s0]:
            s0 += 1
        idx1_pilot, idx2_pilot = Psub.select_pair(s0, chunk_type, return_as='idx')

        # find weight factors for optimal rescaling
        _, rescaling_wts = self.estimate_search_scales(
            idx1_pilot, idx2_pilot, 0.5 * dist_max, pair_type='omega')
        
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
        
        # Build worker argument list 
        worker_args = []
        
        for pid in range(len(pairs_list)):
            if not valid_pairs[pid]:
                continue
            idx1, idx2 = self.PeakSubsets.select_pair(pid, subset_type = chunk_type, return_as='idx')

            worker_args.append((
                pid, idx1, idx2,
                'omega',
                filter_mode,
                weights_eff,
                self.tol_gv, self.tol_eta, self.tol_logI,
                self.n_steps,
            ))

        # Resolve number of workers
        n_cores    = len(os.sched_getaffinity(os.getpid())) - 1
        if n_workers == -1:
            n_workers = n_cores
        else:
            n_workers = max(1, min(n_workers, n_cores))
        
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
        # unique pair identifier in omega_pair_id. -1 for unpaired peaks
        self.merge_outputs(pair_type='omega',
                           drop_broken=True,
                           drop_unpaired=drop_unpaired)

        if doplot:
            _ = plot_pair_distances(self.cf, pair_type = 'omega', bins=50, log_scale=False)
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
    def checks_before_pairing(self, chunk_type='scans'):
        """ Initialization and checks for match_omega_pairs_by_chunks(): """
        
        assert chunk_type in ['frames','scans'], "chunk_type must be 'frames' or 'scans' " 
        self.logger.info('Chunk type for pairing: %s', chunk_type)
              
        if self.PeakSubsets is None:
            self.set_peak_subsets()

        Psub = self.PeakSubsets
        LUT         = getattr(Psub, chunk_type+'_LUT', None)
        pairs_list  = getattr(Psub, chunk_type+'_subsets', None)
        valid_pairs = getattr(Psub, 'valid_'+chunk_type+'_subsets', None)
        if LUT is None:
            raise RuntimeError(
                'Lookup table not found in PeakSubsets.'
                'Run self.PeakSubsets.sort_by_sinogram(), or check input name.')
        if pairs_list is None:
            raise RuntimeError(
                'Pairs list not found in PeakSubsets')
        if valid_pairs is None:
            raise RuntimeError(
                'No list of valid subset in PeakSubsets.'
                'Run self.PeakSubsets.find_valid_subsets()') 
        if (chunk_type=='frames' and self.cf.sortedby != 'dty_omega') or (chunk_type=='scans' and not self.cf.sortedby.startswith('dty')):
            raise RuntimeError(
                'columnfile is not properly sorted. Run PeakSubsets.sort_by_sinogram().')
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
            self.logger.debug('  pilot search: cutoff=%.4f  pairs found: %d',
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
        d_gv_med  = np.median(d_gv)
        d_eta_med = np.median(d_eta)
        d_logI_med = np.median(d_logI)

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

    def run_pairing(self, idx1, idx2, pair_type='omega',
                     filter_mode = 'strict', weights=None):
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
                                verbose     = verbose)
        
        self.outputs.append(results)

        
    def merge_outputs(self, cf_target=None, pair_type='omega',
                      drop_broken=False, drop_unpaired = False):
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
            # compute the id offset for each chunk: cumulative sum of chunk sizes
            chunk_sizes = np.array([o['pair_id'].max() + 1 for o in valid_outputs])
            offsets     = np.concatenate([[0], np.cumsum(chunk_sizes)[:-1]])

            # shift each chunk's pair_ids to the global range and concatenate
            all_idx     = np.concatenate([o['idx'] for o in valid_outputs])
            all_pair_id = np.concatenate([o['pair_id'] + off
                                  for o, off in zip(valid_outputs, offsets)])

            label_col[all_idx] = all_pair_id
            next_id = int(chunk_sizes.sum())

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

def _compute_dist_matrix(space_ref, space_partner, dist_cutoff):
    """
    Build sparse distance matrix between reference and partner search spaces.
    """
    tree_ref     = scipy.spatial.cKDTree(space_ref)
    tree_partner = scipy.spatial.cKDTree(space_partner)
    dij = csr_matrix(tree_ref.sparse_distance_matrix(tree_partner, dist_cutoff))
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

    for i, dist_cutoff in enumerate(dist_steps):
        if free1.sum() < 1 or free2.sum() < 1:
            _print('  No free peaks remaining — stopping early.')
            break
        # indices into global cf
        cur_idx1 = idx1[free1]
        cur_idx2 = idx2[free2]
            
        space_1, _ = _search_space(cf, cur_idx1, _w)
        space_2, _ = _search_space(cf, cur_idx2, _w, flip=pair_type)

        # compute distance matrix and clean it
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
                
        # strict filter_mode: apply tolerance filters within loop
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

        # map identified pairs into global cf
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
    if filter_mode == 'relaxed':
        tol_mask = _filter_pairs_physical(
                cf, idx1[paired_mask], partner1[paired_mask],
                pair_type, tol_gv, tol_eta, tol_logI)
        
        wrong_om = _wrong_omegas(cf, idx1[paired_mask], partner1[paired_mask],
                                 tol=1e-4)
        tol_mask &= ~wrong_om
        
        n_pairs = tol_mask.sum()
        idx1_paired = idx1[paired_mask][tol_mask]
        idx2_paired = partner1[paired_mask][tol_mask]
            
        _print('Applying tolerance filters')
        _print('Total pairs after filtering : {} / {}'.format(n_pairs, n_total))
        _print(_progress_bar(n_pairs, n_total))
        
    else:
        idx1_paired = idx1[paired_mask]
        idx2_paired = partner1[paired_mask]

    if n_pairs < 1:
        _print('  WARNING: No pairs found — returning empty output.')
        result = {'idx'      : np.array([], dtype=int),
                  'pair_id'  : np.array([], dtype=int)}
        return result

    local_ids   = np.arange(n_pairs)
    # interleave [ref0, partner0, ref1, partner1, ...]
    idx_out    = np.empty(2 * n_pairs, dtype=int)
    pid_out    = np.empty(2 * n_pairs, dtype=int)
    idx_out[0::2]    = idx1_paired;  idx_out[1::2]    = idx2_paired
    pid_out[0::2]    = local_ids;    pid_out[1::2]    = local_ids

    _print('\n{}\n'.format('-' * 50))

    result = {'idx'    : idx_out,
              'pair_id': pid_out}
    return result

def _worker_run_pairing(args):
    """ Module-level worker function for parallelised chunk pairing """
    pid, idx1, idx2, pair_type, filter_mode, weights, \
        tol_gv, tol_eta, tol_logI, n_steps = args

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
        verbose=False)
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
    pair_id_col = cf.getcolumn(pair_type + '_pair_id')
    paired = np.nonzero(pair_id_col > -1)[0]
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


def locate_eta_pairs(cf, pairs, y0=0.):
    """
    Fit the centre of mass position of eta-pairs and write results
    back into cf as new columns 'xs' and 'ys'.

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

    xs_pairs, ys_pairs = np.linalg.solve(R, y).T

    # ── write into cf. Overwrite pre-existing xs,ys
    xs = np.full(cf.nrows, np.nan)
    ys = np.full(cf.nrows, np.nan)
    xs[i1] = xs_pairs;  xs[i2] = xs_pairs
    ys[i1] = ys_pairs;  ys[i2] = ys_pairs

    for name, col in (('xs', xs), ('ys', ys)):
        cf.addcolumn(col, name)

    return xs, ys


def locate_omega_pairs(cf, pairs, ds=None, y0=0):
    """
    Fit the centre of mass position of omega-pairs and write results
    back into cf as new columns 'xs' and 'ys'.

    For omega-pairs the linear system used in 'locate_eta_pairs' is near-singular
    (paired peaks are ~180 degrees apart in omega), so the position is recovered from
    two-theta and dty:

    dx  : x-shift from rotation centre, derived from tth asymmetry
    ylab: y-shift in lab frame, taken as dty - y0
    Then rotate (dx, ylab) back into sample coordinates using omega.

    xlab, ylab and omega are averaged to make sure both peaks in eahc pair have the same (xs,ys) coordinates  
    """
    i1, i2 = pairs
    L = cf.parameters.get('distance')
    if ds is not None and 'frelon' in ds.detector:
        L /= 1000  # Frelon dty in mm. divide by 1000 to match units

    xlab = np.empty(cf.nrows, dtype=float)
    ylab = np.empty(cf.nrows, dtype=float)
    xs   = np.empty(cf.nrows, dtype=float)
    ys   = np.empty(cf.nrows, dtype=float)

    # shift from rot center along x-axis: same absolute distance dx but reversed sign 
    tantth1  = np.tan(np.radians(cf.tth[i1]))
    tantth2  = np.tan(np.radians(cf.tth[i2]))
    dx       = L * (tantth1-tantth2) / (tantth1+tantth2)  
    xlab[i1] = dx
    xlab[i2] = -dx
    
    # shift along y-axis given by motor position : use the mean as the best estimate of the true ylab
    dy   = (cf.dty[i1] - cf.dty[i2]) / 2 
    ylab[i1] = y0 + dy
    ylab[i2] = y0 - dy

    # rotate peak coords in sample reference frame (xs,ys).
    # use averaged omega values to make sure paired peaks are relocated to the same spot 
    xy_lab  = np.transpose((xlab,ylab))
    xy_s    = np.transpose((xs,ys))
    
    r = np.radians(cf.omega)
    co,so = np.cos(r), np.sin(r)
    co_av = (co[i1]-co[i2]) / 2 
    so_av = (so[i1]-so[i2]) / 2    

    R = np.transpose( ((co_av , so_av),
                       (-so_av, co_av)), axes=(2,0,1))

    xy_s[i1] = np.einsum('ijk,ik->ij',R,xy_lab[i1])
    xy_s[i2] = np.einsum('ijk,ik->ij',-R,xy_lab[i2])
    xs = xy_s[:,0]
    ys = xy_s[:,1]
    
    for name, col in (('xs', xs), ('ys', ys)):
        cf.addcolumn(col, name)
    return xs, ys


