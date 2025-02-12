# coding: utf-8

from __future__ import print_function, division
from unitcell import unitcell
from columnfile import columnfile
"""
Various utility functions for selecting peaks within columnfiles
"""

import numpy as np
from ImageD11.cImageD11 import array_bin
import ImageD11.unitcell
import ImageD11.refinegrains

def assign_ring_histo(cf, dsmax, hbins, cell):
    """
    Assigns rings using a histogram of peaks in dstar
    Uses the local maxima
    
    cf = columnfile
    dsmax = max dstar to consider
    hbins = number of bins from 0 to dsmax
    cell = ImageD11.unitcell object
    
    adds columns "ring" and "dstarbin"
    
    returns (bincen, histogram, ring_assign )
    """
    cell.makerings( dsmax )
    # This is something like transform.tth_histo
    bedge = np.linspace(0,dsmax,hbins+1) # hbins gaps
    bcen  = (bedge[1:] + bedge[:-1])/2
    istep = hbins/dsmax
    # Top bin allowed is hbins-1
    # dsb = np.floor( cf.ds * istep ).astype(int).clip(0,hbins-1)
    dsb = array_bin( cf.ds, istep, hbins )
    # This is histogram
    h = np.bincount( dsb, weights = cf.sum_intensity, minlength=hbins )
    # Local maxima and minima
    maxes = (h[1:-1]>h[2:])&(h[1:-1]>=h[:-2])
    mins =  (h[1:-1]<=h[2:])&(h[1:-1]<h[:-2])
    bc = bcen[1:-1]
    maxpos = bc[maxes]
    minpos = bc[mins]
    # Mask of use vs do-not-use for the histogram
    ra = np.full(h.shape, -1, dtype=int)
    for i,d in enumerate(cell.ringds):
        idx = np.argmin( abs(maxpos - d) )
        if minpos[idx] > maxpos[idx]:
            lo, hi = ( np.searchsorted( bcen, minpos[idx-1], side='left'), 
                      np.searchsorted(bcen, minpos[idx], side='right') )
            ra[lo:hi]=i
        else:
            lo, hi = ( np.searchsorted( bcen, minpos[idx], side='left'), 
                      np.searchsorted(bcen, minpos[idx+1], side='right' ) )
            ra[lo:hi]=i
    cf.addcolumn( ra[dsb], 'ring' )
    cf.addcolumn( dsb, 'dstarbin' )
    return (bcen,h,ra)


def rings_mask(cf, dstol, dsmax, cell=None):
    """Create a boolean mask for a columnfile
       Based on peak distance in dstar (within dstol) to HKL ring values
       And a maximum dstar value (dsmax)
       Optionally provide an ImageD11.unitcell instance
       Or it will generate one from the parameters in the columnfile"""
    if cell is None:
        cell = ImageD11.unitcell.unitcell_from_parameters(cf.parameters)
    cell.makerings(dsmax)
    m = np.zeros(cf.nrows, bool)
    for v in cell.ringds:
        if v < dsmax:
            m |= (abs(cf.ds - v) < dstol)
    return m

def filter_peaks_by_phase(
        cf:columnfile, 
        dstol:float, 
        dsmax:float, 
        cell:unitcell=None
):
    """
    It filter the given columnfile 
    by the provided unicell, with dstar tolerance and dstar max value
    returns filtered columnfile
    """
    mask = rings_mask(cf=cf, dstol=dstol, dsmax=dsmax, cell=cell)
    cf = cf.copyrows(mask)
    return cf

def sorted_peak_intensity_mask(colf, uself=True, frac=0.995, B=0.2, doplot=None):
    """
    Create a boolean mask for a columnfile based on peaks sorted by fractional intensity.

    Args:
        colf: Input column file object.
        uself: Apply Lorentz factor correction (default: True).
        frac: Fraction of peaks to keep based on intensity.
        B: Thermal factor for intensity correction.
        doplot: Optional plotting argument.

    Returns:
        mask: Boolean mask for selected peaks.
    """
    mask, cums = sorted_peak_intensity_mask_and_cumsum(colf = colf, uself=uself, frac=frac, B=B,)    
    if doplot is not None:
        plot_sorted_peak_intensity(cums, mask, frac, doplot)

    return mask

def plot_sorted_peak_intensity(cums, mask, frac, doplot):
    """
    Helper function to plot sorted peak intensity data.

    Args:
        cums: Cumulative sums of sorted intensities.
        mask: Boolean mask for selected peaks.
        frac: Fraction of peaks to keep.
        doplot: Plot zooming range.
    """
    from matplotlib import pyplot as plt

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    axs[0].plot(cums / cums[-1], ',')
    axs[0].set(xlabel='npks', ylabel='fractional intensity')
    axs[0].plot([mask.sum(), ], [frac, ], "o")

    axs[1].plot(cums / cums[-1], ',')
    axs[1].set(xlabel='npks logscale', ylabel='fractional intensity',
               xscale='log', ylim=(doplot, 1.),
               xlim=(np.searchsorted(cums, doplot), len(cums)))
    axs[1].plot([mask.sum(), ], [frac, ], "o")

    plt.show()


def sorted_peak_intensity_mask_and_cumsum(colf, uself=True, frac=0.995, B=0.2,):
    """
    Create a boolean mask for a columnfile based on peaks sorted by fractional intensity.

    Args:
        colf: Input column file object.
        uself: Apply Lorentz factor correction (default: True).
        frac: Fraction of peaks to keep based on intensity.
        B: Thermal factor for intensity correction.

    Returns:
        mask: Boolean mask for selected peaks.
        cums: Cumulative sums of the (sorted) peak intensities
    """
    # Correct intensities for structure factor (decreases with 2theta)
    cor_intensity = colf.sum_intensity * (np.exp(colf.ds * colf.ds * B))
    if uself:
        lf = ImageD11.refinegrains.lf(colf.tth, colf.eta)
        cor_intensity *= lf
    order = np.argsort(cor_intensity)[::-1]  # Sort peaks by intensity
    sortedpks = cor_intensity[order]
    cums = np.cumsum(sortedpks)
    cums /= cums[-1]

    # TODO
    # anything above this should be calculated once and added to the CF
    # check if the column exists before calculating
    # slider for frac?
    # all above only needs calculating once
    enough = np.searchsorted(cums, frac)
    # Aim is to select the strongest peaks for indexing.
    cutoff = sortedpks[enough]
    mask = cor_intensity > cutoff

    return mask, cums


def select_ring_peaks_by_intensity(cf, dstol=0.005, dsmax=None, frac=0.99, B=0.2, doplot=None,):
    """
    Select peaks based on ring intensity.

    Args:
        cf: Input columnfile + unit cell parameters.
        dstol: Difference in d* for assigning peaks to rings.
        dsmax: High angle cutoff for removing peaks.
        frac: Fractional normalised intensity to keep (removes weak peaks)
        B: Thermal factor to downweight low angle vs high angle peaks for normalised intensity
        doplot: Whether to draw a plot. (float number, range for zoomed plot)
    Returns:
        cfc: Columnfile with selected peaks.
    """
    if dsmax is None:
        dsmax = cf.ds.max()
        cfd = cf
    else:
        cfd = cf.copyrows( cf.ds <= dsmax )
        
    cfc = filter_peaks_by_phase(cf=cfd, dstol=dstol, dsmax=dsmax)
    ms = sorted_peak_intensity_mask(cfc, frac=frac, B=B, doplot=doplot)
    cfc.filter(ms)
    return cfc

