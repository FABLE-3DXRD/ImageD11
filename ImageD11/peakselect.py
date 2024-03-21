# coding: utf-8

from __future__ import print_function, division

"""
Various utility functions for selecting peaks within columnfiles
"""

import numpy as np
from ImageD11.cImageD11 import array_bin

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
