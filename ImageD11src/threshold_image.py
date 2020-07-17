
from __future__ import print_function, division

import numpy as np
from ImageD11 import cImageD11
# Estimate a noise level and determine a threshold

def sigma_cut( img, nsigma=3., cycles=3 ):
    """
    Assume the image has been background corrected and 
    follows a Gaussian noise distribution (e.g. CCD)
    """
    imc = np.ascontiguousarray( img, dtype=np.float32 )
    mean, variance = cImageD11.array_mean_var_cut( imgc.ravel(), n=cycles, cut=nsigma, verbose=0 )
    return cuts[-1], means[-1]


def fix_frelon_lines(img, bkg = None, block=256, nsigma=2.):
    """
    Reduce readout artifacts from frelon fast readout
    """
    cImageD11.frelon_lines
    cut, val = sigma_cut( img, nsigma=2 )
    assert img.shape[1]%block == 0
    blocked = np.reshape( img,
                          (img.shape[0]*img.shape[1]//block, block) )
    bgm = blocked<cut
    drift = (blocked*bgm).sum( axis=1 ) / bgm.sum( axis=1 )
    corrected = np.reshape(img, blocked.shape) - drift[:, np.newaxis]
    return np.reshape( corrected, img.shape )


def binary_clean( mask, out = None ):
    """
    Remove the one pixel stuff
    """
    if out is None:
        out = np.zeros(mask.shape, np.int8)
    npx = cImageD11.clean_mask( mask.astype( np.int8 ), out )
    return npx, out
    

def threshold_image(im, nsigma=2, fix_lines=False, mskmem=None, block=256,
                   verbose=0):
    """
    Process one single image (numpy array)
    """
    if fix_lines:
        im = fix_frelon_lines(im, block=block, nsigma=nsigma)
        if verbose: print("fixed lines")
    cut, val = sigma_cut(im, nsigma)
    if verbose:
        print("cut,val: ",cut, val)
    mask = im > cut
    if verbose:
        print("mask raw:",mask.sum())
    npx, clean_mask = binary_clean( mask, out=mskmem )
    i = np.empty( npx, np.uint16)
    j = np.empty( npx, np.uint16)
    tmp = np.empty( im.shape[0], np.int32 )
    ok = cImageD11.mask_to_coo( clean_mask.astype(np.int8), i, j, tmp )
    if verbose:
        print("mask clean",clean_mask.sum())
    assert ok == 0
    print(i.shape,j.shape,clean_mask.sum())
    frm = sparse_frame( i, j, im.shape, itype=np.uint16, 
                       pixels = {'data': im[clean_mask.astype(np.bool)]} )
    return frm

def threshold_series( imseries, nsigma=2, fix_lines=False):
    mskmem = None
    for im in imseries:
        if mskmem is None:
            mskmem = np.zeros( im.shape, np.bool )
        yield threshold_image( im, nsigma=nsigma, 
                              fix_lines = fix_lines, mskmem=mskmem )
        

    
    
        

