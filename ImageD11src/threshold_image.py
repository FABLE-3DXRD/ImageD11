
import numpy as np
from ImageD11 import cImageD11
# Estimate a noise level and determine a threshold

def sigma_cut( img, nsigma=3. ):
    """
    Assume the image has been background corrected and 
    follows a Gaussian noise distribution (e.g. CCD)
    """
    means=[img.mean(),]
    stds =[img.std(),]
    cuts = [means[-1] + stds[-1]*nsigma,]
    for i in range(3):
        imc = img[img < cuts[-1]]
        means.append( imc.mean() )
        stds.append(  imc.std() )
        cuts.append( means[-1] + stds[-1]*nsigma )
    # return means, stds, cuts
    return cuts[-1]


def fix_frelon_lines(img, cut=None, block=256):
    """
    Reduce readout artifacts from frelon fast readout
    """
    if cut is None:
        cut = sigma_cut( img, nsigma=2 )
    assert img.shape[1]%block == 0
    blocked = np.reshape( np.where(img<cut, img, 0),
                          (img.shape[0]*img.shape[1]/block, block) )
    drift = blocked.mean( axis=1 )
    corrected = np.reshape(img, blocked.shape) - drift[:, np.newaxis]
    return np.reshape( corrected, img.shape )


def binary_clean( mask ):
    """
    Remove the one pixel stuff
    """
    ret = np.zeros(mask.shape, np.int8)
    npx = cImageD11.clean_mask( mask.astype( np.int8 ), ret )
    return npx, ret
    
