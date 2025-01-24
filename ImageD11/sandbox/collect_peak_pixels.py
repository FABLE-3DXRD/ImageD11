
from __future__ import print_function, division

import multiprocessing, time
import numpy as np
import scipy.sparse, h5py
import fabio
from ImageD11 import cImageD11
#from .
import sparse_peak_searching
import sys, os, glob

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
    c = cuts[-1].copy()
    return c


def fix_frelon_lines(img, cut=None, block=256):
    """
    Reduce readout artifacts from frelon fast readout
    """
    if cut is None:
        cut = sigma_cut( img, nsigma=2 )
    assert img.shape[1]%block == 0
    blocked = np.reshape( np.where(img<cut, img, 0),
                          (img.shape[0]*img.shape[1]//block, block) )
    drift = blocked.mean( axis=1 )
    corrected = np.reshape(img, blocked.shape) - drift[:, np.newaxis]
    return np.reshape( corrected, img.shape )


def binary_clean( mask ):
    """
    Remove the one pixel stuff
    for pixel at i,j:
        if p[i,j] and neighbor i+/-1, j+-1
           keep
    TODO : faster in C ?
    """
    m = (mask > 0).astype(int)
    np.add( m[1:], mask[:-1], m[1:] )
    np.add( m[:-1], mask[1:], m[:-1] )
    np.add( m[:,1:], mask[:,:-1], m[:,1:] )
    np.add( m[:,:-1], mask[:,1:], m[:,:-1] )
    return (mask > 0) & (m > 1)


def to_coo( data, mask, name ):
    """
    Convert an image (2D, probably float) and mask (0 versus !0)
    to a sparse array (could be more efficient...)
    """
    m = mask.astype(bool)
    j = np.outer(np.ones(  data.shape[0], dtype=np.uint16),
                 np.arange(data.shape[1], dtype=np.uint16)  )
    ind0, ind1 = j.T[m], j[m]
    # round to the nearest 0.125
    values = np.round( data[m] * 8 ) / 8
    obj = sparse_peak_searching.sparse_frame( ind0, ind1, values,
                                              shape = data.shape,
                                              name = name)
    return obj

    
def segment( fname, bg ):
    """ Does a segmentation of a single file/frame """
    imo = fabio.open(fname)
    imb = imo.data.astype( np.float32 )  - bg
    #print (type(imb))
    imc = fix_frelon_lines( imb, block=128 )
    s = sigma_cut(imc)
    m = imc > s
    cleanmask = m.astype(np.int8)
    cImageD11.clean_mask( m.astype(np.int8), cleanmask )
    spar = to_coo( imc, cleanmask, fname )
    spar.sigma_cut = s
    spar.con_labels()
    spar.max_labels()
    print(fname,spar.nmlabel)
    return spar

def write_frames_to_hdf( frames,
                         hdfname,
                         hdffolder,
                         opts ):
    h = h5py.File( hdfname, "w" )
    for f in frames:
        f.tohdf( h )
    h.close()



def segment_folder( folder , nimages=900
                    , extn=".edf.gz" ):
    print("Starting", folder)
    bg = fabio.open("pks/bkg.edf").data
    start = time.time()
    fnames = [ os.path.join( folder, folder + "%04d%s"%(i, extn)) for i
               in range(nimages) ]    
    spars = [segment(fname, bg) for fname in fnames ]
    return (spars, folder)


def dofolder(f):
    start = time.time()
    spars, folder = segment_folder( f )
    print("Segmented", folder, time.time()-start)
    hname = "hdfs3/%s.hdf"%(folder)
    opts =  { "compression" : "gzip",
              "compression_opts" : 4 }
    sparse_peak_searching.sparse_frames_to_hdf( hname, f, spars )
    #write_frames_to_hdf( spars, hname, folder, opts )
    print("Wrote", folder, time.time()-start)

if __name__=="__main__":
    folders =[ f for f in sorted( glob.glob("Au6_s0*") ) if os.path.isdir(f) ]
    if not os.path.exists("hdfs3"):
        os.mkdir("hdfs3")
    if 1:
        import multiprocessing
        p = multiprocessing.Pool( processes = multiprocessing.cpu_count() )
        for x in p.map(dofolder, folders):
            pass
    else:
        for x in folders:
            dofolder(x)





    



