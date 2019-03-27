
from __future__ import print_function, division

import multiprocessing, time
import numpy as np
import scipy.sparse, h5py
import fabio
from ImageD11 import cImageD11

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
                          (img.shape[0]*img.shape[1]/block, block) )
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


def to_scipy_csr( data, mask ):
    """
    Convert an image (2D, probably float) and mask (0 versus !0)
    to a sparse csr array (could be more efficient...)
    """
    m = mask.astype(np.bool)
    j = np.outer(np.ones(  data.shape[0], dtype=np.int32),
                 np.arange(data.shape[1], dtype=np.int32)  )
    indices = j[m]
    indptr = np.concatenate( ((0,), np.cumsum( m.sum( axis=1 ))))
    values = data[m]
    obj = scipy.sparse.csr_matrix( ( values, indices, indptr ),
                                   shape = data.shape )
    return obj


def write_csr_to_hdf( vals, fname, pname, opts ):
    """
    vals = sparse values to be stored
    fname = hdf file
    pname = hdf path

    opts = compression options
         { shuffle : [True|False],
           compression : ["gzip" | "lzf" | ...  ],
           compression_opts : [ 0-9 | None] }

    Saves :
       data/indices/indptr from scipy.sparse.csr_matrix
       .attrs["h5sparse_shape"]=shape
       .attrs["h5sparse_format"]='csr'

    TODO?: labels for peak labelling

    See also:
       https://github.com/appier/h5sparse
       https://falexwolf.de/blog/171212_sparse_matrices_with_h5py/
       https://anndata.readthedocs.io/en/latest/
       ... Jerome Kieffer mentioned NXsparse for Nexus/SILX ...
    """
    h = h5py.File( fname )
    g = h.require_group( pname )
    g.require_dataset( "data",  
                       vals.data.shape,
                       vals.data.dtype,
                       data=vals.data,
                       **opts )
    g.require_dataset( "indices",
                       vals.indices.shape,
                       vals.indices.dtype,
                       data=vals.indices,                      
                       **opts )
    g.require_dataset( "indptr",
                       vals.indptr.shape,
                       vals.indptr.dtype,
                       data=vals.indptr,
                       **opts )
    g.attrs['h5sparse_format']="csr"
    g.attrs['h5sparse_shape']=vals.shape
    h.close()
    
                       
def read_csr_from_hdf( fname, pname ):
    """
    Reads a sparse matrix from a hdf file(fname)/group(pname)
    """
    h = h5py.File( fname )
    p = h[pname]
    shape = p.attrs['h5sparse_shape']
    assert p.attrs['h5sparse_format'] == 'csr'
    vals = p['data'][:]
    indices = p['indices'][:]
    indptr = p['indptr'][:]
    obj = scipy.sparse.csr_matrix( ( vals, indices, indptr ), shape=shape )
    return obj
    
    
def check_eq_csr( o1, o2 ):
    """ For debugging : compare two sparse images as equal or not """
    if o1.nnz != o2.nnz:
        return False
    if ~(o1.indptr == o2.indptr).all():
        return False
    if ~(o1.indices == o2.indices).all():
        return False
    return np.allclose(o1.data, o2.data)


def csr_join( args ):
    """ 
    Concatenate a bunch of csr together
    """
    shapes = [a.shape for a in args]
    s0 = np.sum( [ a.shape[0] for a in args] ) # nrows
    s1 = a[0].shape[1]                         # ncols
    for a in args:
        assert a.shape[1] == s1
    nnz = np.sum( [ a.nnz for a in args ] )
    ip  = np.zeros( s0+1, np.int32 )
    i0=1
    for a in args:
        ip[i0:i0+s1] = a.indptr[1:] + ip[i0-1]
        i0 += s1
    data    = np.concatenate( [ a.data    for a in args] )
    indices = np.concatenate( [ a.indices for a in args] )
    obj = scipy.sparse.csr_matrix( (data, indices, ip), shape = (s0, s1) )
    return obj
    
def segment( fname , bg):
    """ Does a segmentation of a single file/frame """
    imb = fabio.open(fname).data - bg
    imc=fix_frelon_lines( imb, block=128 )
    s=sigma_cut(imc)
    m = imc > s
    cleanmask=m.astype(np.int8)
    cImageD11.clean_mask( m.astype(np.int8), cleanmask )
    spar = to_scipy_csr( imc, cleanmask )
    return spar

def segment_folder( folder , nimages=900, extn=".edf.gz",
                    opts =  { "compression" : "lzf",
                              "shuffle" : True }):
    print("Starting", folder)
    start = time.time()
    fnames = [ os.path.join( folder, folder + "%04d%s"%(i, extn)) for i
               in range(nimages) ]
    bg = fabio.open("pks/bkg.edf").data
    spars = [segment(fname, bg) for fname in fnames ]
    print("Segmented", folder, time.time()-start)
    blob3d = csr_join( spars )
    hname = "hdfs/%s.hdf"%(folder)
    write_csr_to_hdf( blob3d, hname, folder, opts )
    print("Wrote", folder, time.time()-start)
    return hname






folders = sorted( glob.glob("Au6_s0*") )
import multiprocessing
p = multiprocessing.Pool( multiprocessing.cpu_count()//2 )
for x in p.imap( segment_folder, folders ):
    print("Done",x)

    



