

from __future__ import print_function
#import multiprocessing, os
#cpu = str(multiprocessing.cpu_count())-1
#os.environ['OPENBLAS_NUM_THREADS']=cpu
#os.environ['NUMBA_NUM_THREADS']=cpu
#os.environ['OMP_NUM_THREADS']=cpu

import timeit
from ImageD11 import cImageD11, sparseframe
import numpy as np, fabio, h5py, os
import pylab as pl
import numba
import scipy.ndimage

""" Generator function based image pipeline """

timer = timeit.default_timer
LASTTIME = timer()
def bench(name):
    def outer(func):
        def inner(*args, **kwds):
            global LASTTIME
            gen = func(*args, **kwds)
            for value in gen:
                end = timer()
                print(name,"%.2f ms"%(1e3*(end-LASTTIME)),end=" ")
                LASTTIME=end
                yield value
        return inner
    return outer
        


def filename_series( stem, first, last, fmt="%s%04d.edf"):
    """ Generates the series of filenames to process """
    i = first
    while i <= last:
        fname = fmt%(stem,i)
        print("\n%s"%(os.path.split(fname)[-1]),end=" ")
        yield fname
        i += 1

@bench("open+dark")
def open_sub_dark( names, darkfilename, DARK_N=1 ):
    """ Opens each file in names and subtracts the dark
    yields data with dark subtracted
    """
    drk = fabio.open(darkfilename).data.astype(np.float32)
    drk /= DARK_N
    for name in names:
        img = fabio.open(name)
        img.data = img.data.astype(np.float32) - drk
        yield img

@bench("mask")
def do_mask( imgs, maskfilename ):
    """ applies a mask to the data """
    msk = fabio.open(maskfilename).data.astype(bool)
    for img in imgs:
        img.data[msk] = 0
        yield img

                      
@numba.njit(fastmath=True,parallel=True)
def bgest( buffer, nsigma, out):
    """ guesses the background using a clipped sigma """
    for i in numba.prange(buffer.shape[1]): # slow image
        for j in range(buffer.shape[2]): # along lines
            # first pass with no mask
            # out[i,j] = np.median( buffer[:,i,j] )
            # continue
            s0 = buffer[0,i,j]
            sm = buffer[0,i,j]*0
            sm2 = buffer[0,i,j]*0
            npx = buffer.shape[0]
            for k in range(1,buffer.shape[0]): # one from each frame
                t = buffer[k,i,j]  - s0
                sm += t
                sm2 += t*t
            avg = sm/npx + s0
            if (sm2 - sm*sm/npx) > 0:
                std = np.sqrt( (sm2 - sm*sm/npx)/npx )
            else:
                std = 1.0
            for _ in range(2):
                npx = 0
                s0 = avg
                cut = std * nsigma
                sm *= 0
                sm2*= 0
                for k in range(buffer.shape[0]): # one from each frame
                    t = buffer[k,i,j]  - s0
                    if t < cut:
                        sm += t
                        sm2 += t*t
                        npx += 1
                avg = sm/npx + s0
                if (sm2 - sm*sm/npx) > 0:
                    std = np.sqrt( (sm2 - sm*sm/npx)/npx )
                else:
                    std = 1.0
            out[i,j] = avg

    
            
@bench('bgsub')
def buffered_bgsub( series, nsigma, bufsize=11):
    """
    Sets up a buffers of a series of frames
    """
    buffer = None
    bufo = []
    # init
    nread = 0
    # import pdb; pdb.set_trace()
    for i in range(bufsize):
        frame = next(series)
        if buffer is None:
            buffer = np.zeros( (bufsize, frame.data.shape[0], frame.data.shape[1]),
                               np.float32 )
            bg = np.zeros( frame.data.shape, np.float32 )
            var = np.zeros( frame.data.shape, np.float32 )
        buffer[ i ] = frame.data
        bufo.append(frame)
        nread += 1
    bgest( buffer, nsigma, bg )
    for i in range(bufsize//2+1):
        o = bufo[i%bufsize]
        o.bg = bg
        yield o
    i = i + 1
    for frame in series:
        # mask using old estimae
        bufo[nread%bufsize] = frame
        buffer[nread%bufsize] = frame.data
        bgest( buffer, nsigma, bg )
        o.bg = bg
        yield o
        nread += 1
        i += 1
    while i < nread:
        o = bufo[i%bufsize]
        o.bg = bg
        yield o
        i += 1


@bench('gauss filter')
def filter_gauss( series, sigma ):
    """ applies a 2D gaussian filter prior the thresholding """
    smoothed= None
    for img in series:
        signal = img.data - img.bg
        if smoothed is None:
            smoothed = np.empty_like(img.data)
        scipy.ndimage.gaussian_filter( signal,
                                       sigma,
                                       output = smoothed )
        img.signal= signal
        img.smoothed = smoothed
        yield img

@bench("threshold")
def apply_threshold( series, nsigma ):
    """ converts the data to sparse array """
    msk = None
    for img in series:
        if msk is None:
            msk = np.empty( img.data.shape, np.uint8)
            cln = np.empty( img.data.shape, np.uint8)
            tmp = np.empty( img.data.shape[0], 'i' )
            row = np.empty( img.data.size, np.uint16)
            col = np.empty( img.data.size, np.uint16)        
        a,s = cImageD11.array_mean_var_cut( img.smoothed, nsigma )
        nnz = cImageD11.make_clean_mask( img.smoothed, a+nsigma*s, msk, cln )
        cImageD11.mask_to_coo( cln, row[:nnz], col[:nnz], tmp )
        header = {"filename":img.filename}
        img.mask = cln
        intensity = sparseframe.Container( img.signal[cln>0], **header )
        img.sparseframe = sparseframe.sparse_frame( row[:nnz], col[:nnz],
                                                    img.data.shape,
                                                    itype=np.uint16,
                                                    pixels = {"intensity":
                                                              intensity } )
        yield img
        

def grab_folder():
    import pylab as pl
    ROOT = "/data/id11/wolfgang/new_segmentation/erik/Al4Cu_dct1_"
    os.chdir(ROOT)
    # dark_n is 21
    nsigma = 8.
    stem =  os.path.join(ROOT, "0_rawdata", "Orig", "Al4Cu_dct1_")
    first = 0
    last = 3599
    names = filename_series( stem, first, last )
    minus_dark = open_sub_dark( names, "0_rawdata/Orig/darkend0000.edf", DARK_N=21 )
    masked = do_mask( minus_dark, "median_image-mask.edf" )
    #print(masked, type(masked))
    nsigma = 3
    rolling_bg = buffered_bgsub( masked, nsigma, bufsize=21 )
    smoothed = filter_gauss( rolling_bg, 2. )
    nsigma = 5 
    thresholded = apply_threshold( smoothed, nsigma )

    with h5py.File(os.path.join(ROOT, "0_rawdata","sparse", "peaksmed.hdf"),"w") as h5:
        gr = h5.require_group("peaks") 
        for j, frm in enumerate(thresholded):
            g = gr.require_group(str(j))
            frm.sparseframe.to_hdf_group( g )
            print("nnz",frm.sparseframe.nnz, end=" ")
            if 0:
                pl.figure(1)
                pl.imshow( frm.signal )
                pl.figure(2)
                pl.imshow( frm.mask )
                pl.show()

                
if __name__=="__main__":
    grab_folder()



    
