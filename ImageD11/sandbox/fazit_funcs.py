
import numpy as np
import hdf5plugin, fabio
from ImageD11 import blobcorrector, parameters, transform
import numba



def compute_tth_eta_lut( splinefile, pars, dims):
    """
    Computes look up values of tth, eta for each pixel
    """
    c = blobcorrector.correctorclass( splinefile )
    p = parameters.read_par_file( pars )
    xp, yp = c.make_pixel_lut( dims )
    t, e = transform.compute_tth_eta( (xp.ravel(), yp.ravel()),
                                      **p.parameters )
    t.shape=dims
    e.shape=dims
    return t, e


def tth_guess_bins( tth, npx = 1 ):
    one_px = (tth[1:] - tth[:-1]).max()
    binsize = one_px * npx
    itth = np.round(tth / binsize)
    return itth


def make_lut( itth, eta ):
    ie = 1./361
    sorter = (itth + (ie*eta)).ravel()
    lut = np.argsort( sorter ).astype( np.uint32 ) # read from here
    adrout = np.argsort( lut ).astype( np.uint32 ) # write to here
    return lut, adrout


def forward( data, lut, adrout ):
    out = np.empty_like( data )
    out.flat = data.flat[lut]
    return out

def back( data, lut, adrout ):
    out = np.empty_like( data )
    out.flat = data.flat[adrout]
    return out

@numba.njit
def forward_filt( data ):
    bg = np.empty( data.shape, np.float32 )
    b = data[0]
    bg[0] = b
    gain = 0.5
    sigmap = 0.2
    sigmat = 25.
    msk = np.empty( data.shape, np.uint8 )
    for i in range( 1, len(data) ):
        diff = data[i] - b
        t = sigmap * b + sigmat
        if diff > t:
            diff = t * diff / abs(diff) / 16
            msk[i] = 1
        elif diff < -t:
            diff = t * diff / abs(diff) / 4
            msk[i] = 1
        else:
            msk[i] = 0
        b += diff * gain
        bg[i] = b
    b = data[len(data)-1]
    for i in range( len(data)-2, -2, -1 ):
        diff = data[i] - b
        t = sigmap * b + sigmat
        if diff > t:
            diff = t * diff / abs(diff) / 16
            msk[i] += 1
        elif diff < -t:
            diff = t * diff / abs(diff) / 4
            msk[i] += 1
        else:
            if msk[i] == 1:
                bg[i] = b
        if msk[i] == 0 or msk[i] == 2:
            bg[i] = (bg[i]+b)/2
        b += diff * gain
    return bg, msk
        
    


if __name__=="__main__":
    import sys
    splinefile = sys.argv[1]
    pars = sys.argv[2]
    # dims = int(sys.argv[3]), int(sys.argv[4])
    fname = sys.argv[3]
    try:
        mask = fabio.open(sys.argv[4]).data > 0
    except:
        mask = None
    data  = fabio.open(fname).data
    dims = data.shape
    tth, eta = compute_tth_eta_lut( splinefile, pars, dims)
    print("tth: %.3f - %.3f   eta: %.3f - %.3f"%(
        tth.min(), tth.max(), eta.min(), eta.max()) )
    itth = tth_guess_bins( tth, npx = 0.7)
    print(itth.min(), itth.max())
    end = itth.max() + 1
    if mask is not None:
        itth[mask] += end # throw the masked to the end
    lut, adrout = make_lut( itth, eta )
    print(lut[0])

    trans = forward( data, lut, adrout )
    check = back( trans, lut, adrout )

    assert (check == data).all()

    bg, mr =[x.reshape(trans.shape) for x in forward_filt( trans.ravel() )]
    
    bgim = back( bg, lut, adrout )
    mr2 = back( mr, lut, adrout ) 
    bg1d, m2 =[x.reshape(data.shape) for x in forward_filt( data.ravel() )]

    print(bg, bg.min(), bg.max())
    from matplotlib.colors import LogNorm
    import pylab as pl
    f, ((ax1, ax2),(ax3,ax4)) = pl.subplots( 2,2, figsize=(12,3),
                                             sharex=True, sharey=True )
    ax1.imshow( data , vmax=1000, norm = LogNorm() )
    ax2.imshow( mr2, vmax=2,)
    ax3.imshow( data-bg1d, vmin=-10, vmax=200 )
    ax4.imshow( data-bgim, vmin=-10, vmax=200 )
    
    pl.show()
    
    
