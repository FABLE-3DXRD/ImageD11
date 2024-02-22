import numpy as np
import numba

from tqdm.contrib.concurrent import process_map
import ImageD11.unitcell
import ImageD11.refinegrains
import ImageD11.blobcorrector
from ImageD11.blobcorrector import eiger_spatial


from matplotlib import pyplot as plt


def correct_pixel(pixel, spline_file):
    sr, fr = pixel
    sc, fc = ImageD11.blobcorrector.correctorclass(spline_file).correct(sr, fr)
    return (sc, fc)


def apply_spatial(cf, spline_file):
    # sc = np.zeros(cf.nrows)
    # fc = np.zeros(cf.nrows)
    
    print("Spatial correction...")
    
    raw_pixels = np.vstack((cf['s_raw'], cf['f_raw'])).T
    
    corrected_pixels = process_map(correct_pixel, raw_pixels, [spline_file] * len(raw_pixels), max_workers=63, chunksize=len(raw_pixels)//63)
    
    sc, fc = [list(t) for t in zip(*corrected_pixels)]
        
    cf.addcolumn(sc, "sc")
    cf.addcolumn(fc, "fc")
    
    return cf

def grain_to_rgb(g, ax=(0,0,1)):
    return hkl_to_color_cubic(crystal_direction_cubic(g.ubi, ax))

def crystal_direction_cubic(ubi, axis):
    hkl = np.dot(ubi, axis)
    # cubic symmetry implies:
    #      24 permutations of h,k,l
    #      one has abs(h) <= abs(k) <= abs(l)
    hkl= abs(hkl)
    hkl.sort()
    return hkl

def hkl_to_color_cubic(hkl):
    """
    https://mathematica.stackexchange.com/questions/47492/how-to-create-an-inverse-pole-figure-color-map
        [x,y,z]=u⋅[0,0,1]+v⋅[0,1,1]+w⋅[1,1,1].
            These are:
                u=z−y, v=y−x, w=x
                This triple is used to assign each direction inside the standard triangle
                
    makeColor[{x_, y_, z_}] := 
         RGBColor @@ ({z - y, y - x, x}/Max@{z - y, y - x, x})                
    """
    x,y,z = hkl
    assert x<=y<=z
    assert z>=0
    u,v,w = z-y, y-x, x
    m = max( u, v, w )
    r,g,b = u/m, v/m, w/m
    return (r,g,b)

def hkl_to_pf_cubic(hkl):
    x,y,z = hkl
    assert x<=y<=z
    assert z>=0
    m = np.sqrt((hkl**2).sum())
    return x/(z+m), y/(z+m)

def triangle():
    """ compute a series of point on the edge of the triangle """
    xy = [ np.array(v) for v in ( (0,1,1), (0,0,1), (1,1,1)) ]
    xy += [ xy[2]*(1-t) + xy[0]*t for t in np.linspace(0.1,1,5)]
    return np.array( [hkl_to_pf_cubic( np.array(p) ) for p in xy] )


def calcy(cos_omega, sin_omega, sol):
    return sol[0] + cos_omega*sol[1] + sin_omega*sol[2]

def fity(y, cos_omega, sin_omega, wt=1):
    """
    Fit a sinogram to get a grain centroid
    # calc = d0 + x*co + y*so
    # dc/dpar : d0 = 1
    #         :  x = co
    #         :  y = so
    # gradients
    # What method is being used here???????????
    """
    g = [wt*np.ones(y.shape, float),  wt*cos_omega, wt*sin_omega]
    nv = len(g)
    m = np.zeros((nv,nv),float)
    r = np.zeros( nv, float )
    for i in range(nv):
        r[i] = np.dot( g[i], wt * y )
        for j in range(i,nv):
            m[i,j] = np.dot( g[i], g[j] )
            m[j,i] = m[i,j]
    sol = np.dot(np.linalg.inv(m), r)
    return sol

def fity_robust(dty, co, so, nsigma=5, doplot=False):
    # NEEDS COMMENTING
    cen, dx, dy = fity(dty, co, so)
    calc2 = calc1 = calcy(co, so, (cen, dx, dy))
    # mask for columnfile, we're selecting specific 4D peaks
    # that come from the right place in y, I think?
    selected = np.ones(co.shape, bool)
    for i in range(3):
        err = dty - calc2
        estd = max( err[selected].std(), 1.0 ) # 1 micron
        #print(i,estd)
        es = estd*nsigma
        selected = abs(err) < es
        cen, dx, dy = fity( dty, co, so, selected.astype(float) )
        calc2 = calcy(co, so, (cen, dx, dy))
    # bad peaks are > 5 sigma
    if doplot:
        f, a = plt.subplots(1,2)
        theta = np.arctan2( so, co )
        a[0].plot(theta, calc1, ',')
        a[0].plot(theta, calc2, ',')
        a[0].plot(theta[selected], dty[selected], "o")
        a[0].plot(theta[~selected], dty[~selected], 'x')
        a[1].plot(theta[selected], (calc2 - dty)[selected], 'o')
        a[1].plot(theta[~selected], (calc2 - dty)[~selected], 'x')
        a[1].set(ylim = (-es, es))
        plt.show()
    return selected, cen, dx, dy

def graincen(gid, colf, doplot=True):
    # Get peaks beloging to this grain ID
    m = colf.grain_id == gid
    # Get omega values of peaks in radians
    romega = np.radians(colf.omega[m])
    # Calculate cos and sin of omega
    co = np.cos(romega)
    so = np.sin(romega)
    # Get dty values of peaks
    dty = colf.dty[m]
    selected, cen, dx, dy = fity_robust(dty, co, so, doplot=doplot)
    return selected, cen, dx, dy




@numba.njit(parallel=True)
def pmax(ary):
    """ Find the min/max of an array in parallel """
    mx = ary.flat[0]
    mn = ary.flat[0]
    for i in numba.prange(1,ary.size):
        mx = max( ary.flat[i], mx )
        mn = min( ary.flat[i], mn )
    return mn, mx

@numba.njit(parallel=True)
def palloc(shape, dtype):
    """ Allocate and fill an array with zeros in parallel """
    ary = np.empty(shape, dtype=dtype)
    for i in numba.prange( ary.size ):
        ary.flat[i] = 0
    return ary

# counting sort by grain_id
@numba.njit
def counting_sort(ary, maxval=None, minval=None):
    """ Radix sort for integer array. Single threaded. O(n)
    Numpy should be doing this...
    """
    if maxval is None:
        assert minval is None
        minval, maxval = pmax( ary ) # find with a first pass
    maxval = int(maxval)
    minval = int(minval)
    histogram = palloc( (maxval - minval + 1,), np.int64 )
    indices = palloc( (maxval - minval + 2,), np.int64 )
    result = palloc( ary.shape, np.int64 )
    for gid in ary:
        histogram[gid - minval] += 1
    indices[0] = 0
    for i in range(len(histogram)):
        indices[ i + 1 ] = indices[i] + histogram[i]
    i = 0
    for gid in ary:
        j = gid - minval
        result[indices[j]] = i
        indices[j] += 1
        i += 1
    return result, histogram


@numba.njit(parallel=True)
def find_grain_id(spot3d_id, grain_id, spot2d_label, grain_label, order, nthreads=20):
    """
    Assignment grain labels into the peaks 2d array
    spot3d_id = the 3d spot labels that are merged and indexed
    grain_id = the grains assigned to the 3D merged peaks
    spot2d_label = the 3d label for each 2d peak
    grain_label => output, which grain is this peak
    order = the order to traverse spot2d_label sorted
    """
    assert spot3d_id.shape == grain_id.shape
    assert spot2d_label.shape == grain_label.shape
    assert spot2d_label.shape == order.shape
    T = nthreads
    print("Using",T,"threads")
    for tid in numba.prange( T ):
        pcf = 0 # thread local I hope?
        for i in order[tid::T]:
            grain_label[i] = -1
            pkid = spot2d_label[i]
            while spot3d_id[pcf] < pkid:
                pcf += 1
            if spot3d_id[pcf] == pkid:
                grain_label[i] = grain_id[pcf]
                

def tocolf(pkd, parfile, dxfile, dyfile):
    """ Converts a dictionary of peaks into an ImageD11 columnfile
    adds on the geometric computations (tth, eta, gvector, etc) """
    spat = eiger_spatial(dxfile=dxfile, dyfile=dyfile)
    cf = ImageD11.columnfile.colfile_from_dict(spat(pkd))
    cf.parameters.loadparameters(parfile)
    cf.updateGeometry()
    return cf









def unitcell_peaks_mask(cf, dstol, dsmax):
    cell = ImageD11.unitcell.unitcell_from_parameters(cf.parameters)
    cell.makerings(dsmax)
    m = np.zeros(cf.nrows, bool)
    for v in cell.ringds:
        if v < dsmax:
            m |= (abs(cf.ds - v) < dstol)

    return m

def strongest_peaks(colf, uself=True, frac=0.995, B=0.2, doplot=None):
    # correct intensities for structure factor (decreases with 2theta)
    cor_intensity = colf.sum_intensity * (np.exp(colf.ds*colf.ds*B))
    if uself:
        lf = ImageD11.refinegrains.lf(colf.tth, colf.eta)
        cor_intensity *= lf
    order = np.argsort( cor_intensity )[::-1] # sort the peaks by intensity
    sortedpks = cor_intensity[order]
    cums =  np.cumsum(sortedpks)
    cums /= cums[-1]
    enough = np.searchsorted(cums, frac)
    # Aim is to select the strongest peaks for indexing.
    cutoff = sortedpks[enough]
    mask = cor_intensity > cutoff
    if doplot is not None:
        fig, axs = plt.subplots(1,2,figsize=(10,5))
        axs[0].plot(cums/cums[-1], ',')
        axs[0].set(xlabel='npks',ylabel='fractional intensity')
        axs[0].plot([mask.sum(),], [frac,], "o" )
        axs[1].plot(cums/cums[-1], ',')
        axs[1].set(xlabel='npks logscale',ylabel='fractional intensity', xscale='log', ylim=(doplot,1.), 
                 xlim=(np.searchsorted(cums, doplot), len(cums)))
        axs[1].plot( [mask.sum(),], [frac,], "o" )
        plt.show()
    return mask

def selectpeaks(cf, dstol=0.005, dsmax=10, frac=0.99, doplot=None):
    m = unitcell_peaks_mask(cf, dstol=dstol, dsmax=dsmax)
    cfc = cf.copy()
    cfc.filter(m)
    ms = strongest_peaks(cfc, frac=frac, doplot=doplot)
    cfc.filter(ms)
    return cfc


