
from __future__ import print_function

import sys, os
import timeit
import multiprocessing
import h5py, hdf5plugin, fabio
import functools
import numpy as np
import tqdm
import ImageD11.cImageD11
import scipy.sparse, scipy.spatial, scipy.ndimage

nthreads = int( os.environ[ 'SLURM_CPUS_PER_TASK'] )

def do3dmerge( honame, scan ):
    # read the 2d peak search results
    t0 = timeit.default_timer()
    print('Making 3D merge', end=' ')
    sys.stdout.flush()
    with h5py.File( honame, 'r') as hin:
        g = hin[scan]['peaks2d']
        h5input = g.attrs['h5input']
        h5scan = g.attrs['h5scan']
        detector_path = g.attrs['detector_path']
        omega_path = g.attrs['omega_path']
        omega = g['omega'][:]
        s = g['s_raw'][:]
        f = g['f_raw'][:]
        o = g['o_raw'][:] 
        I = g['s_I'][:]
        n = g['npks'][:]
        s1 = g['s_1'][:]
    # pointers to frames in s,f,o,I,n,s1
    p= np.cumsum( np.concatenate(([0,], n ) ) )
    # make a KDTree for each frame (wastes a bit of memory, but easier for sorting later)
    trees = [ scipy.spatial.cKDTree( np.transpose( (s[p[i]:p[i+1]], f[p[i]:p[i+1]]) ) )
              for i in range(len(n)) ]
    print('made trees',end=' ')
    sys.stdout.flush()
    # because interlaced might not be in order
    order = np.argsort( omega % 360 )
    # peaks that overlap, k : 0 -> npks == len(s|f|o)
    # diagonal
    krow = list(range(len(o)))
    kcol = list(range(len(o)))
    for i in range(1,len(n)):  # match these to previous
        flo = order[i-1]
        fhi = order[i]
        tlo = trees[flo]
        thi = trees[fhi]
        # 1.6 is how close centers should be to overlap
        lol = trees[flo].query_ball_tree( trees[fhi], r = 1.6 ) 
        for srcpk, destpks in enumerate(lol):  # dest is strictly higher than src
            for destpk in destpks:
                krow.append( srcpk + p[flo] )
                kcol.append( destpk + p[fhi] )
    csr = scipy.sparse.csr_matrix( ( np.ones(len(krow), dtype=bool), 
                                    (kcol, krow) ), shape=(len(o),len(o)) )                
    # connected components == find all overlapping peaks
    ncomp, labels= scipy.sparse.csgraph.connected_components( csr, 
                                         directed=False, return_labels=True )
    print('connected components', end=' ')
    sys.stdout.flush()
    # Now merge the properties
    npkmerged = np.bincount( labels, minlength = ncomp )               # number of peaks that were merged
    s3d1      = np.bincount( labels, minlength = ncomp, weights=s1 )   # s_1
    s3dI      = np.bincount( labels, minlength = ncomp, weights=I )    # s_I
    ssI       = np.bincount( labels, minlength = ncomp, weights=I*s )  # s_sI
    sfI       = np.bincount( labels, minlength = ncomp, weights=I*f )  # s_sI
    soI       = np.bincount( labels, minlength = ncomp, weights=I*o )  # s_sI
    s3d = ssI / s3dI
    f3d = sfI / s3dI
    o3d = soI / s3dI
    with h5py.File( honame, 'a') as hin:
        g = hin[scan].require_group('peaks3d')
        g['s_raw'] = s3d
        g['f_raw'] = f3d
        g['omega'] = o3d
        g['sum_intensity'] = s3dI
        g['number_of_pixels'] = s3d1
    t1 = timeit.default_timer()
    print('wrote %.3f/s'%(t1-t0))


class worker: 
    """ subtracts background, custom for ma47050 """
    def __init__(self, bgfile, offset=100, ratio=120):
        self.bg = fabio.open(bgfile).data.astype(np.float32)
        self.threshold = 50    # ADU to zero out image
        self.smoothsigma = 1.  # sigma for Gaussian before labelleing
        self.bgc = 0.9         # fractional part of bg per peak to remove
        self.minpx = 3
        
        self.m_offset = self.bg < offset
        self.mbg = np.mean( self.bg[self.m_offset] )
        self.m_ratio  = self.bg > ratio
        self.invbg = 1 / self.bg[self.m_ratio]
        self.bins = b = np.linspace(0,2,256)
        self.bc = (b[1:]+b[:-1])/2

        self.wrk = np.empty( self.bg.shape, 'b')
        self.labels = np.empty( self.bg.shape, 'i' )
        
    def bgsub(self, img):
        img = img.astype( np.float32 )
        offset = np.mean( img[self.m_offset] ) - self.mbg
        np.subtract( img, offset, img)
        ratio =  img[self.m_ratio] * self.invbg
        h,b = np.histogram( ratio, bins=self.bins )
        htrim = np.where(h < h.max()*0.05, 0, h)
        r = (htrim * self.bc).sum()/htrim.sum()
        np.multiply( img, 1/r, img )
        np.subtract( img, self.bg, img )
        self.offset = offset
        self.scale = r
        return img

    def peaksearch(self, img, omega = 0):
        self.cor = self.bgsub( img )
        # smooth the image for labelling (removes noise maxima)
        self.smoothed = scipy.ndimage.gaussian_filter(self.cor , self.smoothsigma)
        assert self.smoothed.dtype == np.float32
        # zero out the background
        self.mt = self.smoothed < self.threshold
        self.smoothed[self.mt] = 0
        # label on smoothed image
        self.npks = ImageD11.cImageD11.localmaxlabel(self.smoothed, self.labels, self.wrk)
        self.labels[self.mt] = 0
        # now find the borders of each blob : first dilate
        l3 = scipy.ndimage.uniform_filter(self.labels*7, 3)
        self.borders = (self.labels*7) != l3
        # border properties - use the real data or the smoothed? Real data. Might be noisier
        self.blobsouter = ImageD11.cImageD11.blobproperties( self.cor, self.labels*self.borders, self.npks )
        # Computed background per peak
        self.per_peak_bg = np.concatenate( ([0,], self.blobsouter[:,ImageD11.cImageD11.mx_I]))
        self.bgcalc = self.per_peak_bg[self.labels]
        self.m_top = self.cor > self.bgcalc * self.bgc
        self.forprops = self.cor*self.m_top
        self.blobs = ImageD11.cImageD11.blobproperties( self.forprops, 
                                                        self.labels * self.m_top, 
                                                        self.npks,
                                                        omega = omega)
        ImageD11.cImageD11.blob_moments(self.blobs)
        self.enoughpx = self.blobs[:,ImageD11.cImageD11.s_1] >= self.minpx
        self.goodpeaks = self.blobs[self.enoughpx]
        return self.goodpeaks
    
    def plots(self, ax):
        import pylab as pl
        ax.pcolormesh( abs(self.cor), norm=pl.matplotlib.colors.LogNorm() )
        ax.scatter( self.blobs[:, ImageD11.cImageD11.f_raw],
                    self.blobs[:, ImageD11.cImageD11.s_raw], 
                    s = np.sqrt(self.blobs[:, ImageD11.cImageD11.s_1])/np.pi,
                    edgecolors='m',
                    facecolors='None',) 

    
@functools.lru_cache(maxsize=1)
def get_dset(h5name, dsetname):
    """This avoids to re-read the dataset many times"""
    dset = h5py.File(h5name, "r")[dsetname]
    return dset

def pps( arg ):
    hname, dsetname, num, omega, bgfile = arg
    if pps.worker is None:
        pps.worker = worker(bgfile)
    frm = get_dset( hname, dsetname)[num]
    pks = pps.worker.peaksearch( frm, omega=omega )
    return num, pks
    
pps.worker = None

PKSAVE = 's_raw f_raw o_raw s_1 s_I'.split()
PKCOL = [getattr( ImageD11.cImageD11, p) for p in PKSAVE]
#    s_1, s_I, s_I2,\
#    s_fI, s_ffI, s_sI, s_ssI, s_sfI, s_oI, s_ooI, s_foI, s_soI, \
#    bb_mn_f, bb_mn_s, bb_mx_f, bb_mx_s, bb_mn_o, bb_mx_o, \
#    mx_I, mx_I_f, mx_I_s, mx_I_o, dety, detz, \
#    avg_i, f_raw, s_raw, o_raw, f_cen, s_cen, \
#    m_ss, m_ff, m_oo, m_sf, m_so, m_fo


def process( hname, scan, detector, omeganame, bgfile, houtname, mypool ):
    
    with h5py.File( hname, 'r') as hin:
        shp = hin[scan][detector].shape
        omega = hin[scan][omeganame][()]
        
    assert len(omega) == shp[0]
    detname = "/".join((scan, detector))
    args = [ (hname, detname, i, omega[i], bgfile) for i in range(shp[0]) ]
    t0 = timeit.default_timer()
    
    with h5py.File( houtname, "a") as hout:
        g = hout.require_group( scan + '/peaks2d')
        g.attrs['h5input']=hname
        g.attrs['h5scan']=scan
        g.attrs['detector_path']=detector
        g.attrs['omega_path']=omeganame
        g.attrs['bgfile']=bgfile
        g.attrs['script_source'] = open(__file__, 'rb').read()
        g['omega'] = omega
        g['npks'] = np.zeros( len(omega), int ) 
        for name in PKSAVE:
            g.create_dataset( name, 
                              dtype = np.float32, 
                              shape = (1,), 
                              chunks = (4096,),
                              maxshape = (None, ) )
        
        npks = 0
        i = -1
        for num, pks in tqdm.tqdm( mypool.imap( pps, args ), total=shp[0] ):
            i += 1
            if len(pks) == 0:
                continue          
            npk  = len(pks)
            npks += npk
            g['npks'][i] = npk
            for name, col in zip(PKSAVE, PKCOL):
                g[name].resize( (npks,) )
                g[name][-npk:] = pks[:, col]
    # closes, flushes, now merge
    do3dmerge( houtname, scan )
        
if __name__=="__main__":
    houtname = sys.argv[1]
    hname = sys.argv[2]
    detector = sys.argv[3]
    omega = sys.argv[4]
    bgfile = sys.argv[5]
    # '/data/visitor/ma4750/id11/analysis_jw/march2022/background_big_beam.edf'
    scans = sys.argv[6:]
    with multiprocessing.Pool(nthreads) as mypool:
        for scan in scans:
            process( hname, scan, detector, omega, bgfile, houtname, mypool )


        
        
