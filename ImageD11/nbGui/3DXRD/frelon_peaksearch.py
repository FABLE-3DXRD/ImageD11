import h5py
import numpy as np
import scipy

import fabio
import functools
import ImageD11.cImageD11

import multiprocessing
from ImageD11.columnfile import columnfile
from tqdm.contrib.concurrent import process_map

from tqdm.notebook import tqdm

# we will now define our 2D segmenter and 3D merging functions
# this will eventually make its way into ImageD11 properly
# there are some important parameters that we need to ensure are correctly selected for good segmentation
# these are the following lines:
# self.m_offset = self.bg < 80
# self.m_ratio = self.bg > 200
# both in the "bgsub" function of the "worker" class
# these determine which part of the detector image is treated as the background, and which is treated as peaks
# these parameters may need adjusting when going from undeformed to deformed data
# to adjust them, change the numbers, then run the three following cells to plot the results of the segmentation on a single frame
# a good choice of numbers will mean we exclude regions between rings and powder background on the rings
# therefore only peaks created by grains should survive

def do3dmerge(cf_2d_dict, n, omega):   
    s = cf_2d_dict["s_raw"]
    f = cf_2d_dict["f_raw"]
    o = cf_2d_dict["o_raw"]
    I = cf_2d_dict["s_I"]
    s1 = cf_2d_dict["s_1"]
    
    print('Making 3D merge')

    # pointers to frames in s,f,o,I,n,s1
    p = np.cumsum(np.concatenate(([0, ], n)))
    # make a KDTree for each frame (wastes a bit of memory, but easier for sorting later)
    trees = [scipy.spatial.cKDTree(np.transpose((s[p[i]:p[i + 1]], f[p[i]:p[i + 1]])))
             for i in range(len(n))]
    print('made trees')
    # sys.stdout.flush()
    # because interlaced might not be in order
    order = np.argsort(omega % 360)
    # peaks that overlap, k : 0 -> npks == len(s|f|o)
    # diagonal
    krow = list(range(len(o)))
    kcol = list(range(len(o)))
    for i in range(1, len(n)):  # match these to previous
        flo = order[i - 1]
        fhi = order[i]
        tlo = trees[flo]
        thi = trees[fhi]
        # 1.6 is how close centers should be to overlap
        lol = trees[flo].query_ball_tree(trees[fhi], r=1.6)
        for srcpk, destpks in enumerate(lol):  # dest is strictly higher than src
            for destpk in destpks:
                krow.append(srcpk + p[flo])
                kcol.append(destpk + p[fhi])
    csr = scipy.sparse.csr_matrix((np.ones(len(krow), dtype=bool),
                                   (kcol, krow)), shape=(len(o), len(o)))
    # connected components == find all overlapping peaks
    ncomp, labels = scipy.sparse.csgraph.connected_components(csr,
                                                              directed=False, return_labels=True)
    
    print('connected components')
    # sys.stdout.flush()
    # Now merge the properties
    npkmerged = np.bincount(labels, minlength=ncomp)  # number of peaks that were merged
    s3d1 = np.bincount(labels, minlength=ncomp, weights=s1)  # s_1
    s3dI = np.bincount(labels, minlength=ncomp, weights=I)  # s_I
    ssI = np.bincount(labels, minlength=ncomp, weights=I * s)  # s_sI
    sfI = np.bincount(labels, minlength=ncomp, weights=I * f)  # s_sI
    soI = np.bincount(labels, minlength=ncomp, weights=I * o)  # s_sI
    s3d = ssI / s3dI
    f3d = sfI / s3dI
    o3d = soI / s3dI
    
    cf_3d_dict = {
        "s_raw": s3d,
        "f_raw": f3d,
        "omega": o3d,
        "sum_intensity": s3dI,
        "Number_of_pixels": s3d1
    }
    
    spot3d_id = labels
    
    return cf_3d_dict, spot3d_id

class worker:
    """ subtracts background, custom for ma4750 """

    def __init__(self, bgfile):
        self.bg = fabio.open(bgfile).data
        self.threshold = 50  # was 50 # ADU to zero out image
        self.smoothsigma = 1.  # sigma for Gaussian before labelleing
        self.bgc = 0.9  # fractional part of bg per peak to remove
        self.minpx = 3

        self.m_offset = self.bg < 80
        
        self.mbg = np.mean(self.bg[self.m_offset])
        # self.m_ratio = self.bg > 135  # best for undeformed data
        
        self.m_ratio = self.bg > 200
        
        self.bg -= self.mbg  # remove dark
        self.invbg = 1 / self.bg[self.m_ratio]
        self.bins = b = np.linspace(0, 2, 256)
        self.bc = (b[1:] + b[:-1]) / 2

        self.wrk = np.empty(self.bg.shape, 'b')
        self.labels = np.empty(self.bg.shape, 'i')

    def bgsub(self, img):
        img = img.astype(np.float32)
        offset = np.mean(img[self.m_offset])  # remove dark
        np.subtract(img, offset, img)
        ratio = img[self.m_ratio] * self.invbg
        h, b = np.histogram(ratio, bins=self.bins)
        htrim = np.where(h < h.max() * 0.05, 0, h)
        r = (htrim * self.bc).sum() / htrim.sum()
        # Better to scale background to data rather than data to background?
        # np.multiply(img, 1 / r, img)
        np.subtract(img, self.bg * r, img)
        self.offset = offset
        self.scale = r
        return img

    def peaksearch(self, img, omega=0):
        self.cor = self.bgsub(img)
        # smooth the image for labelling (removes noise maxima)
        self.smoothed = scipy.ndimage.gaussian_filter(self.cor, self.smoothsigma)
        assert self.smoothed.dtype == np.float32
        # zero out the background
        self.mt = self.smoothed < self.threshold
        self.smoothed[self.mt] = 0
        # label on smoothed image
        self.npks = ImageD11.cImageD11.localmaxlabel(self.smoothed, self.labels, self.wrk)
        self.labels[self.mt] = 0
        # now find the borders of each blob : first dilate
        l3 = scipy.ndimage.uniform_filter(self.labels * 7, 3)
        self.borders = (self.labels * 7) != l3
        # border properties - use the real data or the smoothed? Real data. Might be noisier
        self.blobsouter = ImageD11.cImageD11.blobproperties(self.cor, self.labels * self.borders, self.npks)
        # Computed background per peak
        self.per_peak_bg = np.concatenate(([0, ], self.blobsouter[:, ImageD11.cImageD11.mx_I]))
        self.bgcalc = self.per_peak_bg[self.labels]
        self.m_top = self.cor > self.bgcalc * self.bgc
        self.forprops = self.cor * self.m_top
        self.blobs = ImageD11.cImageD11.blobproperties(self.forprops,
                                                       self.labels * self.m_top,
                                                       self.npks,
                                                       omega=omega)
        ImageD11.cImageD11.blob_moments(self.blobs)
        self.enoughpx = self.blobs[:, ImageD11.cImageD11.s_1] >= self.minpx
        self.goodpeaks = self.blobs[self.enoughpx]
        return self.goodpeaks

    def plots(self):
        pass
    
    
@functools.lru_cache(maxsize=1)
def get_dset(h5name, dsetname):
    """This avoids to re-read the dataset many times"""
    dset = h5py.File(h5name, "r")[dsetname]
    return dset


def pps(arg):
    hname, dsetname, num, omega, bgfile = arg
    if pps.worker is None:
        pps.worker = worker(bgfile)
    frm = get_dset(hname, dsetname)[num]
    pks = pps.worker.peaksearch(frm, omega=omega)
    return num, pks


pps.worker = None

# PKSAVE = 's_raw f_raw o_raw s_1 s_I m_ss m_ff m_oo m_sf m_so m_fo'.split()
PKSAVE = ["s_1", "s_I", "s_I2", "s_fI", "s_ffI", "s_sI", "s_ssI", "s_sfI", "s_oI", "s_ooI", "s_soI", "s_foI", "mx_I", "mx_I_f", "mx_I_s", "mx_I_o", "bb_mx_f", "bb_mx_s", "bb_mx_o", "bb_mn_f", "bb_mn_s", "bb_mn_o", "avg_i", "f_raw", "s_raw", "o_raw", "m_ss", "m_ff", "m_oo", "m_sf", "m_so", "m_fo"]
PKCOL = [getattr(ImageD11.cImageD11, p) for p in PKSAVE]

def process(ds, bgfile, ncpu):      
    hname = ds.masterfile
    scan_name = ds.scans[0]
    frames_dset = scan_name + "/measurement/" + ds.detector
    omega = ds.omega[0, :]
    
    n_frames = omega.shape[0]

    args = [(hname, frames_dset, i, omega[i], bgfile) for i in range(n_frames)]
    
    all_peaks = process_map(pps, args, chunksize=1)
    
    # make a dict to hold our results in
    cf_2d_dict = {}
    
    # populate it with empty lists
    for title in PKSAVE:
        cf_2d_dict[title] = []
    
    npks = []

    for peak_props in all_peaks:
        this_index, props = peak_props
        
        for title, prop_index in zip(PKSAVE, PKCOL):
            cf_2d_dict[title].extend(props[:, prop_index])
            
        npks.append(props.shape[0])
    
    # convert to numpy arrays:
    for key, value in cf_2d_dict.items():
        cf_2d_dict[key] = np.array(value)
    
    cf_3d_dict, spot3d_id = do3dmerge(cf_2d_dict, npks, omega)
    
    cf_2d = ImageD11.columnfile.colfile_from_dict(cf_2d_dict)
    
    cf_2d.addcolumn(spot3d_id, "spot3d_id")
    cf_2d.addcolumn(cf_2d.o_raw, "omega")
    
    cf_3d = ImageD11.columnfile.colfile_from_dict(cf_3d_dict)
    
    cf_3d.addcolumn(np.arange(cf_3d.nrows), "index")
    
    return cf_2d, cf_3d
