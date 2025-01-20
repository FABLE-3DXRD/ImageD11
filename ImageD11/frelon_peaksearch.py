import h5py
import numpy as np
import scipy

import fabio
import functools
import ImageD11.cImageD11

import multiprocessing
from ImageD11.columnfile import columnfile
from tqdm.contrib.concurrent import process_map

from tqdm.auto import tqdm

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


def do3dmerge(cf_2d_dict, n, omega, radius=1.6):
    """
    Merges the peaks in cf_2d_dict if the center is within radius pixels
    n = the number of peaks on each frame in cf_2d_dict
    omega = the rotation motor position for each frame

    returns cf_3d_dict = merged peaks, as well as the input
    """
    s = cf_2d_dict["s_raw"]
    f = cf_2d_dict["f_raw"]
    o = cf_2d_dict["o_raw"]
    I = cf_2d_dict["s_I"]
    s1 = cf_2d_dict["s_1"]

    print("Making 3D merge")

    # pointers to frames in s,f,o,I,n,s1
    p = np.cumsum(
        np.concatenate(
            (
                [
                    0,
                ],
                n,
            )
        )
    )
    # make a KDTree for each frame (wastes a bit of memory, but easier for sorting later)
    trees = [
        scipy.spatial.cKDTree(np.transpose((s[p[i] : p[i + 1]], f[p[i] : p[i + 1]])))
        for i in range(len(n))
    ]
    print("made trees")
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
        lol = trees[flo].query_ball_tree(trees[fhi], r=radius)
        for srcpk, destpks in enumerate(lol):  # dest is strictly higher than src
            for destpk in destpks:
                krow.append(srcpk + p[flo])
                kcol.append(destpk + p[fhi])
    csr = scipy.sparse.csr_matrix(
        (np.ones(len(krow), dtype=bool), (kcol, krow)), shape=(len(o), len(o))
    )
    # connected components == find all overlapping peaks
    ncomp, labels = scipy.sparse.csgraph.connected_components(
        csr, directed=False, return_labels=True
    )
    print("connected components")
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
        "Number_of_pixels": s3d1,
    }

    cf_2d_dict["spot3d_id"] = labels
    return cf_3d_dict


class worker:
    """subtracts background, custom for ma4750"""

    def __init__(
        self,
        bgfile,
        maskfile,
        threshold=50,
        smoothsigma=1.0,
        bgc=0.9,
        minpx=3,
        m_offset_thresh=100,
        m_ratio_thresh=135,
        flatfile=None,
        darkfile=None,
    ):
        """
        bgfile = image containing the background
        maskfile = detector mask

        threshold = ADU below which to zero out image
        smoothsigma =  sigma for Gaussian filter before labelleing
        bgc = fractional part of bg per peak to remove
        minpx = minimum number of pixels in a peak
        m_offset_thresh = bgfile less than this is constant
        m_ratio_thresh  = bgfile greather than this is scaled?

        flatfile = detector sensitivity image.
        darkfile = detector offset image
        """
        self.flat = None
        self.dark = None
        self.bg = None
        self.mask = None

        if flatfile is not None:
            if isinstance(flatfile, np.ndarray):
                self.flat = flatfile
            else:
                self.flat = fabio.open(flatfile).data

        if darkfile is not None:
            if isinstance(flatfile, np.ndarray):
                self.flat = flatfile
            else:
                self.flat = fabio.open(flatfile).data

        if bgfile is not None:
            if isinstance(bgfile, np.ndarray):
                self.bg = bgfile
            else:
                self.bg = self.correct(fabio.open(bgfile).data)

        if maskfile is not None:
            if isinstance(maskfile, np.ndarray):
                self.mask = maskfile
            else:
                self.mask = fabio.open(maskfile).data

        self.threshold = threshold  # was 50 # ADU to zero out image
        self.smoothsigma = smoothsigma  # sigma for Gaussian before labelleing
        self.bgc = bgc  # fractional part of bg per peak to remove
        self.minpx = minpx

        if self.bg is not None:
            self.m_offset = self.bg < m_offset_thresh
            self.m_ratio = self.bg > m_ratio_thresh
            self.mbg = np.mean(self.bg[self.m_offset])
            self.bg -= self.mbg  # remove dark
            self.invbg = 1 / self.bg[self.m_ratio]

        self.bins = b = np.linspace(0, 2, 256)
        self.bc = (b[1:] + b[:-1]) / 2  # bin cens
        self.wrk = None
        self.labels = None

    def correct(self, img):
        """
        Applies dark + flat if given
        """
        cor = img.astype(np.float32)
        if self.dark is not None:
            np.subtract(img, self.dark, img)
        if self.flat is not None:
            np.divide(img, self.flat, img)
        return cor

    def bgsub(self, img):
        """
        This attempts to remove a background by rescaling the self.bg image.

        """
        img = self.correct(img)
        if self.bg is None:  # use image border
            self.scale = 1
            self.offset = np.median(
                (img[0].min(), img[-1].min(), img[:, 0].min(), img[:, -1].min())
            )
            return img - self.offset
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

    def masksub(self, img):
        # active pixels are 0, masked pixels are 1
        img_masked = img.copy()
        if self.mask is not None:
            img_masked[self.mask == 1] = 0
        return img_masked

    def peaksearch(self, img, omega=0):
        if self.wrk is None:
            self.wrk = np.empty(img.shape, "b")
            self.labels = np.empty(img.shape, "i")

        self.cor = self.bgsub(img)
        self.cor = self.masksub(self.cor)
        # smooth the image for labelling (removes noise maxima)
        self.smoothed = scipy.ndimage.gaussian_filter(self.cor, self.smoothsigma)
        assert self.smoothed.dtype == np.float32
        # zero out the background
        self.mt = self.smoothed < self.threshold
        self.smoothed[self.mt] = 0
        # label on smoothed image
        self.npks = ImageD11.cImageD11.localmaxlabel(
            self.smoothed, self.labels, self.wrk
        )
        self.labels[self.mt] = 0
        # now find the borders of each blob : first dilate
        l3 = scipy.ndimage.uniform_filter(self.labels * 7, 3)
        self.borders = (self.labels * 7) != l3
        # border properties - use the real data or the smoothed? Real data. Might be noisier
        self.blobsouter = ImageD11.cImageD11.blobproperties(
            self.cor, self.labels * self.borders, self.npks
        )
        # Computed background per peak
        self.per_peak_bg = np.concatenate(
            (
                [
                    0,
                ],
                self.blobsouter[:, ImageD11.cImageD11.mx_I],
            )
        )
        self.bgcalc = self.per_peak_bg[self.labels]
        self.m_top = self.cor > self.bgcalc * self.bgc
        self.forprops = self.cor * self.m_top
        self.blobs = ImageD11.cImageD11.blobproperties(
            self.forprops, self.labels * self.m_top, self.npks, omega=omega
        )
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
    hname, dsetname, num, omega, worker_args = arg
    if pps.worker is None:
        pps.worker = worker(**worker_args)
    frm = get_dset(hname, dsetname)[num]
    pks = pps.worker.peaksearch(frm, omega=omega)
    return num, pks


def guess_bg(ds, scan_number=0, start=0, step=9, n=None):
    """
    ds = dataset object
    scan_number = index into ds.scans
    start = 0 = first image to look at
    step is the step for going through frames
    n = number of frames to look at (None == all)

    returns median(  [ frames[j::step,:,:].min() for j in range(step) ] )
    e.g. take the minimum of every step images
         then return the median of these
    """
    with h5py.File(ds.masterfile, "r") as hin:
        dset = hin[ds.scans[scan_number]]["measurement"][ds.detector]
        bgj = []
        end = len(dset)
        if n is not None:
            end = min(end, n + start)
        for j in range(step):  # loop over step offset
            bg = dset[start + j]
            for i in range(start + step, end, step):  # min for this series
                frm = dset[i + j]
                bg = np.where(bg > frm, frm, bg)
            bgj.append(bg)
    return np.median(bgj, axis=0)


pps.worker = None

# PKSAVE = 's_raw f_raw o_raw s_1 s_I m_ss m_ff m_oo m_sf m_so m_fo'.split()
PKSAVE = [
    "s_1",
    "s_I",
    "s_I2",
    "s_fI",
    "s_ffI",
    "s_sI",
    "s_ssI",
    "s_sfI",
    "s_oI",
    "s_ooI",
    "s_soI",
    "s_foI",
    "mx_I",
    "mx_I_f",
    "mx_I_s",
    "mx_I_o",
    "bb_mx_f",
    "bb_mx_s",
    "bb_mx_o",
    "bb_mn_f",
    "bb_mn_s",
    "bb_mn_o",
    "avg_i",
    "f_raw",
    "s_raw",
    "o_raw",
    "m_ss",
    "m_ff",
    "m_oo",
    "m_sf",
    "m_so",
    "m_fo",
]
PKCOL = [getattr(ImageD11.cImageD11, p) for p in PKSAVE]


def process(ds, worker_args, ncpu=None, **kwargs):
    """
    Runs over the first scan in a dataset in parallel

    ds = ImageD11.sinograms.dataset
    work_args = see worker.__init__.__doc__
    ncpu = cores to use (None = guess)

    returns cf_2d, cf_3d
    """
    if ncpu is None:
        nthreads = max(1, ImageD11.cImageD11.cores_available() - 1)
    else:
        nthreads = int(ncpu)
    hname = ds.masterfile
    scan_name = ds.scans[0]
    frames_dset = scan_name + "/measurement/" + ds.detector
    omega = ds.omega[0, :]

    n_frames = omega.shape[0]

    args = [(hname, frames_dset, i, omega[i], worker_args) for i in range(n_frames)]

    all_peaks = process_map(pps, args, chunksize=1, max_workers=nthreads, **kwargs)

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

    cf_3d_dict = do3dmerge(cf_2d_dict, npks, omega)

    # conversion to a dictionary is a dataset method
    cf_2d_dict["omega"] = cf_2d_dict["o_raw"]

    cf_2d_dict.pop("o_raw")  # removes duplicate
    cf_2d = ds.get_colfile_from_peaks_dict(cf_2d_dict)  # does the spatial

    cf_3d_dict["spot3d_id"] = np.arange(len(cf_3d_dict["s_raw"]))
    cf_3d = ds.get_colfile_from_peaks_dict(cf_3d_dict)

    return cf_2d, cf_3d
