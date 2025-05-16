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
    sfI = np.bincount(labels, minlength=ncomp, weights=I * f)  # s_fI
    soI = np.bincount(labels, minlength=ncomp, weights=I * o)  # s_oI
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


def do4dmerge(cf_2d_dict, n, omega, dty, ostep, ystep, radius=1.6):
    """
    Merges the peaks in cf_2d_dict if the center is within radius pixels
    n = the number of peaks on each frame in cf_2d_dict
    omega = the rotation motor position for each frame
    dty = the dty motor position for each frame

    returns cf_4d_dict = merged peaks, as well as the input
    """
    s = cf_2d_dict["s_raw"]
    f = cf_2d_dict["f_raw"]
    o = cf_2d_dict["o_raw"]
    d = cf_2d_dict["dty"]
    I = cf_2d_dict["s_I"]
    s1 = cf_2d_dict["s_1"]

    print("Making 4D merge")

    # pointers to frames in s,f,o,d,I,n,s1
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
    
    # make tree of (omega, dty)
    
    frames_dty = dty.ravel()
    frames_omega = omega.ravel() % 360
    om_dty_tree = scipy.spatial.cKDTree(np.transpose((frames_dty, frames_omega)))
    
    # peaks that overlap, k : 0 -> npks == len(s|f|o|d)
    # diagonal
    krow = list(range(len(o)))
    kcol = list(range(len(o)))
    # iterate over flat list of frame numbers starting from 1 (not 0)
    for i in range(1, len(n)):
        # find the tree point with 
        dist, idx_omlow = om_dty_tree.query((frames_dty[i], frames_omega[i] - ostep))
        if dist < ostep:
            flo_om = idx_omlow
            fhi_om = i
            # 1.6 is how close centers should be to overlap
            lol = trees[flo_om].query_ball_tree(trees[fhi_om], r=radius)
            for srcpk, destpks in enumerate(lol):  # dest is strictly higher than src
                for destpk in destpks:
                    krow.append(srcpk + p[flo_om])
                    kcol.append(destpk + p[fhi_om])
        # over dty
        dist, idx_dtylow = om_dty_tree.query((frames_dty[i] - ystep, frames_omega[i]))
        
        if dist < ystep:
            flo_dty = idx_dtylow
            fhi_dty = i
            lol = trees[flo_dty].query_ball_tree(trees[fhi_dty], r=radius)
            for srcpk, destpks in enumerate(lol):  # dest is strictly higher than src
                for destpk in destpks:
                    krow.append(srcpk + p[flo_dty])
                    kcol.append(destpk + p[fhi_dty])
    
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
    s4d1 = np.bincount(labels, minlength=ncomp, weights=s1)  # s_1
    s4dI = np.bincount(labels, minlength=ncomp, weights=I)  # s_I
    ssI = np.bincount(labels, minlength=ncomp, weights=I * s)  # s_sI
    sfI = np.bincount(labels, minlength=ncomp, weights=I * f)  # s_sI
    soI = np.bincount(labels, minlength=ncomp, weights=I * o)  # s_sI
    sdI = np.bincount(labels, minlength=ncomp, weights=I * d)  # s_sI
    s4d = ssI / s4dI
    f4d = sfI / s4dI
    o4d = soI / s4dI
    d4d = sdI / s4dI

    cf_4d_dict = {
        "s_raw": s4d,
        "f_raw": f4d,
        "omega": o4d,
        "dty": d4d,
        "sum_intensity": s4dI,
        "Number_of_pixels": s4d1,
    }

    cf_2d_dict["spot4d_id"] = labels
    return cf_4d_dict


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
                self.flat = fabio.open(flatfile).data.astype(np.float32)
            # normalise flat
            self.flat /= self.flat.mean()

        if darkfile is not None:
            if isinstance(darkfile, np.ndarray):
                self.dark = darkfile
            else:
                self.dark = fabio.open(darkfile).data.astype(np.float32)

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
            np.subtract(cor, self.dark, cor)
        if self.flat is not None:
            np.divide(cor, self.flat, cor)
        return cor

    def bgsub(self, img, scale_factor=None):
        """
        This attempts to remove a background by rescaling the self.bg image.

        """
        img = self.correct(img)
        if scale_factor is not None:
            img = img * scale_factor
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

    def peaksearch(self, img, omega=0, scale_factor=None):
        if self.wrk is None:
            self.wrk = np.empty(img.shape, "b")
            self.labels = np.empty(img.shape, "i")
        
        self.cor = self.bgsub(img, scale_factor=scale_factor)
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
    hname, dsetname, num, omega, worker_args, scale_factor = arg
    if pps.worker is None:
        pps.worker = worker(**worker_args)
    frm = get_dset(hname, dsetname)[num]
    pks = pps.worker.peaksearch(frm, omega=omega, scale_factor=scale_factor)
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

def segment_master_file(
    master_file,
    frames_dataset,
    omega_angles,
    worker_args,
    num_cpus=None,
    scale_factor=None,
    **process_map_kwargs
):
    """
    Collects peaks for all the frames in the first scan in the dataset 
    using parallel processing.
    """
    num_threads = num_cpus or max(1, ImageD11.cImageD11.cores_available() - 1)
    num_frames = omega_angles.shape[0]
    if scale_factor is None:
        args = [
            (master_file, frames_dataset, i, omega_angles[i], worker_args, None) 
            for i in range(num_frames)
        ]
    else:
        assert len(omega_angles) == len(scale_factor), "The shape of omega_angles and scale_factor should agree"
        args = [
            (master_file, frames_dataset, i, omega_angles[i], worker_args, scale_factor[i]) 
            for i in range(num_frames)
        ] 
    all_frames_peaks_list = process_map(
        pps, 
        args, 
        chunksize=1, 
        max_workers=num_threads, 
        **process_map_kwargs
    )
    return all_frames_peaks_list

def peaks_list_to_dict(all_frames_peaks_list):
    """
    Merges collected peak properties into structured dictionaries
    """
    peaks_2d_dict = {title: [] for title in PKSAVE}
    num_peaks = []

    for peak_index, properties in all_frames_peaks_list:
        for title, column_index in zip(PKSAVE, PKCOL):
            peaks_2d_dict[title].extend(properties[:, column_index])

        num_peaks.append(properties.shape[0])
    
    # Convert lists to numpy arrays
    for key in peaks_2d_dict:
        peaks_2d_dict[key] = np.array(peaks_2d_dict[key])
    
    return peaks_2d_dict, num_peaks

def segment_dataset(
        dataset, 
        worker_args, 
        num_cpus=None, 
        scan_number=0,
        monitor_name=None,
        monitor_ref_func=np.mean,
        **process_map_kwargs
):
    """
    Process a dataset by collecting all 2d peaks (for frames),
    merge them,
    and performing spatial correction (detector plane correction)
    Performs 4D merge if you give more than one scan number
    Returns columnfile_2d, columnfile_3d (and columnfile_4d if relevant)
    
    dataset: ImageD11.sinograms.dataset.Dataset object
    worker_args: arguments to the peaksearch worker
    num_cpus: number of CPUs to multiprocess over
    scan_number: which scan number (or numbers) in ds.scans to segment.
    monitor_name: name of a monitor to look for via ds.set_monitor(name)
    Intensities will be scaled during semgnetation if a name is provided
    monitor_ref_func: function to apply to ds.monitor to generate a reference value
    when we normalise, we multiply by ref_value_func(ds.monitor)/ds.monitor
    we suggest np.mean as an example...
    hint: if you want to return a constant, use this:
    ref_value_func=lambda x: 1e5
    """
    # Step 1: collect all peaks
    if monitor_name is not None:
        dataset.set_monitor(name=monitor_name, ref_value_func=monitor_ref_func)
        dataset.save()
        scale_factor = dataset.monitor_ref/dataset.monitor
    else:
        scale_factor=None
    # see if scan_number is an iterable
    do_4d = False
    try:
        iterator = iter(scan_number)
    except TypeError:
        # scan_number not iterable
        if scale_factor is None:
            all_frames_peaks_list = segment_master_file(
                dataset.masterfile,
                "%s/measurement/%s"%(dataset.scans[scan_number], dataset.detector),
                dataset.omega_for_bins[scan_number],
                worker_args,
                num_cpus,
                None,
                **process_map_kwargs
            )
        else:
            all_frames_peaks_list = segment_master_file(
                dataset.masterfile,
                "%s/measurement/%s"%(dataset.scans[scan_number], dataset.detector),
                dataset.omega_for_bins[scan_number],
                worker_args,
                num_cpus,
                scale_factor[scan_number],
                **process_map_kwargs
            )
        
        # Step 2: merge collected peaks data into dict
        peaks_2d_dict, num_peaks = peaks_list_to_dict(
            all_frames_peaks_list
        )
        
    else:
        # scan number is iterable
        do_4d = True
        # for loop over scan numbers
        # concatenate dicts together
        from collections import defaultdict
        merged_dict = defaultdict(list)
        
        num_peaks_all = []
        
        import time
        for sn in scan_number:
            # get all_frames_peaks_list
            
            if scale_factor is None:
                all_frames_peaks_list = segment_master_file(
                    dataset.masterfile,
                    "%s/measurement/%s"%(dataset.scans[sn], dataset.detector),
                    dataset.omega_for_bins[sn],
                    worker_args,
                    num_cpus,
                    None,
                    **process_map_kwargs
                )
            else:
                all_frames_peaks_list = segment_master_file(
                    dataset.masterfile,
                    "%s/measurement/%s"%(dataset.scans[sn], dataset.detector),
                    dataset.omega_for_bins[sn],
                    worker_args,
                    num_cpus,
                    scale_factor[sn],
                    **process_map_kwargs
                )
            
            peaks_2d_dict, num_peaks = peaks_list_to_dict(all_frames_peaks_list)

            peaks_2d_dict["dty"] = np.full_like(peaks_2d_dict["s_raw"], dataset.dty[sn][0])
            for key, arr in peaks_2d_dict.items():
                merged_dict[key].append(arr)
                
            num_peaks_all.extend(num_peaks)

        peaks_2d_dict = {key: np.concatenate(arrs) for key, arrs in merged_dict.items()}
        num_peaks = num_peaks_all

    # Step 3: Perform 3d merge
    if not do_4d:
        omega_angles = dataset.omega_for_bins[scan_number]
        peak_3d_dict = do3dmerge(peaks_2d_dict, num_peaks, omega_angles)
    
    # Step 4: Optionally perform 4d merge
    if do_4d:
        omega_angles = dataset.omega_for_bins[scan_number]
        dty_values = dataset.dty[scan_number]
        peak_4d_dict = do4dmerge(peaks_2d_dict, num_peaks, omega_angles, dty_values, dataset.ostep, dataset.ystep)

    # Step 4: Spatial Correction
    peaks_2d_dict["omega"] = peaks_2d_dict["o_raw"]
    peaks_2d_dict.pop("o_raw")
    peaks_2d_dict["sum_intensity"] = peaks_2d_dict["s_I"].copy()
    peaks_2d_dict["Number_of_pixels"] = peaks_2d_dict["s_1"].copy()
    columnfile_2d = dataset.get_colfile_from_peaks_dict(peaks_2d_dict)
    
    if do_4d:
        peak_4d_dict["spot4d_id"] = np.arange(len(peak_4d_dict["s_raw"]))
        columnfile_4d = dataset.get_colfile_from_peaks_dict(peak_4d_dict)
        return columnfile_2d, columnfile_4d
    else:
        peak_3d_dict["spot3d_id"] = np.arange(len(peak_3d_dict["s_raw"]))
        columnfile_3d = dataset.get_colfile_from_peaks_dict(peak_3d_dict)
        return columnfile_2d, columnfile_3d
