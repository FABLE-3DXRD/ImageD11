from __future__ import print_function, division

""" Do segmentation of lima/eiger files with no notion of metadata
Blocking will come via lima saving, so about 1000 frames per file

Make the code parallel over lima files ...
    minimal interprocess communication intended
 - read an input file which has the source/destination file names for this job
"""


import sys
import time
import os
import math
import logging
import numpy as np
import h5py
import hdf5plugin
import fabio
import numba
import warnings

from ImageD11 import sparseframe
from ImageD11.sinograms import dataset
import ImageD11.cImageD11


try:
    from bslz4_to_sparse import chunk2sparse
except ImportError:
    chunk2sparse = None

# Code to clean the 2D image and reduce it to a sparse array:
# things we might edit
class SegmenterOptions:
    # These are the stuff that belong to us in the hdf5 file
    # (in our group: lima_segmenter)
    jobnames = (
        "cut",
        "howmany",
        "pixels_in_spot",
        "maskfile",
        "bgfile",
        "cores_per_job",
        "files_per_core",
    )

    # There are things that DO NOT belong to us
    datasetnames = ("limapath", "analysispath", "datapath", "imagefiles", "sparsefiles")

    def __init__(
        self,
        cut=1,  # keep values abuve cut in first look at image
        howmany=100000,  # max pixels per frame to keep
        pixels_in_spot=3,
        maskfile="",
        bgfile="",
        cores_per_job=1,  # avoid buggy multiprocessing. Can't fix.
        files_per_core=8,
    ):
        self.cut = cut
        self.howmany = howmany
        self.pixels_in_spot = pixels_in_spot
        self.maskfile = maskfile
        self.mask = None
        self.bgfile = bgfile
        self.bg = None
        self.files_per_core = files_per_core
        self.cores_per_job = cores_per_job

    def __repr__(self):
        return "\n".join(
            [
                "%s:%s" % (name, getattr(self, name, None))
                for name in self.jobnames + self.datasetnames
            ]
        )

    def setup(self):
        self.thresholds = tuple([self.cut * pow(2, i) for i in range(6)])
        # validate input
        if len(self.maskfile):
            # The mask must have:
            #   0 == active pixel
            #   1 == masked pixel
            self.mask = 1 - fabio.open(self.maskfile).data.astype(np.uint8)
            assert self.mask.min() < 2
            assert self.mask.max() >= 0
        if len(self.bgfile):
            self.bg = fabio.open(self.bgfile).data

    def load(self, h5name, h5group):

        with h5py.File(h5name, "r") as hin:
            grp = hin[h5group]
            pgrp = grp.parent
            for name in self.jobnames:
                if name in grp.attrs:
                    setattr(self, name, grp.attrs.get(name))
            for name in self.datasetnames:
                #     datasetnames = ( 'limapath', 'analysispath', 'datapath',
                #                      'imagefiles', 'sparsefiles' )
                if name in pgrp.attrs:
                    data = pgrp.attrs.get(name)
                    setattr(self, name, data)
                elif name in pgrp:
                    data = pgrp[name][()]
                    if name.endswith("s"):  # plural
                        if isinstance(data, np.ndarray):
                            data = list(data)
                        if isinstance(data[0], np.ndarray) or isinstance(
                            data[0], bytes
                        ):
                            data = [x.decode() for x in data]
                    else:
                        data = str(data)
                    setattr(self, name, data)
                else:
                    logging.warning("Missing " + name)
        self.setup()

    def save(self, h5name, h5group):
        logging.info("saving to " + h5name + "::" + h5group)
        with h5py.File(h5name, "a") as hout:
            grp = hout.require_group(h5group)
            for name in self.jobnames:
                value = getattr(self, name, None)
                print(name, value)
                if value is not None:
                    grp.attrs[name] = value


@numba.njit
def select(img, mask, row, col, val, cut):
    # TODO: This is in now cImageD11.tosparse_{u16|f32}
    # Choose the pixels that are > cut and put into sparse arrays
    k = 0
    for s in range(img.shape[0]):
        for f in range(img.shape[1]):
            if img[s, f] * mask[s, f] > cut:
                row[k] = s
                col[k] = f
                val[k] = img[s, f]
                k += 1
    return k


@numba.njit
def top_pixels(nnz, row, col, val, howmany, thresholds):
    """
    selects the strongest pixels from a sparse collection
    - THRESHOLDS should be a sorted array of potential cutoff values to try
    that are higher than the original cutoff used to select data
    - howmany is the maximum number of pixels to return
    """
    # quick return if there are already few enough pixels
    if nnz <= howmany:
        return nnz
    # histogram of how many pixels are above each threshold
    h = np.zeros(len(thresholds), dtype=np.uint32)
    for k in range(nnz):
        for i, t in enumerate(thresholds):
            if val[k] > t:
                h[i] += 1
            else:
                break
    # choose the one to use. This is the first that is lower than howmany
    tcut = thresholds[-1]
    for n, t in zip(h, thresholds):
        if n < howmany:
            tcut = t
            break
    # now we filter the pixels
    n = 0
    for k in range(nnz):
        if val[k] > tcut:
            row[n] = row[k]
            col[n] = col[k]
            val[n] = val[k]
            n += 1
            if n >= howmany:
                break
    return n


OPTIONS = None  # global. Nasty.


class frmtosparse:
    def __init__(self, mask, dtype):
        # cache the mallocs on this function. Should be one per process
        self.row = np.empty(mask.size, np.uint16)
        self.col = np.empty(mask.size, np.uint16)
        self.val = np.empty(mask.size, dtype)
        self.mask = mask

    def __call__(self, frm, cut):
        nnz = select(frm, self.mask, self.row, self.col, self.val, cut)
        return nnz, self.row[:nnz], self.col[:nnz], self.val[:nnz]


def clean(nnz, row, col, val):
    #flake8: global OPTIONS
    if nnz == 0:
        return None
    if nnz > OPTIONS.howmany:
        nnz = top_pixels(nnz, row, col, val, OPTIONS.howmany, OPTIONS.thresholds)
        # Now get rid of the single pixel 'peaks'
        #   (for the mallocs, data is copied here)
        s = sparseframe.sparse_frame(
            row[:nnz].copy(), col[:nnz].copy(), OPTIONS.mask.shape
        )
        s.set_pixels("intensity", val[:nnz].copy())
    else:
        s = sparseframe.sparse_frame(row, col, OPTIONS.mask.shape)
        s.set_pixels("intensity", val)
    if OPTIONS.pixels_in_spot <= 1:
        return s
    # label them according to the connected objects
    s.set_pixels("f32", s.pixels["intensity"].astype(np.float32))
    npk = sparseframe.sparse_connected_pixels(
        s, threshold=0, data_name="f32", label_name="cp"
    )
    # only keep spots with more than 3 pixels ...
    #    mom = sparseframe.sparse_moments( s,
    #                                     intensity_name="f32",
    #                                     labels_name="cp" )
    #    npx = mom[:, cImageD11.s2D_1]
    npx = np.bincount(s.pixels["cp"], minlength=npk)
    pxcounts = npx[s.pixels["cp"]]
    pxmsk = pxcounts >= OPTIONS.pixels_in_spot
    if pxmsk.sum() == 0:
        return None
    sf = s.mask(pxmsk)
    return sf


def reader(frms, mask, cut, start=0):
    """
    iterator to read chunks or frames and segment them
    returns sparseframes
    """
    assert start < len(frms)
    if (
        (chunk2sparse is not None)
        and ("32008" in frms._filters)
        and (not frms.is_virtual)
        and (OPTIONS.bg is None)
    ):
        print("# reading compressed chunks")
        fun = chunk2sparse(mask, dtype=frms.dtype)
        for i in range(start, frms.shape[0]):
            filters, chunk = frms.id.read_direct_chunk((i, 0, 0))
            npx, row, col, val = fun.coo(chunk, cut)
            spf = clean(npx, row, col, val)
            yield spf
    else:
        fun = frmtosparse(mask, frms.dtype)
        for i in range(start, frms.shape[0]):
            frm = frms[i]
            if OPTIONS.bg is not None:
                frm = frm.astype(np.float32) - OPTIONS.bg
            npx, row, col, val = fun(frm, cut)
            spf = clean(npx, row, col, val)
            yield spf


def segment_lima(args):
    """Does segmentation on a single hdf5
    srcname,
    destname,
    dataset
    """
    srcname, destname, dataset = args
    # saving compression style:
    opts = {
        "chunks": (10000,),
        "maxshape": (None,),
        "compression": "lzf",
        "shuffle": True,
    }
    start = time.time()
    with h5py.File(destname, "a") as hout:
        with h5py.File(srcname, "r") as hin:
            if dataset not in hin:
                print("Missing", dataset, "in", srcname)
                return
            # TODO/fixme - copy some headers over
            print("# ", srcname, destname, dataset)
            print("# time now", time.ctime(), "\n#", end=" ")
            frms = hin[dataset]
            g = hout.require_group(dataset)
            row = g.create_dataset("row", (0,), dtype=np.uint16, **opts)
            col = g.create_dataset("col", (0,), dtype=np.uint16, **opts)
            # can go over 65535 frames in a scan
            # num = g.create_dataset("frame", (1,), dtype=np.uint32, **opts)
            sig = g.create_dataset("intensity", (0,), dtype=frms.dtype, **opts)
            nnz = g.create_dataset("nnz", (frms.shape[0],), dtype=np.uint32)
            g.attrs["itype"] = np.dtype(np.uint16).name
            g.attrs["nframes"] = frms.shape[0]
            g.attrs["shape0"] = frms.shape[1]
            g.attrs["shape1"] = frms.shape[2]
            npx = 0
            nframes = frms.shape[0]
            if OPTIONS.mask is None:
                # put in a dummy now that we have the frame shape
                OPTIONS.mask = np.ones((frms.shape[1], frms.shape[2]), dtype=np.uint8)
            for i, spf in enumerate(reader(frms, OPTIONS.mask, OPTIONS.cut)):
                if i % 100 == 0:
                    if spf is None:
                        print("%4d 0" % (i), end=",")
                    else:
                        print("%4d %d" % (i, spf.nnz), end=",")
                    sys.stdout.flush()
                if spf is None:
                    nnz[i] = 0
                    continue
                if spf.nnz + npx > len(row):
                    row.resize(spf.nnz + npx, axis=0)
                    col.resize(spf.nnz + npx, axis=0)
                    sig.resize(spf.nnz + npx, axis=0)
                row[npx:] = spf.row[:]
                col[npx:] = spf.col[:]
                sig[npx:] = spf.pixels["intensity"]
                # num[npx:] = i
                nnz[i] = spf.nnz
                npx += spf.nnz
            g.attrs["npx"] = npx
    end = time.time()
    print("\n# Done", nframes, "frames", npx, "pixels  fps", nframes / (end - start))
    return destname

    # the output file should be flushed and closed when this returns


OPTIONS = None


def initOptions(h5name, jobid):
    global OPTIONS
    OPTIONS = SegmenterOptions()
    OPTIONS.load(h5name, "lima_segmenter")
    OPTIONS.jobid = jobid


def main(h5name, jobid):
    #flake8: global OPTIONS
    initOptions(h5name, jobid)
    options = OPTIONS
    assert options is not None
    assert OPTIONS is not None
    args = []
    files_per_job = options.cores_per_job * options.files_per_core  # 64 files per job
    start = options.jobid * files_per_job
    end = min((options.jobid + 1) * files_per_job, len(options.imagefiles))
    for i in range(start, end):
        args.append(
            (
                os.path.join(options.datapath, options.imagefiles[i]),  # src
                os.path.join(options.analysispath, options.sparsefiles[i]),  # dest
                options.limapath,
            )
        )
    if options.cores_per_job > 1:
        ImageD11.cImageD11.check_multiprocessing(patch=True)  # avoid fork
        import multiprocessing

        with multiprocessing.Pool(
            processes=options.cores_per_job,
            initializer=initOptions,
            initargs=(h5name, jobid),
        ) as mypool:
            donefile = sys.stdout
            for fname in mypool.map(segment_lima, args, chunksize=1):
                donefile.write(fname + "\n")
                donefile.flush()
    else:
        for arg in args:
            fname = segment_lima(arg)
            print(fname)
            sys.stdout.flush()

    print("All done")


def setup_slurm_array(dsname, dsgroup="/", pythonpath=None):
    """Send the tasks to slurm"""
    dso = dataset.load(dsname, dsgroup)
    nfiles = len(dso.sparsefiles)
    dstlima = [os.path.join(dso.analysispath, name) for name in dso.sparsefiles]
    done = 0
    for d in dstlima:
        if os.path.exists(d):
            done += 1
    print("total files to process", nfiles, "done", done)
    if done == nfiles:
        return None
    sdir = os.path.join(dso.analysispath, "slurm")
    if not os.path.exists(sdir):
        os.makedirs(sdir)

    sparsefilesdir = os.path.split(dstlima[0])[0]
    if not os.path.exists(sparsefilesdir):
        os.makedirs(sparsefilesdir)

    options = SegmenterOptions()
    options.load(dsname, dsgroup + "/lima_segmenter")
    options.setup()
    if options.mask is not None:
        print("# Opened mask",
              options.maskfile,
              " %.2f %% pixels are active" % (100 * options.mask.mean()),
             )
    files_per_job = options.files_per_core * options.cores_per_job
    jobs_needed = math.ceil(nfiles / files_per_job)
    sbat = os.path.join(sdir, "lima_segmenter_slurm.sh")
    if pythonpath is None:
        cmd = sys.executable
    else:
        cmd = "PYTHONPATH=%s %s" % (pythonpath, sys.executable)
    command = (
        "%s -m ImageD11.sinograms.lima_segmenter segment %s $SLURM_ARRAY_TASK_ID"
        % (cmd, dsname)
    )
    with open(sbat, "w") as fout:
        fout.write(
            """#!/bin/bash
#SBATCH --job-name=array-lima_segmenter
#SBATCH --output=%s/lima_segmenter_%%A_%%a.out
#SBATCH --error=%s/lima_segmenter_%%A_%%a.err
#SBATCH --array=0-%d
#SBATCH --time=02:00:00
# define memory needs and number of tasks for each array job
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=%d
#
date
echo Running on $HOSTNAME : %s
OMP_NUM_THREADS=1 %s > %s/lima_segmenter_$SLURM_ARRAY_TASK_ID.log 2>&1
date
"""
            % (sdir, sdir, jobs_needed, options.cores_per_job, command, command, sdir)
        )
    logging.info("wrote " + sbat)
    return sbat

def sbatchlocal(fname, cores=None):
    """ 
    Execute a grid batch job on the local machine
    Loops over the array submission for you
    """
    import concurrent.futures
    if cores is None:
        cores = ImageD11.cImageD11.cores_available()
    lines = open(fname,'r').readlines()
    for line in lines:
        if line.find('--array=')>=0:
            start, end = line.split('=')[-1].split('-')
    commands = [ 'SLURM_ARRAY_TASK_ID=%d bash %s > %s_%d.log'%(i, fname, fname, i )
                for i in range(int(start), int(end)) ]
    with concurrent.futures.ThreadPoolExecutor(max_workers = cores) as pool:
        for _ in pool.map( os.system, commands ):
            pass

def setup(
    dsname,
    cut=None,  # keep values abuve cut in first look at image
    howmany=100000,  # max pixels per frame to keep
    pixels_in_spot=3,
    maskfile="",
    bgfile="",
    cores_per_job=1,
    files_per_core=8,
    pythonpath=None,
):
    """
    Writes options into the dataset file
    cut=None is replaced by 1 for eiger, 25 otherwise
    pythonpath -> point to a non-default install
    """
    dso = dataset.load(dsname)
    if cut is None:
        if "eiger" in dso.limapath:
            cut = 1
        else:
            cut = 25  # 5 sigma
    # Use mask from dataset if no mask provided and dataset has one
    if len(maskfile) == 0 and hasattr(dso, "maskfile"):
        maskfile = dso.maskfile
    if len(maskfile) == 0 and ("eiger" in dso.limapath):
        warnings.warn("Eiger detector needs a maskfile that is missing")
    # Use background file from dataset if no bgfile provided and dataset has one
    if len(bgfile) == 0 and hasattr(dso, "bgfile"):
        bgfile = dso.bgfile
    options = SegmenterOptions(
        cut=cut,
        howmany=howmany,
        pixels_in_spot=pixels_in_spot,
        maskfile=maskfile,
        bgfile=bgfile,
        cores_per_job=cores_per_job,
        files_per_core=files_per_core,
    )
    options.save(dsname, "lima_segmenter")
    return setup_slurm_array(dsname, pythonpath=pythonpath)


def segment():
    # Uses a hdffile from dataset.py to steer the processing
    # everything is passing via this file.
    #
    h5name = sys.argv[2]
    jobid = int(sys.argv[3])
    main(h5name, jobid)


if __name__ == "__main__":

    if sys.argv[1] == "setup":
        setup(sys.argv[2])

    if sys.argv[1] == "segment":
        segment()
