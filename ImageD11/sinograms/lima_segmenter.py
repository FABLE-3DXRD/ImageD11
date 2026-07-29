from __future__ import print_function, division

"""Do segmentation of lima/eiger files with no notion of metadata
Blocking will come via lima saving, so about 1000 frames per file

Make the code parallel over lima files ...
    minimal interprocess communication intended
 - read an input file which has the source/destination file names for this job
"""

import os

# because multiprocessing
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import logging
import math
import sys
import time
import warnings

import fabio
import h5py
import hdf5plugin
import numba
import numpy as np

import ImageD11.cImageD11
from ImageD11 import sparseframe
from ImageD11.sinograms import dataset

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
        cores_per_job=16,
        files_per_core=4,
    ):
        self.cut = cut
        self.howmany = howmany
        self.pixels_in_spot = pixels_in_spot
        self.maskfile = maskfile
        self.mask = None
        self.bgfile = bgfile
        self.bg = None
        self.bg_is_stack = False
        self.bg_dataset = None
        self.bgname = None
        self.files_per_core = files_per_core
        self.cores_per_job = cores_per_job
        print(self.bgfile)

    def __repr__(self):
        return "\n".join(
            [
                "%s:%s" % (name, getattr(self, name, None))
                for name in self.jobnames + self.datasetnames
            ]
        )

    def setup(self):
        self.thresholds = tuple([self.cut * pow(2, i) for i in range(6)])
        self.mask = None
        self.bg = None
        self.bg_is_stack = False
        self.bg_dataset = None
        self.bgname = None
        if len(self.maskfile):
            self.mask = 1 - fabio.open(self.maskfile).data.astype(np.uint8)
            assert self.mask.min() < 2
            assert self.mask.max() >= 0
        if len(self.bgfile):
            bgname = self.bgfile
            if isinstance(bgname, bytes):
                bgname = bgname.decode()
            if "::" in bgname:
                bgname, bgdataset = bgname.split("::", 1)
            else:
                bgdataset = getattr(self, "limapath", None)
            if h5py.is_hdf5(bgname):
                self.bg_is_stack = True
                self.bg_dataset = bgdataset
                self.bgname = bgname
            else:
                self.bg = fabio.open(bgname).data

    def load(self, h5name, h5group):
        with h5py.File(h5name, "r") as hin:
            grp = hin[h5group]
            pgrp = grp.parent
            for name in self.jobnames:
                if name in grp.attrs:
                    setattr(self, name, grp.attrs.get(name))
            for name in self.datasetnames:
                if name in pgrp.attrs:
                    data = pgrp.attrs.get(name)
                    setattr(self, name, data)
                elif name in pgrp:
                    data = pgrp[name][()]
                    if name.endswith("s"):
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
    if nnz <= howmany:
        return nnz
    h = np.zeros(len(thresholds), dtype=np.uint32)
    for k in range(nnz):
        for i, t in enumerate(thresholds):
            if val[k] > t:
                h[i] += 1
            else:
                break
    tcut = thresholds[-1]
    for n, t in zip(h, thresholds):
        if n < howmany:
            tcut = t
            break
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


OPTIONS = None


class frmtosparse:
    def __init__(self, mask, dtype):
        self.row = np.empty(mask.size, np.uint16)
        self.col = np.empty(mask.size, np.uint16)
        self.val = np.empty(mask.size, dtype)
        self.mask = mask

    def __call__(self, frm, cut):
        nnz = select(frm, self.mask, self.row, self.col, self.val, cut)
        return nnz, self.row[:nnz], self.col[:nnz], self.val[:nnz]


def clean(nnz, row, col, val, config_options=None):
    if config_options is None:
        options = OPTIONS
    else:
        options = config_options

    if nnz == 0:
        return None
    if nnz > options.howmany:
        nnz = top_pixels(nnz, row, col, val, options.howmany, options.thresholds)
        s = sparseframe.sparse_frame(
            row[:nnz].copy(), col[:nnz].copy(), options.mask.shape
        )
        s.set_pixels("intensity", val[:nnz].copy())
    else:
        s = sparseframe.sparse_frame(row, col, options.mask.shape)
        s.set_pixels("intensity", val)
    if options.pixels_in_spot <= 1:
        return s
    s.set_pixels("f32", s.pixels["intensity"].astype(np.float32))
    npk = sparseframe.sparse_connected_pixels(
        s, threshold=0, data_name="f32", label_name="cp"
    )
    npx = np.bincount(s.pixels["cp"], minlength=npk)
    pxcounts = npx[s.pixels["cp"]]
    pxmsk = pxcounts >= options.pixels_in_spot
    if pxmsk.sum() == 0:
        return None
    sf = s.mask(pxmsk)
    return sf


def reader(frms, mask, cut, start=0, bgfrms=None):
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
        and (bgfrms is None)
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
            if bgfrms is not None:
                frm = frm.astype(np.float32) - bgfrms[i]
            elif OPTIONS.bg is not None:
                frm = frm.astype(np.float32) - OPTIONS.bg
            npx, row, col, val = fun(frm, cut)
            spf = clean(npx, row, col, val)
            yield spf


def segment_lima(args):
    srcname, destname, dataset = args
    opts = {
        "chunks": (10000,),
        "maxshape": (None,),
        "compression": "lzf",
        "shuffle": True,
    }
    start = time.time()
    with h5py.File(destname, "a") as hout, h5py.File(srcname, "r") as hin:
        if dataset not in hin:
            print("Missing", dataset, "in", srcname)
            return
        print("# ", srcname, destname, dataset)
        print("# time now", time.ctime(), "\n#", end=" ")
        frms = hin[dataset]

        g = hout.require_group(dataset)
        row = g.create_dataset("row", (0,), dtype=np.uint16, **opts)
        col = g.create_dataset("col", (0,), dtype=np.uint16, **opts)
        sig = g.create_dataset("intensity", (0,), dtype=frms.dtype, **opts)
        nnz = g.create_dataset("nnz", (frms.shape[0],), dtype=np.uint32)
        g.attrs["itype"] = np.dtype(np.uint16).name
        g.attrs["nframes"] = frms.shape[0]
        g.attrs["shape0"] = frms.shape[1]
        g.attrs["shape1"] = frms.shape[2]
        npx = 0
        nframes = frms.shape[0]

        if OPTIONS.mask is None:
            OPTIONS.mask = np.ones((frms.shape[1], frms.shape[2]), dtype=np.uint8)

        if OPTIONS.bg_is_stack:
            with h5py.File(OPTIONS.bgname, "r") as bg_h5:
                bgdataset = OPTIONS.bg_dataset or dataset
                if bgdataset not in bg_h5:
                    raise KeyError(
                        "Missing background dataset %s in %s"
                        % (bgdataset, OPTIONS.bgname)
                    )
                bgfrms = bg_h5[bgdataset]
                if bgfrms.shape != frms.shape:
                    raise ValueError(
                        "Background stack shape %s does not match frame stack shape %s"
                        % (bgfrms.shape, frms.shape)
                    )

                for i, spf in enumerate(
                    reader(frms, OPTIONS.mask, OPTIONS.cut, bgfrms=bgfrms)
                ):
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
                    nnz[i] = spf.nnz
                    npx += spf.nnz
        else:
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
                nnz[i] = spf.nnz
                npx += spf.nnz

        g.attrs["npx"] = npx

    end = time.time()
    print("\n# Done", nframes, "frames", npx, "pixels  fps", nframes / (end - start))
    return destname


OPTIONS = None


def initOptions(h5name, jobid):
    global OPTIONS
    OPTIONS = SegmenterOptions()
    OPTIONS.load(h5name, "lima_segmenter")
    OPTIONS.jobid = jobid


def main(h5name, jobid):
    initOptions(h5name, jobid)
    options = OPTIONS
    assert options is not None
    assert OPTIONS is not None
    args = []
    files_per_job = options.cores_per_job * options.files_per_core
    start = options.jobid * files_per_job
    end = min((options.jobid + 1) * files_per_job, len(options.imagefiles))
    for i in range(start, end):
        args.append(
            (
                os.path.join(options.datapath, options.imagefiles[i]),
                os.path.join(options.analysispath, options.sparsefiles[i]),
                options.limapath,
            )
        )
    if options.cores_per_job > 1:
        import multiprocessing

        ctx_type = "spawn"
        if "linux" in sys.platform:
            ctx_type = "fork"
        ctx = multiprocessing.get_context(ctx_type)
        num_processes = min(options.cores_per_job, len(args))
        print("# Starting pool with", num_processes, " workers...", flush=True)

        mypool = ctx.Pool(
            processes=num_processes,
            initializer=initOptions,
            initargs=(h5name, jobid),
        )
        try:
            donefile = sys.stdout
            for fname in mypool.imap_unordered(segment_lima, args, chunksize=1):
                donefile.write(str(fname) + "\n")
                donefile.flush()
        except Exception as e:
            print("Main loop error: ", e, file=sys.stderr)
            mypool.terminate()
        finally:
            print("# Closing pool...", flush=True)
            mypool.close()
            mypool.join()
            print("# Pool closed.", flush=True)
    else:
        for arg in args:
            fname = segment_lima(arg)
            print(fname)
            sys.stdout.flush()
    print("# All done")


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
        print(
            "# Opened mask",
            options.maskfile,
            " %.2f %% pixels are active" % (100 * options.mask.mean()),
        )
    files_per_job = options.files_per_core * options.cores_per_job
    jobs_needed = math.ceil(nfiles / files_per_job) - 1
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
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=%d
#
date
echo Running on $HOSTNAME : %s
export TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo $TMPDIR
mkdir $TMPDIR
OMP_NUM_THREADS=1 %s > %s/lima_segmenter_$SLURM_ARRAY_TASK_ID.log 2>&1
rm -rf $TMPDIR
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
    lines = open(fname, "r").readlines()
    for line in lines:
        if line.find("--array=") >= 0:
            start, end = line.split("=")[-1].split("-")
    commands = [
        "SLURM_ARRAY_TASK_ID=%d bash %s > %s_%d.log" % (i, fname, fname, i)
        for i in range(int(start), int(end))
    ]
    with concurrent.futures.ThreadPoolExecutor(max_workers=cores) as pool:
        for _ in pool.map(os.system, commands):
            pass


def setup(
    dsname,
    cut=None,
    howmany=100000,
    pixels_in_spot=3,
    maskfile="",
    bgfile="",
    cores_per_job=16,
    files_per_core=4,
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
            cut = 25
    if len(maskfile) == 0 and hasattr(dso, "maskfile"):
        maskfile = dso.maskfile
    if len(maskfile) == 0 and ("eiger" in dso.limapath):
        warnings.warn("Eiger detector needs a maskfile that is missing")
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
    h5name = sys.argv[2]
    jobid = int(sys.argv[3])
    main(h5name, jobid)


if __name__ == "__main__":
    if sys.argv[1] == "setup":
        setup(sys.argv[2])

    if sys.argv[1] == "segment":
        segment()
