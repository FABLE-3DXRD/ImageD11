from __future__ import print_function, division

import os, h5py, numpy as np
import fast_histogram
import logging

import ImageD11.grain
import ImageD11.unitcell
import ImageD11.sinograms.properties
from ImageD11.blobcorrector import correct_cf_with_dxdyfiles, correct_cf_with_spline, correct_cf_with_h5files
from ImageD11.columnfile import colfile_from_dict

"""
TO DO: 

- debug / convert this to work on f2scan and/or fscan2d data
- treat 360 scans starting near the center (monitor/hits)
- send pyFAI style integration jobs
- mca harvesting for fluo tomo
- scalar reconstructions (counters, roi, etc)
- 1d reconstructions (fluoct xrdct)
- run peaksearch/segmentations
- peak sinograms
- linking or merging peaks across dty/omega
- indexing and reconstructing grains

"""


# POSSIBLE_DETECTOR_NAMES = ("frelon3", "eiger")


def guess_chunks(name, shape):
    if name == "omega":
        return (shape[0], 1)
    if name == "dty":
        return (1, shape[1])
    return shape


def guess_omega_step( omega, rptcut=0.02 ):
    """
    Estimates the step in omega angle data for multi-turn

    Case 1: more turns try to interleave intentionally (irrationals)
        - step should reduce smoothly as more frames are added
        - example: step = pi or sqrt(2) or 180 * (3 - np.sqrt(5)) or any irrational
        - data improve resolution 'forever'
    Case 2: more turns repeat the same angles over and over
        - step stays the same when more frames are added
        - example: 0.25, 0.3, 0.5, 1.0, etc
    Case 1&2: Prime rationals. Frames start interleaving, eventually repeat.
        - example: 17 degree step. Needs 17 turns, then starts repeating.
        - example: 1/3607 degree step. Needs a lot of turns, then starts repeating.

    rptcut = tolerance to decide if frames are the same angle.
             A fraction of the largest step found.
    """
    omega = np.asarray(omega)
    v = (omega - omega.min()) % 360    # Values mod 360 to put frames in order
    v.sort()           # Adjacent frames in angle
    dv = v[1:]-v[:-1]  # Step from one frame to the next
    # If we have (many) intentional repeats the pattern is
    # 0,0,0,step,0,0,0,step,0,0,0,step,...
    # The max is the largest step we might want to use
    # The min may be zero
    guess = dv[dv>(dv.max()*rptcut)].mean()
    # print('mx, avg',dv.max(), dv.mean(), guess)
    return guess

class DataSet:
    """One DataSet instance per detector!"""

    # simple strings or ints
    ATTRNAMES = (
        "dataroot",
        "analysisroot",
        "sample",
        "dset",
        "shape",
        "dsname",
        "datapath",
        "analysispath",
        "masterfile",
        "limapath",
        "detector",
        "omegamotor",
        "dtymotor",
        "monitorname",
        "monitor_ref",
        "pksfile",
        "sparsefile",
        "parfile",
        "e2dxfile",
        "e2dyfile",
        "detectorh5",
        "splinefile",
        "maskfile",
        "bgfile",
        "pksfile",
        "col4dfile",
        "col3dfile",
        "col2dfile",
        "grainsfile",
        "sparsefile",
        "icolfile",
        "pbpfile",
    )
    STRINGLISTS = ("scans", "imagefiles", "sparsefiles")
    # sinograms
    NDNAMES = (
        "omega",
        "omega_for_bins",
        "dty",
        "nnz",
        "frames_per_file",
        "nlm",
        "frames_per_scan",
        "monitor",
        "ybinedges", "ybincens",
        "obinedges", "obincens"
    )

    def __init__(
        self,
        dataroot=".",
        analysisroot=".",
        sample="sample",
        dset="dataset",
        detector="eiger",
        omegamotor="rot_center",
        dtymotor="dty",
        filename=None,
        analysispath=None,
    ):
        """The things we need to know to process data"""

        # defaults to eiger and nanoscope station, can be overwritten with init parameters detector, omegamotor and dtymotor

        self.detector = detector  # frelon3
        self.limapath = None  # where is the data in the Lima files

        self.omegamotor = omegamotor  # diffrz
        self.dtymotor = dtymotor  # diffty

        self.dataroot = dataroot  # folder to find {sample}/{sample}_{dset}
        self.analysisroot = analysisroot  # where to write or find sparse data
        self.sample = sample  # from bliss path
        self.dset = dset  # from bliss path

        self.dsname = "_".join((self.sample, self.dset))

        # paths for raw data

        self.datapath = os.path.join(self.dataroot, self.sample, self.dsname)
        self.masterfile = os.path.join(self.datapath, self.dsname + ".h5")

        # These are in order ! The order of the lists is important - all things should match.
        self.scans = None  # read from master or load from analysis
        self.frames_per_scan = (
            None  # how many frames (and motor positions) in each scan row.
        )
        self.imagefiles = None  # List of strings. w.r.t self.datapath
        self.frames_per_file = None  # how many frames in this file (Lima files)
        self.sparsefiles = None  # maps sparse files to self.imagefiles

        self.shape = (0, 0)
        self.omega = None
        self.dty = None
        self.monitor = None
        self.monitorname = None
        self.monitor_ref = None
        self.ybinedges = None
        self.ybincens = None
        self.y0 = None
        self.obinedges = None
        self.obincens = None

        self._peaks_table = None
        self._pk2d = None
        self._pk4d = None

        self.dsfile = None

        # paths for processed data
        self.analysispath = None # default
        # Loaded
        if filename is not None:
            self.dsfile = filename
            self.load(filename)
        # Supplied overwrites
        if analysispath is not None: 
            # Can be loaded with the dataset 
            self.analysispath = analysispath

        self.update_paths()

    def update_paths(self, force=False, verbose = False):
        # paths for processed data
        # root of analysis for this dataset for this sample:
        self.analysispath_default = os.path.join(
            self.analysisroot, self.sample, self.dsname
        )
        if self.analysispath is None:
            self.analysispath = self.analysispath_default

        self.dsfile_default = os.path.join(
            self.analysispath, self.dsname + "_dataset.h5"
        )
        # at the moment, set self.dsfile to be the default
        # if save or load is ever called, this will be replaced
        if self.dsfile is None:
            self.dsfile = self.dsfile_default
        # They should be saved / loaded with the dataset.
        for name, extn in [
            ("pksfile", "_peaks_table.h5"),
            ("col4dfile", "_peaks_4d.h5"),
            ("col3dfile", "_peaks_3d.h5"),
            ("col2dfile", "_peaks_2d.h5"),
            ("grainsfile", "_grains.h5"),
            ("sparsefile", "_sparse.h5"),
            ("icolfile", "_icolf.h5"),  # subset peaks selected for indexing (pbp)
            ("pbpfile", "_pbp.txt"),  # point by point raw output
            ("refmanfile", "_refine_manager.h5"),  # PBPRefine object for refinement
            ("refpeaksfile", "_refine_peaks.h5"),  # icolf for refinement
            ("refmapfile", "_refine_map_in.h5"),  # input pbp map for refinement
            ("refoutfile", "_refine_map_out.h5"),  # output pbp map from refinement
        ]:
            # If the user has got a different name (via loading or edit), we keep that
            if (getattr(self, name, None) is None) or force:
                # Otherwise, these are the defaults.
                setattr(self, name, os.path.join(self.analysispath, self.dsname + extn))
                if verbose:
                    print('updated', getattr( self, name, None ) )
            else:
                if verbose:
                    print('not updated', getattr( self, name, None ) )

    def __repr__(self):
        r = []
        for name in "dataroot analysisroot sample dset".split():
            r.append('%s = "%s"' % (name, getattr(self, name)))
        r.append("shape = ( %d, %d)" % tuple(self.shape))
        if self.scans is not None:
            r.append(
                "# scans %d from %s to %s"
                % (len(self.scans), self.scans[0], self.scans[-1])
            )
        return "\n".join(r)

    def compare(self, other):
        """Try to see if the load/save is working"""
        from types import FunctionType

        sattrs = set([name for name in vars(self) if name[0] != "_"])
        oattrs = set([name for name in vars(self) if name[0] != "_"])
        if sattrs != oattrs:
            logging.info("Attribute mismatch " + str(sattrs) + " != " + str(oattrs))
            return False
        for a in sattrs:
            s = getattr(self, a)
            if isinstance(s, FunctionType):
                continue
            o = getattr(other, a)
            t = type(s)
            if type(o) != type(s):
                logging.info("Type mismatch %s %s" % (str(t), str(a)))
                return False
            if t == np.ndarray:
                if s.shape != o.shape:
                    logging.info("Shape mismatch %s %s" % (str(s.shape), str(o.shape)))
                    return False
                if (s != o).all():
                    logging.info("Data mismatch " + str(a))
                    return False
            else:
                if s != o:
                    logging.info("Data mismatch ")
                    return False
        logging.info("Dataset objects seem to match!")
        return True

    def report(self):
        print(self)
        print("# Collected %d missing %d" % (self.check_images()))
        print("# Segmented %d missing %d" % (self.check_sparse()))

    def import_all(self, scans=None, shape=None):
        # collect the data
        self.import_scans(scans=scans)
        # lima frames
        self.import_imagefiles()
        # motor positions
        self.import_motors_from_master()
        self.guess_shape()
        self.guessbins()
        # pixels per frame
        try:
            self.import_nnz()
        except:
            logging.info("nnz not available. Segmentation done?")

    def import_from_sparse(self, hname, scans=None, shape=None):
        """
        hname = hdf5 file containing sparse pixels (and motors)
        dataset = a dataset instance to import into
        scans = defaults to reading all "%d.1" scans in the file
                give a list to read in some other order or a subset
        """
        self.sparsefile = hname
        if scans is None:
            with h5py.File(hname, "r") as hin:
                # Read all in numerical order
                scans = list(hin["/"])
                order = np.argsort([float(v) for v in scans if v.endswith(".1")])
                self.scans = [scans[i] for i in order]
        else:
            self.scans = scans
        self.masterfile = hname  # hacky, motors come from the sparsefile
        self.import_nnz_from_sparse()  # must exist
        self.import_motors_from_master()
        # self.guess_shape() # fails with sparse
        if shape is not None:
            self.shape = shape
            self.nnz = np.array(self.nnz).reshape(self.shape)
        else:
            self.shape = self.nnz.shape
        self.omega = np.array(self.omega).reshape(self.shape)
        self.dty = np.array(self.dty).reshape(self.shape)
        if len(scans) == 1 and self.shape[0]>1:
            file_nums = np.arange(self.shape[0]*self.shape[1]).reshape(self.shape)
            self.scans = [
                "%s::[%d:%d]" % (self.scans[0], row[0], row[-1] + 1)
                for row in file_nums
            ]

        self.guessbins()

    def import_scans(self, scans=None, hname=None):
        """Reads in the scans from the bliss master file"""
        # we need to work out what detector we have at this point
        # self.guess_detector()
        if hname is None:
            hname = self.masterfile
        frames_per_scan = []
        with h5py.File(hname, "r") as hin:
            if scans is None:
                scans = [
                    scan
                    for scan in list(hin["/"])
                    if (
                        scan.endswith(".1")
                        and ("measurement" in hin[scan])
                        and (self.detector in hin[scan]["measurement"])
                        and (self.omegamotor in hin[scan]["measurement"])
                    )
                ]
            goodscans = []
            for scan in scans:
                # Make sure that this scan has a measurement from our detector
                if self.detector not in hin[scan]["measurement"]:
                    print("Bad scan", scan)
                else:
                    try:
                        frames = hin[scan]["measurement"][self.detector]
                    except KeyError as e:  # Thrown by h5py
                        print("Bad scan", scan, ", h5py error follows:")
                        print(e)
                        continue
                    if len(frames.shape) == 3:  # need 1D series of frames
                        goodscans.append(scan)
                        frames_per_scan.append(frames.shape[0])
                    else:
                        print("Bad scan", scan)

        self.scans = goodscans
        self.frames_per_scan = frames_per_scan

        logging.info("imported %d scans from %s" % (len(self.scans), hname))
        return self.scans

    def import_imagefiles(self):
        """Get the Lima file names from the bliss master file, also scan_npoints"""
        # self.import_scans() should always be called before this function, so we know the detector
        self.imagefiles = []
        self.frames_per_file = []
        with h5py.File(self.masterfile, "r") as hin:
            bad = []
            for i, scan in enumerate(self.scans):
                if ("measurement" not in hin[scan]) or (
                    self.detector not in hin[scan]["measurement"]
                ):
                    print("Bad scan", scan)
                    bad.append(scan)
                    continue
                frames = hin[scan]["measurement"][self.detector]
                self.imageshape = frames.shape[1:]
                for vsrc in frames.virtual_sources():
                    self.imagefiles.append(vsrc.file_name)
                    self.frames_per_file.append(
                        vsrc.src_space.shape[0]
                    )  # not sure about this
                    # check limapath
                    if self.limapath is None:
                        self.limapath = vsrc.dset_name
                    assert self.limapath == vsrc.dset_name
        self.frames_per_file = np.array(self.frames_per_file, int)
        self.sparsefiles = [
            os.path.join(
                "sparsefiles", name.replace("/", "_").replace(".h5", "_sparse.h5")
            )
            for name in self.imagefiles
        ]
        logging.info("imported %d lima filenames" % (np.sum(self.frames_per_file)))

    def import_motors_from_master(self):
        """read the motors from the lima file
        you need to import the imagefiles first
        these will be the motor positions to accompany the images
        # could also get these from sparse files if saved
        """
        # self.guess_motornames()
        self.omega = [
            None,
        ] * len(self.scans)
        self.dty = [
            None,
        ] * len(self.scans)
        with h5py.File(self.masterfile, "r") as hin:
            bad = []
            for i, scan in enumerate(self.scans):
                # Should always be there, if not, filter scans before you get to here
                om = hin[scan]["measurement"][self.omegamotor][()]
                if len(om) == self.frames_per_scan[i]:
                    self.omega[i] = om
                else:  # hope the first point was good ? Probably corrupted MUSST data.
                    self.omega[i] = [
                        om[0],
                    ]
                    bad.append(i)
                # this can be an array or a scalar
                # read from h5:
                dty = hin[scan]["instrument/positioners"][self.dtymotor]
                if len(dty.shape) == 0:
                    self.dty[i] = np.full(self.frames_per_scan[i], dty[()])
                elif dty.shape[0] == self.frames_per_scan[i]:
                    self.dty[i] = dty[:]
                else:
                    # corrupted MUSST?
                    self.dty[i] = np.full(self.frames_per_scan[i], dty[0])
        for b in bad:
            dom = [
                (abs(self.omega[i][0] - self.omega[b])[0], i)  # always length-1 arrays, take first element
                for i in range(len(self.scans))
                if i not in bad
            ]
            # dom is a list of tuples of (first omega value, i)
            # make it into an array
            dom = np.array(dom)
            
            if len(dom) > 0:
                j = int(dom[np.argmin(dom[:,0])][1])  # get argmin of omega column of dom, go there, then take the corresponding i
                self.omega[b] = self.omega[j]  # best match
                print(
                    "replace bad scan omega", b, self.scans[b], "with", j, self.scans[j]
                )
        logging.info("imported omega/dty")

    def guess_shape(self):
        npts = np.sum( self.frames_per_scan )
        if os.path.exists(self.masterfile):
            # strip [i::j] from self.scans if already there:
            seen = set()
            scans = []
            snames = [ s.split('::')[0] for s in self.scans ]
            for s in snames:
                if s not in seen:
                    seen.add( s )
                    scans.append( s )
            # number of turns
            rotations = []
            for i, scan in enumerate(scans):
                with h5py.File(self.masterfile, "r") as hin:
                    s = hin[scan]
                    title = s["title"].asstr()[()]
                    # print("Scan title", title)
                    if title.split()[0] == "fscan2d":
                        s0 = s["instrument/fscan_parameters/slow_npoints"][()]
                        s1 = s["instrument/fscan_parameters/fast_npoints"][()]
                        if s0 > 1:
                            file_nums = np.arange(s0 * s1).reshape((s0, s1))
                            # slice notation means last frame+1 to be inclusive
                            rotations += [
                                "%s::[%d:%d]" % (scan, row[0], row[-1] + 1)
                                for row in file_nums
                                ]
                        else:
                            rotations += [ scan, ]
                    elif title.split()[0] == "f2scan":
                        # good luck ? Assuming rotation was the inner loop here:
                        step = s["instrument/fscan_parameters/step_size"][()]
                        s1 = int(np.round(360 / step))
                        s0 = self.frames_per_scan[i] // s1
                        # logging.warning("Dataset might need to be reshaped")
                        if s0 > 1:
                            file_nums = np.arange( s0 * s1 ).reshape((s0, s1))
                            if (s1 * s0) != self.frames_per_scan[i]:
                                logging.warning( 'scan %s problem in guessing f2scan shape s1 = %d s0 = %s nframes = %d'%(
                                    scan, s1, s0, self.frames_per_scan[i] ) )
                            rotations += [
                                "%s::[%d:%d]" % (scan, row[0], row[-1] + 1)
                                for row in file_nums
                                ]
                        else:
                            rotations += [ scan, ]
                    else:
                        s0 = 1
                        s1 = npts
                        rotations.append( scan )
            self.scans = rotations
        if len(self.scans) >= 1:
            s0 = len(self.scans)
            s1 = npts // s0
        else:
            # no scans
            s0 = 0
            s1 = 0
        self.shape = s0, s1
        if np.prod(self.shape) != npts:
                print("Warning: irregular scan - might be bugs in here")
                print(npts, len(self.scans))
        self.omega = np.array(self.omega).reshape(self.shape)
        self.dty = np.array(self.dty).reshape(self.shape)
        logging.info(
                "sinogram shape = ( %d , %d ) imageshape = ( %d , %d)"
                % (self.shape[0], self.shape[1], self.imageshape[0], self.imageshape[1])
            )

    def guessbins(self):
        """
        Attempts to estimate the step size in the data by looking at the numbers
        in self.omega and self.dty that should already have self.shape reflecting
        the length of the individual scans.

        Perhaps this is the wrong approach. But we don't have it in the bliss data.

        The data might not be on a regular grid.
        """
        ny, nomega = self.shape
        if self.obincens is None:
            self.omin = self.omega.min()
            self.omax = self.omega.max()
            if (self.omax - self.omin) > 360:
                # multi-turn scan...
                self.omin = 0.0
                self.omax = 360.0
                self.omega_for_bins = self.omega % 360
                # assume the first scan is representative
                # if you have different steps in different scans ... that is bad
                self.ostep = guess_omega_step( self.omega[0] )
                #                                               include endpoint
                self.obincens = np.arange(self.omin, self.omax + self.ostep*0.1, self.ostep)
                if self.obinedges is None:
                    self.obinedges = np.arange(
                       self.omin - self.ostep / 2, self.omax + self.ostep / 1.9, self.ostep
                    )
            else:
                self.omega_for_bins = self.omega
                self.ostep = (self.omax - self.omin) / (nomega - 1)
                self.obincens = np.linspace(self.omin, self.omax, nomega)
        else: # self.obincens was loaded
            self.omin = self.obincens[0]
            self.omax = self.obincens[-1]
            self.ostep = np.mean(self.obincens[1:] - self.obincens[:-1])
            if (self.omax - self.omin)>=360:
                self.omega_for_bins = self.omega % 360
            else:
                self.omega_for_bins = self.omega
        if self.obinedges is None: # catches last 3 else here.
            self.obinedges = np.linspace(
               self.omin - self.ostep / 2, self.omax + self.ostep / 2, nomega + 1
            )
        # values 0, 1, 2
        # shape = 3
        # step = 1
        if self.ybincens is not None:
            self.ymin = self.ybincens[0]
            self.ymax = self.ybincens[-1]
        else:
            self.ymin = self.dty.min()
            self.ymax = self.dty.max()
            self.ybincens = np.linspace(self.ymin, self.ymax, ny)
        if ny > 1:
            self.ystep = (self.ymax - self.ymin) / (ny - 1)
        else:
            self.ystep = 1
        if self.ybinedges is None:
            self.ybinedges = np.linspace(
                self.ymin - self.ystep / 2, self.ymax + self.ystep / 2, ny + 1
            )
        # try to estimate y0 from ybins. Maybe not the best approach?
        # Rotation centre position should be somewhere in the master hdf5 file
        if self.y0 is None:  
            self.y0 = np.round((self.ymin + self.ymax) / 2, 6)
                  
    def correct_bins_for_half_scan(self, y0 = 0.0):
        """
        Pad self.ybincens / self.ybinedges around the bin nearest to y0
        so that the dataset becomes symmetric. 
        The original measured bins are never moved; only virtual bins are added
        on whichever side is shorter. A boolean mask records which bins are real.
 
        Sets self.ybin_real_mask : bool array on self.ybincens, False on virtual bins
        """
        ystep = self.ystep
        # Recover original measured bins
        if hasattr(self, 'ybin_real_mask'):
            yc_orig = np.asarray(self.ybincens, dtype=float)[self.ybin_real_mask]
        else:
            yc_orig = np.asarray(self.ybincens, dtype=float).copy()
        # get the bin closest to y0
        central_bin    = int(np.argmin(np.abs(yc_orig - y0)))
        central_val    = yc_orig[central_bin]
        # Distance from central_val to each end of the measured range
        lo_dist = central_val - yc_orig[0]
        hi_dist = yc_orig[-1] - central_val
        half    = np.ceil(max(lo_dist, hi_dist) / ystep) * ystep
        # Build symmetric grid centred on central_val, with the same ystep
        n_half    = int(round(half / ystep))
        new_cens  = central_val + np.arange(-n_half, n_half + 1) * ystep
        new_edges = central_val - ystep / 2 + np.arange(-n_half, n_half + 2) * ystep
        # Real-bin mask: True where new_cens matches an orig_cens within ystep/4
        real_mask = np.array([
            np.any(np.abs(yc_orig - cv) < ystep / 4)
            for cv in new_cens])
        # Update bins
        self.ybincens       = new_cens
        self.ybinedges      = new_edges
        self.ybin_real_mask = real_mask
        n_pad = int((~real_mask).sum())
        print(
        "[correct_bins_for_half_scan]  y0 = {:.4f}, central bin value = {:.4f}, central bin id: {}, "
        "halfrange={:.4f}, n_bins={} "
        "({} real + {} padded).".format(y0, central_val, central_bin, half,
                                        len(new_cens), real_mask.sum(), n_pad))

    # ------------------------------------------------------------------
    # Peak-pairing helpers: sort dataset by dty scans / dty-omega frames
    # ------------------------------------------------------------------
    def get_dty_pairs(self):
        """
        Find mirror pairs of dty scans on each side of the central y0:
            scan hi :  dty = y0 + Δy
            scan lo :  dty = y0 - Δy

        Sets self.dty_pairs: list of Pair namedtuples ('Pair', ['yi_hi', 'yi_lo', 'dty_hi', 'dty_lo'])
        where yi_* are indices into self.ybincens.
        """
        from collections import namedtuple
        y0 = self.y0
        if y0 != 0:
            self.correct_bins_for_half_scan(y0)
 
        yc       = self.ybincens.copy()
        n_y      = len(yc)
        yi0      = n_y // 2
 
        Pair = namedtuple('Pair', ['yi_hi', 'yi_lo', 'dty_hi', 'dty_lo'])
 
        def _make_pair(yi_hi, yi_lo):
            return Pair(
                yi_hi  = yi_hi,
                yi_lo  = yi_lo,
                dty_hi = float(yc[yi_hi]),
                dty_lo = float(yc[yi_lo]))
 
        # If odd number of steps (odd n_y), central bin is paired with itself; then mirror pairs outward
        pairs = []
        for yi in range(0, yi0 + 1):
            pairs.append(_make_pair(yi0 + yi, yi0 - yi))
 
        self.dty_pairs = pairs
        print("[get_dty_pairs] {n} pairs of dty scans (y0={y0})".format(
            n=len(pairs), y0=y0))
 
    def get_omega_dty_pairs(self):
        """
        Find mirror pairs of frames (A, B) in omega and dty:
            frame A :  dty = y0 + Δy ,  omega = ω
            frame B :  dty = y0 - Δy ,  omega = (ω + 180) mod 360
 
        Sets self.omega_dty_pairs: list of Pair namedtuples
        ('Pair', ['yi_A','yi_B','omi_A','omi_B', 'dty_A','dty_B','omega_A','omega_B'])
        with indices yi*, omi* into self.ybincens and self.obincens
        """
        from collections import namedtuple
        y0 = self.y0
        if y0 != 0:
            self.correct_bins_for_half_scan(y0)
 
        yc       = self.ybincens.copy()
        oc       = self.obincens.copy()
        n_y      = len(yc)
        n_o      = len(oc)
        k        = n_o // 2          # index shift for a 180-deg omega step
 
        # Sanity checks
        if self.omax - self.omin < 181:
            raise ValueError(
                "Half-rotation acquisition: cannot pair frames in omega.")
        if n_o % 2 != 0:
            raise ValueError(
                "omega_bin_centers must have an even number of bins, "
                "got {n}. Switch to dty-scan-wise pairing.".format(n=n_o))
        if not np.allclose((oc[k:] - oc[:k]) % 360, 180.0, atol=self.ostep / 2):
            raise ValueError(
                "self.obincens[i + k] is not self.obincens[i] + 180 deg for all i. "
                "Check that omega covers exactly 360 deg with uniform spacing.")
 
        Pair = namedtuple('Pair', ['yi_A', 'yi_B', 'omi_A', 'omi_B', 'dty_A', 'dty_B', 'omega_A', 'omega_B'])
 
        def _make_pair(oi_A, yi_A, oi_B, yi_B):
            return Pair(
                yi_A    = yi_A,
                yi_B    = yi_B,
                omi_A   = oi_A,
                omi_B   = oi_B,
                dty_A   = float(yc[yi_A]),
                dty_B   = float(yc[yi_B]),
                omega_A = float(oc[oi_A]),
                omega_B = float(oc[oi_B]))
 
        pairs = []
 
        if n_y % 2 == 1:
            # Odd number of dty steps → central bin exists
            yi0 = n_y // 2
            # Off-centre bins: pair mirror scans (yi0+yi, yi0-yi)
            for yi in range(1, yi0 + 1):
                yi_hi = yi0 + yi
                yi_lo = yi0 - yi
                for oi in range(n_o):
                    oi_lo = (oi + k) % n_o
                    pairs.append(_make_pair(oi, yi_hi, oi_lo, yi_lo))
 
            # Central bin: pair lower-omega half with upper half (same yi)
            for oi in range(k):
                pairs.append(_make_pair(oi, yi0, oi + k, yi0))
 
        else:
            # Even number of dty steps → no central bin
            half = n_y // 2
            for yi_hi in range(half, n_y):
                yi_lo = n_y - 1 - yi_hi
                for oi in range(n_o):
                    oi_lo = (oi + k) % n_o
                    pairs.append(_make_pair(oi, yi_hi, oi_lo, yi_lo))
 
        self.omega_dty_pairs = pairs
 
        n_y0  = sum(1 for p in pairs if p.yi_A == p.yi_B)
        n_off = len(pairs) - n_y0
        print(
            "[get_omega_dty_pairs] {ntot} pairs of 2D frames (y0={y0}): "
            "{ny0} at y0 (omega-only), {noff} off-y0.".format(
                ntot=len(pairs), y0=y0, ny0=n_y0, noff=n_off))
  
    def sort_by_dty(self, cf):
        """
        Sort columnfile cf by dty and build an index table for quick peak selection by dty scan.
 
        Sets ds.dty_scans_index : dict { yi : (slice(lo, hi), n_peaks, total_intensity) }
        yi is a dty-scan index matching those in self.dty_pairs. 
        Scans with zero peaks are absent from the dict.
        """
        yc  = self.ybincens
        ybe = self.ybinedges
 
        # Remove any stale index tables
        for attr in ('dty_scans_index', 'frames_index'):
            if hasattr(self, attr):
                delattr(self, attr)
 
        # Assign each peak to its dty bin and sort
        dtyi  = np.clip(
            np.searchsorted(ybe, cf.dty, side="right") - 1, 0, len(yc) - 1)
        order = np.argsort(dtyi)
        cf.reorder(order)
        cf.sortedby = "dty"
 
        dtyi_s = dtyi[order]
        _, first_idx, counts = np.unique(
            dtyi_s, return_index=True, return_counts=True)
 
        dty_scans_index = {}
        for lo, cnt in zip(first_idx, counts):
            key  = int(dtyi_s[lo])
            slc  = slice(int(lo), int(lo + cnt))
            Itot = cf.sum_intensity[slc].sum()
            dty_scans_index[key] = (slc, cnt, Itot)
 
        self.dty_scans_index = dty_scans_index
 
        n_empty = len(yc) - len(dty_scans_index)
        print(
            "[sort_by_dty]  {npk} peaks sorted; "
            "{nne} non-empty dty scans ({ne} scans have no peaks).".format(
                npk=cf.nrows, nne=len(dty_scans_index), ne=n_empty))

    def sort_by_sinogram(self, cf):
        """
        Sort columnfile cf in sinogram space (omega primary, dty secondary)
        and build a frame-level index table for quick peak selection.

        Sets self.frames_index : dict { (omi, yi) : (slice(lo, hi), n_peaks, total_intensity) }
        Keys match those in self.omega_dty_pairs.
        Frames with zero peaks are absent from the dict.
        """
        oc  = self.obincens
        yc  = self.ybincens
        obe = self.obinedges
        ybe = self.ybinedges
 
        # Remove any stale index tables
        for attr in ('dty_scans_index', 'frames_index'):
            if hasattr(self, attr):
                delattr(self, attr)
 
        # Assign each peak to its (omega, dty) bin
        i_omega = np.clip(
            np.searchsorted(obe, cf.omega, side="right") - 1, 0, len(oc) - 1)
        j_dty   = np.clip(
            np.searchsorted(ybe, cf.dty,   side="right") - 1, 0, len(yc) - 1)
 
        # Sort in-place: omega primary, dty secondary
        order = np.lexsort((j_dty, i_omega))
        cf.reorder(order)
        cf.sortedby = "omega_dty"
 
        i_omega_s = i_omega[order]
        j_dty_s   = j_dty[order]
 
        bin_pairs = np.column_stack([i_omega_s, j_dty_s])
        _, first_idx, counts = np.unique(
            bin_pairs, axis=0, return_index=True, return_counts=True)
 
        frames_index = {}
        for lo, cnt in zip(first_idx, counts):
            key  = (int(i_omega_s[lo]), int(j_dty_s[lo]))
            slc  = slice(int(lo), int(lo + cnt))
            Itot = cf.sum_intensity[slc].sum()
            frames_index[key] = (slc, cnt, Itot)
 
        self.frames_index = frames_index
 
        n_total = len(oc) * len(yc)
        n_empty = n_total - len(frames_index)
        print(
            "[sort_by_sinogram]  {npk} peaks sorted; "
            "{nne} non-empty frames indexed "
            "({ne} frames have no peaks).".format(
                npk=cf.nrows, nne=len(frames_index), ne=n_empty))
        
    def valid_pairs(self):
        """
        Return a boolean array marking pairs where both members contain
        at least one peak in the current index table (frames_index or dty_scans_index)
        Automatically detects which index table is present in self.
        """
        index_table = next(
            (a for a in ('frames_index', 'dty_scans_index')
             if hasattr(self, a)), None)
        if index_table is None:
            raise RuntimeError(
                "No index table found in dataset. "
                "Run sort_by_sinogram() or sort_by_dty() first.")
 
        pair_attr = 'omega_dty_pairs' if 'frames' in index_table else 'dty_pairs'
        pairs     = getattr(self, pair_attr, None)
        idx_dict  = getattr(self, index_table)
 
        valid = np.zeros(len(pairs), dtype=bool)
 
        for i, p in enumerate(pairs):
            if 'omega' in pair_attr:
                key_A = (p.omi_A, p.yi_A)
                key_B = (p.omi_B, p.yi_B)
            else:
                key_A = p.yi_hi
                key_B = p.yi_lo
            if key_A in idx_dict and key_B in idx_dict:
                valid[i] = True
 
        print(
            "[valid_pairs]  {ntot} pairs checked in {pl}\n"
            "  found {nv} pairs where both members have peaks "
            "({pct:.2f}% valid).".format(
                ntot=len(valid), pl=pair_attr,
                nv=valid.sum(),
                pct=valid.sum() / len(valid) * 100))
 
        return valid

    def check_dty_symmetry(self, cf, saveplot=False):
        """
        Plot N_peaks and total_intensity vs. abs(dty) for mirror dty_scans
        (y0+Δy, y0-Δy) and check their correlation.
 
        If alignment is correct, values should match between mirror scans.
        Significant mismatch suggests an incorrect y0, sample movement, or
        a beam issue during scanning.
        """
        import matplotlib.pyplot as plt

        if not hasattr(self, 'dty_pairs'):
            self.get_dty_pairs(self.y0)
        if not hasattr(self, 'dty_scans_index'):
            self.sort_by_dty(cf)
        
        valid = self.valid_pairs()
        
        # paired dty scans on the hi and lo side
        scan_id_hi = np.array( [p.yi_hi for p in self.dty_pairs] )
        scan_id_lo = np.array( [p.yi_lo for p in self.dty_pairs] )

        # x-axis: dty position (flipped for lo side)
        x_hi = self.ybincens[scan_id_hi]
        x_lo = -self.ybincens[scan_id_lo]
        
        # y-axis: n_peaks and log(sum_intensity) + ratios hi/lo 
        npks_hi = [self.dty_scans_index[yi][1] if v else np.nan for yi,v in zip(scan_id_hi, valid)]
        npks_lo = [self.dty_scans_index[yi][1] if v else np.nan for yi,v in zip(scan_id_lo, valid)]
        sumI_hi = [self.dty_scans_index[yi][2] if v else np.nan for yi,v in zip(scan_id_hi, valid)]
        sumI_lo = [self.dty_scans_index[yi][2] if v else np.nan for yi,v in zip(scan_id_lo, valid)]
        log_sumI_hi = np.log10(sumI_hi)
        log_sumI_lo = np.log10(sumI_lo)
        ratio_npks  = np.array(npks_hi) / np.array(npks_lo)
        ratio_sumI  = log_sumI_hi / log_sumI_lo
        # correlation coeff between hi and lo sides 
        corr_npks = np.corrcoef(np.array(npks_hi)[valid],
                                np.array(npks_lo)[valid])[0, 1]
        corr_sumI = np.corrcoef(log_sumI_hi[valid],
                                log_sumI_lo[valid])[0, 1]

        # plot figure
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
        axes = axes.flatten()
 
        # Panel 1: n_peaks per scan
        axes[0].plot(x_hi[1:], npks_hi[1:], '.-', label='hi-side (dty > 0)')
        axes[0].plot(x_lo[1:], npks_lo[1:], '.-', label='lo-side (dty < 0)')
        axes[0].set_ylabel('npks per scan')
        axes[0].set_title('Number of peaks per scan')
        axes[0].legend()
        axes[0].text(0.05, 0.9, 'Corr: {:.3f}'.format(corr_npks),
                     transform=axes[0].transAxes)
 
        # Panel 2: log(total intensity) per scan
        axes[1].plot(x_hi[1:], log_sumI_hi[1:], '.-')
        axes[1].plot(x_lo,     log_sumI_lo,      '.-')
        axes[1].set_title('Total intensity per scan')
        axes[1].set_ylabel('log10(sumI) per scan')
        axes[1].text(0.05, 0.9, 'Corr: {:.3f}'.format(corr_sumI),
                     transform=axes[1].transAxes)
 
        # Panel 3: ratio of n_peaks (hi/lo)
        axes[2].plot(x_hi[1:], ratio_npks[1:], 'k.-')
        axes[2].axhline(1, color='gray', linestyle='--')
        axes[2].set_ylabel('Ratio npks (hi/lo)')
        axes[2].set_xlabel('dty (flipped for lo)')
 
        # Panel 4: ratio of log(sumI) (hi/lo)
        axes[3].plot(x_hi[1:], ratio_sumI[1:], 'k.-')
        axes[3].axhline(1, color='gray', linestyle='--')
        axes[3].set_ylabel('Ratio log10(sumI) (hi/lo)')
        axes[3].set_xlabel('dty (flipped for lo)')
 
        fig.suptitle('dty alignment - ' + str(self.dsname), fontsize=14)
        plt.tight_layout()
 
        if saveplot:
            fname = os.path.join(
                self.analysispath, self.dsname + '_dty_alignment.svg')
            fig.savefig(fname, format='svg')
        return fig
 

    def get_ring_current_per_scan(self):
        """Gets the ring current for each scan (i.e rotation/y-step)
        Stores it inside self.ring_currents_per_scan and a scaled version inside self.ring_currents_per_scan_scaled"""
        if not hasattr(self, "ring_currents_per_scan"):
            ring_currents = []
            with h5py.File(self.masterfile, "r") as h5in:
                for scan in self.scans:
                    ring_current = float(h5in[scan]["instrument/machine/current"][()])
                    ring_currents.append(ring_current)

            self.ring_currents_per_scan = np.array(ring_currents)
            self.ring_currents_per_scan_scaled = np.array(
                ring_currents / np.max(ring_currents)
            )

    def get_monitor(self, name="fpico6"):
        # masterfile or sparsefile
        hname = self.masterfile
        if hasattr(self, "sparsefile") and os.path.exists(self.sparsefile):
            hname = self.sparsefile
        monitor = []
        with h5py.File(hname, "r") as hin:
            for scan in self.scans:
                if scan.find("::") > -1:
                    snum, slc = scan.split("::")
                    lo, hi = [int(v) for v in slc[1:-1].split(":")]
                    mon = hin[snum]["measurement"][name][lo:hi]
                else:
                    mon = hin[scan]["measurement"][name][:]
                monitor.append(mon)

        return np.concatenate(monitor).reshape(self.shape)
    
    def reset_peaks_cache(self):
        """
        Clear cached peaks table - relevant if you set a monitor which will change intensities in columnfiles.
        """
        import warnings
        if self._pk2d is not None:
            # we have an existing 2D peaks table
            warnings.warn("Clearing cached pk2d")
            self._pk2d = None
        
        if self._pk4d is not None:
            # we have an existing 4D peaks table
            warnings.warn("Clearing cached pk4d")
            self._pk4d = None
        
        if os.path.exists(self.col2dfile):
             warnings.warn("I found an existing 2D colfile on disk - you probably want to remake this with ds.get_cf_2d(ignore_existing=True)")
                
        if os.path.exists(self.col4dfile):
             warnings.warn("I found an existing 4D colfile on disk - you probably want to remake this with ds.get_cf_4d(ignore_existing=True)")
    
    def set_monitor(self, name="fpico6", ref_value_func=np.mean):
        """
        Sets self.monitor and self.monitor_ref after calling self.get_monitor()
        Clears cached pk2d and pk4d so they can be re-computed
        
        ref_value_func: function to apply to self.monitor to generate a reference value
        when we normalise, we multiply by ref_value_func(self.monitor)/self.monitor
        we suggest np.mean as an example...
        hint: if you want to return a constant, use this:
        ref_value_func=lambda x: 1e5
        """
        self.monitor = self.get_monitor(name=name)
        self.monitor_ref = ref_value_func(self.monitor)
        
        self.reset_peaks_cache()
        
    
    def get_monitor_pk2d(self, pk2d, name='fpico6'):
        """
        To be used to normalise the peaks 2d
        """
        if self.monitor is None:
            monitor = self.get_monitor(name)
        else:
            monitor = self.monitor
        iy = np.digitize( pk2d['dty'], self.ybinedges ) - 1
        io = np.digitize( pk2d['omega'], self.obinedges ) - 1 
        #pk2d['iy'] = iy  # cache these too ?
        #pk2d['io'] = io
        return monitor[ iy, io ]

    def guess_detector(self):
        """Guess which detector we are using from the masterfile"""

        # open the masterfile
        hname = self.masterfile
        detectors_seen = []
        scan = "1.1"

        with h5py.File(hname, "r") as hin:
            # go through the first scan, and see what detectors we can see
            for measurement in list(hin[scan]["measurement"]):
                if measurement.attrs.get("interpretation") == "image":
                    detectors_seen.append(measurement)

        if len(detectors_seen) != 1:
            raise ValueError(
                "More than one detector seen! Can't work out which one to process."
            )
        else:
            self.detector = detectors_seen[0]

    #     def guess_motornames(self):
    #         '''Guess which station we were using (Nanoscope or 3DXRD) from which motors are in instrument/positioners'''
    #         from ImageD11.sinograms.assemble_label import HEADERMOTORS_NSCOPE, HEADERMOTORS_TDXRD, HEADERMOTORS
    #         # open the masterfile
    #         hname = self.masterfile
    #         motors_seen = []
    #         scan = "1.1"

    #         with h5py.File( hname, 'r' ) as hin:
    #             # go through the first scan, and see what motors we can see
    #             for positioner in list(hin[scan]['instrument/positioners']):
    #                 if positioner in HEADERMOTORS:
    #                     motors_seen.append(positioner)

    #         using_nscope = False
    #         using_tdxrd = False
    #         for motor in motors_seen:
    #             if motor in HEADERMOTORS_NSCOPE:
    #                 using_nscope = True
    #             elif motor in HEADERMOTORS_TDXRD:
    #                 using_tdxrd = True

    #         if using_nscope and using_tdxrd:
    #             raise ValueError("Found both nscope and tdxrd motors in positioners, not sure which one we were using!")

    #         if using_nscope:
    #             self.omegamotor = 'rot_center'
    #             self.dtymotor = 'dty'
    #         elif using_tdxrd:
    #             self.omegamotor = 'diffrz'
    #             self.dtymotor = 'diffty'

    def sinohist(self, weights=None, omega=None, dty=None, method="fast", return_edges=False):
        """Bin some data onto the sinogram histogram"""
        omin = self.omega_for_bins.min()
        omax = self.omega_for_bins.max()
        omin_edge = omin - self.ostep / 2
        omax_edge = omax + self.ostep / 2

        bins = len(self.obincens), len(self.ybincens)
        rng = (
            (omin_edge, omax_edge),
            (self.ybinedges[0], self.ybinedges[-1]),
        )
        if isinstance(weights, np.ndarray):
            wt = weights.ravel()
        else:
            wt = weights
        if omega is None:
            omega = self.omega_for_bins
        if dty is None:
            dty = self.dty
            
        if (self.omega.max() - self.omega.min()) > 360:
            om_mod = omega % 360
        else:
            om_mod = omega
        
        if method == "numpy":
            ret = np.histogram2d(
                om_mod.ravel(), dty.ravel(), weights=wt, bins=bins, range=rng
            )
            histo = ret[0]
        elif method == "fast":
            histo = fast_histogram.histogram2d(
                om_mod.ravel(), dty.ravel(), weights=wt, bins=bins, range=rng
            )
        if return_edges:
            om_edges = np.linspace(rng[0][0], rng[0][1], bins[0])
            dty_edges = np.linspace(rng[1][0], rng[1][1], bins[1])
            return histo, om_edges, dty_edges
        else:
            return histo

    def get_phases_from_disk(self):
        if not hasattr(self, "parfile") or self.parfile is None:
            raise AttributeError("Need self.parfile to load phases!")
        return ImageD11.unitcell.Phases(self.parfile)

    @property
    def peaks_table(self):
        if self._peaks_table is None:
            self._peaks_table = ImageD11.sinograms.properties.pks_table.load(
                self.pksfile
            )
        return self._peaks_table

    @property
    def pk2d(self):
        if self._pk2d is None:
            if self.monitor is not None:
                # we normalise
                scale_factor = self.monitor_ref/self.monitor
                self._pk2d = self.peaks_table.pk2d(self.omega_for_bins, self.dty, scale_factor=scale_factor)
            else:
                # don't normalise
                self._pk2d = self.peaks_table.pk2d(self.omega_for_bins, self.dty)
        return self._pk2d

    @property
    def pk4d(self):
        if self._pk4d is None:
            if self.monitor is not None:
                # we normalise
                scale_factor = self.monitor_ref/self.monitor
                self._pk4d = self.peaks_table.pk2dmerge(self.omega_for_bins, self.dty, scale_factor=scale_factor)
            else:
                # don't normalise
                self._pk4d = self.peaks_table.pk2dmerge(self.omega_for_bins, self.dty)
        return self._pk4d

    def get_colfile_from_peaks_dict(self, peaks_dict=None):
        """Converts a dictionary of peaks (peaks_dict) into an ImageD11 columnfile
        adds on the geometric computations (tth, eta, gvector, etc)
        Uses self.pk2d if no peaks_dict provided"""
        # TODO add optional peaks mask

        if peaks_dict is None:
            peaks_dict = self.pk2d

        cf = colfile_from_dict(peaks_dict)

        # Define spatial correction
        if hasattr(self, "detectorh5") and (self.detectorh5 is not None):
            cf = correct_cf_with_h5files(cf, self.detectorh5, self.detector)
        elif hasattr(self, "e2dxfile") and (self.e2dxfile is not None):
            cf = correct_cf_with_dxdyfiles(cf, self.e2dxfile, self.e2dyfile)
        elif hasattr(self, "splinefile") and (self.splinefile is not None):
            cf = correct_cf_with_spline(cf, self.splinefile)
        else:
            print('No spatial correction files supplied. Will return uncorrected file.')
        # Generate columnfile from peaks table
        return cf

    def update_colfile_pars(self, cf, phase_name=None):
        """Load parameters and update geometry for colfile"""
        cf.parameters.loadparameters(self.parfile, phase_name=phase_name)
        cf.updateGeometry()

    def get_cf_2d(self, ignore_existing=False):
        if os.path.exists(self.col2dfile) and not ignore_existing:
            print("Loading existing colfile from", self.col2dfile)
            return self.get_cf_2d_from_disk()
        return self.get_colfile_from_peaks_dict()

    def get_cf_4d(self, ignore_existing=False):
        if os.path.exists(self.col4dfile) and not ignore_existing:
            print("Loading existing colfile from", self.col4dfile)
            return self.get_cf_4d_from_disk()
        return self.get_colfile_from_peaks_dict(peaks_dict=self.pk4d)

    def get_cf_2d_from_disk(self):
        cf_2d = ImageD11.columnfile.columnfile(self.col2dfile)
        return cf_2d

    def get_cf_3d_from_disk(self):
        cf_3d = ImageD11.columnfile.columnfile(self.col3dfile)
        return cf_3d

    def get_cf_4d_from_disk(self):
        cf_4d = ImageD11.columnfile.columnfile(self.col4dfile)
        return cf_4d

    def get_grains_from_disk(self, phase_name=None):
        group_name = "grains"
        if phase_name is not None:
            group_name = phase_name
        grains = ImageD11.grain.read_grain_file_h5(
            self.grainsfile, group_name=group_name
        )
        if phase_name is not None and hasattr(self, "phases"):
            print("Adding reference unitcells from self.phases")
            for g in grains:
                g.ref_unitcell = self.phases.unitcells[phase_name]
        return grains

    def save_grains_to_disk(self, grains, phase_name=None):
        group_name = "grains"
        if phase_name is not None:
            group_name = phase_name
        ImageD11.grain.write_grain_file_h5(
            self.grainsfile, grains, group_name=group_name
        )

    def import_nnz(self):
        """Read the nnz arrays from the sparsefiles"""
        nnz = []
        for spname in self.sparsefiles:
            with h5py.File(os.path.join(self.analysispath, spname), "r") as hin:
                nnz.append(hin[self.limapath]["nnz"][:])
        self.nnz = np.concatenate(nnz).reshape(self.shape).astype(np.int32)
        logging.info(
            "imported nnz, average %f" % (self.nnz.mean())
        )  # expensive if you are not logging it.

    def import_nnz_from_sparse(self):
        """Read the nnz arrays from the sparsefiles"""
        with h5py.File(self.sparsefile, "r") as hin:
            self.nnz = np.array([hin[scan]["nnz"][:] for scan in self.scans])
        logging.info(
            "imported nnz, average %f" % (self.nnz.mean())
        )  # expensive if you are not logging it.
        self.frames_per_scan = [len(nnz) for nnz in self.nnz]

    #    def compute_pixel_labels(self):
    # this should instead from from the pk2d file generated by sinograms/properties.py
    #        nlm = []
    #        for spname in self.sparsefiles:
    #            n, l = peaklabel.add_localmax_peaklabel( os.path.join( self.analysispath, spname ),
    #                                                     self.limapath )
    #            nlm.append(n)
    #        self.nlm = np.concatenate( nlm ).reshape( self.shape )

    #    def import_nlm(self):
    # this should instead from from the pk2d file generated by sinograms/properties.py
    #        """ Read the Nlmlabels
    #        These are the number of localmax peaks per frame
    #        """
    #        nlm = []
    #        for spname in self.sparsefiles:
    #            with h5py.File( os.path.join( self.analysispath, spname ), "r" ) as hin:
    #                nlm.append( hin[self.limapath]['Nlmlabel'][:] )
    #        self.nlm = np.concatenate( nlm ).reshape( self.shape )
    #        logging.info('imported nlm, max %d'%(self.nlm.max()))

    def check_files(self, path, filenames, verbose=0):
        """See whether files are created or not"""
        # images collected
        done = 0
        missing = 0
        for fname in filenames:
            fullname = os.path.join(path, fname)
            if os.path.exists(fullname):
                done += 1
            else:
                missing += 1
                if verbose > 0:
                    print("missing", fullname)
                    verbose -= 1
        return done, missing

    def check_images(self):
        """Is the experiment finished ?"""
        return self.check_files(self.datapath, self.imagefiles)

    def check_sparse(self):
        """Has the segmentation been done ?"""
        return self.check_files(self.analysispath, self.sparsefiles, verbose=2)

    def save(self, h5name=None, h5group="/"):
        if h5name is None:
            if os.path.exists( self.dsfile ):
                h5name = self.dsfile

            # none supplied, so use default path
            h5name = self.dsfile_default
            # make sure parent directories exist
            # ensure that the analysis path exists
            dsfile_folder = os.path.dirname(self.dsfile_default)
            if not os.path.exists(dsfile_folder):
                os.makedirs(dsfile_folder)

        ZIP = {"compression": "gzip"}

        with h5py.File(h5name, "a") as hout:
            grp = hout[h5group]
            # Simple small objects
            for name in self.ATTRNAMES:
                data = getattr(self, name, None)
                if data is not None:
                    grp.attrs[name] = data
                    # The string lists
            for name in self.STRINGLISTS:
                data = getattr(self, name, None)
                if data is not None and len(data):
                    sdata = np.array(data, "S")
                    ds = grp.require_dataset(
                        name,
                        shape=sdata.shape,
                        chunks=sdata.shape,
                        dtype=h5py.string_dtype(),
                        **ZIP
                    )
                    ds[:] = sdata
            #
            for name in self.NDNAMES:
                data = getattr(self, name, None)
                if data is not None:
                    data = np.asarray(data)
                    try:
                        chunks = guess_chunks(name, data.shape)
                        ds = grp.require_dataset(
                            name,
                            shape=data.shape,
                            chunks=chunks,
                            dtype=data.dtype,
                            **ZIP
                        )
                        ds[:] = data
                    except:
                        print(name)
                        print(len(data))
                        print(data.shape)
                        print(chunks)
                        raise

        # if we got here, we saved the file successfully
        self.dsfile = h5name

    def load(self, h5name=None, h5group="/"):

        if h5name is None:
            if os.path.exists( self.dsfile ):
                h5name = self.dsfile
            elif os.path.exists( self.dsfile_default ):
                # none supplied, so use default path
                h5name = self.dsfile_default
            else:
                raise Exception( "Filename for dataset not found")

        """ Recover this from a hdf5 file """
        with h5py.File(h5name, "r") as hin:
            grp = hin[h5group]
            for name in self.ATTRNAMES:
                if name in grp.attrs:
                    setattr(self, name, grp.attrs.get(name))
            self.shape = tuple(self.shape)  # hum
            for name in self.NDNAMES:
                if name in grp:
                    data = grp[name][()]
                    setattr(self, name, data)
            for name in self.STRINGLISTS:
                if name in grp:
                    stringlist = list(grp[name][()])
                    if hasattr(stringlist[0], "decode") or isinstance(
                        stringlist[0], np.ndarray
                    ):
                        data = [s.decode() for s in stringlist]
                    else:
                        data = stringlist
                    setattr(self, name, data)
        self.guessbins()

        # analysis paths can only be calculated once
        self.update_paths()
        # if we got here, we loaded the file successfully
        self.dsfile = h5name

        return self


def load(h5name, h5group="/"):
    ds_obj = DataSet(filename=h5name)
    return ds_obj


# Example
#    s = dataset(
#        dataroot = "/data/visitor/ma5415/id11/20221027",
#        analysisroot = "/data/visitor/ma5415/id11/20221027/analysis/SparsePixels",
#        sample = "NSCOPE_SiMo1000_6",
#        dset = "DT600nm_Z1" )
#    s.import_all()


def check(dataroot, analysisroot, sample, dset, destination, scans=None):
    h5o = DataSet(
        dataroot=dataroot, analysisroot=analysisroot, sample=sample, dset=dset
    )
    h5o.import_all(scans=scans)
    h5o.save(destination)

    print("Checking: Read back from hdf5")
    t = load(destination)
    t.report()
    return t.compare(h5o)


if __name__ == "__main__":
    import sys

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

    dataroot = sys.argv[1]
    analysisroot = sys.argv[2]
    sample = sys.argv[3]
    dset = sys.argv[4]
    destination = sys.argv[5]

    check(dataroot, analysisroot, sample, dset, destination)
