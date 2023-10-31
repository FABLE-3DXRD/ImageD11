from __future__ import print_function, division

import os, h5py, numpy as np
import fast_histogram
import logging

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


def guess_chunks(name, shape):
    if name == "omega":
        return (shape[0], 1)
    if name == "dty":
        return (1, shape[1])
    return shape


class DataSet:

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
    )
    STRINGLISTS = ("scans", "imagefiles", "sparsefiles")
    # sinograms
    NDNAMES = ("omega", "dty", "nnz", "frames_per_file", "nlm", "frames_per_scan")

    def __init__(self, dataroot=".", analysisroot=".", sample="sample", dset="dataset"):
        """The things we need to know to process data"""

        self.detector = "eiger"  # frelon3
        self.limapath = None  # where is the data in the Lima files

        self.omegamotor = "rot_center"  # diffrz
        self.dtymotor = "dty"  # diffty

        self.dataroot = dataroot  # folder to find {sample}/{sample}_{dset}
        self.analysisroot = analysisroot  # where to write or find sparse data
        self.sample = sample  # from bliss path
        self.dset = dset  # from bliss path

        self.dsname = "_".join((self.sample, self.dset))

        self.datapath = os.path.join(self.dataroot, self.sample, self.dsname)
        self.analysispath = os.path.join(self.analysisroot, self.sample, self.dsname)
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

    def import_scans(self, scans=None, hname=None):
        """Reads in the scans from the bliss master file"""
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
                    )
                ]
            goodscans = []
            for scan in scans:
                frames = hin[scan]["measurement"][self.detector]
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
        # npts = None
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
            name.replace("/", "_").replace(".h5", "_sparse.h5")
            for name in self.imagefiles
        ]
        logging.info("imported %d lima filenames" % (np.sum(self.frames_per_file)))

    def import_motors_from_master(
        self,
    ):  #  could also get these from sparse files if saved
        """read the motors from the lima file
        you need to import the imagefiles first
        these will be the motor positions to accompany the images
        """
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
                (abs(self.omega[i][0] - self.omega[b]), i)
                for i in range(len(self.scans))
                if i not in bad
            ]
            if len(dom) > 0:
                j = np.argmin(dom[0][1])
                self.omega[b] = self.omega[j]  # best match
                print(
                    "replace bad scan omega", b, self.scans[b], "with", j, self.scans[j]
                )
        logging.info("imported omega/dty")

    def guess_shape(self):
        """Guess the shape if it was not given"""
        npts = np.sum(self.frames_per_scan)
        if len(self.scans) == 1:  # probably fscan2d or f2scan
            with h5py.File(self.masterfile, "r") as hin:
                s = hin[self.scans[0]]
                title = s["title"].asstr()[()]
                print("Scan title", title)
                if title.split()[0] == "fscan2d":
                    s0 = s["instrument/fscan_parameters/slow_npoints"][()]
                    s1 = s["instrument/fscan_parameters/fast_npoints"][()]
                elif title.split()[0] == "f2scan":
                    # good luck ? Assuming rotation was the inner loop here:
                    step = s["instrument/fscan_parameters/step_size"][()]
                    s1 = int(np.round(360 / step))
                    s0 = npts // s1
                else:
                    s0 = 1
                    s1 = npts
        if len(self.scans) > 1:
            s0 = len(self.scans)
            s1 = npts // s0
        self.shape = s0, s1
        if np.prod(self.shape) != npts:
            print("Warning: irregular scan - might be bugs in here")
            print(npts, len(self.scans))
        self.omega = np.array(self.omega)
        self.dty = np.array(self.dty)
        logging.info(
            "sinogram shape = ( %d , %d ) imageshape = ( %d , %d)"
            % (self.shape[0], self.shape[1], self.imageshape[0], self.imageshape[1])
        )

    def guessbins(self):
        ny, nomega = self.shape
        self.omin = self.omega.min()
        self.omax = self.omega.max()
        if (self.omax - self.omin) > 360:
            # multi-turn scan...
            self.omin = 0.0
            self.omax = 360.0
            self.omega_for_bins = self.omega % 360
        else:
            self.omega_for_bins = self.omega
        # values 0, 1, 2
        # shape = 3
        # step = 1
        self.ostep = (self.omax - self.omin) / (nomega - 1)
        self.ymin = self.dty.min()
        self.ymax = self.dty.max()
        if ny > 1:
            self.ystep = (self.ymax - self.ymin) / (ny - 1)
        else:
            self.ystep = 1
        self.obincens = np.linspace(self.omin, self.omax, nomega)
        self.ybincens = np.linspace(self.ymin, self.ymax, ny)
        self.obinedges = np.linspace(
            self.omin - self.ostep / 2, self.omax + self.ostep / 2, nomega + 1
        )
        self.ybinedges = np.linspace(
            self.ymin - self.ystep / 2, self.ymax + self.ystep / 2, ny + 1
        )

    def sinohist(self, weights=None, omega=None, dty=None, method="fast"):
        """Bin some data onto the sinogram histogram"""
        bins = len(self.obincens), len(self.ybincens)
        rng = (
            (self.obinedges[0], self.obinedges[-1]),
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
        if method == "numpy":
            ret = np.histogram2d(
                omega.ravel(), dty.ravel(), weights=wt, bins=bins, range=rng
            )
            histo = ret[0]
        elif method == "fast":
            histo = fast_histogram.histogram2d(
                omega.ravel(), dty.ravel(), weights=wt, bins=bins, range=rng
            )
        return histo

    def import_nnz(self):
        """Read the nnz arrays from the scans"""
        nnz = []
        for spname in self.sparsefiles:
            with h5py.File(os.path.join(self.analysispath, spname), "r") as hin:
                nnz.append(hin[self.limapath]["nnz"][:])
        self.nnz = np.concatenate(nnz).reshape(self.shape).astype(np.int32)
        logging.info(
            "imported nnz, average %f" % (self.nnz.mean())
        )  # expensive if you are not logging it.

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

    def save(self, h5name, h5group="/"):

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

    def load(self, h5name, h5group="/"):
        """Recover this from a hdf5 file"""
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
        return self


def load(h5name, h5group="/"):
    return DataSet().load(h5name, h5group)


def import_from_sparse(
    hname,
    omegamotor="instrument/positioners/rot",
    dtymotor="instrument/positioners/dty",
):
    ds = DataSet()
    with h5py.File(hname, "r") as hin:
        scans = list(hin["/"])
        order = np.argsort([float(v) for v in scans if v.endswith(".1")])
        scans = [scans[i] for i in order]
        dty = [hin[scan][dtymotor][()] for scan in scans]
        omega = [hin[scan][omegamotor][()] for scan in scans]
        nnz = [hin[scan]["nnz"][()] for scan in scans]
    #        nlm = [hin[scan]['Nlmlabel'][()] for scan in scans]
    ds.scans = scans
    ds.nnz = nnz
    ds.nnz = np.array(nnz)
    ds.shape = ds.nnz.shape
    ds.omega = np.zeros(ds.nnz.shape, float)
    for i, o in enumerate(omega):
        if isinstance(o, float) or (len(o) == len(ds.nnz[i])):
            ds.omega[i] = o
        if len(o) > len(ds.nnz[i]):
            ds.omega[i] = ds.omega[i - 2]  # guess zig zag
            # warning here

    ds.dty = np.zeros(ds.nnz.shape, float)
    for i, o in enumerate(dty):
        if isinstance(o, float) or (len(o) == len(ds.nnz[i])):
            ds.dty[i] = o
        else:
            raise Exception("Cannot read %d dty %s %s" % (i, str(o), str(o.shape)))
    #    assert ds.nlm.shape == ds.shape
    try:
        ds.guess_scans()
    except:
        print("warning, guess scans failed")
    try:
        ds.guessbins()
    except:
        print("warning, guessbins failed")
    return ds


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
