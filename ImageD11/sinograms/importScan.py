import logging
from pathlib import Path
import h5py
import numpy as np

logger = logging.getLogger(__name__)


class ImportScan:
    """
    Handles importing scan and keep the state of imported scan
    """

    def __init__(
        self,
        masterfile: str = None,
        detector: str = "eiger",
        omegamotor: str = "rot_center",
        dtymotor: str = "dty",
        load_import_scan_file: str = None,
    ):
        """
        Initialize the scan importer.

        Args:
            masterfile (str): Path to the master .h5 file.
            detector (str): Type of detector (e.g., "eiger", "frelon3").
            omegamotor (str): Motor for omega rotation (default: "rot_center").
            dtymotor (str): Motor for y-axis translation (default: "dty").
            load_import_scan_file (str):
                Path to load previously stored scan data h5file.
        """
        self.masterfile = masterfile
        self.detector = detector
        self.omegamotor = omegamotor
        self.dtymotor = dtymotor

        if load_import_scan_file:
            self.load_my_scan(load_import_scan_file)
            logger.info(f"Loaded scan data from {load_import_scan_file}")
        else:
            if not self.masterfile or not Path(self.masterfile).exists():
                logger.error(f"Master file {self.masterfile} does not exist.")
                raise FileNotFoundError(
                    f"Provided master file does not exist: {self.masterfile}"
                )

            self._initialize()

    def _initialize(self):
        """Initialize scan attributes."""
        self.scans = None  # list of good scans
        self.frames_per_scan = None  # list of frames per good scan
        self.imagefiles = None  # list[str], filenames from scans
        self.frames_per_file = (
            None  # numpy array, have number of frames per file (scan)
        )
        self.imageshape = 0  # shape of frame,
        self.limapath = None  # str, path
        self.omega = None  # list of omega values for the number of scans
        self.dty = None  # list of translation values for the number of scans
        self.shape = None
        # floats, omeaga min, max and step
        self.omin, self.omax, self.ostep = None, None, None
        self.omega_for_bins = None  # list[int]
        self.obincens = None  # numpy array
        self.obinedges = None  # numpy array

        self.ymin, self.ymax, self.ystep = None, None, None
        self.ybincens = None  # numpy array
        self.ybinedges = None  # numpy array

    def import_scans(self, scans: list[str] = None):

        with h5py.File(self.masterfile, "r") as hin:
            available_scans = list(hin.keys())

            if scans is None:
                scans = [
                    scan
                    for scan in available_scans
                    if scan.endswith(".1")
                    and "measurement" in hin[scan]
                    and self.detector in hin[scan]["measurement"]
                    and self.omegamotor in hin[scan]["measurement"]
                ]

            good_scans = []
            frames_per_scan = []

            for scan in scans:
                try:
                    frames = hin[scan]["measurement"][self.detector]
                    if len(frames.shape) == 3:
                        good_scans.append(scan)
                        frames_per_scan.append(frames.shape[0])
                    else:
                        logger.warning(f"Bad scan: {scan} (Invalid shape)")
                except KeyError as e:
                    logger.warning(f"Bad scan: {scan}, h5py error: {e}")

        self.scans = good_scans
        self.frames_per_scan = frames_per_scan

        logger.info(f"Imported {len(self.scans)} scans from {self.masterfile}")
        return self.scans

    def import_imagefiles(self):
        """Extracts image file names and frame counts from the master file."""
        if not self.scans:
            logger.error("Call import_scans() before importing image files.")
            raise AttributeError("Scans must be imported first.")

        self.imagefiles = []
        self.frames_per_file = []

        with h5py.File(self.masterfile, "r") as hin:
            for scan in self.scans:
                try:
                    frames = hin[scan]["measurement"][self.detector]
                    self.imageshape = frames.shape[1:]
                    for vsrc in frames.virtual_sources():
                        self.imagefiles.append(vsrc.file_name)
                        self.frames_per_file.append(vsrc.src_space.shape[0])
                        if self.limapath is None:
                            self.limapath = vsrc.dset_name
                        assert self.limapath == vsrc.dset_name
                except KeyError:
                    logger.warning(f"Skipping scan {scan} due to missing data.")

        self.frames_per_file = np.array(self.frames_per_file, dtype=int)
        logger.info(f"Imported {np.sum(self.frames_per_file)} frames.")

    def import_motors_from_master(self):
        """Extracts omega and dty motor positions from the master file."""
        if not self.scans:
            logger.error("Call import_scans() and import_imagefiles() first.")
            raise AttributeError("Scans and images must be imported first.")

        self.omega = [None] * len(self.scans)
        self.dty = [None] * len(self.scans)

        with h5py.File(self.masterfile, "r") as hin:
            for i, scan in enumerate(self.scans):
                om = hin[scan]["measurement"].get(self.omegamotor, None)
                dty = hin[scan]["instrument/positioners"].get(self.dtymotor, None)

                if om is not None:
                    self.omega[i] = (
                        om[()] if len(om) == self.frames_per_scan[i] else [om[0]]
                    )
                if dty is not None:
                    self.dty[i] = (
                        np.full(self.frames_per_scan[i], dty[()])
                        if dty.shape == ()
                        else dty[:]
                    )

        logger.info("Imported omega/dty motor positions.")

    def guess_shape(self):
        """Guesses the shape of the dataset."""
        if not self.frames_per_scan:
            logger.error("Call import_scans() first.")
            raise AttributeError("Scans must be imported first.")

        total_frames = np.sum(self.frames_per_scan)
        num_scans = len(self.scans)
        if num_scans == 1:
            with h5py.File(self.masterfile, "r") as hin:
                s = hin[self.scans[0]]
                title = s["title"].asstr()[()]
                print("Scan title", title)
                if title.split()[0] == "fscan2d":
                    s0 = s["instrument/fscan_parameters/slow_npoints"][()]
                    s1 = s["instrument/fscan_parameters/fast_npoints"][()]
                    file_nums = np.arange(s0 * s1).reshape((s0, s1))
                    # slice notation means last frame+1 to be inclusive
                    self.scans = [
                        "%s::[%d:%d]" % (self.scans[0], row[0], row[-1] + 1)
                        for row in file_nums
                    ]
                elif title.split()[0] == "f2scan":
                    # good luck ? Assuming rotation was the inner loop here:
                    step = s["instrument/fscan_parameters/step_size"][()]
                    s1 = int(np.round(360 / step))
                    s0 = total_frames // s1
                    logger.warning("Dataset might need to be reshaped")
                else:
                    s0 = 1
                    s1 = total_frames
            self.shape = (s0, s1)
        else:
            self.shape = (num_scans, total_frames // num_scans)

        logger.info(f"Dataset shape guessed: {self.shape}")
        self.omega = np.array(self.omega).reshape(self.shape)
        self.dty = np.array(self.dty).reshape(self.shape)
        if np.prod(self.shape) != total_frames:
            logger.warning("irregular scan")

    def guessbins(self):
        """Determines binning for omega and y positions."""
        if not self.shape:
            logger.error("Call guess_shape() first.")
            raise AttributeError("Shape must be determined first.")

        n_y, n_omega = self.shape
        self.omin, self.omax = np.min(self.omega), np.max(self.omega)
        if (self.omax - self.omin) > 360:
            # multi-turn scan...
            self.omin = 0.0
            self.omax = 360.0
            self.omega_for_bins = self.omega % 360
        else:
            self.omega_for_bins = self.omega
        self.ostep = (self.omax - self.omin) / (n_omega - 1)

        self.ymin, self.ymax = self.dty.min(), self.dty.max()
        self.ystep = (self.ymax - self.ymin) / (n_y - 1) if n_y > 1 else 1

        self.obincens = np.linspace(self.omin, self.omax, n_omega)
        self.ybincens = np.linspace(self.ymin, self.ymax, n_y)
        self.obinedges = np.linspace(
            self.omin - self.ostep / 2, self.omax + self.ostep / 2, n_omega + 1
        )
        self.ybinedges = np.linspace(
            self.ymin - self.ystep / 2, self.ymax + self.ystep / 2, n_y + 1
        )
        logger.info("Binning parameters estimated.")

    def import_all(self, scans: list[str] = None):
        """Imports all necessary scan data."""
        self.import_scans(scans)
        self.import_imagefiles()
        self.import_motors_from_master()
        self.guess_shape()
        self.guessbins()
        logger.info("The dataset scan folder are loaded")

    def store_my_scan(self, file_path: str) -> None:
        """Stores computed scan data in an HDF5 file."""
        pt_file = Path(file_path)

        if not pt_file.parent.exists():
            logger.error("Directory does not exist: %s", pt_file.parent)
            raise FileNotFoundError(f"Invalid directory: {pt_file.parent}")

        with h5py.File(pt_file, "w") as h5file:
            # Save basic properties
            attributes = {
                "masterfile": str(self.masterfile),
                "detector": self.detector,
                "omegamotor": self.omegamotor,
                "dtymotor": self.dtymotor,
                "imageshape": self.imageshape,
                "limapath": self.limapath or "None",
                "omin": self.omin,
                "omax": self.omax,
                "ostep": self.ostep,
                "ymin": self.ymin,
                "ymax": self.ymax,
                "ystep": self.ystep,
                "shape": self.shape,
            }
            h5file.attrs.update({k: v for k, v in attributes.items() if v is not None})

            # Save datasets
            datasets = {
                "scans": np.array(self.scans, dtype="S"),
                "frames_per_scan": self.frames_per_scan,
                "imagefiles": np.array(self.imagefiles, dtype="S"),
                "frames_per_file": self.frames_per_file,
                "omega": self.omega,
                "dty": self.dty,
                "obincens": self.obincens,
                "obinedges": self.obinedges,
                "ybincens": self.ybincens,
                "ybinedges": self.ybinedges,
                "omega_for_bins": self.omega_for_bins,
            }
            for key, value in datasets.items():
                if value is not None:
                    h5file.create_dataset(key, data=value)

        logger.info("Scan data stored successfully in '%s'", file_path)

    def load_my_scan(self, file_path: str) -> None:
        """Loads computed scan data from an HDF5 file."""
        pt_file = Path(file_path)

        if not pt_file.exists():
            logger.error("File does not exist: %s", pt_file)
            raise FileNotFoundError(f"Invalid scan file: {pt_file}")

        self._initialize()  # Reset attributes before loading

        with h5py.File(pt_file, "r") as h5file:
            # Load basic properties
            self.masterfile = h5file.attrs["masterfile"]
            self.detector = h5file.attrs["detector"]
            self.omegamotor = h5file.attrs["omegamotor"]
            self.dtymotor = h5file.attrs["dtymotor"]
            self.imageshape = tuple(h5file.attrs["imageshape"])
            self.limapath = h5file.attrs.get("limapath", None)

            # Load datasets
            self.scans = [scan.decode() for scan in h5file["scans"][()]]
            self.frames_per_scan = h5file["frames_per_scan"][()]
            self.imagefiles = [file.decode() for file in h5file["imagefiles"][()]]
            self.frames_per_file = h5file["frames_per_file"][()]

            for attr in [
                "omega",
                "dty",
                "obincens",
                "obinedges",
                "ybincens",
                "ybinedges",
                "omega_for_bins",
            ]:
                if attr in h5file:
                    setattr(self, attr, h5file[attr][()])

            # Load optional attributes
            self.shape = tuple(h5file.attrs.get("shape", ()))
            self.omin = h5file.attrs.get("omin")
            self.omax = h5file.attrs.get("omax")
            self.ostep = h5file.attrs.get("ostep")
            self.ymin = h5file.attrs.get("ymin")
            self.ymax = h5file.attrs.get("ymax")
            self.ystep = h5file.attrs.get("ystep")

        logger.info("Scan data loaded successfully from '%s'", file_path)
