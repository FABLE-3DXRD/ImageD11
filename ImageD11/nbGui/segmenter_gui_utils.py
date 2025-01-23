import h5py
import numpy as np
import os


def guess_ESRF_paths():  # This should be in silx somewhere?
    """ Locates:
    dataroot     holds raw data       in folders dataroot     + {sample}/{sample}_{dataset}
    analysisroot holds output results in folders analysisroot + {sample}/{sample}_{dataset}
    """
    path_items = os.getcwd().split('/')
    if 'visitor' in path_items:
        idx = path_items.index('visitor')
        experiment, session = path_items[idx + 1], path_items[idx + 3]
        path = os.path.join("/data", "visitor", experiment, "id11", session)
        return [os.path.join(path, folder) for folder in
                ("RAW_DATA", "PROCESSED_DATA")]
    return "", ""


def printsamples(dataroot):
    samples = sorted([name for name in os.listdir(dataroot)
                      if os.path.isdir(os.path.join(dataroot, name))])
    print("Samples:\n\t " + "\n\t".join(sorted(samples)))


def printdatasets(dataroot, sample):
    sroot = os.path.join(dataroot, sample)
    print("Datsets:\n\t " + "\n\t".join(sorted(
        [name[len(sample) + 1:] for name in os.listdir(sroot)
         if os.path.isdir(os.path.join(sroot, name))
         and name.startswith(sample)])))


def chooseframe(dset, scan=None, idx=None, counter="_roi1", fetch_raw_image=False):
    """
    Locate a busy frame from the dataset and optionally fetch the raw image.

    Args:
        dset: Dataset object.
        scan: Scan to use (default is middle scan).
        idx: Frame index (None to auto-detect).
        counter: Counter to use for locating the frame.
        fetch_raw_image: Whether to fetch the raw image for the frame.

    Returns:
        scan, idx, raw_image (if fetch_raw_image is True, otherwise None)
    """
    if scan is None:
        scan = dset.scans[len(dset.scans) // 2]

    raw_image = None
    with h5py.File(dset.masterfile, "r") as hin:
        ctr = dset.detector + counter

        if idx is None:
            if scan.find("::") > -1:  # 1.1::[10000:12000]  etc
                lo, hi = [int(v) for v in scan[:-1].split("[")[1].split(":")]
                scan = scan.split("::")[0]
                roi1 = hin[scan]["measurement"][ctr][lo:hi]
                idx = np.argmax(roi1) + lo
            else:  # "1.1"
                roi1 = hin[scan]["measurement"][ctr][:]
                idx = np.argmax(roi1)

        print("Using frame", idx, "from scan", scan)

        # placing the raw image
        if fetch_raw_image:
            raw_image = hin[scan + "/measurement/" + dset.detector][idx]

    return scan, idx, raw_image
