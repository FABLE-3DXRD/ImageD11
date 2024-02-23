from __future__ import print_function
import h5py, os, time, numpy as np
import numba
import fabio

"""WARNING: work in progresss"""

import ImageD11.sinograms.dataset


# given sample + dset -> make a single file with all the sparse pixels in it.
#  ... this is not really needed, but saves redoing the pixel labelling code


SCANMOTORS = (
    "diffrz",
    "diffrz_cen360",
    "diffrz_center",
    "fpico4",
    "fpico3",
    "diffty",
    "diffty_center",
    "rot",
    "rot_cen360",
    "rot_center",
    "fpico6",
    "dty",
    "dty_center",
)
HEADERMOTORS = (
    "diffty",
    "diffrz",
    "samtx",
    "samty",
    "samtz",
    "diffry",
    "samrx",
    "samry",
    "dty",
    "rot",
    "pz",
    "px",
    "py",
    "shtz",
    "shty",
    "shtx",
)


def testready(dsobj):
    """assume they are finish if they exist. Might be still writing ???"""
    done, todo = dsobj.check_sparse()
    return (done > 0) and (todo == 0)


def getsparse(dsobj, num):
    """
    dsobj = DataSet object
    returns the pixels from the sparse segmentation
    """
    with h5py.File(
        os.path.join(dsobj.analysispath, dsobj.sparsefiles[num]), "r"
    ) as hin:
        row = hin[dsobj.limapath]["row"][:]
        col = hin[dsobj.limapath]["col"][:]
        intensity = hin[dsobj.limapath]["intensity"][:]
        nnz = hin[dsobj.limapath]["nnz"][:]
    return row, col, intensity, nnz


@numba.njit(boundscheck=True)
def filterpixels(cutimage, row, col, intensity, nnz):
    """
    cutimage = 2D image of thresholds to use
    row, col, intensity = pixels
    nnz = pixels per frame for a stack of images in row/col/intensity
    """
    # 1: scan all of the pixels to decide on the output space required
    nnz_out = np.zeros_like(nnz)
    ntotal = 0
    ip = 0
    for iframe in range(len(nnz)):
        for j in range(nnz[iframe]):
            if intensity[j + ip] > cutimage[row[j + ip], col[j + ip]]:
                nnz_out[iframe] += 1
                ntotal += 1
        ip += nnz[iframe]
    row_out = np.zeros((ntotal,), dtype=row.dtype)
    col_out = np.zeros((ntotal,), dtype=col.dtype)
    intensity_out = np.zeros((ntotal,), dtype=intensity.dtype)
    nout = 0
    ip = 0
    for iframe in range(len(nnz)):
        for j in range(nnz[iframe]):
            if intensity[j + ip] > cutimage[row[j + ip], col[j + ip]]:
                row_out[nout] = row[j + ip]
                col_out[nout] = col[j + ip]
                intensity_out[nout] = intensity[j + ip]
                nout += 1
        ip += nnz[iframe]
    return row_out, col_out, intensity_out, nnz_out


def harvest_masterfile(
    dset, outname, scanmotors=SCANMOTORS, headermotors=HEADERMOTORS, cutimage=None
):
    """
    dset = ImageD11.sinograms.dataset.DataSet object
    outname = sparse file to write

    cut = image = 2D image of cutoffs for keeping pixels (e.g. keep high angle
                  and use a higher threshold for low angle or noisy)
    """
    opts = {
        "chunks": (10000,),
        "maxshape": (None,),
        "compression": "lzf",
        "shuffle": True,
    }
    with h5py.File(outname, "a") as hout:
        hout.attrs["h5input"] = dset.masterfile
        print("Harvesting", dset.masterfile, end=": ")
        with h5py.File(dset.masterfile, "r") as hin:
            done = []
            for scan in dset.scans:
                if scan.find("::"):
                    scan = scan.split("::")[0]
                if scan in done:
                    continue
                gin = hin[scan]
                bad = False
                for check in ("title", "measurement", "measurement/" + dset.detector):
                    if check not in hin[scan]:
                        print(scan, "missing", check, "skipping")
                        bad = True
                if bad:
                    print("Skipping", scan)
                    continue
                title = hin[scan]["title"][()]
                g = hout.require_group(scan)
                gm = g.require_group("measurement")
                for m in scanmotors:  # vary : many
                    if m in gin["measurement"]:
                        data = data = gin["measurement"][m][:]
                        ds = gm.require_dataset(m, shape=data.shape, dtype=data.dtype)
                        ds[()] = data
                gip = g.require_group("instrument/positioners")
                for m in headermotors:  # fixed : scalar
                    if "instrument/positioners/%s" % (m) in gin:
                        data = gin["instrument/positioners"][m][()]
                        ds = gip.require_dataset(m, shape=data.shape, dtype=data.dtype)
                        ds[()] = data
                try:
                    frms = gin["measurement"][dset.detector]
                except Exception as e:
                    print(e)
                    print(list(gin))
                    print(list(gin["measurement"]))
                    print(dset.detector)
                    raise
                g.attrs["itype"] = frms.dtype.name
                g.attrs["nframes"] = frms.shape[0]
                g.attrs["shape0"] = frms.shape[1]
                g.attrs["shape1"] = frms.shape[2]
                print(scan, end=", ")
                done.append(scan)
            print()

        # Finished with master file. Now harvest the segmented files.
        idx = 0
        print("Loading pixels:", end=" ")
        for scan in done:
            g = hout.require_group(scan)
            for name in "row", "col":
                if name not in g:
                    g.create_dataset(name, shape=(0,), dtype=np.uint16, **opts)
            if "intensity" not in g:
                g.create_dataset(
                    "intensity", shape=(0,), dtype=g.attrs["itype"], **opts
                )
            nfrm = g.attrs["nframes"]
            g.require_dataset("nnz", shape=(nfrm,), dtype=np.uint32)
            nstart = nread = npx = pstart = 0
            while nread < nfrm:
                row, col, intensity, nnz = getsparse(dset, idx)
                if cutimage is not None:
                    row, col, intensity, nnz = filterpixels(
                        cutimage, row, col, intensity, nnz
                    )
                idx += 1  # loop over sparse files in this scan
                nread = nstart + len(nnz)  # number of frames in this limafile
                g["nnz"][nstart:nread] = nnz
                nstart = nread
                pread = pstart + len(row)  # number of pixels in this limafile
                g["row"].resize((pread,))
                g["row"][pstart:pread] = row
                g["col"].resize((pread,))
                g["col"][pstart:pread] = col
                g["intensity"].resize((pread,))
                g["intensity"][pstart:pread] = intensity
                pstart = pread
            print(scan, end=", ")
        print()
    return outname


def main(dsname, outname=None, cutimage=None):
    """
    dsname = Dataset describing the masterfile + segmentation etc
    outname = sparse pixels file to write. Defaults to dset.sparsefile
    """
    if isinstance(dsname, ImageD11.sinograms.dataset.DataSet):
        dset = dsname
    else:
        dset = ImageD11.sinograms.dataset.load(dsname)
    if outname is None:
        outname = dset.sparsefile
    if cutimage is not None:
        cutimage = fabio.open(cutimage).data
    harvest_masterfile(dset, outname, cutimage=cutimage)


if __name__ == "__main__":
    import sys

    main(sys.argv[1], sys.argv[2])
