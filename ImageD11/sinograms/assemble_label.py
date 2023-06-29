
from __future__ import print_function
import h5py, os, time, numpy as np

"""WARNING: work in progresss"""

import ImageD11.sinograms.dataset


# given sample + dset -> make a single file with all the sparse pixels in it.
#  ... this is not really needed, but saves redoing the pixel labelling code


SCANMOTORS = ("diffrz", "diffrz_cen360", "diffrz_center", "fpico4", "fpico3", "diffty", "diffty_center",
              "rot", "rot_cen360", "rot_center", "fpico6", "dty", "dty_center")
HEADERMOTORS = ("diffty", "diffrz", "samtx", "samty", "samtz", "diffry", "samrx", "samry",
                "dty", "rot", "pz", "px", "py", "shtz", "shty", "shtx")


def testready(dsobj):
    """ assume they are finish if they exist. Might be still writing ??? """
    done, todo = dsobj.check_sparse()
    return (done>0) and (todo == 0)
                                                       

def getsparse( dsobj, num, titles = ('row','col','intensity','nnz') ):
    """ 
    dsobj = DataSet object
    returns the pixels from the sparse segmentation
    """
    with h5py.File( os.path.join( dsobj.analysispath, dsobj.sparsefiles[num]) , "r" ) as hin:
        pixels = { name : hin[dsobj.limapath][name][:] for name in titles }
    return pixels


def harvest_masterfile( dset, outname,
                        scanmotors=SCANMOTORS,
                        headermotors=HEADERMOTORS, ):
    """
    dset = ImageD11.sinograms.dataset.DataSet object
    outname = sparse file to write
    """    
    opts = {
        "chunks": (10000,),
        "maxshape": (None,),
        "compression": "lzf",
        "shuffle": True,
    }
    with h5py.File(outname, "a") as hout:
        hout.attrs["h5input"] = dset.masterfile
        print("Harvesting",dset.masterfile,end=": ")
        with h5py.File(dset.masterfile, "r") as hin:
            for scan in dset.scans:
                gin = hin[scan]
                bad = False
                for check in ('title','measurement','measurement/'+dset.detector):
                    if check not in hin[scan]:
                        print(scan,"missing",check,'skipping')
                        bad = True
                if bad:
                    print("Skipping", scan)
                    continue
                title = hin[scan]["title"][()]
                g = hout.require_group(scan)
                gm = g.require_group("measurement")
                for m in scanmotors:  # vary : many
                    if m in gin["measurement"]:
                        data = data=gin["measurement"][m][:]
                        ds = gm.require_dataset(m, shape=data.shape, dtype = data.dtype )
                        ds[()] = data
                gip = g.require_group("instrument/positioners")
                for m in headermotors:  # fixed : scalar
                    if "instrument/positioners/%s"%(m) in gin:
                        data=gin["instrument/positioners"][m][()]
                        ds  = gip.require_dataset(m, shape = data.shape, dtype = data.dtype )
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
                print(scan, end=', ')
            print()
                        
        # Finished with master file. Now harvest the segmented files.
        idx = 0
        titles = ('row','col','intensity','nnz')
        print('Loading pixels:',end=' ')
        for scan in dset.scans:
            g = hout.require_group( scan )
            for name in 'row', 'col':
                if name not in g:
                    g.create_dataset(name, shape = (0,), dtype=np.uint16, **opts)
            if "intensity" not in g:
                g.create_dataset("intensity", shape = (0,), dtype=g.attrs['itype'], **opts)
            nfrm = g.attrs['nframes']
            g.require_dataset("nnz", shape = (nfrm,), dtype=np.uint32)
            nstart = nread = npx = pstart = 0
            while nread < nfrm:
                pixels = getsparse( dset, idx, titles )
                idx += 1 # loop over sparse files in this scan
                nread = nstart + len(pixels['nnz']) # number of frames in this limafile
                g['nnz'][nstart : nread] = pixels['nnz']
                nstart = nread
                pread = pstart + len(pixels['row']) # number of pixels in this limafile
                for name in 'row', 'col', 'intensity':
                    g[name].resize( (pread, ) )
                    g[name][pstart:pread] = pixels[name]
                pstart = pread
            print(scan,end=', ')
        print()
    return outname

    
def main( dsname, outname ):
    dset = ImageD11.sinograms.dataset.load( dsname )
    harvest_masterfile( dset, outname )
    
if __name__=="__main__":
    import sys
    main( sys.argv[1], sys.argv[2] )