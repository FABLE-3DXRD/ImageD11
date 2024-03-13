#!/usr/bin/env python

from __future__ import print_function
import os, time
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import h5py, hdf5plugin
import numpy as np, fabio, numba
from ImageD11 import sparseframe, cImageD11


# Code to clean the 2D image and reduce it to a sparse array:

@numba.njit
def select( img, msk, row, col, val, cut=1 ):
    # Choose the pixels that are > cut and put into sparse arrays
    k = 0
    for so in (0,550,1100,1650):
        for s in range(so,so+512):
            for fo in (0,515,1040,1555):
                for f in range(fo,fo+512):
                    if img[s,f]>cut:
                        if msk[s,f]: # skip masked
                            continue
                        row[k] = s
                        col[k] = f
                        val[k] = img[s,f]
                        k += 1
    return k

class chooser( object ):
    def __init__(self, msk, pixels_in_spot = 3, cut = 1):
        """ 
        msk = mask to say which pixels to ignore
        pixels_in_spot = how many pixels touching to keep the spot
        cut = intensity threshold
        """
        self.msk = msk
        self.pixels_in_spot=pixels_in_spot
        self.cut=cut
        self.row = np.zeros(msk.shape, np.uint16).ravel()
        self.col = np.zeros(msk.shape, np.uint16).ravel()
        self.val = np.zeros(msk.shape, np.float32).ravel()
        
    def __call__(self, img):
        # choose the pixels that are != 0 and not masked
        nnz = select( img, self.msk, self.row, self.col, self.val,
                      cut=self.cut )
        # copy these into a "sparse" format
        s = sparseframe.sparse_frame( self.row[:nnz], self.col[:nnz],
                                      img.shape )
        s.set_pixels("intensity", self.val[:nnz])
        # label them according to the connected objects
        sparseframe.sparse_connected_pixels( s, threshold = self.cut,
                                             data_name='intensity',
                                             label_name='cp')
        # only keep spots with more than 3 pixels ...
        mom = sparseframe.sparse_moments( s, intensity_name='intensity',
                                          labels_name='cp' )
        npx = mom[:,cImageD11.s2D_1]
        pxcounts = npx[s.pixels['cp']-1]
        msk = pxcounts > self.pixels_in_spot
        if msk.sum() == 0:
            return None
        r = s.mask( msk )
        sparseframe.sparse_connected_pixels( r, threshold = self.cut, 
                                            label_name='connectedpixels',
                                             data_name='intensity' )
        return r

def segment_scans(fname, scans, outname, sparsify_func,
                  scanmotors = ['rot','rot_center','fpico6'],
                  headermotors = ['dty', 'rot'],
                  detector = 'eiger', ):
    """ Does segmentation on a series of scans in hdf files:
    fname : input hdf file
    scans : input groups in the hdf files [h[s] for s in scans] will be treated
    outname : output hdf file to put the results into
    sparsify_func : function to select which pixels to keep
    headermotors : things to copy from input instrument/positioners to output
    scanmotors : datasets to copy from input measurement/xxx to output
    """
    opts = { 'chunks' : (10000,), 'maxshape' : (None,),
             'compression':'lzf', 'shuffle': True }
    with h5py.File( outname, "a") as hout:
        with h5py.File(fname,"r") as hin:
            for scan in scans:
                if scan.endswith(".2"): # for fscans
                    continue
                print(fname,outname,scan,time.ctime())
                gin = hin[scan]
                bad = 0
                for m in scanmotors:
                    if ('measurement' in gin) and m in gin['measurement']:
                        pass
                    else:
                        bad += 1
                if bad:
                    print("SKIPPING BAD SCAN",scan)
                    continue
                hdr = {}
                g = hout.create_group( scan )
                gm = g.create_group('measurement')
                for m in scanmotors: # vary : many
                    gm.create_dataset(m, data = gin['measurement'][m][:] )
                gip = g.create_group('instrument/positioners')
                for m in headermotors: # fixed : scalar
                    gip.create_dataset(m, data = gin['instrument/positioners'][m][()] )
                    
                frms = gin['measurement'][detector]
                row = g.create_dataset( 'row', (1,), dtype=np.uint16, **opts )
                col = g.create_dataset( 'col', (1,), dtype=np.uint16, **opts )
                # can go over 65535 frames in a scan 
                num = g.create_dataset( 'frame', (1,), dtype=np.uint32,  **opts )
                sig = g.create_dataset( 'intensity', (1,), dtype=frms.dtype, **opts )
                nnz = g.create_dataset( 'nnz', (frms.shape[0],), dtype=np.uint32)
                g.attrs['itype']  = np.dtype(np.uint16).name
                g.attrs['nframes']= frms.shape[0]
                g.attrs['shape0'] = frms.shape[1]
                g.attrs['shape1'] = frms.shape[2]
                npx = 0
                for i,frm in enumerate(frms):
                    spf = sparsify_func( frm )
                    if spf is None:
                        continue
                    if spf.nnz + npx > len(row):
                        row.resize( spf.nnz + npx, axis=0 )
                        col.resize( spf.nnz + npx, axis=0 )
                        sig.resize( spf.nnz + npx, axis=0 )
                        num.resize( spf.nnz + npx, axis=0 )
                    row[npx:] = spf.row
                    col[npx:] = spf.col
                    sig[npx:] = spf.pixels['intensity']
                    num[npx:] = i
                    nnz[i] = spf.nnz
                    npx += spf.nnz
                    if i%100==0:
                        print("%4d %d"%(i, spf.nnz), end=" ")
                        sys.stdout.flush()


def main( hname, msk, outname):
    msk = fabio.open(mskname).data.astype(np.uint8)
    print("Masked pixels",msk.sum(),"should be 269126")
    print("hname",hname)
    print("outname",outname)
    with h5py.File(hname,"r") as h:
        scans = [snum for _,snum in sorted(
            [(float(num),num) for num in list(h)])]
        print(scans)
        print(scans[0])
        try:
            print(list(h[scans[0]]['measurement']))
        except:
            print(list(h[scans[0]]))
            raise
    chf = chooser(msk)
    segment_scans(hname, scans, outname, chf )



if __name__ =="__main__":
    import sys
    hname = os.path.abspath(sys.argv[1])
    mskname = "/data/id11/nanoscope/Eiger/20200925.msk"
    path, fname = os.path.split(hname)
    outpath = path.replace("/id11/", "/id11/analysis/") # unix !!
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    outname = os.path.join( outpath, fname.replace(".h5","_sparse.h5"))
    if os.path.exists( outname ):
        print("Output already exists, exiting")
        print(outname)
        sys.exit()
    print( hname, outname)
    main( hname, mskname, outname )
        


