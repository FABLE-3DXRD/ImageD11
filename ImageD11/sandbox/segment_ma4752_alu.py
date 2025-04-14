
from __future__ import print_function
# Code to clean the 2D image and reduce it to a sparse array:
# things we might edit

CUT = 2           # keep values abuve cut in first look at image
howmany = 10000   # max pixels per frame to keep
# additional thresholds to test if you have too many pixels (tuple)
thresholds = (4,8,16,32,64,128,256)   
# Remove single pixel spots. How many 8-connected pixels needed in a spot to keep it?
pixels_in_spot = 3
# measurements that might be involved in scanning to transfer to output
SCANMOTORS = ['rot', 'rot_cen360', 'rot_center', 'fpico6', 'dty', 'dty_center']
# scalars from scan/instrument/positioners/XXXX to transfer to output
HEADERMOTORS = ['dty', 'rot', 'pz', 'px', 'py','shtz','shty','shtx' ]
SCANTYPES = ['fscan','f2scan','fscan2d']
DETECTORNAME = 'eiger'
MASKFILE = "/data/id11/nanoscope/Eiger/mask_20210428.edf"
CORES = None # guess the number of cores to use

import sys, os
HNAME = sys.argv[1]


###################### should not need to change much below here

# validate input

assert CUT >= 0
def checkth():
    for i in range(1,len(thresholds)):
        assert thresholds[i] > thresholds[i-1]
checkth()

# we will use multiple processes, each with 1 core
# all sub-processes better do this!
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
os.environ['OMP_NUM_THREADS']   = "1"    
os.environ['NUMBA_NUM_THREADS'] = "1"
os.environ['MKL_NUM_THREADS']   = "1"    

import time
import h5py, hdf5plugin, numpy as np, numba, fabio
        
# pip install ImageD11 --no-deps # if you do not have it yet:
from ImageD11 import sparseframe, cImageD11




@numba.njit
def select( img, msk, row, col, val ):
    #flake8: global CUT
    cut = CUT
    # Choose the pixels that are > cut and put into sparse arrays
    k = 0
    for s in range(img.shape[0]):
        for f in range(img.shape[1]):
            if img[s,f]>cut:
                if msk[s,f]: # skip masked
                    continue
                row[k] = s
                col[k] = f
                val[k] = img[s,f]
                k += 1
    return k

@numba.njit
def top_pixels( nnz, row, col, val, howmany ):
    """
    selects the strongest pixels from a sparse collection
    - thresholds should be a sorted array of potential cutoff values to try
    that are higher than the original cutoff used to select data
    - howmany is the maximum number of pixels to return
    """
    # quick return if there are already few enough pixels
    #flake8: global thresholds
    if nnz <= howmany:
        return nnz
    # histogram of how many pixels are above each threshold
    h = np.zeros(len(thresholds), dtype=np.uint32)
    for k in range(nnz):
        for i,t in enumerate(thresholds):
            if val[k] > t:
                h[i] += 1
            else:
                break
    # choose the one to use. This is the first that is lower than howmany
    tcut = thresholds[-1]
    for n,t in zip( h, thresholds ):
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

    
def choose_parallel( args ):
    """ reads a frame and sends back a sparse frame """
    h5name, address, frame_num = args
    with h5py.File( h5name, "r" ) as h:
        frm = h[address][frame_num]
    row = np.empty(msk.size, np.uint16)
    col = np.empty(msk.size, np.uint16)
    val = np.empty(msk.size, frm.dtype)
    nnz = select( frm, msk, row, col, val)
    #flake8: global howmany
    if nnz == 0:
        sf = None
    else:
        if nnz > howmany:
            nnz = top_pixels( nnz, row, col, val, howmany )
        # Now get rid of the single pixel 'peaks'
        s = sparseframe.sparse_frame( row[:nnz].copy(), col[:nnz].copy(),
                                      frm.shape )
        s.set_pixels("intensity", val[:nnz].copy())
        #flake8: global pixels_in_spot
        if pixels_in_spot <= 1:
            sf = s
        else: 
            # label them according to the connected objects
            sparseframe.sparse_connected_pixels( s, threshold = CUT,
                                                 data_name='intensity',
                                                 label_name='cp')
            # only keep spots with more than 3 pixels ...
            mom = sparseframe.sparse_moments( s, intensity_name='intensity',
                                                 labels_name='cp' )
            npx = mom[:,cImageD11.s2D_1]
            pxcounts = npx[s.pixels['cp']-1]
            pxmsk = pxcounts >= pixels_in_spot
            if pxmsk.sum() == 0:
                sf = None
            else:
                sf = s.mask( pxmsk )
    return frame_num, sf

def segment_scans(fname, scans, outname, mypool,
                  scanmotors = SCANMOTORS,
                  headermotors = HEADERMOTORS,
                  detector = DETECTORNAME):
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
    ndone = 0
    with h5py.File( outname, "a") as hout:
        for scan in scans:
            if scan.endswith(".2"): # for fscans
                continue
            with h5py.File(fname, "r") as hin:
                hout.attrs['h5input'] = fname
                gin = hin[scan]
                bad = 1
                title = hin[scan]['title'][()]
                print('# ',scan,title)
                print('# time now',time.ctime(),"\n#", end =' ')
                for scantype in SCANTYPES:
                    if title.find(scantype) >= 0:
                        bad = 0
                if ('measurement' in gin) and (detector not in gin['measurement']):
                    bad += 1
                if bad:
                    print("# SKIPPING BAD SCAN",scan)
                    continue
                hdr = {}
                g = hout.create_group( scan )
                gm = g.create_group('measurement')
                for m in scanmotors: # vary : many
                    if m in gin['measurement']:
                        gm.create_dataset(m, data = gin['measurement'][m][:] )
                    # else:
                        # the thing requested is not there - ignore the problem for now
                gip = g.create_group('instrument/positioners')
                for m in headermotors: # fixed : scalar
                    gip.create_dataset(m, data = gin['instrument/positioners'][m][()] )
                try:
                    frms = gin['measurement'][detector]
                except:
                    print(list(gin))
                    print(list(gin['measurement']))
                    print(detector)
                    raise
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
                address = scan+"/measurement/"+detector
                nimg = frms.shape[0]
                args = [ (fname, address, i) for i in range(nimg) ]
            # CLOSE the input file here
            #flake8: global CORES
            # approx 1440 frames, 40 cores
            chunksize = max(1,len(args)//CORES//8)
            for i, spf in mypool.map( choose_parallel, args, 
                                         chunksize=chunksize, timeout=60):
                if spf is None:
                    nnz[i] = 0
                    continue
                if i % 500 == 0:
                    print("%4d %d"%(i, spf.nnz), end=",")
                if spf.nnz + npx > len(row):
                    row.resize( spf.nnz + npx, axis=0 )
                    col.resize( spf.nnz + npx, axis=0 )
                    sig.resize( spf.nnz + npx, axis=0 )
                    num.resize( spf.nnz + npx, axis=0 )
                row[npx:] = spf.row[:]
                col[npx:] = spf.col[:]
                sig[npx:] = spf.pixels['intensity']
                num[npx:] = i
                nnz[i] = spf.nnz
                npx += spf.nnz
            ndone += nimg
            print("\n# Done",scan,nimg,ndone)
            sys.stdout.flush()
    return ndone
                
    # the output file should be flushed and closed when this returns
    
def main( hname, outname, mypool ):
    # read all the scans
    print("hname = ",repr(hname))
    print("outname = ",repr(outname))
    with h5py.File(hname,"r") as h:
        scans = list(h['/'])
        print('scans = ', repr(scans))
    print('#',time.ctime())
    start = time.time()
    nfrm = segment_scans(hname, scans, outname, mypool )
    end = time.time()
    print('#',time.ctime())
    print("# Elapsed",end-start,"/s,   f.p.s %.2f"%(nfrm/(end-start)))


if __name__=="__main__":
    
    import concurrent.futures
    
    
    outh5 = os.path.split(HNAME)[-1].replace(".h5","_sparsefull.h5")
    outpath= "/data/visitor/ma4752/id11/segmented"
    OUTNAME = os.path.join( outpath, outh5 )

    # each subprocess gets their own, no need if we fork
    msk = 1 - fabio.open(MASKFILE).data  
    
    if CORES is None:
        try:
            CORES = int(os.environ['SLURM_JOB_CPUS_PER_NODE'])-1
        except:
            CORES = os.cpu_count() - 1
    CORES = max(1,CORES)
    print("# Aiming to use",CORES,"processes")
    with concurrent.futures.ProcessPoolExecutor( max_workers = CORES ) as mypool:
        main( HNAME, OUTNAME, mypool )
    
