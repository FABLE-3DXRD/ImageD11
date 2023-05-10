
"""WARNING: work in progresss"""
# FIXME!!!! import ImageD11.sinograms.dataset !!!

from __future__ import print_function
import h5py, os, time, numpy as np

# given sample + dset -> make a single file with all the sparse pixels in it.
#  ... this is not really needed, but saves redoing the pixel labelling code



SCANMOTORS = ("diffrz", "diffrz_cen360", "diffrz_center", "fpico4", "fpico3", "diffty", "diffty_center")
HEADERMOTORS = ("diffty", "diffrz", "samtx", "samty", "samtz", "diffry", "samrx", "samry")

SCANMOTORS = ("rot", "rot_cen360", "rot_center", "fpico6", "dty", "dty_center")
HEADERMOTORS = ("dty", "rot", "pz", "px", "py", "shtz", "shty", "shtx")
DETECTORNAME = 'eiger'

def get_sparse_file_name( eigerfile, rawpath ):
    """
    Map the eiger file name onto the sparse peaks filename stored in rawpath
    """
    return os.path.join(rawpath, eigerfile.replace("/eiger","_eiger").replace('.h5','_sparse.h5') )
#    return os.path.join(rawpath, eigerfile.replace(".h5","_sparse.h5"))


def getsparsepath( sample, dsname ):
    return os.path.join( os.getcwd() + '/sparse', sample, dsname )

def getoutpath( sample, dsname ):
    return os.path.join( os.getcwd() + '/sparse', dsname + "_sparse.h5" )

def testready(experiment, sample, dsname):
    ok = []
    with h5py.File( os.path.join(experiment, sample, dsname, dsname + '.h5'), 'r' ) as hin:
        scans = list(hin['/'])
        for scan in scans:
            if scan.endswith('.2'):
                continue
            if ('measurement' not in hin[scan]) or ('eiger' not in hin[scan]['measurement']):
                print("Missing measurement/eiger", sample, dsname, scan)
                return False
                                   
            dso = hin[scan]['measurement/eiger']
            
            for vs in dso.virtual_sources():
                ef = get_sparse_file_name( vs.file_name, getsparsepath(sample, dsname)) 
                done = os.path.exists(ef)
                if not done:
                    print("?",sample, dsname, vs.file_name, ef)
                    return False
    return True

#assert testready("KPLNO","KPLNO_z50")
                                                   

def getpixels( dset, samplename, dsetname ):
    """ 
    dset = h5py DataSet object for the eiger detector
    sparsefilename = converter function that takes eiger filename and gives sparse name
    
    """
    pixels = {}
    titles = 'row','col','intensity','nnz'
    for t in titles:
        pixels[t] = []
    frame0 = 0
    for vs in dset.virtual_sources():
        fsparse = ef = get_sparse_file_name( vs.file_name, getsparsepath(samplename, dsetname)) 
        with h5py.File(  fsparse, "r" ) as hin:
            fds = vs.dset_name
            for name in 'row','col','intensity','nnz':
                pixels[name].append( hin[vs.dset_name][name][:] )# read the data
#            pixels['frame'].append( hin[vs.dset_name]['frame'][:] + frame0 ) # fir
            frame0 += vs.src_space.shape[0]
    for t in titles:
        pixels[t] = np.concatenate(pixels[t])
    return pixels


def segment_scans_assemble_raw( experiment, sample, dsetname,
                   scanmotors=SCANMOTORS,
                   headermotors=HEADERMOTORS,
                   detector=DETECTORNAME ):
    rawpath = getsparsepath( sample, dsetname )
    fname = os.path.join( experiment, sample, dsetname, dsetname+".h5" )
    outname = getoutpath( sample, dsetname )
    if os.path.exists( outname ):
        print("Already done", outname)
        return outname
    print("# output", outname)
    print("# input", fname, rawpath)
    
    opts = {
        "chunks": (10000,),
        "maxshape": (None,),
        "compression": "lzf",
        "shuffle": True,
    }
    ndone = 0
    with h5py.File(outname, "a") as hout:
        hout.attrs["h5input"] = fname
        with h5py.File(fname, "r") as hin:
            scans = list(hin['/'])
            for scan in scans:
                if scan.endswith(".2"):  # for fscans
                    continue
                gin = hin[scan]
                bad = False
                for check in ('title','measurement','measurement/'+detector):
                    if check not in hin[scan]:
                        print(scan,"missing",check,'skipping')
                        bad = True
                if bad:
                    print("Skipping", scan)
                    continue
                title = hin[scan]["title"][()]
                print(scan, end=',')
                g = hout.create_group(scan)
                gm = g.create_group("measurement")
                for m in scanmotors:  # vary : many
                    if m in gin["measurement"]:
                        gm.create_dataset(m, data=gin["measurement"][m][:])
                gip = g.create_group("instrument/positioners")
                for m in headermotors:  # fixed : scalar
                    if "instrument/positioners/%s"%(m) in gin:
                        gip.create_dataset(m, data=gin["instrument/positioners"][m][()])
                try:
                    frms = gin["measurement"][detector]
                    pixels = getpixels( frms, sample, dsetname )
                except:
                    print(list(gin))
                    print(list(gin["measurement"]))
                    print(detector)
                    raise
                # FIXME : can't we do this using external links?
                opts['chunks'] = ( len(pixels['nnz']), )
                g.create_dataset("row", data=pixels['row'], dtype=np.uint16, **opts)
                g.create_dataset("col", data=pixels['col'], dtype=np.uint16, **opts)
                # can go over 65535 frames in a scan
#                g.create_dataset("frame", data=pixels['frame'], dtype=np.uint32, **opts)
                g.create_dataset("intensity", data=pixels['intensity'], dtype=frms.dtype, **opts)
                g.create_dataset("nnz", data=pixels['nnz'], dtype=np.uint32)
                g.attrs["itype"] = np.dtype(np.uint16).name
                g.attrs["nframes"] = frms.shape[0]
                g.attrs["shape0"] = frms.shape[1]
                g.attrs["shape1"] = frms.shape[2]
            print("\n# Done", scan)
    return outname


def main( experiment, samples_dsets ):
        # experiment= '/data/visitor/hc5185/id11/20230505/RAW_DATA'
    for sample , dset in samples_dsets:
        outname = segment_scans_assemble_raw(experiment, sample, "{sample}_{dset}".format(**locals()) )
        
        
        
    