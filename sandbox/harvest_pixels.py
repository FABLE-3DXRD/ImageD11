
from __future__ import print_function, division

import glob, os, multiprocessing, time
import numpy as np
import fabio, h5py
from ImageD11 import sparseframe, cImageD11

def ftomonames( scanpars, folder):
    """ nimages, extn, interlaced, iflip  in scanpars
    Generate the ftomo filename series """
    nums = range(scanpars['nimages'])
    if scanpars['interlaced']:
        stems = [ os.path.join( folder, folder+"0_%04d"+"%s"%(scanpars['extn']) ),
                  os.path.join( folder, folder+"1_%04d"+"%s"%(scanpars['extn']) ) ]
        if scanpars['iflip']:
            s1 = [stems[0]%num for num in nums]
            s2 = [stems[1]%num for num in nums[::-1] ]
        else:
            s1,s2 = [[stem%num for num in nums] for stem in stems]
        # interleave two lists of same length
        return [val for pair in zip(s1,s2) for val in pair]
    else:
        stem = os.path.join( folder, folder+"%04d"+"%s"%(scanpars['extn']) )
        return [stem%num for num in nums]
    
   
def segment( fname, bg, datamem, datapath, block=128, nsigma=3. ):
    """ Does a segmentation of a single file/frame """
    imo = fabio.open(os.path.join( datapath, fname) ) # data
#    print(fname)
    oname = os.path.split( os.path.splitext(fname)[0] )[-1]
    # subtract dark and promote to float32 (without overflow on uint16 I hope)
    cImageD11.uint16_to_float_darksub( datamem.ravel(), bg.ravel(), imo.data.ravel())
    # check offset and noise
    avg,sig = cImageD11.array_mean_var_cut( datamem, cut=nsigma )
    # remove the frelon lines
    datamem.shape = imo.data.shape[0]*(imo.data.shape[1]//block) , block
    cImageD11.frelon_lines( datamem, avg + sig*nsigma )
    datamem.shape = imo.data.shape
    # overwrites datamem (should copy before if you want it
    avg,sig = cImageD11.array_mean_var_cut( datamem, cut=nsigma )
    threshold = avg + sig * nsigma
    # select the pixels to keep:
    rawmask = datamem > threshold
    cleanmask = rawmask.astype(np.int8)
    cImageD11.clean_mask( rawmask.astype(np.int8), cleanmask )
    # save some metadata:
    header = {"filename": fname }
    for k in ("SRCUR","time_of_day"):
        if k in imo.header:
            header[k] = imo.header[k]
    spar = sparseframe.from_data_mask( cleanmask, datamem, header )
    return spar


def write_frames_to_hdf( frames,
                         datapath,
                         fnames,
                         hdfname,
                         hdffolder ):
    """ saves a 1D scan (single ftomo) to a single hdf """
    h = h5py.File( hdfname, "a" )
    g = h.require_group( hdffolder )
    g.attrs['datapath'] = datapath
    for f,name in zip(frames,fnames):
        gname = os.path.split( os.path.splitext(name)[0])[-1]
        f.to_hdf_group( g.require_group( gname ) )
    h.close()


def dofolder( args ):
    """ Processes the files in one folder (single ftomo scan) """
    folder, scanpars = args
    start = time.time()
    fnames = ftomonames( scanpars, folder ) # , nimages, extn, interlaced, iflip )
    hname = os.path.join("hdfs", folder+".hdf")
    bg = fabio.open(
        os.path.join( scanpars['bgfolder'],
                      folder,folder+"_median.edf")
        ).data.astype(np.float32)
    datamem = np.zeros_like( bg )
    spars = [segment(fname, bg, datamem, scanpars['datapath'], block=128)
             for fname in fnames ]
    print("Segmented", folder, time.time()-start)
    write_frames_to_hdf( spars, scanpars['datapath'], fnames, hname, "/" )
    ret = "Wrote %s %.4f /s"%(folder, time.time()-start)
    return ret
 

def harvest_pixels(folders):
    if not os.path.exists("hdfs"):
        os.mkdir( "hdfs" )
    import multiprocessing
    nproc = multiprocessing.cpu_count()//4 # io
    with multiprocessing.Pool( processes = nproc ) as p:
        for x in p.imap_unordered( dofolder, folders ):
            print("Done",x)
        

            
if __name__=="__main__":


    scanpars = {
        "extn" : '.edf',
        "nimages" : 360,
        "interlaced" : True,
        "iflip" : True,
        "datapath" : '/data/id11/nanoscope/Commissioning/2018Apr/difftomo_al_quick',
        "bgfolder" : "/data/id11/jon/difftomo_al_quick/oldmethod",
    }

    folders = [("Al_y%03d_"%(i),scanpars) for i in range(-60,61)]
    
    harvest_pixels( folders )
    

