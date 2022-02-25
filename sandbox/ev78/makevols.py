
from __future__ import print_function

import fabio, numpy as np, h5py

# read in .txt file and .dat files

def readtxt(fname):
    """
    Reads txt file written by spec and fills in table of:
      filenumber, diffty, diffrz, pico3, current
    """
    cnext = cindex = None
    h = {}
    vals = []
    ctrs = "diffty pico3 notsrcu".split()
    h["pico3"]=1
    h["notsrcu"]=1
    for line in open(fname):
        if line[0] == "#":
            if cnext is not None:
                h[cnext] = float(line.split()[cindex])
                cnext = None
            for c in ctrs:
                if line.find(c)>0:
                    cnext=c
                    cindex=line.split().index(c)
                    break
            continue
        items = line.split()
        if len(items) == 0:
            continue
        fname = items[0]
        angle = 0.5*(float(items[4]) + float(items[5]))
        #print fname, fabio.getnum(fname), diffty, angle, pico3
        vals.append([fabio.getnum(fname), angle] + [h[c] for c in ctrs])
        # print vals[-1]
    return np.array(vals)
                
def decidegrid(vals):
    """
    Figure out a position and angle step from vals array
      number, angle, ypos
    """
    y0 = vals[0,2] # y of first image
    for i,y in enumerate(vals[:,2]):
        if y!=y0: break
    blocksize = i
    if len(vals) % blocksize != 0:
        print("Blocksize seems to be",blocksize, end=' ')
        print("you have",len(vals),"images. Problem.")
    vals.shape = vals.shape[0]/blocksize,blocksize,vals.shape[1]
    # Find the best angular grid mapping onto this
    a0=vals[0,:,1]
    a1=vals[1,:,1]
    offseterr=[abs(a0[:-i]-a1[:-i][::-1]).sum() for i in range(1,20)]
    return vals, np.argmin(offseterr)

def myloadtxt(fname, header=17):
    return np.array([float(l.split()[1]) for l in open(fname).readlines()[header:]])

def readdat(stem,vp):
    fname = "%s%04d.dat"%(stem,vp[0,0,0])
    dat = np.loadtxt( fname )
    xaxis = dat[:,0]
    fullarray = np.zeros((vp.shape[0],vp.shape[1],len(dat)), np.float32)
    n=vp.shape[0]
    # Normalise to srcur here:
    srcurmean = vp[:,:,-1].mean()
    for i in range(n):
        sys.stdout.write("%4d / %4d \r"%(i,n))
        sys.stdout.flush()
        for j in range(vp.shape[1]):
            fname="%s%04d.dat"%(stem,vp[i,j,0])
            scale = srcurmean/vp[i,j,-1]
            fullarray[i,j,:] = myloadtxt( fname ) * scale
    return xaxis, fullarray

def savearray(ar,name,grp):
    if name not in grp:
        dataset = grp.create_dataset(
            name=name,
            shape=ar.shape,
            dtype="float32")
        dataset[:] = ar[:]

def writehdf(fname, blob, tth, vp):
    hdf5path = "DiffTomo/NXdata/sinogram"
    dirname = "DiffTomo/NXdata"
    basename= "sinogram"
    h = h5py.File(fname)
    group = h.require_group(dirname)
    savearray(tth, "zaxis", group)
    savearray(vp[:,:,2].mean(axis=1), "yaxis", group)
    savearray(vp[:,:,1].mean(axis=0), "xaxis", group)
    savearray(vp, "valsarray", group)
    savearray(blob, basename, group)
    h.close()




if __name__=="__main__":
    import sys, os
    HOME = "/data/visitor/ev78/id11"
    os.chdir(HOME)
    stem = sys.argv[1]
    print(stem)
    txtfile = "%s/%s.txt"%(stem,stem)
    # file with motor positions and image filenames
    vals = readtxt(txtfile)
    # figure array dimensions
    vals, offset = decidegrid(vals)
    print(vals.shape)
    # Read 1D integrated data (pyFAI was already run on it)
    xaxis, blob = readdat( "%s/%s"%(stem,stem), vals)

    # Clipping the extra patterns due to scan problem (too many images)
    b=blob[:,:-offset].copy()
    # Throw out the last offset images (from decidegrid code)
    b[0::2]=blob[0::2,:-offset]
    # Flip every other one as we scan 0->180, 180->0
    b[1::2]=blob[1::2,:-offset][:,::-1]
    v = vals[:,:-offset].copy()
    v[0::2] = vals[0::2,:-offset]
    v[1::2] = vals[1::2,:-offset][:,::-1]
    outname = os.path.join("/data/visitor/ev78/id11/voldata",stem)+".hdf"
    writehdf(outname,  b, xaxis, v)
