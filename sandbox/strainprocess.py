
from __future__ import print_function


import pyFAI, fabio, glob, os, sys, numpy as np, pylab as pl


import fitstrain

HOME="/data/visitor/ma3114/id11/4B_mechanical_test"
poni="/data/visitor/ma3114/id11/Calibration/CeO2_calibrant_scan2.poni"
# spline="/data/visitor/ma3114/id11/Calibration/frelon21_mar16.spline"
mask = fabio.open("/data/visitor/ma3114/id11/fit2d.msk").data

NTTH = 3000
NETA = 3600

ai = pyFAI.load(poni)


todo = glob.glob("scan*")
#todo = ["scan_20161124_160720_",]
todo.sort()



#I,tth,eta = ai.integrate2d( np.ones(mask.shape), NTTH, NETA, mask=mask )
#I.sort(axis=0)
#inds = [I.shape[0]-np.searchsorted(0.1,col)/2 for col in I.T]

# better to compute the tth vals from a unit cell and find them...
PEAKS = [ (1015, 150),
          (1141, 70 ),
          (1437, 150),
          (1617, 70 ),
          (1760, 100),
          (1896, 70 ),
          (2033, 70 ) ]



def fitpeaks( I, tth, eta):
    fits=[]
    for x0,xw in PEAKS:
        cen, sIx, sI, lo, hi = fitstrain.fitstrain( I, x0, xw )
        th = tth[lo:hi]
        mytth = np.interp( cen, np.arange(len(th)), th )
        sol, calc = fitstrain.fitcen( mytth, np.sqrt( sI ), eta )
        err = np.abs(mytth - calc)
        # remove outliers
        cut = np.median(err)*3
        w = np.where( err > cut, 0, np.sqrt(sI) )
        sol, calc = fitstrain.fitcen( mytth, w, eta )
        if 0:
            #print cen, mytth[lo:hi]
            pl.plot(mytth,"+")
            pl.plot(calc)
            pl.show()
        fits.append(sol)
    return fits

def hdr(o):
    cntrs = o.header['counter_mne'].split()
    cntrs.index('admet1')
    cntvs = o.header['counter_pos'].split()
    return "  ".join(( cntvs[cntrs.index('admet1')], 
                       cntvs[cntrs.index('admet2')] ))

for direc in todo:
    os.chdir(os.path.join(HOME,direc))
    edfs = glob.glob("*.edf")
    edfs.sort()
    sumdata = np.zeros((2048,2048),np.float32)
    print("# Doing folder",direc, end=' ')
    for e in edfs:
        np.add( sumdata, fabio.open(e).data, sumdata)
        print(e, end=' ')
        sys.stdout.flush()
    print()
    I, tth , eta = ai.integrate2d( sumdata, NTTH, NETA,
                                   mask = mask,
                                   polarization_factor=0.0,
                                   unit="2th_deg",
                            )
    fits = fitpeaks( I , tth, eta )
    print(direc,hdr(fabio.open(e)), end=' ')
    for f in fits:
        print((len(f)*"  %g")%tuple(f),"    ", end=' ')
    print()
#            I.sort( axis = 0 )
            
