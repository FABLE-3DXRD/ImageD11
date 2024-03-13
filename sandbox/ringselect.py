
from __future__ import print_function

import sys
from ImageD11.columnfile import columnfile
from ImageD11.unitcell import unitcell_from_parameters
from ImageD11.parameters import read_par_file
from ImageD11.transform import compute_tth_eta
import numpy as np
import pylab as pl

# put parfile, npeaks and tol and this would be general purpose!
c = columnfile(sys.argv[1])
p = read_par_file(sys.argv[2])
tol = float(sys.argv[3])
tthmax = float(sys.argv[4])
outfile = sys.argv[5]
npx = int(sys.argv[6])
top = int(sys.argv[7])
c.filter( c.Number_of_pixels > npx )

u = unitcell_from_parameters( p )
w = p.get("wavelength")

tth, eta = compute_tth_eta( (c.sc, c.fc), **p.parameters)

dsmax = 2*np.sin(1.03*tthmax*np.pi/360)/w
u.makerings(dsmax)

mask = np.zeros( c.nrows, dtype=bool )
for d in u.ringds:
    tthc = np.arcsin(w*d/2)*360/np.pi
    M = len(u.ringhkls[d])
    print ("hkl",u.ringhkls[d][-1],M, end=" " )
    sel =  abs( tth - tthc ) < tol
    nring = sel.sum()
    print (nring, end=" " )
    if nring == 0:
        print()
        continue
    intensities = c.IMax_int[sel]
    idxsel = np.argsort( intensities )
    npkstokeep = top * M
    idxkeep = np.arange(c.nrows)[sel] 
    print( npkstokeep, len(idxkeep), end=" ")
    if npkstokeep < nring:
        idxs = np.take( idxkeep, idxsel[-npkstokeep:] )
    else:
        idxs = idxkeep 
    print( len(idxkeep),len(idxsel[-npkstokeep:] ), end=" ")
    np.put( mask, idxs, True)
    print( mask.sum() )

pl.plot( c.tth, c.IMax_int, "r." , ms = 1)
    
c.filter( mask )
pl.plot( c.tth, c.IMax_int, "g.", ms=1 )
pl.show()
c.writefile( sys.argv[5] )

