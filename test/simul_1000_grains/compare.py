
from __future__ import print_function
from ImageD11.grain import read_grain_file, write_grain_file
from ImageD11.columnfile import columnfile
import numpy as np

ideal = read_grain_file("ideal.map")

fnames = "f2.map allgrid_fittrans.map allgrid_indexer.map allgrid_makemap.map allgrid_teo.map".split()
#fnames = "teo.map",

def check(gl1, gl2):
    dTs = 0
    dVs = 0
    mxT = 0
    mxTi = 0
    mxV = 0
    mxVi = 0
    for i, (g1, g2) in enumerate(zip( gl1, gl2)):
        dT = g1.translation - g2.translation
        dvol = np.linalg.det(g1.ubi) - np.linalg.det(g2.ubi)
        t = abs(dT).max()
        dTs += t
        dVs += abs(dvol)
        if t > mxT:
            mxTi, mxT = i, t
        if abs(dvol) > mxV:
            mxVi, mxV = i, abs(dvol)
    print( dTs/i, dVs/i, mxT, mxTi, mxV, mxVi )
    


print( "fname avg_pos_err avg_vol_err max_pos_err grain_max_pos_err max_vol_err grain_max_vol_err")
for f in fnames:
    print (f,end=" ")
    check(read_grain_file(f), ideal)



gbad = read_grain_file("shaken.map")[159]
write_grain_file("gbad.ubi",[gbad,])
