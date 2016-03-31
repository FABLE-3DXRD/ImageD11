
import os, sys, numpy as np
from ImageD11.columnfile import columnfile

flippers = {
 "_f0" : [[1,0],[0,1]],
 "_f1" : [[-1,0],[0,1]],
 "_f2" : [[1,0],[0,-1]],
 "_f3" : [[-1,0],[0,-1]],
 "_f4" : [[0,1],[1,0]],
 "_f5" : [[0,-1],[1,0]],
 "_f6" : [[0,1],[-1,0]],
 "_f7" : [[0,-1],[-1,0]],
}
c = columnfile( sys.argv[1] )
splinefile = sys.argv[2]
sf = c.s_raw.copy(), c.f_raw.copy()
for name in flippers.keys():
    newname = sys.argv[1].replace(".flt",name)+".flt"
    newsf = np.dot( flippers[name], sf )
    c.s_raw[:], c.f_raw[:] = newsf
    c.writefile( newname )
    cmd = "fix_spline.py %s %s %s"%(
        newname, newname.replace(".flt","_s.flt"), splinefile ) 
    os.system( cmd )
