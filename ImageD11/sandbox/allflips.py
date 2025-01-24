
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
dims = (2048,2048)
sf = c.s_raw.copy(), c.f_raw.copy()
for name in list(flippers.keys()):
    newname = sys.argv[1].replace(".flt",name)+".flt"
    newsf = np.dot( flippers[name], sf )
    origin = np.dot( flippers[name], dims ).clip(-np.inf, 0 )
    c.s_raw[:], c.f_raw[:] = (newsf.T + origin).T
    c.writefile( newname )

    cmd = "fix_spline.py %s %s %s"%(
        newname, newname.replace(".flt","_s.flt"), splinefile ) 
    os.system( cmd )

    # Now try to reverse flip
    d = columnfile( newname.replace(".flt","_s.flt") )
    d.s_raw[:] = sf[0] # copy back original
    d.f_raw[:] = sf[1]
    sfc = np.array( (d.sc.copy(), d.fc.copy()) ) 
    newsfc = np.dot( np.linalg.inv( flippers[name] ) , (sfc.T - origin).T ) 
    d.sc[:], d.fc[:] = newsfc
    d.writefile( newname.replace(".flt","_sfl.flt") )
