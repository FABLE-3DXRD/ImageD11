"""
compute_xlylzl - Function signature:
  compute_xlylzl(s,f,p,r,dist,xlylzl,[n])
Required arguments:
  s : input rank-1 array('f') with bounds (n)
  f : input rank-1 array('f') with bounds (n)
  p : input rank-1 array('d') with bounds (4)
  r : input rank-2 array('d') with bounds (3,3)
  dist : input rank-1 array('d') with bounds (3)
  xlylzl : in/output rank-2 array('d') with bounds (3,n)
Optional arguments:
  n := len(s) input int
"""
import sys
from ImageD11.columnfile import columnfile
from ImageD11.parameters import parameters
from ImageD11.transform import compute_xyz_lab, detector_rotation_matrix
from ImageD11 import fImageD11

import numpy as np, time

start=time.time()
c = columnfile(sys.argv[1])
p = parameters()
p.loadparameters(sys.argv[2])

try:
    pks = [c.sc, c.fc]
except:
    pks = [c.xc, c.yc]

testtilts = [-0.5, 0., 0.24]


tilts = [(x,y,z) for x in testtilts for y in testtilts for z in testtilts]
pars = np.array( ( p.get("z_center"),
                   p.get("y_center"),
                   p.get("z_size"),
                   p.get("y_size")))
dist = np.array( (p.get("distance"),0.,0.))

for o11,o12,o21,o22 in [ [ 1,0,0,1],
                         [-1,0,0,1],
                         [-1,0,0,-1],
                         [ 1,0,0,-1],
                         [ 0,1,1,0],
                         [ 0,-1,1,0],
                         [ 0,-1,-1,0],
                         [ 0,1,-1,0]]:
    for tx,ty,tz in tilts:
        p.parameters['tilt_x']=tx
        p.parameters['tilt_y']=ty
        p.parameters['tilt_z']=tz
        p.parameters['o11']=o11
        p.parameters['o12']=o12
        p.parameters['o21']=o21
        p.parameters['o22']=o22

        xlylzl = compute_xyz_lab( pks,
                          **p.parameters )

        dmat = detector_rotation_matrix(
            p.get('tilt_x'), p.get('tilt_y'), p.get('tilt_z'))
        fmat = np.array( [[ 1,   0,   0],
                          [ 0, p.get('o22'), p.get('o21')],
                          [ 0, p.get('o12'), p.get('o11')]])
        r = np.dot(dmat, fmat)

        outxyz = np.zeros( (c.nrows, 3), np.float )
        fImageD11.compute_xlylzl( pks[0].copy(), 
                                  pks[1].copy(),
                                  pars, r.T, dist,
                                  outxyz.T)
        error = np.abs(outxyz - xlylzl.T)

#        print tx,ty,ty,
#        print o11,o12,o21,o22,error.mean()

        assert error.mean() < 1e-6



# f2py -m fImageD11 -c fImageD11.f90 --f90flags="-fopenmp" -lgomp -lpthread  --fcompiler=gnu95
