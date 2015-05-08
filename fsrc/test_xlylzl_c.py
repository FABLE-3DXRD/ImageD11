

import sys
from ImageD11.columnfile import columnfile
from ImageD11.parameters import parameters
from ImageD11.transform import compute_xyz_lab, detector_rotation_matrix
import cImageD11_wrap

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

        outxyz = np.zeros( (c.nrows, 3), np.float32 )
        cImageD11_wrap.compute_xlylzl_wrap( c.sc.astype(np.float32).copy(), 
                                       c.fc.astype(np.float32).copy(),
                                       pars, r.ravel(), dist,
                                       outxyz)
        error = np.abs(outxyz - xlylzl.T)/np.abs(outxyz+1)
        if  error.mean() > 1e-5:
            print tx,ty,ty
            print outxyz[0,:]
            print xlylzl[:,0]
            print o11,o12,o21,o22,error.mean()

        assert error.mean() < 1e-5

print "ok"


# f2py -m fImageD11 -c fImageD11.f90 --f90flags="-fopenmp" -lgomp -lpthread  --fcompiler=gnu95
