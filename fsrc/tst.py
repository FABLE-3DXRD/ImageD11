"""
 fImageD11.compute_gv?
Type:       fortran
String Form:<fortran object>
Docstring:
compute_gv - Function signature:
  gv = compute_gv(xlylzl,omega,wvln,wedge,chi,t,[n])
Required arguments:
  xlylzl : input rank-2 array('f') with bounds (3,n)
  omega : input rank-1 array('f') with bounds (n)
  wvln : input float
  wedge : input float
  chi : input float
  t : input rank-1 array('f') with bounds (3)
Optional arguments:
  n := shape(xlylzl,1) input int
Return objects:
  gv : rank-2 array('f') with bounds (3,n)
"""
import sys
from ImageD11.columnfile import columnfile
from ImageD11.parameters import parameters
from ImageD11.transform import \
    compute_xyz_lab, \
    compute_tth_eta_from_xyz, \
    compute_g_vectors
import ImageD11.transform
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

wvln = float(p.get("wavelength"))
wedge = float(p.get("wedge"))
chi =  float(p.get("chi"))
t = np.array( [float(p.get("t_x")),
               float(p.get("t_y")),
               float(p.get("t_z"))], np.float32)
print "Reading %.4f"%(time.time()-start),
start = time.time()
for i in range(10):
    xlylzl = compute_xyz_lab( pks,
                          **p.parameters )
    tth, eta = compute_tth_eta_from_xyz( xlylzl,
                                     c.omega,
                                     **p.parameters)
    gv = compute_g_vectors(tth, eta, c.omega*p.get('omegasign'),
                       wvln,
                       wedge,
                       chi)

d1 = time.time()-start   
osi = p.get('omegasign')
gvf = np.zeros(xlylzl.shape,np.float, order='F')
start = time.time()
for i in range(100):
    fImageD11.compute_gv( xlylzl, c.omega, osi, wvln, wedge, chi, t, gvf )
d2 = time.time()-start
print "%.4f %.4f"%(d1,d2),"Ratio %.2f"%(10*d1/d2)
err = gv-gvf
scor =abs(err).mean(axis=1)/abs(gv).mean(axis=1)
assert (scor < 1e-6).all(), "Mismatch"

