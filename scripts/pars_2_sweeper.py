#!/usr/bin/env python

from __future__ import print_function

"""
Adapted to print detector line required for grainsweeper.3D
Jette Oddershede, DTU Physics, jeto@fysik.dtu.dk, April 2013
"""

import sys

import numpy as np

from ImageD11 import transform, parameters


try:
    parfile = sys.argv[1]
except:
    print("Usage: %s parfile v h mm")
    print("generates parameters for grainsweeper from parameters for ImageD11")
    print("v and h are the detector dimension in number of pixels")
    print("mm should be 1 if par in mm and 1e-3 if par in microns")
    sys.exit()
try:
    pixv=float(sys.argv[2])
    pixh=float(sys.argv[3])
    mm=float(sys.argv[4])
except:
    pixv=2048.
    pixh=2048.
    mm=1
v = (pixv-1)/2.
h = (pixh-1)/2.

print("Getting parameters from",parfile)

p = parameters.parameters()
p.loadparameters(parfile)

pars = p.get_parameters()

xyz = np.transpose(
    transform.compute_xyz_lab([[v],[h]],
                          **pars))

ks = list(pars.keys())
ks.sort()
for k in ks:
    print(k,pars[k])
print()
print("central pixel position")
print("detector vertical",v,"horizontal",h)
print("real space x, y, z = ", xyz)

print("\n")
print("detector %i %i %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f stem format %i %i %i %i" % \
    (pixv,
     pixh,
     mm*pars['y_size'],
     mm*pars['z_size'],
     mm*xyz[0,0],
     mm*xyz[0,1],
     mm*xyz[0,2],
     180*pars['tilt_x']/np.pi,
     180*pars['tilt_y']/np.pi,
     180*pars['tilt_z']/np.pi,
     pars['o11'],
     pars['o12'],
     pars['o21'],
     pars['o22']))
