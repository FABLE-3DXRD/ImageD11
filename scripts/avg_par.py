#!/usr/bin/env python
from __future__ import print_function

import sys

from ImageD11.parameters import *
allp = []

for i in range(int(sys.argv[2])):
    p = parameters()
    p.loadparameters("%d.par"%(i))
    allp.append(p)

g = parameters()
g.loadparameters("0.par")

import numpy

for pname in "distance tilt_x tilt_y tilt_z wedge y_center z_center".split():
    d = numpy.array( [float(p.get(pname)) for p in allp] )
    order=d.argsort()
    d.sort()
    print(pname,d[1:-1].mean(),d.std(), (d-d.mean())/d.std(),order[0],order[-1])
    g.set(pname, d[1:-1].mean())

for pname in "t_x t_y t_z".split():
    g.set(pname, 0.0)
print("write to",sys.argv[1])
g.saveparameters(sys.argv[1])


