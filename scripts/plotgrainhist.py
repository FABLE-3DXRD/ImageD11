#!/usr/bin/env python

# ImageD11_v1.0 Software for beamline ID11
# Copyright (C) 2008  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  0211-1307  USA





from ImageD11.refinegrains import *

import sys

try:
    flt = sys.argv[1]
    par = sys.argv[2]
    ubi = sys.argv[3]
    tol = float(sys.argv[4])
    if len(sys.argv)>5:
        nbins = int(sys.argv[5])
except:
    print "Usage: %s flt par ubi tol [nbins=30]"%(sys.argv[0])
    sys.exit()

if len(sys.argv)>6:
    o=refinegrains(OmFloat=True, OmSlop=float(sys.argv[6]))
else:
    o=refinegrains(OmFloat=False)
                  

o.loadparameters(par)
o.readubis(ubi)
o.loadfiltered(flt)
o.tolerance=tol
o.generate_grains()
o.assignlabels()
from matplotlib.pylab import *
d = o.scandata[flt]
bins = array(arange(0.0,tol*11./10,tol/nbins))
ng = int(maximum.reduce(d.labels))+1
drl = [ compress(d.labels==i, d.drlv2) for i in range(ng)]
dp5 = [sqrt(di) for di in drl]
hl = [ hist(dpi, bins)[0] for dpi in dp5]
print "hl0:",hl[0].shape, bins.shape
if bins.shape[0] != hl[0].shape[0]:
    bins = (bins[1:] + bins[:-1])/2
cla()

for i in range(ng): 
    plot(bins,hl[i],label=str(i))


    

print " "*10,
for j in range(ng):
    print "%5d"%(j),
print
for i in range(len(bins)):
    print "%10.6f"%(bins[i]),
    for j in range(ng):
        print "%5d"%(hl[j][i]),
    print


show()

