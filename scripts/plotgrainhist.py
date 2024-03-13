#!/usr/bin/env python

from __future__ import print_function

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





from ImageD11.refinegrains import refinegrains
import numpy as np

import sys

try:
    flt = sys.argv[1]
    par = sys.argv[2]
    ubi = sys.argv[3]
    tol = float(sys.argv[4])
    if len(sys.argv)>5:
        nbins = int(sys.argv[5])
    else:
        nbins = 30
except:
    print("Usage: %s flt par ubi tol [nbins=30] [omega_slop]"%(sys.argv[0]))
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

import matplotlib.pylab as pl
# indexed peaks only
d = o.scandata[flt]
d.filter(d.labels >= 0)

drlv_bins = np.linspace( 0, tol, nbins )
ng = int(d.labels.max())+1
drlv = np.sqrt( d.drlv2 )
dp5 = [ drlv[d.labels==i] for i in range(ng)]
hl = [ np.histogram(dpi, drlv_bins)[0] for dpi in dp5 ]
print("hl0:",hl[0].shape, drlv_bins.shape)
if drlv_bins.shape[0] != hl[0].shape[0]:
    plotbins = (drlv_bins[1:] + drlv_bins[:-1])/2
    
pl.subplot(211)
for i in range(ng): 
    pl.plot(plotbins,hl[i],label=str(i))

if 1:
    print(" "*10, end=' ')
    for j in range(ng):
        print("%5d"%(j), end=' ')
    print()
    for i in range(len(drlv_bins)-1):
        print("%10.6f"%(drlv_bins[i]), end=' ')
        for j in range(ng):
            print("%5d"%(hl[j][i]), end=' ')
        print()
pl.subplot(212)
pl.hist2d( drlv, d.labels, (drlv_bins, np.arange(-0.5,ng,1.)),vmin=0.5)
pl.colorbar()
pl.ylabel("Grain")
pl.xlabel("drlv")
pl.show()

