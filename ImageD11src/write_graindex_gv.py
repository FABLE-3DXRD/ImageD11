
from __future__ import print_function

# ImageD11_v0.4 Software for beamline ID11
# Copyright (C) 2005  Jon Wright and Soren Schmidt
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


"""
Conversion of Soren's matlab script for making the file final.log
which is suitable for input into graindex
"""

import numpy as np

from math import sqrt, pi

def make_ds_list(cell,limit=2.):
    """
    Generates a list of d-spacings
    """
    print("Generating hkls with unit cell:", cell)
    cell.makerings(limit)
    ds_list=[]
    keys = list(cell.ringhkls.keys()) # d-spacings
    keys.sort()
    ptype=0
    for ky in keys:
        ptype=ptype+1
        ds = ky
        hklmax = -1e10
        hmax = -1e10
        kmax = -1e10
        hkls_in_ring = cell.ringhkls[ky]
        hkls_in_ring.sort()
        hkls_in_ring.reverse()
        print(hkls_in_ring)
        for h,k,l in hkls_in_ring:
            if h+k+l >= hklmax and h >= hmax and k>kmax:
                hklmax = h+k+l
                hmax = h
                kmax = k
                ds_string="(%d%d%d)"%(h,k,l)
        ds_list.append([ds,ptype,ds_string])
    return ds_list

def get_ds_string(g,ds_list):
    """
    Attempt to emulate the graindex (hkl) syntax.
    Obviously you would not want a (10,0,0) peak!
    """
    length = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2])
    # this is 1/d
    min_diff = abs(length - ds_list[0][0])
    ptype = ds_list[0][1]
    ds_string = ds_list[0][2]
    for item in ds_list[1:]:
        diff = abs(length - item[0])
        if diff < min_diff:
            min_diff=diff
            ptype = item[1]
            ds_string = item[2]
    return ds_string, ptype


def write_graindex_gv(outfilename,gv,tth,eta,omega,intensity,unitcell):
    """
    call with array of gvectors, tth, eta and omega vals

    Using an indexing object to get the ring assignments
    """
    outputfile = open(outfilename,"w")
    ds_list=make_ds_list(unitcell) #
    print(ds_list)
    order = np.argsort(tth)
    nr=0
    for i in order:
        nr=nr+1
        ds_string,ptype = get_ds_string(gv[:,i],ds_list)
        outputfile.write("%i %f %f %f %s %i 0 0 0 %f %f %f 0 0 1 %.2f\n"%(
           nr,           # line number
           omega[i],    # Omega angle
           eta[i],      # eta angle
           tth[i],      # two_theta angle
           ds_string,   # hkl in format (111) or (222) etc
           ptype,       # assignment to hkl peaks
           gv[0,i]*2*pi, # The g-vector
           gv[1,i]*2*pi,
           gv[2,i]*2*pi,
           intensity[i]
           ))
    outputfile.close()
