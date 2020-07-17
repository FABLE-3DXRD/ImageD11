#!/usr/bin/env python

from __future__ import print_function




# ImageD11_v1.0 Software for beamline ID11
# Copyright (C) 2005-2007  Jon Wright
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
Convert continuous rings to diffraction spots by artificially
setting rows and columns to zero

Note: Trivial script
"""

from fabio import openimage
import numpy as np
import sys

if __name__=="__main__":

    def usage():
        print("Usage: %s infile outfile stepsize"%(sys.argv[0]))
        print("Puts rows/cols of zeros into image spaced by stepsize")
        print("Optionally do radial cuts of ~1 degree:")
        print("Usage: %s infile outfile1 outfile2 ci cj"%(sys.argv[0]))
        sys.exit()
    
    try:
        inputdata = openimage.openimage(sys.argv[1])
    except:
        usage()

    if len(sys.argv)==4:
        try:
            step = int(sys.argv[3])
        except:
            usage()
        inputdata.data[::step,:]=0
        inputdata.data[:,::step]=0
        inputdata.write(sys.argv[2], force_type=inputdata.data.dtype)
        outname = sys.argv[2]
        
    if len(sys.argv)==6:
        ci , cj = float(sys.argv[4]),float(sys.argv[5])

        si, sj = inputdata.data.shape
        i, j = np.mgrid[0:si,0:sj]
        phi = np.arctan2( i - ci, j - cj ) * 361 / np.pi

        m1 = (phi%2).astype(int) == 0
        
        m2 = 1 - m1

        d1 = inputdata.data*m1
        d2 = inputdata.data*m2
        inputdata.data  = d1.astype(inputdata.data.dtype)
        inputdata.write( sys.argv[2] )
        inputdata.data  = d2.astype(inputdata.data.dtype)
        inputdata.write( sys.argv[3] )
