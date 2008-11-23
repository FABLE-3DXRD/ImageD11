#!/usr/bin/python




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
import sys

if __name__=="__main__":
    try:
        inputdata = openimage.openimage(sys.argv[1])
        step = int(sys.argv[3])
    except:
        print "Usage: %s infile outfile stepsize"%(sys.argv[0])
        print "Puts rows/cols of zeros into image spaced by stepsize"

    for i in range(0,inputdata.data.shape[0],step):
        inputdata.data[i,:]=0
    for i in range(0,inputdata.data.shape[1],step):
        inputdata.data[:,i]=0


    inputdata.write(sys.argv[2], force_type=inputdata.data.dtype)

