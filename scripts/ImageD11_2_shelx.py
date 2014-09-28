#!/usr/bin/python


#!/bliss/users/blissadm/python/bliss_python/suse82/bin/python


# ImageD11_v01.0 Software for beamline ID11
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
Script for getting integrated intensities from .flt + .ubi file
from the command line
"""


from ImageD11 import columnfile
from math import cos, sin, radians

def lf(pk,tth,eta):
    """
    3dxrd lorentz factor
    """
    return pk/abs(sin(radians(eta)))/abs(cos(radians(tth)))


if __name__=="__main__":
    import sys
    try:
        c = columnfile.columnfile(sys.argv[1])
        assert "h" in c.titles
        assert "sum_intensity" in c.titles
        assert "tth" in c.titles
        assert "eta" in c.titles
        for i in range(c.nrows):
            print ("%4d"*3)%(c.h[i],c.k[i],c.l[i])
            print lf(c.sum_intensity[i], c.tth[i], c.eta[i])
            # Should print the error here, if you know it
    except:
        print "Usage: %s fltfile"
        raise
