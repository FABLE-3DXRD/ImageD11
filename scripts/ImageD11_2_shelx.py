#!/usr/bin/env python

from __future__ import print_function


#
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


import sys
from math import cos, sin, radians
from ImageD11 import columnfile


def lorfac(peak, tth, eta):
    """
    3dxrd lorentz factor
    """
    return peak/abs(sin(radians(eta)))/abs(cos(radians(tth)))


if __name__ == "__main__":
    try:
        CF = columnfile.columnfile(sys.argv[1])
        assert "h" in CF.titles
        assert "sum_intensity" in CF.titles
        assert "tth" in CF.titles
        assert "eta" in CF.titles
        for i in range(CF.nrows):
            print(("%4d"*3) % (CF.h[i], CF.k[i], CF.l[i]))
            print(lorfac(CF.sum_intensity[i], CF.tth[i], CF.eta[i]))
            # Should print the error here, if you know it
    except:
        print("Usage: %s fltfile")
        raise
