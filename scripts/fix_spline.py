#!/usr/bin/env python

from __future__ import print_function

# 
# ImageD11_v1.0 Software for beamline ID11
# Copyright (C) 2005-2009  Jon Wright
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



from six.moves import input
import sys, os
    

"""
Script for repairing use of incorrect spline file
"""
from ImageD11 import columnfile, blobcorrector


def fix_flt( inname, outname , cor ):
    inc = columnfile.columnfile( inname )
    for i in range( inc.nrows ):
        inc.sc[i], inc.fc[i] = cor.correct( inc.s_raw[i], inc.f_raw[i] )
    inc.writefile( outname )

def help():
    print(sys.argv[0] + " columnfile_in columnfile_out splinefile")
    print("use splinefile name perfect to remove distortion")
    print("replaces fc and sc columns with corrected numbers")
        

def mymain():
    # If we are running from a command line:
    inname = sys.argv[1]
    if not os.path.exists(inname) or len(sys.argv) < 4:
        help()
        sys.exit()
    outname = sys.argv[2]
    if os.path.exists(outname):
        if not input("Sure you want to overwrite %s ?"%(outname)
                         )[0] in ['y','Y']:
            sys.exit()
    splinename = sys.argv[3]
    if splinename == 'perfect':
        cor = blobcorrector.perfect()
    else:
        cor = blobcorrector.correctorclass( splinename )

    fix_flt( inname, outname, cor )

if __name__=="__main__":
    mymain()
