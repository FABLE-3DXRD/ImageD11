#!/usr/bin/env python

from __future__ import print_function

# Command line driver only - calls ImageD11.peaksearcher
# ... which can also be called by gui

## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py



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
# MERCHANTABILITY or FITESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  0211-1307  USA



"""
Script for peaksearching images from the command line

Uses the connectedpixels extension for finding blobs above a threshold
and the blobcorrector(+splines) for correcting them for spatial distortion

Defines one function (peaksearch) which might be reused
"""
import h5py
import hdf5plugin # first!
import time

# For benchmarking
reallystart = time.time()

if __name__=="__main__":
    # If we are running from a command line:
    myparser = None
    try:
        from argparse import ArgumentParser
        from ImageD11 import peaksearcher
        parser = ArgumentParser()
        myparser = peaksearcher.get_options(parser)
        options , args = myparser.parse_known_args()
        peaksearcher.peaksearch_driver(options, args)
    except:
        if myparser is not None:
            myparser.print_help()
        print("\n\n And here is the problem:\n")
        raise
end = time.time()
t = end-reallystart
print("Total time = %f /s" % ( t ))

