#!/usr/bin/env python

from __future__ import print_function
from six.moves import input
"""
Applies a new spatial correction based on spd style distortion files
"""
from ImageD11 import columnfile
import fabio
import sys, os
import numpy as np

def spatialfix( cf, dx, dy, flip=""):
    if "V" in flip:
        cf.s_raw = 2048 - cf.s_raw
        print("did a vertical flip")
    if "H" in flip:
        cf.f_raw = 2048 - cf.f_raw
        print("did a horizontal flip")
    # Closest pixel indices
    si = (cf.s_raw+0.5).astype(int).clip(0,2047)
    fi = (cf.f_raw+0.5).astype(int).clip(0,2047)
    # Add shifts
    cf.sc[:] = cf.s_raw + dy[ si, fi ]
    cf.fc[:] = cf.f_raw + dx[ si, fi ]
    return cf

try:
    dx = fabio.open(sys.argv[1]).data
    dy = fabio.open(sys.argv[2]).data
    cf = columnfile.columnfile(sys.argv[3])
except:
    print("Usage dxfile dyfile oldcolfile newcolfile")
    raise

if os.path.exists( sys.argv[4] ):
    print("About to overwrite",sys.argv[4], end=' ')
    response = input("OK? ")
    if response[0] not in "yY":
        print("So I will exit then")
        sys.exit()

if len(sys.argv)>5:
    flip = sys.argv[5]

cf = spatialfix( cf, dx, dy )

cf.writefile( sys.argv[4] )

