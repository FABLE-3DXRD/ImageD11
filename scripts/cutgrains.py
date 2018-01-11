#!/usr/bin/env python
from __future__ import print_function


"""
Removes grains from a grain file if they have less than
npks peaks

Usage : cutgrains.py   ubi_in  ubi_out  npks
"""

from ImageD11.grain import read_grain_file, write_grain_file
import sys
try:
    GRAINS = read_grain_file(sys.argv[1])
    NPKS = int(sys.argv[3])
    KEEP = []
    for g in GRAINS:
        if int(g.npks) >= NPKS:
            KEEP.append(g)
    write_grain_file( sys.argv[2], KEEP)
except:
    print(__doc__)
    raise
            
