

"""
Applies a new spatial correction based on spd style distortion files
"""
from ImageD11 import columnfile
import fabio
import sys, os
import numpy as np

try:
    dx = fabio.open(sys.argv[1]).data
    dy = fabio.open(sys.argv[2]).data
    cf = columnfile.columnfile(sys.argv[3])
except:
    print "Usage dxfile dyfile oldcolfile newcolfile"

if os.path.exists( sys.argv[4] ):
    print "About to overwrite",sys.argv[4],
    if raw_input("OK? ")[0] not in 'yY':
        print "So I will exit then"
        sys.exit()

# Closest pixel indices
si = (cf.s_raw+0.5).astype(np.int)
fi = (cf.f_raw+0.5).astype(np.int)

# Add shifts
cf.sc[:] = cf.s_raw + dy[ si, fi ]
cf.fc[:] = cf.f_raw + dx[ si, fi ]

cf.writefile( sys.argv[4] )

