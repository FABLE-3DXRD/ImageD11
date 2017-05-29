

from ImageD11.grain import read_grain_file , write_grain_file
import sys, numpy as np
gl = read_grain_file( sys.argv[1] )

s = 10.
for g in gl:
    g.translation += (np.random.random(3)-0.5)*s
    g.ubi += (np.random.random((3,3))-0.5)*0.002

write_grain_file( sys.argv[2], gl )
