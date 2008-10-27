#!/sware/exp/fable/standalone/redhate4-a64/bin/python

from ImageD11.grain import read_grain_file, write_grain_file
import sys
grains = read_grain_file(sys.argv[1])
np = int(sys.argv[3])
keep = []
for g in grains:
    if int(g.npks)>=np:
        keep.append(g)
write_grain_file( sys.argv[2], keep)
