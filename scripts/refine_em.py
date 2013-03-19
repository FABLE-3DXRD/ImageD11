#!/usr/bin/python

from ImageD11.columnfile  import columnfile
from ImageD11.grain import read_grain_file, write_grain_file
import sys, os

c = columnfile( sys.argv[1] )
g = read_grain_file( sys.argv[2] )

for i in range(len(g)):
    #g[i].translation[2] = 0.0
    write_grain_file("%d.ubi"%(i),[g[i]])
    d = c.copy()
    d.filter( d.labels == i )
    d.writefile("%d.flt"%(i))

    os.system( "fitgrain.py -p %s -u %d.ubi -U %d.ubi -P %d.par -f %d.flt -x t_z"%(sys.argv[3],i,i,i,i))

#    for j in range(3):
#            os.system( "fitgrain.py -u %d.ubi -U %d.ubi -p %d.par -P %d.par -f %d.flt"%(i,i,i,i,i))
