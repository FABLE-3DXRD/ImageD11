#!/usr/bin/env python
from __future__ import print_function

from ImageD11.columnfile  import columnfile
from ImageD11.grain import read_grain_file, write_grain_file
import sys, os, multiprocessing

try:
    c = columnfile( sys.argv[1] )
    g = read_grain_file( sys.argv[2] )
    parfile = sys.argv[3]
    cmds = []
except:
    print( "Usage: %s colfile.flt.new grains.map parameters.par  --omega_slop=1 etc"%(sys.argv[0]))
    sys.exit()
    

for i in range(len(g)):
    #g[i].translation[2] = 0.0
    write_grain_file("%d.ubi"%(i),[g[i]])
    d = c.copy()
    d.filter( d.labels == i )
    d.writefile("%d.flt"%(i))
    cmd = "fitgrain.py -p %s -u %d.ubi -U %d.ubi -P %d.par -f %d.flt -x t_z"%(
        parfile,i,i,i,i)
    for extra_arg in sys.argv[4:]:
        cmd += " "+extra_arg
    cmds.append( cmd )

p = multiprocessing.Pool( multiprocessing.cpu_count() )
p.map( os.system, cmds )

