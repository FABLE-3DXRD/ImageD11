
from __future__ import print_function


import sys
from matplotlib.pylab import *

f = open(sys.argv[1],"r")

for line in f.readlines():
    vals = line.split()
    x, y = float(vals[0]), float(vals[1])
    npks = [ int(v) for v in vals[2:] ]
    best = argmax(npks)
    print(x,y,best+1,npks[best])
