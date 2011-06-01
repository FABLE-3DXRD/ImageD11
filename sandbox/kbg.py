

import numpy, sys, os
from fabio.openimage import openimage
# Switch to conventional optparse and fabio fileseries
stem,first,last,step = sys.argv[1:5]
first = int(first)
last  = int(last)
step  = int(step)

if len(sys.argv) >= 6:
    print "Using start image",sys.argv[5]
    startimage = openimage(sys.argv[5]).data
else:
    startimage = 100.0

# FIXME!!!
# initial estimate should be the dark image
#    this takes care of hot pixels
#
# Noise level should be use defined, or estimated from dark or
# images or historical
# ... probably 16 is only OK for low count levels with no real background
#     it would actually be dominated by shot noise in some cases

class kbg(object):
    def __init__(self, x0, p0, Q=0.01):
        self.x = x0
        self.p = p0
        self.p0= p0
        self.Q = Q

    def update(self, z):
        self.p = self.p + self.Q
        # R = noise estimate
        err = z - self.x
        R = numpy.where( err > self.p0, err, self.p0)
        K = self.p/(self.p + R)
        self.x = self.x + K * err
        self.p = ( 1 - K ) * self.p

# These numbers are for a frelon, mean 100, var = 4*4
obj = kbg( startimage, 16.)

images = range(first,last+1,step)
import random
random.seed(42)
random.shuffle( images )
for i in images:
    fname = "%s%04d.edf"%(stem,i)
    try:
        d = openimage(fname).data
    except:
        print "BAD",fname
        continue
    obj.update( d.astype(numpy.float32) )
#    print fname,"\r",
#    sys.stdout.flush()
#    for i in [0,771294,2075034,1044+2048*221]:
#        print "%d %6d %7.2f %7.2f "%(i,d.flat[i],obj.x.flat[i],obj.p.flat[i]),
#    print


out = openimage("%s%04d.edf"%(stem,first))
testdata = out.data.copy()
out.data = obj.x
out.write("kbg.edf",force_type=numpy.float32)
out.data = obj.p
out.write("kbgv.edf",force_type=numpy.float32)

#from matplotlib.pylab import *
#h,b = numpy.histogram( testdata - obj.x , bins = arange(-100,100))
#plot((b[:-1]+b[1:])*.5,h,"-")
#show()
