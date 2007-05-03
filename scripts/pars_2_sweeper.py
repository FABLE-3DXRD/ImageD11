#!/bliss/users/blissadm/python/bliss_python/suse82/bin/python

from ImageD11 import transform, parameters

import sys, Numeric

try:
    parfile = sys.argv[1]
except:
    print "Usage: %s parfile [v h]"
    print "generates parameters for grainsweeper from parameters for ImageD11"
    sys.exit()
try:
    v=float(sys.argv[2])
    h=float(sys.argv[3])
except:
    v=512.
    h=1536/2.

print "Getting parameters from",parfile

p = parameters.parameters()
p.loadparameters(parfile)

pars = p.get_parameters()

xyz = Numeric.transpose(
    transform.compute_tth_eta([[v],[h]],
                          return_pixel_xyz=True,
                          **pars))

ks = pars.keys()
ks.sort()
for k in ks:
    print k,pars[k]
print
print "central pixel position"
print "detector vertical",v,"horizontal",h

print "real space x, y, z = ", [pars['distance'],0,0] + xyz

