
from __future__ import print_function

import sys
from ImageD11.columnfile import columnfile
from ImageD11.unitcell import unitcell_from_parameters
from ImageD11.parameters import read_par_file
from ImageD11.transform import compute_tth_eta
import numpy as np

# put parfile, npeaks and tol and this would be general purpose!
c = columnfile(sys.argv[1])
p = read_par_file(sys.argv[2])
tol = float(sys.argv[3])
u = unitcell_from_parameters( p )
w = p.get("wavelength")

try:
    tth = c.tth
except:
    tth, eta = compute_tth_eta( (c.sc, c.fc), **p.parameters)

dsmax = 2*np.sin(tth.max()*np.pi/360)/w
u.makerings(dsmax)


tthc = [np.arcsin(w*d/2)*360/np.pi for d in u.ringds]

tthe = np.abs( np.outer( tth, np.ones( len(tthc) )) - tthc)
        

tthe_i = np.argmin( tthe, axis=1 )
dtth = tthe.min(axis=1)
print(dtth, dtth.shape)
if 0:
    import pylab
    pylab.plot( tth, eta, "+", ms=1)
    pylab.plot( tthc, np.zeros(len(tthc)), "|", ms=20)
    pylab.hist(dtth, bins=1000)
    pylab.show()

m = dtth < tol
c.filter(m & (tth<16))
c.writefile(sys.argv[4])
