import sys
from ImageD11.columnfile import *
from ImageD11.unitcell import unitcell
# put parfile, npeaks and tol and this would be general purpose!
c = columnfile(sys.argv[1])
u = unitcell((3.66,3.66,3.66,90,90,90),"F")
u.makerings(1.)
print u.ringds
tth = [arcsin(0.144088*d/2)*360/pi for d in u.ringds]
m = abs(c.tth - tth[0])< 0.1
for t in tth:
    m = m | ( abs(c.tth - t)<0.1)
c.filter(m)
c.writefile(sys.argv[2])
