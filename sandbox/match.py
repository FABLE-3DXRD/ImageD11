

from ImageD11 import sym_u
from ImageD11 import grain
import sys, numpy as np, xfab.tools

print "# Usage: %s grains1.map grains2.map group angletol  distancetol"%(sys.argv[0])
print "# Matches grains using symmetry group from:  "
print "#  [cubic|hexagonal|trigonal|tetragonal|orthorhombic|monoclinic_[a|b|c]|triclinic]"


g1l = grain.read_grain_file(sys.argv[1])
g2l = grain.read_grain_file(sys.argv[2])

try:
    h = getattr( sym_u, sys.argv[3] )()
except:
    print "# No group!"
    print "#Using cubic"
    h = sym_u.cubic()
print "# Symmetry considered"
for o in h.group:
    print o

try:
    tolangle = float(sys.argv[4])
except:
    tolangle = 1.0

try:
    toldist = float(sys.argv[5])
except:
    toldist = 100.0



da=np.zeros( (len(g1l),len(g2l)), np.float)
dt2=np.zeros( (len(g1l),len(g2l)), np.float)

for g in g1l + g2l:
    g.u = xfab.tools.ubi_to_u_b(g.ubi)[0]
    assert (abs(np.dot(g.u, g.u.T) - np.eye(3)).ravel()).sum() < 1e-6
    

for i,g1 in enumerate(g1l):
    minangle = 180.0
    best = None
    symubis = [np.dot(o, g1.ubi) for o in h.group]
    symus   = [xfab.tools.ubi_to_u_b(ubi)[0] for ubi in symubis]
    printed = False
    for j,g2 in enumerate(g2l):
        angle = 180.
        sg = None
        #print "g2.ubi",g2.ubi
        for k,symu in enumerate(symus):
            umis = np.dot( symu.T, g2.u )
            if (abs(np.dot(umis, umis.T)- np.eye(3)).ravel()).sum() > 1e-6:
                print i,j,k
                print symubis[k]
                print umis
                1/0
            symangle = np.degrees(np.arccos( min(1, (np.trace( umis ) - 1.0)/2.0)))
            if symangle < angle:
                angle = symangle
                sg = k
        dt = g1.translation - g2.translation
        da[i,j] = angle
        dt2[i,j] = np.sqrt((dt*dt).sum())
        if angle < tolangle and np.sqrt(dt2[i,j]) < toldist:
            #print "# Found a match!", h.group[sg].ravel()
            print i,j,"angle  %.4f"%(angle),"distance %.3f"%(dt2[i,j])
            #print "# grain ubi,t"
            #print g1.ubi,g1.translation
            #print g2.ubi,g2.translation
            printed = True
    if not printed:
        ja = np.argmin(da[i])
        jd = np.argmin(dt2[i])
        print "# not matched %d min angle %d  %.4f  %.3f"%(i,ja,da[i,ja],dt2[i,ja])
        print "# not matched %d min dist  %d  %.4f  %.3f"%(i,ja,da[i,jd],dt2[i,jd])


np.save( "da",da)
np.save( "dt2",dt2)
