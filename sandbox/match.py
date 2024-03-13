
from __future__ import print_function


from ImageD11 import sym_u
from ImageD11 import grain
import sys, numpy as np, xfab.tools

print("# Usage: %s grains1.map grains2.map group angletol  distancetol"%(sys.argv[0]))
print("# Matches grains using symmetry group from:  ")
print("#  [cubic|hexagonal|trigonal|rhombohedralP|tetragonal|orthorhombic|monoclinic_[a|b|c]|triclinic]")


g1l = grain.read_grain_file(sys.argv[1])
g2l = grain.read_grain_file(sys.argv[2])
for g in g1l + g2l:
    if g.translation is None:
        g.translation = np.zeros(3)
try:
    h = getattr( sym_u, sys.argv[3] )()
except:
    print("# No group!")
    print("#Using cubic")
    h = sym_u.cubic()
print("# Symmetry considered")
for o in h.group:
    print(o)

try:
    tolangle = float(sys.argv[4])
except:
    tolangle = 1.0

try:
    toldist = float(sys.argv[5])
except:
    toldist = 100.0

dt0 = np.array([0.,0.,0.])
if len(sys.argv)>=9:
    try:
        dt0 = np.array([float(x) for x in sys.argv[6:9]])
    except:
        raise



da=np.zeros( (len(g1l),len(g2l)), float)
dt2=np.zeros( (len(g1l),len(g2l)), float)

for g in g1l + g2l:
    # g.u = xfab.tools.ubi_to_u_b(g.ubi)[0]
    assert (abs(np.dot(g.u, g.u.T) - np.eye(3)).ravel()).sum() < 1e-6
    
dtsum = np.zeros(3)
ndt = 0
for i,g1 in enumerate(g1l):
    minangle = 180.0
    best = None
    symubis = [np.dot(o, g1.ubi) for o in h.group]
#    symus   = [xfab.tools.ubi_to_u_b(ubi)[0] for ubi in symubis]
    asymusT   = [xfab.tools.ubi_to_u_b(ubi)[0].T for ubi in symubis]
    printed = False
    trace = np.trace
    pi = np.pi
    for j,g2 in enumerate(g2l):
        sg = None
        aumis = np.dot(asymusT, g2.u)
        # print aumis.shape
        # trace(aumis,axis1=1,axis2=2)
        arg = (aumis[:,0,0]+aumis[:,1,1]+aumis[:,2,2] - 1. )/2.
        asymangle = np.arccos(np.clip(arg, -1, 1))
        sg = np.argmin(asymangle)
        angle = asymangle[sg]*180.0/pi
        # translation
        dt = g1.translation - g2.translation + dt0
        da[i,j] = angle
        dt2[i,j] = np.sqrt((dt*dt).sum())
        if angle < tolangle and dt2[i,j] < toldist:
            #print "# Found a match!", h.group[sg].ravel()
            dtsum += dt
            ndt += 1
            print(i,j,"angle  %.4f"%(angle),"distance %.3f"%(dt2[i,j]),"\t", end=' ')
            print(3*"%.1f "%tuple(dt),"\t", end=' ')
            print(3*"%.1f "%tuple(dtsum/ndt))
            #print "# grain ubi,t"
            #print g1.ubi,g1.translation
            #print g2.ubi,g2.translation
            printed = True
    if not printed:
        ja = np.argmin(da[i])
        jd = np.argmin(dt2[i])
        print( "# not matched %d min angle %d  %.4f  %.3f"%(i,ja,da[i,ja],dt2[i,ja]))
        print( "# not matched %d min dist  %d  %.4f  %.3f"%(i,jd,da[i,jd],dt2[i,jd]))



np.save( "da",da)
np.save( "dt2",dt2)
