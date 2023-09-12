from __future__ import print_function, division

from ImageD11 import sym_u
from ImageD11 import grain
import sys, numpy as np, xfab.tools

print ("# Usage: %s grains1.map grains2.map group angletol  distancetol"%(sys.argv[0]))
print ("# Matches grains using symmetry group from:  ")
print ("#  [cubic|hexagonal|trigonal|rhombohedralP|tetragonal|orthorhombic|monoclinic_[a|b|c]|triclinic]")


g1l = grain.read_grain_file(sys.argv[1])
g2l = grain.read_grain_file(sys.argv[2])
for g in g1l + g2l:
    if g.translation is None:
        g.translation = np.zeros(3)

print( "read in the grains")

try:
    h = getattr( sym_u, sys.argv[3] )()
except:
    print( "# No group!")
    print( "#Using cubic")
    h = sym_u.cubic()
    print ("Made a cubic group")

print ( "# Symmetry considered")
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
    #g.u = xfab.tools.ubi_to_u_b(g.ubi)[0]
    assert (abs(np.dot(g.U, g.U.T) - np.eye(3)).ravel()).sum() < 1e-6
    
dtsum = np.zeros(3)
ndt = 0
toldist2 = toldist*toldist

g2t = np.array( [g2.translation for g2 in g2l] )

for i,g1 in enumerate(g1l):
    minangle = 180.0
    best = None
    symubis = [np.dot(o, g1.ubi) for o in h.group]
#    symus   = [xfab.tools.ubi_to_u_b(ubi)[0] for ubi in symubis]
    asymusT   = np.array([xfab.tools.ubi_to_u_b(ubi)[0].T for ubi in symubis])
    printed = False
    trace = np.trace
    pi = np.pi

    dT = g1.translation - g2t
    dT2 = (dT*dT).sum(axis=1)

    inds = np.arange( len(g2t) )[dT2 < toldist2]

    for j in inds:
        g2 = g2l[j]
        
        dt = g1.translation - g2.translation 
        #dt2 = (dt*dt).sum()
        #if dt2 > toldist2:
        #    continue
        
        sg = None
        aumis = np.dot(asymusT, g2.U)
        arg = (aumis[:,0,0]+aumis[:,1,1]+aumis[:,2,2] - 1. )/2.
        asymangle = np.arccos(np.clip(arg, -1, 1))
        sg = np.argmin(asymangle)
        angle = asymangle[sg]*180.0/pi        
        # translation
        if angle < tolangle :
            #print "# Found a match!", h.group[sg].ravel()
            dtsum += dt
            ndt += 1
            print( i,j,"angle  %.4f"%(angle),"distance %.3f"%(np.sqrt(dT2[j])),"\t", end="")
            print( 3*"%.1f "%tuple(dt),"\t",end="")
            print( 3*"%.1f "%tuple(dtsum/ndt))
            printed = True
    if not printed:
        print ("# not matched %d "%(i))



