
from __future__ import print_function

import glob, sys
import numpy as np
from ImageD11 import grain, sym_u
import xfab.symmetry


# open up a set of ubi files

ubifl = sys.argv[1]
symmetry = sys.argv[2]
angtol = float(sys.argv[3])


symnum = {
    "triclinic": 1,
    "monoclinic":2,
    "orthorhombic":3, 
    "tetragonal":4,
    "trigonal":5,
    "hexagonal":6,
    "cubic":7 } [ symmetry ]


h = getattr(sym_u, symmetry )() # cubic()

gl = grain.read_grain_file(ubifl)

pks = list(set([ tuple(np.dot(o, hkl)) for o in h.group for hkl in 
                 [(1,0,0),(1,1,0),(1,1,1),
                  (2,1,0),(2,1,1),(3,2,1),
              ] ]))
frids = []
for h,k,l in pks:
    if (h<=0) and (k<=0) and (l<=0) and (-h,-k,-l) in pks:
        pks.remove( (h,k,l) )
            
dsu = [ (h*h+k*k+l*l, (h,k,l)) for h,k,l in pks ]
dsu.sort()

hkls = np.transpose( [p[1] for p in dsu] )
npks = len(pks)

# pre-declare arrays to try to make it a bit quicker when running

gv1  = np.dot( gl[0].ub,  hkls )
gv2  = np.dot( gl[0].ubi, hkls )
hklr  = np.dot( gl[0].ub, gv1 )
hkli  = hklr.copy()
rmodg2 = 1/(gv1*gv1).sum(axis=0)
angerr2 = np.zeros( npks, float )



tol = 2.35*angtol*np.pi/180.0
tol2 = tol*tol
print( "#", end = ' ')
for i,g1 in enumerate(gl):
    np.dot( g1.ub, hkls, out=gv1 )
    for j in range(i+1, len(gl)):
        g2 = gl[j]
        # hkl = ubi.gv
        np.dot( g2.ubi, gv1, out=hklr )
        # hkl(int)
        np.round( hklr, out=hkli )
        # gcalc = ub.hkli
        np.dot( g2.ub, hkli, out=gv2 )
        # error in gcalc
        np.subtract( gv2, gv1, out=gv2 ) # overwrite in place
        # squared
        np.multiply( gv2, gv2, out=gv2 ) 
        # |error|
        np.sum( gv2, axis=0, out=angerr2 )
        # divided by |g|*|g|
        np.multiply( angerr2, rmodg2, angerr2 )
        # Angle squared in radians
        scor  = np.exp(-angerr2/tol2)
        tscor = scor.mean()
        if tscor > 0.1:
            dist = np.sqrt( ((g1.translation-g2.translation)**2).sum() )
            if dist > 500:
                continue
            print("\n# grains %d %d  frac %.3f scor %.3f"%(i,j,tscor,scor.sum()),end=' ')
            rotangs = xfab.symmetry.Umis( g1.u, g2.u, 7 )
            angs = [a[1] for a in rotangs]
            dist = np.sqrt( ((g1.translation-g2.translation)**2).sum() )
            print(" misori %.3f  dist %.1f"%(np.min(angs), dist))
            print("# g1.t % 7.3f % 7.3f % 7.3f"%tuple(g1.translation))
            print("# g2.t % 7.3f % 7.3f % 7.3f"%tuple(g2.translation))
            print("# i (h k l) -> (h k l)  Angle  Scor ")
            ip = 0
            order = np.argsort(scor)
            for k in order[::-1]:
                if scor[k] > 0.001:
                    h1,k1,l1 = hkls.T[k]
                    h2,k2,l2 = hklr.T[k]
                    ip += 1
                    print( "%d ( % d % d % d ) -> ( % .3f % .3f % .3f ) %.3f %.5f"%( ip,
                        h1,k1,l1,  h2,k2,l2,  np.degrees(np.sqrt(angerr2[k])),scor[k]))
            print("#",end=' ')
                                                             
