
import glob, sys
import numpy as np
from ImageD11 import grain, sym_u
import xfab.symmetry


# open up a set of ubi files

ubifl = sys.argv[1]
angtol = float(sys.argv[2])
symmetry = sys.argv[3]

symnum = {
    "triclinic": 1,
    "monoclinic":2,
    "orthorhombic":3, 
    "tetragonal":4,
    "trigonal":5,
    "hexagonal":6,
    "cubic":7 } [ symmetry ]


h = getattr(sym_u, symmetry )() # cubic()

#1: Triclinic
#   

gl = grain.read_grain_file(ubifl)

pks = list(set([ tuple(np.dot(o, hkl)) for o in h.group for hkl in 
                 [(1,0,0),(1,1,0),(1,1,1),(1,1,3),(1,2,3)] ] ))
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
angerr2 = np.zeros( npks, np.float )



tol = 2.35*angtol*np.pi/180.0
tol2 = tol*tol
print "#",
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
        tscor = scor.sum()
        if tscor > 4:
            print "\n# grains %d %d  frac %.3f scor %.3f"%(i,j,tscor/npks,tscor),
            rotangs = xfab.symmetry.Umis( g1.u, g2.u, 7 )
            angs = [a[1] for a in rotangs]
            print " misori %.3f"%(np.min(angs))
            ip = 0
            for k in range(npks):
                if scor[k] > 0.25:
                    h1,k1,l1 = hkls.T[k]
                    h2,k2,l2 = hkli.T[k]
                    ip += 1
                    print "%d ( % d % d % d ) -> ( % d % d % d ) %.3f"%( ip,
                        h1,k1,l1,  h2,k2,l2,  np.degrees(np.sqrt(angerr2[k])))
            print "#",
                                                             
