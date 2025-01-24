# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# Some exploration of the intexing in Joel's code. Orient along a scattering vector and then do a 1D search rotating around that vector. The idea is to get to something which searches all of the places something can be found in orientation space. But avoiding the 3D problem for small datasets.

# %matplotlib inline 
#notebook
import numpy as np, pylab as pl
from ImageD11 import indexing, unitcell, grain, transform, cImageD11
from scipy.spatial.transform import Rotation

# Generate some peaks
h,k,l = [x.ravel() for x in np.mgrid[-4:5,-4:5,-4:5]]

uc = unitcell.unitcell( [7.89, 8.910, 9.1011, 92, 97, 99], "P")

orient = Rotation.from_euler("XYZ",(10,20,31)).as_matrix()
ub = np.dot( orient, uc.B )
ubi = np.linalg.inv( ub )

gcalc = np.dot( ub, (h,k,l) )

modg = np.sqrt(gcalc*gcalc).sum(axis=0)

tth, eta, omega = transform.uncompute_g_vectors( gcalc, 0.3 )

pl.plot(tth, eta[0],".")
pl.plot(tth, eta[1],".")

selected_peak = sp = np.argmin(abs(tth - 4.17))
print(tth[sp], eta[0][sp], omega[0][sp], eta[1][sp], omega[1][sp], gcalc[:,sp], h[sp], k[sp], l[sp])

# Now we get to the interesting part. We have simulated some data for a triclinic crystal and picked out some specific peak that we want to work with. The idea is to come up with a "new" way of doing indexing. It is mostly based on Joel Berniers fibre texture thing. 
# - We assume that we have assigned this peak h,k,l indices based on the two theta.
# - We find a rotation which takes unitcell B matrix to put the gvector parallel to this gvector.
# - We then apply rotations around this g-vector to generate new orientations.

gx,gy,gz = gobs = gcalc[:,sp].copy()
h0,k0,l0 = -1,2,0
g0 = np.dot( uc.B, (h0,k0,l0))
# angle between observed reflection and this computed peak:
cosa = np.dot( g0, gobs )/ np.sqrt( np.dot(g0,g0) * np.dot(gobs,gobs ))
ang  = np.arccos(cosa)
vec  = np.cross( g0, gobs )
print(vec, ang)
direc = vec / np.sqrt( np.dot(vec,vec))

u0 = Rotation.from_rotvec( ang * direc).as_matrix()
ub0 = np.dot( u0, uc.B)
gtest = np.dot( ub0, (h0,k0,l0))
print(gtest, gobs)

# Now we want to rotate around gtest. This is a vector in reciprocal space. Ideally we want to apply the rotation in real space, such that dot(UBI, g0) == h0 == ( B-1. U-1 ) is invariant. The vector to rotate around is therefore the g-vector which is in real space(?)
#
#

ubi0 = np.linalg.inv( ub0 )
hkl0 = np.array((h0,k0,l0),float)
unitH = hkl0 / np.sqrt( np.dot(hkl0,hkl0)  )
unitG = gtest / np.sqrt( np.dot(gtest,gtest)  )
angles_to_test = np.linspace( -np.pi, np.pi, 360 )

testmats = [ np.dot( ubi0, Rotation.from_rotvec( unitG*x ).as_matrix() ) for x in angles_to_test ]
rmats    = [ Rotation.from_rotvec( unitG*x ).as_matrix() for x in angles_to_test ]
print(np.dot(testmats[0],gtest))

gobs = gcalc.T.copy()
import timeit
start = timeit.default_timer()
npk = [cImageD11.score( m, gobs, 0.05) for m in testmats] 
print("Scoring takes", timeit.default_timer()-start)

pl.figure()
pl.plot( angles_to_test, npk,"-")

# Conclusion: This appears to work. It takes a while to test all the potential orientations however. 
# Perhaps we can expand the mathematics to try to figure out some contributions a bit quicker?
#
# Our equation is:
# hkl = integer =  UBI0 . rotation_on_g0 . gvobs 
#
# Each g-vector has a component along the g0 axis which will be invariant.
#
# The other components could be found in an axis/angle expansion.
# - gaxis + gs.sin(angle) + gc.cos(angle)
# - h,k,l = ubi0.gaxis + ubi0.gs.sin(angle) + ubi0.gc.cos(angle)
# We can try to find h, k, l in terms of the rotation angle as solutions to this equation?
#
# sin/cos have range -1,1 : gives ranges for hmax, hmin, kmax, kmin, lmax, lmin
#

ux,uy,uz = unitG
ubi0_sin = np.dot(ubi0, np.array( ((0,-uz, uy),(uz, 0, -ux),(-uy, ux, 0))))
ubi0_cos = np.dot(ubi0,np.outer( unitG, unitG ))

tcos = np.dot(ubi0 - ubi0_cos, gobs.T) 
tsin = np.dot(ubi0_sin      , gobs.T) 
tconst = np.dot(ubi0_cos, gobs.T )

pl.figure()
j = 567
pl.plot( angles_to_test, tconst[0,j] + np.cos(angles_to_test)*tcos[0,j] + np.sin(angles_to_test)*tsin[0,j],"-",label='h')
pl.plot( angles_to_test, tconst[1,j] + np.cos(angles_to_test)*tcos[1,j] + np.sin(angles_to_test)*tsin[1,j],"-",label='k')
pl.plot( angles_to_test, tconst[2,j] + np.cos(angles_to_test)*tcos[2,j] + np.sin(angles_to_test)*tsin[2,j],"-",label='l')
pl.legend()

# +
# hvalues
hmin, hmax = tconst[0,j] - np.sqrt( tsin[0,j]**2 + tcos[0,j]**2 ), tconst[0,j] + np.sqrt( tsin[0,j]**2 + tcos[0,j]**2 )
kmin, kmax = tconst[1,j] - np.sqrt( tsin[1,j]**2 + tcos[1,j]**2 ), tconst[1,j] + np.sqrt( tsin[1,j]**2 + tcos[1,j]**2 )
lmin, lmax = tconst[2,j] - np.sqrt( tsin[2,j]**2 + tcos[2,j]**2 ), tconst[2,j] + np.sqrt( tsin[2,j]**2 + tcos[2,j]**2 )
    # http://en.wikipedia.org/wiki/List_of_trigonometric_identities#Linear_combinations
    # a.cos(x) + b.sin(x) = c cos(x+p)
    # with p = atan(-b/a)    
def aclip(x):
    return np.arctan2( np.sin(x),np.cos(x))

assert np.allclose(aclip(angles_to_test), angles_to_test)
    
p = np.arctan2( -tsin[0,j], tcos[0,j] )
sab = np.sign(tcos[0,j])*np.sqrt( tcos[0,j]**2 + tsin[0,j]**2  )
#np.fix rounds towards 0
hr = np.arange(np.fix(hmin), np.fix(hmax)+1)
rhs = hr - tconst[0,j]
sth = np.clip(rhs/sab,-2,2)
a1 = aclip( np.arccos(sth) - p )
a2 = aclip(-np.arccos(sth) - p )
# sin(arccos(x)) = sqrt(1-x**2) 
# sin(arctan(x)) = x/sqrt(1+x**2)
# cos(arccos(x)) = x
# cos(arctan(x)) = 1/sqrt(1+x**2)
# sin(A - B) = sin A cos B - cos A sin B
# cos(A - B) = cos A cos B + sin A sin B
# ... this leads directly to sin/cos without finding angle
#print(h,a1,a2)
pl.plot(a1,hr,"o")
pl.plot(a2,hr,"+")
p = np.arctan2( -tsin[1,j], tcos[1,j] )
sab = np.sign(tcos[1,j])*np.sqrt( tcos[1,j]**2 + tsin[1,j]**2  )
for k in range(int(np.round(kmin)), int(np.round(kmax))+1):
    rhs = float(k) - tconst[1,j]
    sth = rhs/sab
    if sth > -1 and sth < 1:
        a1 = aclip( np.arccos( sth ) - p )
        a2 = aclip(-np.arccos( sth ) - p )
#        print(k,a1,a2)
        pl.plot([a1,],[k,],"o")
        pl.plot([a2,],[k,],"+")
p = np.arctan2( -tsin[2,j], tcos[2,j] )
sab = np.sign(tcos[2,j])*np.sqrt( tcos[2,j]**2 + tsin[2,j]**2  )
for l in range(int(np.round(lmin)), int(np.round(lmax))+1):
    rhs = float(l) - tconst[2,j]
    sth = rhs/sab
    if sth > -1 and sth < 1:
        a1 = aclip(np.arccos( sth ) - p )
        a2 = aclip(- np.arccos( sth ) - p )
        #print(l,a1,a2)
        pl.plot([a1,],[l,],"o")
        pl.plot([a2,],[l,],"+")

# -


