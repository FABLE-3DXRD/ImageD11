
from __future__ import print_function

# Generate a series of test cases for index_unknown.py

# Things to vary:
#   Number of grains
#   Symmetry
#   Unit cells all the same - phase pure
#   Mixtures of different grains
#   Data range in omega, two theta




# Use range in terms of cube root of cell volume
# value 1 == first hkl peak should appear
# value 10 == tenth hkl == (10,0,0) should appear
dataranges = [ (0,5),
        (0,10),
        (2,10),
        (5,20) ]

omegaranges = [ 1, 10, 30, 60, 90, 180 ]
from ImageD11 import transform
import random
from math import sqrt, sin
from numpy import dot, array, pi

gen = random.Random( 42 ) 

def omat( gen ):
    ca = gen.random() 
    sa = sqrt(1 - ca*ca)
    rz = array( [ [ ca , sa , 0 ], [ -sa, ca, 0] , [0,0,1]] )
    ca = gen.random() 
    sa = sqrt(1 - ca*ca)
    ry = array( [ [ ca , 0, sa ], [0,1,0 ], [ -sa, 0, ca] ] )
    ca = gen.random() 
    sa = sqrt(1 - ca*ca)
    rx = array( [ [1,0,0], [ 0, ca , sa ], [ 0, -sa, ca]] )
    return dot( rx, dot( ry, rz ) )


# Go straight to a mixture
# Use cases of very large unit cell - fft
#              rather small unit cell - real space

def rcell(gen):
    cell = []
    for i in range(3):
        cell.append( gen.random()*30.0+20 )
    cell.append( 90 )
    if gen.choice([0,1]) == 0:
        cell.append( gen.random()*30.+90 )
    else:
        cell.append(90)
    cell.append( 90 )
    lattice = gen.choice( [ "P", 'C', "P", 'I', "P" ] ) 
    return cell, lattice

from ImageD11 import unitcell
wavelength = None
dsmax = None
dsmin = None
def makedata():
    global wavelength
    global dsmax
    global dsmin
    u = unitcell.unitcell( *rcell(gen) )
    o = omat(gen)
    #print dot(o,o.T)
    ub = dot( o , u.B )

    print (u.tostring())
    meancell = pow(u.B[0,0]*u.B[1,1]*u.B[2,2], -1./3.)
    if dsmax is None:
        dsmax = 20./meancell
        dsmin = 5./meancell

    hkls = u.gethkls(dsmax)
    # print "DSMAX",dsmax, meancell,len(hkls)

    myh = []
    for ds,h in hkls:
        if ds > dsmin:
            g =  dot( ub, h)
     #       print h,ds, g, sqrt(dot( g,g))
            myh.append( h )

    myh = array(myh)
    gve = dot( ub, myh.T )
    # set wavelength to make the largest angle be 45 degrees
    if wavelength is None:
        print ("set wavelength")
        wavelength = 2*sin(pi/8.0)/dsmax
    tth, eta, omega = transform.uncompute_g_vectors( gve, wavelength ) 
    e0,e1 = eta
    o0,o1 = omega

    gfinal = []

    omax = 5 

    for i in range(len(tth)):
        if o0[i] > 0 and o0[i] < omax:
            gfinal.append( gve.T[i] )
        if o1[i] > 0 and o1[i] < omax:
            gfinal.append( gve.T[i] )
    f = open("ideal.ubi","a") 
    f.write("cell = '" + u.tostring() + "'\n")
    f.write('ub = ' +  repr( ub ) )
    f.write("\nnpks = %d\n"%(len(gfinal)))
    f.close()
    return gfinal

def write_gve( gvecs, name):
    f=open(name, "w")
    f.write("1 2 3 4 5 6 P\n")
    f.write("#  gx  gy  gz   xc  omega\n")
    for g in gvecs:
        f.write("%f  %f  %f 0 0 0 0 0\n"%tuple(g))
    f.close()

gv = makedata() + makedata() + makedata() + makedata() + makedata()

write_gve( gv, "test.vecs")
import os, sys
cmd = sys.executable + " " + os.path.join("..","..","scripts","index_unknown.py")
print( cmd)
os.system(cmd + " -g test.vecs" + 
 " -t 0.1 -f 0.05  --fft -s 20 -v 100  -m 40 -r 1. -k 10 -n 256")

