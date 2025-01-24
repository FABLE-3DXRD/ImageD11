
from __future__ import print_function
import sys
from math import pi, cos, sqrt, exp, log, asin
from numpy import dot
from numpy.linalg import inv

def rad(x): return x*pi/180.

a,b,c,alp,bet,gam = 11.8819,11.8419,16.7686,rad(90.0),rad(90.2079),rad(90.)


mt = [  [ a*a, a*b*cos(gam), a*c*cos(bet)  ],
        [ a*b*cos(gam), b*b, b*c*cos(alp) ],  
        [ a*c*cos(bet), b*c*cos(alp), c*c ]]
t0 = 0.671
wvln = 0.16653
rmt = inv(mt)
for line in open(sys.argv[1]).readlines():
    h = int(line[:4])
    k = int(line[4:8])
    l = int(line[8:12])
    ds = dot((h,k,l), dot(rmt, (h,k,l)))
    i = float(line[12:20])
    s = float(line[20:28])
    tth = 2*asin( sqrt(ds)*wvln/2 )
    fac = ( 1 - exp( log(t0) /cos(tth) ) )/(1 - t0)
    print("%4d%4d%4d%8.2f%8.2f"%(h,k,l,i/fac,s/fac))

