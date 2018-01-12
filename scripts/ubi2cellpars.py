#!/usr/bin/env python

from __future__ import print_function

"""
Print cell parameters from a ubi file
"""

from ImageD11 import indexing
import sys, numpy, math

def norm2(v): return math.sqrt(numpy.dot(v,v))

def try_to_reduce(u):
    # print "reducing",u
    u = u.copy()
    for i in range(3):
        v = u[i]
        best = norm2(v), v
        # print "i,best",i,best
        for j in range(i+1,3):
            for sign in [1,-1]:
                tv = v + u[j]*sign
                ts = norm2(tv),tv
                # print ts, ts[0],best[0]
                while ts[0] < best[0]:
                    # print "got better"
                    best = ts
                    # print "new best",best
                    tv = ts[1] + u[j]*sign
                    ts = norm2(tv),tv
        u[i] = best[1]
    #print "returning",u
    return u
            


for ubifile in sys.argv[1:]:
    print("Reading",ubifile, end=' ')
    try:
        ul = indexing.readubis(ubifile)
    except:
        print("error")
        continue
    print() 
    print("i    a       b       c     alpha    beta     gamma    volume")
    cell_all = []
    volume = []
    for u,i  in zip(ul,list(range(len(ul)))):
        print(i, end=' ') 
        cell = indexing.ubitocellpars(u)
        vol = cell[0]*cell[1]*cell[2]*numpy.sin(cell[3]*numpy.pi/180.)*numpy.sin(cell[4]*numpy.pi/180.)*numpy.sin(cell[5]*numpy.pi/180.)
        print("%.5f "*6 % cell + "%.5f" %vol)
        cell_all.append(cell)
        volume.append(vol)
        continue
        #print u
        red = try_to_reduce(u)
        #print red
        #print u
        if (red != u).any():
            print("Reduced ubi")
            print(red)
            print("%.5f "*6 % indexing.ubitocellpars(red))

    cell_all = numpy.asarray(cell_all)
    print("average")
    (a,b,c,alpha,beta,gamma) = numpy.mean(cell_all,axis=0)#,keepdims=False)
    print("%.5f "*7 %(a,b,c,alpha,beta,gamma,numpy.average(volume)))

    
