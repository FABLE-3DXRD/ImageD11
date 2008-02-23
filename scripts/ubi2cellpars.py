#!/usr/bin/python

"""
Print cell parameters from a ubi file
"""

from ImageD11 import indexing
import sys

for ubifile in sys.argv[1:]:
    print "Reading",ubifile,
    try:
        ul = indexing.readubis(ubifile)
    except:
        print "error"
        continue
    print 
    for u,i  in zip(ul,range(len(ul))):
        print i, 
        print "%.5f "*6 % indexing.ubitocellpars(u)
    
