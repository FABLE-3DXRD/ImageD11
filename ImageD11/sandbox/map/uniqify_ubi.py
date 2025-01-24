
from __future__ import print_function
import glob, sys, os

import numpy as np

from ImageD11 import indexing, sym_u, closest


# open up a set of ubi files


ubifl = glob.glob(sys.argv[1]+"/*.ubi")

ubi_all = {}




#print "Make cubic point group"
c = sym_u.cubic()
c.makegroup()

for fname in ubifl:
    #print fname

    # Flip to standard setting
    ubi_all[fname] = [ sym_u.find_uniq_u( ubi, c ) for ubi in
                       indexing.readubis(fname) ]


# Now try to make these into a unique set.
#
# For each ubi:
#     compute the reciprocal peak positions of central lattice
#     see if these match a previous one.

pks = [

        ( 5, 3,-1), ( 5,-3, 1), (-5, 3, 1),
        (-5,-3, 1), (-5, 3,-1), ( 5,-3,-1),
        (-5,-3,-1), ( 5, 3, 1),
        ( 5, 1,-3), ( 5,-1, 3), (-5, 1, 3),
        (-5,-1, 3), (-5, 1,-3), ( 5,-1,-3),
        (-5,-1,-3), ( 5, 1, 3),
        ( 1, 3,-5), ( 1,-3, 5), (-1, 3, 5),
        (-1,-3, 5), (-1, 3,-5), ( 1,-3,-5),
        (-1,-3,-5), ( 1, 3, 5),
        ( 3, 5,-1), ( 3,-5, 1), (-3, 5, 1),
        (-3,-5, 1), (-3, 5,-1), ( 3,-5,-1),
        (-3,-5,-1), ( 3, 5, 1),
        ( 1, 5,-3), ( 1,-5, 3), (-1, 5, 3),
        (-1,-5, 3), (-1, 5,-3), ( 1,-5,-3),
        (-1,-5,-3), ( 1, 5, 3),
        ( 3, 1,-5), ( 3,-1, 5), (-3, 1, 5),
        (-3,-1, 5), (-3, 1,-5), ( 3,-1,-5),
        (-3,-1,-5), ( 3, 1, 5),

        ]
#print len(pks)
hkl = np.array(pks,float).T

uniq_ubis = []
names = list(ubi_all.keys())

tol = 0.25
for name in names:
    this_list = ubi_all[name]
    for ubi, i in zip(this_list,list(range(len(this_list)))):
        gv = np.dot(np.linalg.inv(ubi), hkl)
        seen = 0
        for j in range(len(uniq_ubis)):
            u = uniq_ubis[j][0]
            npk = closest.score(u,gv.T,tol)
            if npk == len(pks):
                # seen you before
                uniq_ubis[j][1].append((name, i))
                uniq_ubis[j][2].append(ubi)
                uniq_ubis[j][0] = np.add.reduce( uniq_ubis[j][2] )/ \
                                            len(uniq_ubis[j][2])
                seen += 1
            #if npk > 50 and npk < len(pks):
            #    print "close...",npk,
        if seen == 0:
            uniq_ubis.append([ubi,[(name,i)], [ubi] ])
        if seen > 1:
            print("uniq seen more than once",ubi,i,name)

print("# Found",len(uniq_ubis),"unique ubi matrices")

dsu = [ (len(uo[1]),(uo[1],uo[0])) for uo in uniq_ubis]
dsu.sort()

for l,uo in dsu[::-1]:
    ubi = uo[1]
    # print "\n# found",len(uo[0]),"times"
    # Find average ubi matrix for each time...
    usum = np.zeros((3,3),float)
    j = 0
    for name, i in uo[0]:
        usum = usum + ubi_all[name][i]
        j += 1
    u = usum/j
    print("%f %f %f\n"  %(u[0,0],u[0,1],u[0,2]), end=' ')
    print("%f %f %f\n"  %(u[1,0],u[1,1],u[1,2]), end=' ')
    print("%f %f %f\n\n"%(u[2,0],u[2,1],u[2,2]), end=' ')

j = 0
for l, uo in dsu[::-1]:
    print("# ", uo[1][0])
    print("# ", uo[1][1])
    print("# ", uo[1][2])
    j += 1
    for name, i in uo[0]:
        n = os.path.split(name)[1]
        x = int(n.split("_")[0][1:])
        y = int(n.split("_")[1].split(".")[0][1:])
        print("%d %d %d"%(x,y,j))

