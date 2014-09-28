import glob, sys

import numpy as np

from ImageD11 import closest, indexing, transformer


# open up a set of ubi files

parfilename = sys.argv[1]

ubifl = glob.glob(sys.argv[2])

#for a in glob.glob(sys.argv[3:]):
#   ubifl += glob.glob(a)

ubi_all = {}
for fname in ubifl:
    print fname
    ubi_all[fname] = indexing.readubis(fname)
print "OK"




# Now try to make these into a unique set.
#
# For each ubi:
#     compute the reciprocal peak positions of central lattice
#     see if these match a previous one.

pks = [
        ( 0, 0, 1), ( 0, 1, 0), ( 1, 0, 0),
        ( 0, 0,-1), ( 0,-1, 0), (-1, 0, 0),
        ( 0, 1, 1), ( 1, 1, 0), ( 1, 0, 1),
        ( 0,-1, 1), ( 1,-1, 0), (-1, 0, 1),
        ( 0, 1,-1), (-1, 1, 0), ( 1, 0,-1),
        ( 0,-1,-1), (-1,-1, 0), (-1, 0,-1),
        ( 1, 1, 1), (-1, 1, 1), ( 1,-1, 1),
        ( 1, 1,-1), (-1,-1, 1), (-1, 1,-1),
        ( 1,-1,-1), (-1,-1,-1)
        ]

hkl = np.array(pks,np.float).T

uniq_ubis = []
names = ubi_all.keys()
#dsu = [ (int( np.split(".")[0].split("_")[-1] ), n) for n in names ]
#dsu.sort()
#names = [d[1] for d in dsu]

tol = 0.05
for name in names[15:]:
    this_list = ubi_all[name]
    for ubi, i in zip(this_list,range(len(this_list))):
        gv = np.dot(np.linal.inv(ubi),hkl)
        seen = 0
        for j in range(len(uniq_ubis)):
            u = uniq_ubis[j][0]
            npk=closest.score(u,gv.T,tol)
            if npk == len(pks):
                # seen you before
                uniq_ubis[j][1].append((name, i))
                seen += 1
            if npk > 12 and npk < len(pks):
                print "close...",npk,
        if seen == 0:
            uniq_ubis.append([ubi,[(name,i)]])
        if seen > 1:
            print "uniq seen more than once",ubi,i,name

print "Found",len(uniq_ubis),"unique ubi matrices"


z={}
for line in open("job_z.txt","r").readlines():
    j,zh = line.split()

    z [ int(j)-1 ] = float(zh)

def get_z(name):
    j=int(name.split(".")[0].split("_")[-1])
    return j,z[j]


for entry in uniq_ubis:
    print "\n\n"
    for i in range(3):
        for j in range(3):
            print "# ubi[%d,%d] = %f"%(i,j,ubi[i,j])
    ubi = entry[0]

    for (name,i) in entry[1]:
        j,zh = get_z(name)
        print j,"%7.5f"%(zh),name,i,

        # we have the ubi matrix here.
        # we want to refine and score this ubi against some data
        try:
            t = transformer.transformer()
            t.loadfileparameters(parfilename)
            t.loadfiltered(name.replace(".ubi",""))
            t.computegv()
        except:
            print name
            raise

        # bit of a hack - between 10 and 11 degrees
        h = np.matrixmultiply(ubi, t.gv)
        hint = np.floor(h + 0.5).astype(np.int) # rounds down
        diff = h - hint
        drlv2 = np.sum(diff * diff,0)
        ind = np.compress( drlv2 < tol*tol, np.arange(t.twotheta.shape[0]))
        avg = 4
        npix = 3
        intensity = np.sum((t.finalpeaks[avg,:]*t.finalpeaks[avg,:])[ind])
        print closest.score(entry[0],t.gv.T,tol), intensity


# now re
