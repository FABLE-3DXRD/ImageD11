

import specfile, sys, numpy
from fabio.edfimage import edfimage
fname =  sys.argv[1] 
s = specfile.Specfile( fname )

outname = sys.argv[2]

try:
    first, last = [int(v) for v in  s.list().split(":")]
    scans = [ (i,"%d:%d"%(i,i) )for i in range(first, last+1) ]
except:
    first =1
    last = len(s.list().split(","))
    scans = [ (i, "1.%d"%(i)) for i in range(first,last+1)]

x = None
A = -1.007402
B = 0.010251

np = len(s.select( scans[0][1] ).mca(1))

outar = numpy.zeros( (len(scans), np ), numpy.int32)
j=0
for i,label in scans:
    print i,
    outar[ j ] =  s.select(label).mca(1)
    j += 1

e = edfimage( data = outar, header = { "specfile":"fname" } )

e.write( outname, force_type = numpy.int32 )
    
    
