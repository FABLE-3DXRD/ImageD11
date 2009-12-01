

import specfile, sys
fname =  sys.argv[1] 
s = specfile.Specfile( fname )
outname = fname.replace(".mca","")

print s.list()
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
for i,label in scans: 
    y = s.select(label).mca(1)
    outfname = "%s_%04d.dat"%( outname , i - 1)
    print outfname
    output= open( outfname, "w")
    for j in range(len(y)):
        output.write("%f %d\n"%(A+B*j, y[j]))
    output.close()
    
    
    
