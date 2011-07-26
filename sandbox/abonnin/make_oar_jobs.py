
import os
n = 0
i = -12
j = -12
nline = 0
while i < 12.1:
    f = open("%d.sh"%(nline),"w")
    f.write("#!/bin/bash\n")
    while j < 12.1:
        f.write( "./idx_point.sh %f %f %d\n"%(i,j,n))
        n += 1
        j += 0.5
    i += 0.5
    
    j = -12
    f.close()
    os.system("chmod a+x %d.sh"%(nline))
    os.system("oarsub ./%d.sh"%(nline))
    nline += 1
