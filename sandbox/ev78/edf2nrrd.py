

import fabio, sys, numpy as np
stem = sys.argv[1]
if 1:
    out=open("%s.bin"%(stem),"wb")
    for i in range(1024):
        im=fabio.open("%s/%s_%04d.edf"%(stem,stem,i)).data
        im = (im+10)*10
        out.write(im.clip(0,255).astype(np.uint8))

im=fabio.open("%s/%s_%04d.edf"%(stem,stem,0))

open("%s.nrrd"%(stem),"w").write("""NRRD0001
# my first nrrd
type: uchar
dimension: 3
sizes: 1024 1024 1024
encoding: raw
data file: %s.bin
"""%(stem))
