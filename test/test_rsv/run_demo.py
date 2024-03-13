

import fabio, sys, os, gzip, numpy as np

z = np.zeros( (512,512), np.uint16)
fabio.edfimage.edfimage(z).write( gzip.GzipFile("dark.edf.gz","wb") )
fabio.edfimage.edfimage(z+1).write( gzip.GzipFile("flat.edf.gz","wb") )

# Generates simulated data

os.system("PolyXSim.py -i Fe3O4.inp")
# Reconstruct the data using the given ubi
stem=os.path.join("rsvdemo","rsv")
os.system("rsv_mapper.py -n " + stem + " -f 0 -l 356 -F .edf.gz -O flat.edf.gz -d dark.edf.gz" + 
            " -p " + stem + ".par -u " + stem +  ".ubi -x 4 -b 10 -c 3 -o Fe3O4.h5" )


