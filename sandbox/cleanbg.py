
from __future__ import print_function, division
#!/usr/bin/python2.7

import numpy as np, fabio


# Read in an ftomo interlaced scan
#
# Remove the background with an image-by-image filtering
#
# Use a local minimum filter with low value clipping
#
# Handle the odd/even passes separately

def ftomo_series_names( stem, nimages ):
    for i in range(nimages):
    	yield "%s0_%04d.edf" % ( stem, i)
    	yield "%s1_%04d.edf" % ( stem, i)




def cleanbg( stem, nimages, newstem ):
   # How many images to filter over to get bg
   bufsize = 11
   omega = 0.
   omegastep = 180.0/1440

   buffy = np.zeros( (bufsize, 2048, 2048), np.uint16 )
   buffo = []
   # Fill initial buffer to do make initbg
   for i, name in enumerate(ftomo_series_names( stem, nimages )):
       if i == bufsize: 
           break
       o = fabio.open( name )
       o.header["Omega"] = omega
       o.header["OmegaStep"] = omegastep
       omega += omegastep 
       o.name = name
       buffy[i] = o.data.copy()
       buffo.append( o )
       print(i,name)
   bg = np.clip( buffy.min( axis = 0 ), 90, 65535 ) - 90
   # Deal with images up to bufsize/2
   for i in range(bufsize//2):
       outname = newstem + "%04d.edf"%(i)
       print("write",outname, buffo[i].name)
       np.subtract( buffo[i].data , bg, buffo[i].data )
       buffo[i].write( outname )

   # Now go through the bulk of the scan
   for i, name in enumerate(ftomo_series_names( stem, nimages )):
       if i < bufsize  :
           continue
       bp = i%bufsize
       print("open",name,bp)
       o = fabio.open( name )
       o.header["Omega"] = omega
       o.header["OmegaStep"] = omegastep
       omega += omegastep 
       o.name = name
       buffo[bp] = o
       buffy[bp] = o.data.copy()
       bg = np.clip( buffy.min( axis = 0 ), 90, 65535 ) - 90
       # Out image is i - bufsize/2
       op = (i-bufsize//2-1)%bufsize
       np.subtract( buffo[op].data , bg, buffo[op].data )
       outname = newstem + "%04d.edf"%(i-bufsize/2-1)
       print(i,outname, buffo[op].name)
       buffo[op].write( outname )
   n=i
   for i in range(bufsize//2):
       outname = newstem + "%04d.edf"%(n+i-bufsize/2)
       bp = (n+i-bufsize//2)%bufsize
       print("write",outname, buffo[bp].name)
       np.subtract( buffo[bp].data , bg, buffo[bp].data )
       buffo[bp].write( outname )


if __name__ == "__main__":
   import sys, glob, os
#   stem = sys.argv[1]
#   nimages = int(sys.argv[2])
   nimages = 720
#   newstem = sys.argv[3]
   

   for indir in glob.glob( "rubydiff*"):
#   dirs =[ sys.argv[1] , ]
#   for indir in dirs:
       print(indir)
       if not os.path.isdir(indir):

           continue

       outdir = os.path.join("diffproc", indir)
       if not os.path.exists(outdir):
           print("mkdir",outdir)
           os.mkdir(outdir)
       
       # check if last image exists for input:
       lastin = os.path.join(indir,indir+"1_0720.edf")
       lastout = os.path.join(outdir,indir+"1438.edf")
       if not os.path.exists( lastin ):
           print("Missing image",lastin,"skip")
           continue
       if not os.path.exists( lastout ):
           print("Missing result",lastout)
           cleanbg( os.path.join(indir, indir), 
                    720, 
                    os.path.join(outdir,indir) )

       
