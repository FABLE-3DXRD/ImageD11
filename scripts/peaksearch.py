


# ImageD11_v0.4 Software for beamline ID11
# Copyright (C) 2005  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


import time
reallystart=time.time()
from math import sqrt
import sys,glob,os.path
from ImageD11 import blobcorrector
from ImageD11 import opendata
from ImageD11 import connectedpixels
import Numeric

def peaksearch(filename, outputfile, corrector, blobim , thresholds):
   """
   Searches for peaks
   """
   t0=time.time()
   f=open(outputfile,"aq")
   data_object = opendata.openedf(filename)
   picture = data_object.data
   print "%s"%(filename),
   f.write("\n\n# File %s\n"%(filename))
   f.write("# Processed on %s\n"%(time.asctime()))
   f.write("# Spatial correction from %s\n"%(corrector.splinefile))
   f.write("# SPLINE X-PIXEL-SIZE %f\n"%(corrector.xsize))
   f.write("# SPLINE Y-PIXEL-SIZE %f\n"%(corrector.ysize))
   for item in data_object.header.keys():
      try:
         f.write("# %s = %s\n"%(item,data_object.header[item]))
      except KeyError:
         pass
   if blobim.shape != data_object.data.shape:
      raise "Incompatible blobimage buffer for file %s"%(filename)
   for threshold in thresholds:
     npks=0
     np=connectedpixels.connectedpixels(picture,blobim,threshold,verbose=0)
     npi,sum,sumsq,com0,com1,com00,com01,com11=connectedpixels.blobproperties(picture,blobim,np,verbose=0)
     f.write("\n#Threshold level %f\n"%(threshold))
     f.write( "# Number_of_pixels Average_counts    x   y     xc   yc      sig_x sig_y cov_xy\n")

     for  i in range(len(npi)):
       if npi[i]>1:    # Throw out one pixel peaks (div zero)
         npks=npks+1
         n   = npi[i]
         avg = sum[i]/n                             # Average intensity
         si  = sqrt((sumsq[i] - n*avg*avg)/(n-1.))  # Standard dev on intensity
         c0  = com0[i]/sum[i]                       # Centre of mass in index 0
         c1  = com1[i]/sum[i]                       # Centre of mass in index 1
         try:
            c00 = sqrt((com00[i]/sum[i] - c0*c0))
         except:
            c00 = 0. # this means a line of pixels and rounding errors
         try:
            c11 = sqrt((com11[i]/sum[i] - c1*c1))
         except:
            c11 = 0.
         try:
            c01 = (com01[i]/sum[i] - c0*c1)/c00/c11
         except:
            c01=0.
         # Spatial corrections :
         c0c, c1c = corrector.correct(c0,c1)
         s = "%d  %f    %f %f    %f %f    %f %f %f\n"%(n, avg, c0, c1, c0c, c1c,  c00, c11, c01)
         f.write(s)
     print "T=%-5d n=%-5d;"%(int(threshold),npks),
   f.close()
   print " time %f/s"%(time.time()-t0)
   return npks # Number of peaks found


if __name__=="__main__":

   try:
      stem = sys.argv[1]
      outfile = sys.argv[2]
      corrector=blobcorrector.correctorclass(sys.argv[5])
      thresholds = [float(t) for t in sys.argv[6:]]
      first = int(sys.argv[3])
      last = int(sys.argv[4])
      files = ["%s%04d%s"%(stem,i,".edf") for i in range(first,last+1)]
      blobim=Numeric.zeros(opendata.openedf(files[0]).data.shape,Numeric.Int)
      if len(files)==0:
         raise "No files found for stem %s"%(stem)
      files.sort()
      start = time.time()
      print "File being treated in -> out, elapsed time"
      nt=len(thresholds)
      for filein in files:
         peaksearch(filein, outfile, corrector, blobim, thresholds)
   except:
      print "Usage: %s filename  outputfile first last spline threshold [threshold...]"
      raise
      sys.exit()

end=time.time()
t=end-reallystart
print "Total time = %f /s, per image = %f /s"%(t,t/len(files))
