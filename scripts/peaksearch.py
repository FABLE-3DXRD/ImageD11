#!/bliss/users/blissadm/python/bliss_python/suse82/bin/python


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



"""
Script for peaksearching images from the command line

Uses the connectedpixels extension for finding blobs above a threshold
and the blobcorrector(+splines) for correcting them for spatial distortion

Defines one function (peaksearch) which might be reused
"""

import time
# For benchmarking
reallystart=time.time()

from math import sqrt
import sys,glob,os.path
from ImageD11 import blobcorrector
from ImageD11 import opendata
from ImageD11 import connectedpixels
import Numeric

def peaksearch(filename, outputfile, corrector, blobim , thresholds, dark=None, flood=None):
   """
   Searches for peaks, arguments are:
   filename   : image filename to search for peaks in
   outputfile : ascii file to receive results
   corrector  : ImageD11.blobcorrector object to correct peak positions
   blobim     : Numeric.Int array of same shape as data image to hold blob
                assignments
   thresholds : List of threshold values to use (use more than one if you are
                not sure what to use), any iterable will do, not just list
   """
   t0=time.time()
   f=open(outputfile,"aq") # Open the output file for appending
   # Assumes an edf file for now - TODO - should be independent
   data_object = opendata.openedf(filename)
   picture = data_object.data
   # picture is (hopefully) a 2D Numeric array of type UInt16
   if dark != None:
      picture=picture-dark
   if flood != None:
      picture=picture/flood
   print "%s"%(filename), # Progress indicator
   #
   # Transfer header information to output file
   # Also information on spatial correction applied
   f.write("\n\n# File %s\n"%(filename))
   f.write("# Processed on %s\n"%(time.asctime()))
   try:
      f.write("# Spatial correction from %s\n"%(corrector.splinefile))
      f.write("# SPLINE X-PIXEL-SIZE %f\n"%(corrector.xsize))
      f.write("# SPLINE Y-PIXEL-SIZE %f\n"%(corrector.ysize))
   except:
      pass
   for item in data_object.header.keys():
      try:
         f.write("# %s = %s\n"%(item,data_object.header[item]))
      except KeyError:
         pass
   # Check the blobim matches the data for shape
   if blobim.shape != data_object.data.shape:
      raise "Incompatible blobimage buffer for file %s"%(filename)
   #
   # Now peaksearch at each threshold level
   for threshold in thresholds:
      # 
      npks=0
      #
      # Call the c extension to do the peaksearch, on entry:
      #
      # picture = 2D Numeric UInt16 array (of your data)
      # blobim  = 2D Numeric Int array of same shape as data
      # threshold = float - pixels above this number are put into objects
      # verbose = 0/1 - whether or not to print anything
      #
      np=connectedpixels.connectedpixels(picture,blobim,threshold,verbose=0)
      #
      # One return - np is the number of peaks found (number of connnected objects)
      # picture should be unchanged
      # blobim should be overwritten with integers giving the peak number for that pixel
      #   or zero if the pixel was below the threshold
      #
      # Now get the properties of the identified objects, on entry:
      # picture = data, as before
      # blobim  = peak assignments
      # np      = the number of peaks to find
      # verbose = 0/1 - whether to print or not
      #
      npi,sum,sumsq,com0,com1,com00,com01,com11=connectedpixels.blobproperties(picture,blobim,np,verbose=0)
      #
      # This returns a series of 1D Numeric arrays:
      #  npi = Number of pixels in objects        .....  for each blob pixel   sum 1
      #  sum = sum of intensity in blobs TODO rename this                      sum data[i][j]
      #  sumsq = sum of intensity squared                                      sum data[i][j] * data[i][j]
      #  com0  = sum of intensity multiplied by fast index:                    sum i*data[i][j]
      #  com1  = sum of intensity multiplied by slow index:                    sum j*data[i][j]
      #  com00 = sum of intensity multiplied by fast index squared:            sum i*i*data[i][j]
      #  com01 = sum of intensity multiplied by fast index and slow index:     sum i*j*data[i][j]
      #  com00 = sum of intensity multiplied by slow index squared:            sum j*j*data[i][j]
      #
      # Now write results out for this threshold level
      f.write("\n#Threshold level %f\n"%(threshold))
      f.write( "# Number_of_pixels Average_counts    x   y     xc   yc      sig_x sig_y cov_xy\n")

      for  i in range(len(npi)): # Loop over peaks
         if npi[i]>1:    # Throw out one pixel peaks (div zero)
            npks=npks+1
            n   = npi[i]
            avg = sum[i]/n                             # Average intensity
            si  = sqrt((sumsq[i] - n*avg*avg)/(n-1.))  # Standard dev on intensity
            c0  = com0[i]/sum[i]                       # Centre of mass in index 0
            c1  = com1[i]/sum[i]                       # Centre of mass in index 1
            # Covariances - try except to allow for zeros
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
            # Spatial corrections, c0c and c1c are the distortion corrected centre of mass :
            try:
               c0c, c1c = corrector.correct(c0,c1)
            except:
               c0c, c1c =c0,c1
            s = "%d  %f    %f %f    %f %f    %f %f %f\n"%(n, avg, c0, c1, c0c, c1c,  c00, c11, c01)
            f.write(s)
      print "T=%-5d n=%-5d;"%(int(threshold),npks),
   # Close the output file
   f.close()
   # Finish progress indicator for this file
   print " time %f/s"%(time.time()-t0)
   return npks # Number of peaks found


if __name__=="__main__":
   # If we are running from a command line:
   try:
      from optparse import OptionParser
      parser = OptionParser()
      parser.add_option("-n","--namestem",action="store", type="string", dest="stem",
                        help="Name of the files up the digits part, eg mydata in mydata0000.edf" )
      parser.add_option("-f","--first",action="store", type="int", dest="first",default=0,
                        help="Number of first file to process, default=0")
      parser.add_option("-l","--last",action="store", type="int", dest="last",
                        help="Number of last file to process")
      parser.add_option("-o","--outfile",action="store", type="string", dest="outfile",default="peaks.out",
                        help="Output filename, default=peaks.out")
      parser.add_option("-d","--darkfile",action="store", type="string", dest="dark",default=None,
                        help="Dark current filename, to be subtracted, default=None")
      parser.add_option("-D","--darkfileoffset",action="store", type="int", dest="darkoffset",default=100,
                        help="Constant to subtract from dark to avoid overflows, default=100")
      s="/data/opid11/inhouse/Frelon2K/spatial2k.spline"
      parser.add_option("-s","--splinefile",action="store", type="string", dest="spline",default=s,
                        help="Spline file for spatial distortion, default=%s"%(s))
      parser.add_option("-p","--perfect_images",action="store",type="choice",
                        choices=["Y","N"],default="N",dest="perfect",
                        help="Ignore spline Y|N, default=N")
      parser.add_option("-t","--threshold",action="append", type="float", dest="thresholds", default=None,
                        help="Threshold level, you can have several")
      options , args = parser.parse_args()
      stem =        options.stem
      outfile =     options.outfile
      first =       options.first
      last =        options.last
      if options.spline=="N":
         corrector=blobcorrector.correctorclass(options.spline)
      else:
         corrector=None
      # List comprehension - convert remaining args to floats
      thresholds = [float(t) for t in options.thresholds]
      # Generate list of files to proces
      files = ["%s%04d%s"%(stem,i,".edf") for i in range(first,last+1)]
      # Make a blobimage the same size as the first image to process
      if len(files)==0:
         raise "No files found for stem %s"%(stem)
      blobim=Numeric.zeros(opendata.openedf(files[0]).data.shape,Numeric.Int)
      # Not sure why that was there (I think if glob was used)
      # files.sort()
      if options.dark!=None:
         dark=(opendata.openedf(options.dark).data - options.darkoffset).astype(Numeric.UInt16)
      else:
         dark=None
      start = time.time()
      print "File being treated in -> out, elapsed time"
      for filein in files:
         peaksearch(filein, outfile, corrector, blobim, options.thresholds,dark=dark)
   except:
#      print "Usage: %s filename  outputfile first last spline threshold [threshold...]"%(sys.argv[0])
      print
      # Raise the exception to find out in more detail where we went wrong
      raise
      sys.exit()

end=time.time()
t=end-reallystart
print "Total time = %f /s, per image = %f /s"%(t,t/len(files))
