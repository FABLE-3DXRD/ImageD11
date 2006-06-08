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
Script for making a pseudo-dark file as the minimum of a file series

Defines one class (minimum_image) which might perhaps be reused
"""

import time
# For benchmarking
reallystart=time.time()

import sys,glob,os.path
from ImageD11 import opendata
import Numeric

class minimum_image:
   """
   Used for forming the minimum of a series of images
   """
   def __init__(self,filename=None,image=None):
      self.minimum_image=None
      if image is not None:
         # if an image is supplied we take that as the initial minimum
         self.minimum_image=image
      if filename!=None:
         # if a filename is supplied we take the minimum of the image in that
         # file and the one stored
         self.add_file(filename)
         
   def add_file(self,filename):
      """
      Include another file
      """
      data_object = opendata.openedf(filename)
      picture = data_object.data
      if self.minimum_image is None:
         self.minimum_image = picture.copy()
      else:
         # Check dimensions match
         if self.minimum_image.shape == picture.shape:
            self.minimum_image = Numeric.minimum(self.minimum_image, picture)
         else:
            raise Exception("Incompatible image dimensions")

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
      parser.add_option("-o","--outfile",action="store", type="string", dest="outfile",default="bkg.edf",
                        help="Output filename, default=bkg.edf")
      parser.add_option("-s","--step",action="store", type="int", dest="step",default=1,
                        help="step - every nth image")
      options , args = parser.parse_args()
      stem =        options.stem
      outfile =     options.outfile
      first =       options.first
      last =        options.last
      step = options.step
      # Generate list of files to proces
      files = ["%s%04d%s"%(stem,i,".edf") for i in range(first,last+1,step)]
      outputfile = open(outfile,"wb")
      if len(files)==0:
         raise "No files found for stem %s"%(stem)
      start = time.time()
      mi = minimum_image(files[0])
      for filein in files[1:]:
         print filein
         mi.add_file(filein)
      # finally write out the answer
      # model header + data
      f=open(files[0],"rb")
      fh=f.read(1024)
      while fh.find("}\n")<0: 
         fh+=f.read(1024)
      f.close()
      outputfile.write(fh)
      outputfile.write(mi.minimum_image.astype(Numeric.UInt16).tostring())
      outputfile.close()
   except:
      parser.print_help()
      raise
end=time.time()
t=end-reallystart
print "Total time = %f /s, per image = %f /s"%(t,t/len(files))
