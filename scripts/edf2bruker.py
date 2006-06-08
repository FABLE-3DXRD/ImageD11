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
A script to convert edf images into bruker format
"""


from Numeric import *
from ImageD11 import opendata

class darkflood:
    def __init__(self,darkfile=None,darkoffset=None,floodfile=None,floodmultiplier=None):
        self.darkfile=darkfile
        self.darkoffset=darkoffset
        self.floodfile=floodfile
        self.floodmultiplier=None
        #
        self.darkimage = None
        self.floodimage = None
        
    def readdark(self,darkfile):
        try:
            self.darkdata = opendata.opendata(darkfile)
            self.darkfile=darkfile
            self.darkimage = self.darkdata.data
            if self.darkoffset is None:
                self.dataoffset = 100.0
        except:
            print "problem, darkfile was",darkfile
            raise
            
    def readflood(self,floodfile):
        try:
            self.flooddata = opendata.opendata(floodfile)
            self.floodfile=floodfile
            self.floodimage = self.flooddata.data
            if self.floodmultiplier is None:
                centre = self.floodimage[100:-100,100:-100]
                npix = centre.shape[0]*centre.shape[1]
                self.floodmultiplier = sum(ravel(centre).astype(Float32))/npix
        except:
            print "problem, floodfile was",floodfile
            raise
            
    def correct(self,data):
        tin = data.typecode()
        # c0 = data.shape[0]/2
        # c1 = data.shape[1]/2
        # print data[c0,c1]
        if self.darkimage is not None:
            cor = data.astype(Float32) - self.darkimage
            # print cor[c0,c1]
        if self.floodimage is not None:
            cor = cor / self.floodimage
            # print cor[c0,c1]
        if self.floodmultiplier is not None:
            cor = cor * self.floodmultiplier
            # print cor[c0,c1]
        if self.darkoffset is not None:
            cor = cor + self.darkoffset
            # print cor[c0,c1]
        cor =  where(cor>0.1,cor,0.) # truncate zero
        # print cor[c0,c1]
        return cor.astype(tin)
    
            




    
class edf2bruker:

    def __init__(self,dark,flood,template, darkoffset=100,distance=5.0):
        self.distance=distance
        self.darkflood=darkflood(darkoffset=darkoffset)
        self.darkflood.readdark(dark)
        self.darkflood.readflood(flood)
        self.templatefile=template
        self.header=opendata.readbrukerheader(template)
        self.h = self.header['headerstring']

    def putitem(self, TITLE, new):
        if len(new)!=80:
            print new
            raise Exception("Sorry, you are about to corrupt the bruker header, len was %d"%(len(new)))
        h=self.h
        p = h.find(TITLE)
        if p!=-1:
            self.h = h[:p]+new+h[p+80:]
        else:
            raise Exception(TITLE+" not found")
        
        
    def convert(self,filein,fileout):
        # Read input file
        data_in = opendata.opendata(filein)
        # Apply dark and flood
        corrected_image = self.darkflood.correct(data_in.data)
        # make new header
        try:
            om = -float(data_in.header["Omega"])
            oms= -float(data_in.header["OmegaStep"])
        except:
            om = 0.
            oms = 0.
        self.putitem("ANGLES",
                     "ANGLES :%14f%14f%14f%14f"%(0,0,om,90)+" "*(80-14*4-8))
        self.putitem("DISTANC",
                     "DISTANC:%14f"%(self.distance)+" "*(80-14-8))
        self.putitem("RANGE",
                     "RANGE  :     %9f"%( abs(oms) ) + " "*58) 
        self.putitem("INCREME:",
                     "INCREME:     %9f"%( oms ) + " "*58 )
        self.putitem("START",
                     "START  :%14f"%( om )+" "*(80-14-8))
        self.putitem("ENDING",
                     "ENDING :%14f%14f%14f%14f"%(0,0, om + oms ,90)+" "*(80-14*4-8))
        self.putitem("NROWS",
                     "NROWS  :%10d"%(corrected_image.shape[0])+" "*62)
        self.putitem("NCOLS",
                     "NCOLS  :%10d"%(corrected_image.shape[1])+" "*62)
        # beam centre and sample to detector distance and wavelength ??
        # x=hs.find("CENTER")
        outf = open(fileout,"wb")
        outf.write(self.h)
        outf.write(transpose(corrected_image).astype(UInt16).tostring())
        outf.close()
        

if __name__=="__main__":
   import sys
   try:
      from optparse import OptionParser
      parser = OptionParser()
      parser.add_option("-n","--namestem",action="store", type="string", dest="stem",
                        help="Name of the files up the digits part, eg mydata in mydata0000.edf" )
      parser.add_option("-f","--first",action="store", type="int", dest="first",default=0,
                        help="Number of first file to process, default=0")
      parser.add_option("-l","--last",action="store", type="int", dest="last",default=0,
                        help="Number of last file to process, default=0")
      parser.add_option("-d","--dark",action="store", type="string", dest="dark",
                        help="Dark current")
      parser.add_option("-F","--Flood",action="store", type="string", dest="flood",
                        default="/data/opid11/inhouse/Frelon2K/Ags_mask0000.edf",
                        help="Flood field")
      parser.add_option("-D","--distance",action="store",type="float", dest="distance",default=5.0,help="Sample to detector distance")

      parser.add_option("-t","--template",action="store", type="string", dest="template",
                        default = "/data/opid11/inhouse/Frelon2K/brukertemplate.0000")

      options, args = parser.parse_args()


      converter = edf2bruker(options.dark , options.flood , options.template, distance=options.distance)
      
      for i in range(options.first, options.last+1):
          filein = opendata.makename( options.stem, i, ".edf" )
          fileout = opendata.makename( options.stem+"_bruker_0.", i, "" )
          print filein,
          converter.convert(filein,fileout)
          print fileout
          sys.stdout.flush()


   except:
       parser.print_help()

       raise

