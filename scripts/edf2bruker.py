#!/usr/bin/python
## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py

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
    def __init__(self,darkfile=None,darkoffset=None,floodfile=None,floodmultiplier=None,border=None):
        self.darkfile=darkfile
        self.darkoffset=darkoffset
        self.floodfile=floodfile
        self.border=border
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
            print "No dark file"
            self.darkdata = None
            self.darkimage= None
            self.darkfile = None
        if self.darkoffset is None:
            self.dataoffset = 100.0

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
            print "No flood file"
            self.flooddata = None
            self.floodimage= None
            self.floodfile = None
            self.floodmultiplier = None
            
    def correct(self,data):
        tin = data.dtype.char
        # Start by copying
        cor = data.astype(Float32).copy()
        if self.darkimage is not None:
            cor = cor - self.darkimage
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
        if self.border is not None:
            # set the edges to zero
            b=self.border
            cor[:b,:]=0
            cor[:,:b]=0
            cor[-b:,:]=0
            cor[:,-b:]=0
        cor =  where(cor>0.1,cor,0.) # truncate zero
        
        # print cor[c0,c1]
        return cor.astype(tin)







class edf2bruker:

    def __init__(self,dark,flood,template, darkoffset=100,distance=5.0, border=None, wvln=0.5,omegasign=1.):
        self.distance=distance
        self.omegasign=omegasign
        self.wvln=wvln
        self.darkflood=darkflood(darkoffset=darkoffset,border=border)
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
            sgn = self.omegasign
            om =  float(data_in.header["Omega"])*sgn
            oms=  float(data_in.header["OmegaStep"])*sgn
        except:
            om = 0.
            oms = 0.
        self.putitem("WAVELEN",
                     "WAVELEN:%14.7f%14.7f%14.7f"%(self.wvln,0,0)+" "*(80-14*3-8))
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
        parser.add_option("-e","--extn",action="store",type="string",dest="extn",
                          default="edf",
                          help="Filename extension, probably edf, might be cor")
        parser.add_option("-n","--namestem",action="store", type="string", dest="stem",
                          help="Name of the files up the digits part, eg mydata in mydata0000.edf" )
        parser.add_option("-f","--first",action="store", type="int", dest="first",default=0,
                          help="Number of first file to process, default=0")
        parser.add_option("-l","--last",action="store", type="int", dest="last",default=0,
                          help="Number of last file to process, default=0")
        parser.add_option("-d","--dark",action="store", type="string", dest="dark",
                          help="Dark current")
        parser.add_option("-w","--wavelength",action="store", type="float", dest="wvln",
                          default = 0.5,
                          help="Wavelength")
        parser.add_option("-s","--sign of omega",action="store", type="float", dest="omegasign",
                          default = 1.0,
                          help="Sign of ID11 omega rotation, +1 is right handed")
        parser.add_option("-F","--Flood",action="store", type="string", dest="flood",
#                         default="/data/opid11/inhouse/Frelon2K/Ags_mask0000.edf",
                          help="Flood field")
        parser.add_option("-D","--distance",action="store",type="float", dest="distance",
                          default=5.0,help="Sample to detector distance")
        parser.add_option("-b","--border",action="store",type="int", dest="border",
                          default=0,help="Border of image to zero out")

        parser.add_option("-t","--template",action="store", type="string", dest="template",
                          default = "/data/opid11/inhouse/Frelon2K/brukertemplate.0000")

        options, args = parser.parse_args()


        converter = edf2bruker(options.dark , options.flood , options.template,
                               distance=options.distance,
                               border=options.border,
                               wvln=options.wvln,
                               omegasign=options.omegasign)

        for i in range(options.first, options.last+1):
            filein = opendata.makename( options.stem, i, "."+options.extn )
            fileout = opendata.makename( options.stem+"_bruker_0.", i, "" )
            print filein,
            converter.convert(filein,fileout)
            print fileout
            sys.stdout.flush()


    except:
        parser.print_help()

        raise

