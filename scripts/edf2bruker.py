#!/sware/exp/fable/standalone/redhate4-a64/bin/python
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

import time
import numpy
from fabio.openimage import openimage
from fabio.brukerimage import brukerimage
from fabio import file_series

class darkflood:
    """ apply dark and flood corrections """
    def __init__(self,
                 darkfile = None,
                 darkoffset = None,
                 floodfile = None,
                 floodmultiplier = None,
                 border = None,
                 powfac = 1.0,
                 overflow= 65534,
                 detrend = None
                 ):
        self.overflow=overflow
        self.darkfile = darkfile
        self.darkoffset = darkoffset
        self.floodfile = floodfile
        self.border = border
        self.floodmultiplier = None
        #
        self.darkimage = None
        self.floodimage = None
        self.flmult = None
        self.powfac = powfac
        self.detrend = detrend


    def readdark(self,darkfile):
        """ read the dark"""
        try:
            self.darkdata = openimage(darkfile)
            self.darkfile=darkfile
            self.darkimage = self.darkdata.data.astype(numpy.float32)
            if self.detrend is not None:
                self.do_detrend( self.darkimage )
            if self.powfac != 1.0:
                print "apply 0.96 to dark"
                self.darkimage = numpy.power(self.darkimage, 0.96)
        except:
            print "No dark file"
            self.darkdata = None
            self.darkimage= None
            self.darkfile = None

    def readflood(self, floodfile):
        """ read the flood """
        try:
            self.flooddata = openimage(floodfile)
            self.floodfile = floodfile
            self.floodimage = self.flooddata.data.astype(numpy.float32)
            if self.floodmultiplier is None:
                centre = self.floodimage[100:-100,100:-100]
                npix = centre.shape[0]*centre.shape[1]
                self.floodmultiplier = numpy.sum(numpy.ravel(
                        centre).astype(numpy.float32))/npix
                self.flmult = 1 / (self.floodimage * self.floodmultiplier)
        except:
            print "No flood file"
            self.flooddata = None
            self.floodimage= None
            self.floodfile = None
            self.floodmultiplier = None
            
    def correct(self, data, detrend = None):
        """ correct the data """
        tin = data.dtype.char
        # Start by copying
        cor = data.astype(numpy.float32).copy()
        msk = numpy.where( cor > self.overflow,
                           65534,
                           0 )
        if self.powfac != 1.0:
            # print "applying 0.96"
            numpy.power(cor, 0.96, cor)
        if self.darkimage is not None:
            numpy.subtract(cor, self.darkimage,cor)
            # print cor[c0,c1]
        if self.detrend is not None:
            cor = self.do_detrend( cor )
        if self.flmult is not None:
            numpy.multiply(cor , self.flmult, cor)
            # print cor[c0,c1]
        #if self.darkoffset is not None:
        #    cor = cor + self.darkoffset
        if self.border is not None:
            # set the edges to zero
            b=self.border
            cor[:b,:]=0
            cor[:,:b]=0
            cor[-b:,:]=0
            cor[:,-b:]=0
        # Should we bother with this - yes please - noisy pixels overflow
        cor =  numpy.where(cor>0.1, cor, 0.) # truncate zero
        cor =  numpy.where(msk > 1, msk, cor) # mask overflows
        # print cor[c0,c1]
        return cor.astype(tin)



    def do_detrend( self, ar):
        if self.detrend is None:
            return ar
        print "detrending"
        s = ar.copy()
        np = ar.shape[1]/2
        s[:,:np].sort()
        s[:,np:].sort()
        s1 = ar.copy()
        n = self.detrend
        o1 = s[:,:n].sum(axis=1)/n
        o2 = s[:,np:(np+n)].sum(axis=1)/n    
        s1[:,:np] = (s1[:,:np].T - o1 + o1.mean() ).T
        s1[:,np:] = (s1[:,np:].T - o2 + o2.mean() ).T
        return s1



class edf2bruker:

    def __init__(self, 
                 dark, 
                 flood, 
                 template, 
                 darkoffset = 100,
                 distance = 5.0, 
                 border = None, 
                 wvln = 0.5,
                 omegasign = 1.,
                 powfac = 1.0,
                 overflow=65534,
                 detrend = None ):
        """ converts edf (esrf id11 frelon) to bruker (esrf id11 smart6500)"""
        self.distance = distance
        self.overflow = overflow
        self.omegasign = omegasign
        self.wvln = wvln
        self.powfac = powfac
        self.darkflood = darkflood(darkoffset = darkoffset,
                                   border = border,
                                   powfac = self.powfac,
                                   overflow = self.overflow,
                                   detrend = detrend )
        self.darkflood.readdark(dark)
        self.darkflood.readflood(flood)
        self.templatefile = template
        im = brukerimage()
        im.readheader(template)
        self.header = im.header
        self.h = im.__headerstring__
        self.last_time = time.time()

    def putitem(self, TITLE, new):
        """ put something in the header"""
        if len(new)!=80:
            print new
            raise Exception(
       "Sorry, you are about to corrupt the bruker header, len was %d"%(
                    len(new)))
        h = self.h
        p = h.find(TITLE)
        if p!=-1:
            self.h = h[:p]+new+h[p+80:]
        else:
            raise Exception(TITLE+" not found")

    def tlog(self, msg):
        import time
        new = time.time()
        print "%4.2f %s"%(new-self.last_time, msg),
        self.last_time = new
        

    def convert(self,filein,fileout):
        """ convert a file """
        # Read input file
        data_in = openimage(filein)
        #self.tlog("read")
        # Apply dark and flood
        corrected_image = self.darkflood.correct(data_in.data)
        #self.tlog("corr")
        # make new header
        try:
            sgn = self.omegasign
            om =  float(data_in.header["Omega"])*sgn
            oms=  float(data_in.header["OmegaStep"])*sgn
        except:
            om = 0.
            oms = 0.
        self.putitem("WAVELEN",
                     "WAVELEN:%14.7f%14.7f%14.7f"%(self.wvln,0,0)+\
                         " "*(80-14*3-8))
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
                     "ENDING :%14f%14f%14f%14f"%(0,0, om + oms ,90)+\
                         " "*(80-14*4-8))
        self.putitem("NROWS",
                     "NROWS  :%10d"%(corrected_image.shape[0])+" "*62)
        self.putitem("NCOLS",
                     "NCOLS  :%10d"%(corrected_image.shape[1])+" "*62)
        # beam centre and sample to detector distance and wavelength ??
        # x=hs.find("CENTER")
        #self.tlog("header edit")
        outf = open(fileout,"wb")
        outf.write(self.h)
        #self.tlog("write header")
        outf.write(numpy.ravel(
            numpy.transpose(corrected_image.astype(numpy.uint16))).tostring())
        #self.tlog("write")
        outf.close()
        self.tlog("/s")


if __name__=="__main__":
    import sys
    try:
        from optparse import OptionParser
        parser = OptionParser()
        parser.add_option("-e","--extn",action="store",type="string",
                          dest="extn",
                          default=".edf",
                          help="Filename extension, probably edf, might be cor")
        parser.add_option("-n","--namestem",action="store", 
                          type="string", dest="stem",
    help="Name of the files up the digits part, eg mydata in mydata0000.edf" )
        parser.add_option("-f","--first",action="store", 
                          type="int", dest="first",default=0,
                          help="Number of first file to process, default=0")
        parser.add_option("-l","--last",action="store", type="int", 
                          dest="last",default=0,
                          help="Number of last file to process, default=0")
        parser.add_option("-d","--dark",action="store", type="string", 
                          dest="dark",
                          help="Dark current")
        parser.add_option("-w","--wavelength",action="store", 
                          type="float", dest="wvln",
                          default = 0.5,
                          help="Wavelength")
        parser.add_option("-s","--sign of omega",action="store", 
                          type="float", dest="omegasign",
                          default = 1.0,
                   help="Sign of ID11 omega rotation, +1 is right handed")
        parser.add_option("-F","--Flood",action="store", 
                          type="string", dest="flood",
            default="/data/id11/inhouse/Frelon4M/F4M_Sm_July08.edf",
#           default="/data/id11/inhouse/Frelon2K/Ags_mask0000.edf",
                          help="Flood field")
        parser.add_option("-D","--distance",action="store",type="float", 
                          dest="distance",
                          default=5.0,help="Sample to detector distance")
        parser.add_option("-b","--border",action="store",type="int", 
                          dest="border",
                          default=None,
                          help="Border of image to zero out")

        parser.add_option("-p","--powfac",action="store",type="float", 
                          dest="powfac",
                          default=1.0,
                          help="power to raise image to (1.0 or 0.96)")

        parser.add_option("-t","--template",action="store", type="string", 
                          dest="template",
             default = "/data/id11/inhouse/Frelon2K/brukertemplate.0000")


        parser.add_option("-o","--overflow",action="store",type="float", 
                          dest="overflow",
                          default=65534,
                          help="Overflow value, greater than this set to 65534")
        
        parser.add_option("-j","--jthreads",action="store",type="int", 
                          dest="jthreads",
                          default=1,
                          help="Number of threads to use for processing")
        

        parser.add_option("-R","--detRend",action="store",type="int", 
                          dest="detrend",
                          default=None,
                          help="Number of pixels in background filter")
        

        options, args = parser.parse_args()





        from ImageD11.ImageD11_thread import ImageD11_thread
        class threaded_converter(ImageD11_thread):
            def __init__(self, files, c, name="converter"):
                self.files = files
                self.c = c
                ImageD11_thread.__init__(self,myname=name)
            def ImageD11_run(self):
                for fin, fout in self.files:
                    c.convert(fin, fout)
                    print fout
    
        # Make jthreads converters
        converters = [ edf2bruker(options.dark , options.flood ,
                                  options.template,
                                  distance=options.distance,
                                  border=options.border,
                                  wvln=options.wvln,
                                  omegasign=options.omegasign,
                                  powfac = options.powfac,
                                  overflow = options.overflow,
                                  detrend = options.detrend,
                                  )
                       for j in range(options.jthreads)]

        tmp_in = file_series.numbered_file_series(
            options.stem,
            options.first,
            options.last,
            options.extn)

        import os
        files_in = []
        for f in tmp_in:
            if os.path.exists(f):
                files_in.append(f)
            else:
                print "Missing image",f
           
        files_out = file_series.numbered_file_series(
            options.stem+"_bruker_0.",
            options.first,
            options.last,
            "") # extn

        allfiles = zip(files_in, files_out)

        # Divide files over threads
        fl = [ allfiles[j::options.jthreads] for j in range(options.jthreads) ]

        # Create threads
        my_threads = [ threaded_converter( f, c) for c, f in zip(converters, fl) ]

        # Off we go!
        for t in my_threads:
            t.start()

    except:
        parser.print_help()

        raise

