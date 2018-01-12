#!/usr/bin/env python

from __future__ import print_function
## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py



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
from pyFAI import detectors, _distortion



class darkflood:
    """ apply dark and flood corrections """
    def __init__(self,
                 darkfile = None,
                 darkoffset = None,
                 floodfile = None,
                 floodmultiplier = None,
                 splinefile = None,
                 border = None,
                 powfac = 1.0,
                 overflow= 65534,
                 detrend = None,
                 monitorcol = None,
                 monitorval = None,
                 maskfilename = None
                 ):
        self.overflow=overflow
        self.darkfile = darkfile
        self.darkoffset = darkoffset
        self.floodfile = floodfile
        self.splinefile = splinefile
        self.distnorm = None
        self.border = border
        self.floodmultiplier = None
        #
        self.darkimage = None
        self.floodimage = None
        self.flmult = None
        self.powfac = powfac
        self.detrend = detrend
        self.monitorcol = monitorcol
        self.monitorval = monitorval
        if maskfilename is not None:
            self.mask = 1-openimage( maskfilename ).data
        else:
            self.mask = None


    def readdark(self,darkfile):
        """ read the dark"""
        try:
            self.darkdata = openimage(darkfile)
            self.darkfile = darkfile
            self.darkimage = self.darkdata.data.astype(numpy.float32)
            if self.powfac != 1.0:
                print("apply 0.96 to dark")
                self.darkimage = numpy.power(self.darkimage, 0.96)
        except:
            print("No dark file")
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
#                centre = self.floodimage[100:-100,100:-100]
#                npix = centre.shape[0]*centre.shape[1]
#                self.floodmultiplier = numpy.sum(numpy.ravel(
#                        centre).astype(numpy.float32))/npix
                self.flmult = 1 / (self.floodimage)
                print("using flood from",floodfile)
        except:
            print("No flood file")
            self.flooddata = None
            self.floodimage= None
            self.floodfile = None
            self.floodmultiplier = None

    def readspline(self, splinefile):
        """read the spline """
        self.det = detectors.FReLoN(splinefile)
        self.dis = _distortion.Distortion(self.det)
        self.dis.calc_LUT_size()
        self.dis.calc_LUT()
        # remove the pixel size normalisation
        im = numpy.ones((2048,2048),numpy.float32)
        im2 = self.dis.correct(im)
        im2 = numpy.where(im2<0.1,1,im2)
        self.distnorm = 1/im2


    def correct(self, dataobj, detrend = None):
        """ correct the data """
        tin = dataobj.data.dtype
        # Start by copying
        cor = dataobj.data.astype(numpy.float32).copy()
        self.report = ""
        msk = numpy.where( cor > self.overflow,
                           65534,
                           0 )
        if self.powfac != 1.0:
            self.report += "powfac %f;"%(self.powfac)
            numpy.power(cor, 0.96, cor)
        if self.darkimage is not None:
            numpy.subtract(cor, self.darkimage, cor)
            self.report += "dark;"
        if self.detrend is not None:
            assert cor.dtype == numpy.float32, cor.dtype
            cor = self.do_detrend( cor )
            self.report += "detrend;"
        if self.flmult is not None:
            numpy.multiply(cor , self.flmult, cor)
            self.report +="flood;"
        if self.monitorcol is not None:
            scal = self.monitorval / float( dataobj.header[self.monitorcol] )
            numpy.multiply( cor, scal, cor )
            self.report += '%s %.4f;'%(self.monitorcol, scal )
        if self.border is not None:
            # set the edges to zero
            b=self.border
            cor[:b,:]=0
            cor[:,:b]=0
            cor[-b:,:]=0
            cor[:,-b:]=0
            self.report += "border b(%d)=0;"%(self.border)
        if self.splinefile is not None:
            cor = self.dis.correct( cor ) 
            numpy.multiply(cor, self.distnorm, cor)
            self.report += "splinen;"
        if self.darkoffset is not None:
            # print "applying offset of",self.darkoffset,
            numpy.add(cor, self.darkoffset, cor)
            self.report += "+darkoffset %.2f;"%(self.darkoffset)
        if self.mask is not None:
            numpy.multiply(cor, self.mask, cor)
            self.report += "fit2d style mask"
        # Should we bother with this - yes please - noisy pixels overflow
        cor =  numpy.where(cor>0.1, cor, 0.) # truncate zero
        self.report += ">0;"
        cor =  numpy.where(msk != 0 , msk, cor) # mask overflows
        self.report += "msk>%.2f"%(self.overflow)
        ret = cor.astype(tin)
        return ret



    def do_detrend( self, ar):
        if self.detrend is None:
            return ar
        # print "detrending",
        s = ar.copy()
        np = ar.shape[1]/2
        s[:,:np].sort()
        s[:,np:].sort()
        n = self.detrend
        nb = 5 # bad pixels (eg negative outliers)
        o1 = s[:,nb:(n+nb)].sum(axis=1)/n
        o2 = s[:,(np+nb):(np+n+nb)].sum(axis=1)/n
        s1 = ar.copy()
        s1[:,:np] = (ar[:,:np].T - o1 + o1.mean() ).T
        s1[:,np:] = (ar[:,np:].T - o2 + o2.mean() ).T
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
                 detrend = None,
                 monitorval = None,
                 monitorcol = None,
                 maskfilename=None,
                 splinefile = None,
                 omega_zero=0,
                 chi_zero=0):
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
                                   splinefile = splinefile,
                                   detrend = detrend,
                                   monitorval = monitorval,
                                   monitorcol = monitorcol,
                                   maskfilename=maskfilename)
        self.darkflood.readdark(dark)
        self.darkflood.readflood(flood)
        self.darkflood.readspline(splinefile)
        self.templatefile = template
        im = brukerimage()
        im.readheader(template)
        self.header = im.header
        self.h = im.__headerstring__
        self.last_time = time.time()
        self.omega_zero = omega_zero
        self.chi_zero = chi_zero

    def putitem(self, TITLE, new):
        """ put something in the header"""
        if len(new)!=80:
            print(new)
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
        new = time.time()
        print("%4.2f %s"%(new-self.last_time, msg), end=' ')
        self.last_time = new


    def convert(self,filein,fileout):
        """ convert a file """
        # Read input file
        data_in = openimage(filein)
        #self.tlog("read")
        # Apply dark and flood
        corrected_image = self.darkflood.correct(data_in)
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
        # angles are twotheta omega phi chi
        try:
            hphi = float(data_in.header["hphi"])
        except:
            hphi = 0.0
        try:
            hchi = float(data_in.header["hchi"])
        except:
            hchi = -45.0
        if corrected_image.shape[0] != corrected_image.shape[1]:
            # we expand to square and pad the centre
            sz = max(corrected_image.shape[0], corrected_image.shape[1])
            sh = corrected_image.shape
            im = numpy.zeros((sz,sz),corrected_image.dtype)
            numpy.add( im[0:sh[0],0:sh[1]], corrected_image[0:sh[0],0:sh[1]], im[0:sh[0],0:sh[1]] )
            corrected_image = im
        self.putitem("ANGLES",
                     "ANGLES :%14f%14f%14f%14f"%(0,om-self.omega_zero,hphi,hchi-self.chi_zero) +
                     " "*(80-14*4-8))
        self.putitem("DISTANC",
                     "DISTANC:%14f"%(self.distance)+" "*(80-14-8))
        self.putitem("RANGE",
                     "RANGE  :     %9f"%( abs(oms) ) + " "*58)
        self.putitem("INCREME:",
                     "INCREME:     %9f"%( oms ) + " "*58 )
        self.putitem("START",
                     "START  :%14f"%( om )+" "*(80-14-8))
        self.putitem("ENDING",
                     "ENDING :%14f%14f%14f%14f"%(0,om+oms,hphi,hchi)+" "*(80-14*4-8))
        self.putitem("NROWS",
                     "NROWS  :%10d"%(corrected_image.shape[0])+" "*62)
        self.putitem("NCOLS",
                     "NCOLS  :%10d"%(corrected_image.shape[1])+" "*62)
        self.putitem("AXIS",
                     "AXIS   :%10d"%(2)+" "*62) # omega
        # beam centre and sample to detector distance and wavelength ??
        # x=hs.find("CENTER")
        #self.tlog("header edit")
        outf = open(fileout,"wb")
        outf.write(self.h)
        #self.tlog("write header")
        outim = numpy.flipud(corrected_image)
        outf.write( outim.astype(numpy.uint16).tostring() )
        #self.tlog("write")
        outf.close()
        self.tlog("/s "+fileout + " " + self.darkflood.report + "\n")


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
        parser.add_option("-N","--newnamestem",action="store",
                          type="string", dest="newstem", default=None,
                          help="Name of the files up the digits part for output")
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

            default=None,
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
             default = "/data/id11/3dxrd/inhouse/Frelon2K/brukertemplate.0000")


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


        parser.add_option("-O","--darkoffset",action="store",type="float",
                          dest="darkoffset",
                          default=100.0,
                          help="Offset platform to add to image")

        parser.add_option("-r","--run",action="store",type="int",
                          dest="run",
                          default=0,
                          help="Run number for bruker")

        parser.add_option("--monitorcol", action="store", type="string",
                           dest="monitorcol",
                           default = None,
                           help="Header value for incident beam intensity")


        parser.add_option("--monitorval", action="store", type="float",
                          dest="monitorval",
                          default = None,
                          help="Incident beam intensity value to normalise to")

        parser.add_option("--maskfilename",
                          action="store", type="string",
                          dest="maskfilename",
                          default = None,
                          help="Apply a fit2d style mask to image")

        parser.add_option("--splinefile",
                          action="store", type="string",
                          dest="splinefile",
                          default = None,
                          help="Apply a spatial distortion" )

        parser.add_option("--omega_zero",
                          action="store", type="float",
                          dest="omega_zero",
                          default = 0.0,
                          help="Omega zero from smart")

        parser.add_option("--chi_zero",
                          action="store", type="float",
                          dest="chi_zero",
                          default = 0.0,
                          help="chi zero from smart")

        parser.add_option("--ndigits",
                          action="store", type="int",
                          dest="ndigits",
                          default = 4,
                          help="Number of digits in output name")





        options, args = parser.parse_args()

        if options.stem is None:
            print("Fill -n option or use --help")
            sys.exit()


        from ImageD11.ImageD11_thread import ImageD11_thread
        import ImageD11.ImageD11_thread
        class threaded_converter(ImageD11_thread):
            def __init__(self, files, c, name="converter"):
                self.files = files
                self.c = c
                ImageD11_thread.__init__(self,myname=name)
            def ImageD11_run(self):
                for fin, fout in self.files:
                    self.c.convert(fin, fout)
                    if self.ImageD11_stop_now():
                        break


        # Make jthreads converters
        converters = [ edf2bruker(options.dark , options.flood ,
                                  options.template,
                                  darkoffset=options.darkoffset,
                                  splinefile=options.splinefile,
                                  distance=options.distance,
                                  border=options.border,
                                  wvln=options.wvln,
                                  omegasign=options.omegasign,
                                  powfac = options.powfac,
                                  overflow = options.overflow,
                                  detrend = options.detrend,
                                  monitorval = options.monitorval,
                                  monitorcol = options.monitorcol ,
                                  maskfilename = options.maskfilename,
                                  omega_zero = options.omega_zero,
                                  chi_zero = options.chi_zero
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
                print("Missing image",f)


        ostem = options.stem
        if options.newstem is not None:
            ostem = options.newstem

        files_out = file_series.numbered_file_series(
            ostem+"_%d."%(options.run),
            options.first+1,
            options.last+1,
            "", # extn
            options.ndigits )

        allfiles = list(zip(files_in, files_out))

        # Divide files over threads
        fl = [ allfiles[j::options.jthreads] for j in range(options.jthreads) ]

        # Create threads
        my_threads = [ threaded_converter( f, c) for
                       c, f in zip(converters, fl) ]

        # Off we go!
        for t in my_threads:
            t.start()

        nalive = 1
        while nalive > 0:
            try:
                nalive = 0
                for t in my_threads:
                    if t.isAlive():
                        nalive += 1
                time.sleep(1)
            except KeyboardInterrupt:
                ImageD11.ImageD11_thread.stop_now = True
                print("Got your control-c - terminating threads")
                for t in my_threads:
                    if t.isAlive():
                        t.join(timeout=10)


    except:
        parser.print_help()

        raise

