#!/usr/bin/python



## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py

#! /bliss/users/blissadm/python/bliss_python/suse82/bin/python


# ImageD11_v1.0 Software for beamline ID11
# Copyright (C) 2005-2007  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  0211-1307  USA



"""
Script for peaksearching images from the command line

Uses the connectedpixels extension for finding blobs above a threshold
and the blobcorrector(+splines) for correcting them for spatial distortion

Defines one function (peaksearch) which might be reused
"""

import time

# For benchmarking
reallystart = time.time()

from math import sqrt
import sys , glob , os.path

from ImageD11 import blobcorrector
from ImageD11.labelimage import labelimage
import numpy.oldnumeric as Numeric


# Generic file format opener from fabio
from fabio.openimage import openimage
from fabio import file_series



class timer:
    def __init__(self):
        self.start = time.time()
        self.now = self.start
    def tick(self,msg=""):
        now = time.time()
        print msg,"%.2f/s"%(now-self.now),
        self.now = now
    def tock(self,msg=""):
        self.tick(msg)
        print "%.2f/s"% (self.now-self.start)


OMEGA = 0
OMEGASTEP = 1.0
OMEGAOVERRIDE = False

def peaksearch( filename , 
                corrector , 
                thresholds , 
                outputfile , 
                labims, 
                dark = None,
                flood = None,
                darkoffset = 10):
    """
    file_series : fabio file series object, supports iteration
    
    corrector : spatial and dark/flood , linearity etc
    
    thresholds : [ float[1], float[2] etc]
    
    outputfile : the 2d output file, name+".spt"
               : the 3d output files name+"_threshold.flt"
    """
    t = timer()
    f = open(outputfile,"aq") # Open the output file for appending
    # Assumes an edf file for now - TODO - should be independent
    try:
        data_object = openimage(filename)
    except KeyboardInterrupt:
        raise
    except:
        sys.stdout.write(filename+" not found\n")
        return 0

    picture = data_object.data.astype(Numeric.Float32)
    # print picture[512,512]
    # picture is a 2D array
    #t.tick("read")
    if dark != None:
        picture = picture - dark
        #print dark[512,512]
        #print picture[512,512]
    #t.tick("dark")
    if flood != None:
        #print picture[512,512]
        picture = picture/flood
        #print picture[512,512]
    if darkoffset != None:
        picture = picture + darkoffset
    #fd=open("debug.bin","wb")
    #fd.write(picture.tostring())
    #fd.close()
    t.tick(filename+" io/cor") # Progress indicator
    # Transfer header information to output file
    # Also information on spatial correction applied
    f.write("\n\n# File %s\n" % (filename))
    f.write("# Processed on %s\n" % (time.asctime()))
    try:
        f.write("# Spatial correction from %s\n" % (corrector.splinefile))
        f.write("# SPLINE X-PIXEL-SIZE %s\n" % (str(corrector.xsize)))
        f.write("# SPLINE Y-PIXEL-SIZE %s\n" % (str(corrector.ysize)))
    except:
        pass
    for item in data_object.header.keys():
        if item == "headerstring": # skip
            continue
        try:
            f.write("# %s = %s\n" % (item,
                    str(data_object.header[item]).replace("\n"," ")))
        except KeyError:
            pass
    #
    # Now peaksearch at each threshold level
    for threshold in thresholds:
        labelim = labims[threshold]
        if labelim.shape != picture.shape:
            raise "Incompatible blobimage buffer for file %s" %(filename)
        #
        if not data_object.header.has_key("Omega") or OMEGAOVERRIDE:
            # Might have imagenumber or something??
            global OMEGA
            ome = OMEGA
            OMEGA += OMEGASTEP
            # print "Overriding the OMEGA"
        else: # Might have imagenumber or something??
            ome = float(data_object.header["Omega"])
            # print "Reading from header"
        #
        # Do the peaksearch
        f.write("# Omega = %f\n"%(ome))
        labelim.peaksearch(picture, threshold, ome)
        f.write("# Threshold = %f\n"%(threshold))
        f.write("# npks = %d\n"%(labelim.npk))
        #
        if labelim.npk > 0:
            labelim.output2dpeaks(f)
        labelim.mergelast() 
        print "T=%-5d n=%-5d;" % (int(threshold),labelim.npk),
        # Close the output file
    f.close()
    # Finish progress indicator for this file
    t.tock()
    sys.stdout.flush()
    return None 


if __name__=="__main__":
    # If we are running from a command line:
    try:
        from optparse import OptionParser
        parser = OptionParser()
        parser.add_option("-n", "--namestem", action="store",
            dest="stem", type="string",
            help="Name of the files up the digits part  "+\
                 "eg mydata in mydata0000.edf" )
        parser.add_option("-F", "--format", action="store",
            dest="format",default=".edf", type="string",
            help="Image File format, eg edf or bruker" )
        parser.add_option("-f", "--first", action="store",
            dest="first", default=0, type="int",
            help="Number of first file to process, default=0")
        parser.add_option("-l", "--last", action="store",
            dest="last", type="int",default =0,
            help="Number of last file to process")
        parser.add_option("-o", "--outfile", action="store",
            dest="outfile",default="peaks.spt", type="string",
            help="Output filename, default=peaks.spt")
        parser.add_option("-d", "--darkfile", action="store",
            dest="dark", default=None,  type="string",
            help="Dark current filename, to be subtracted, default=None")
        parser.add_option("-D", "--darkfileoffset", action="store",
            dest="darkoffset", default=10, type="int",
            help="Constant to subtract from dark to avoid overflows, default=100")
        s="/data/opid11/inhouse/Frelon2K/spatial2k.spline"
        parser.add_option("-s", "--splinefile", action="store",
            dest="spline", default=s, type="string",
            help="Spline file for spatial distortion, default=%s" % (s))
        parser.add_option("-p", "--perfect_images", action="store",
               type="choice", choices=["Y","N"], default="N", dest="perfect",
                          help="Ignore spline Y|N, default=N")
        parser.add_option("-O", "--flood", action="store", type="string",
                          default=None, dest="flood",
                          help="Flood file, default=None")
        parser.add_option("-t", "--threshold", action="append", type="float",
             dest="thresholds", default=None,
             help="Threshold level, you can have several")
        parser.add_option("--OmegaFromHeader", action="store_false",
                          dest="OMEGAOVERRIDE", default=False, 
                          help="Read Omega values from headers [default]")
        parser.add_option("--OmegaOverRide", action="store_true",
                          dest="OMEGAOVERRIDE", default=False, 
                          help="Override Omega values from headers")
        parser.add_option("-S","--step", action="store",
                          dest="OMEGASTEP", default=1.0, type="float",
                          help="Step size in Omega when you have no header info")
        parser.add_option("-T","--start", action="store",
                          dest="OMEGA", default=0.0, type="float",
                          help="Start position in Omega when you have no header info")

        parser.add_option("--ndigits", action="store", type="int",
                dest = "ndigits", default = 4,
                help = "Number of digits in file numbering [4]")
        parser.add_option("-P", "--padding", action="store",
               type="choice", choices=["Y","N"], default="Y", dest="padding",
                          help="Is the image number to padded Y|N, e.g. should no. 1 be 0001 or just 1 in image name, default=Y")
        options , args = parser.parse_args()


        if options.thresholds is None:
            parser.print_help()
            print "\nNo thresholds supplied [-t 1234]\n"
            sys.exit()

        if len(args) == 0  and options.stem is None:
            parser.print_help()
            print "\nNo files to process?\n"        
            sys.exit()

        # What to do about spatial
 
        if options.perfect=="N" and os.path.exists(options.spline):
            print "Using spatial from",options.spline
            corrfunc = blobcorrector.correctorclass(options.spline)
        else:
            print "Avoiding spatial correction"
            corrfunc = blobcorrector.perfect()



        # Get list of filenames to process
        if len(args) > 0 :
            # We assume unlabelled arguments are filenames 
            file_series_object = file_series.file_series(args)
        else:
            if options.format in ['bruker', 'BRUKER', 'Bruker']:
                extn = ""
                corrfunc.orientation = "bruker"
            elif options.format == 'GE':
                extn = ""
            else:
                extn = options.format
                corrfunc.orientation = "edf"
            file_series_object = file_series.numbered_file_series(
                    options.stem,
                    options.first,
                    options.last,
                    extn,
                    options.ndigits,
                    options.padding )

        # Output files:

        outfile =     options.outfile
        if outfile [-4:] != ".spt":
            outfile = outfile + ".spt"
            print "Your output file must end with .spt, changing to ",outfile

        # Omega overrides

        OMEGA = options.OMEGA
        OMEGASTEP = options.OMEGASTEP
        OMEGAOVERRIDE = options.OMEGAOVERRIDE 
        # Make a blobimage the same size as the first image to process


        # List comprehension - convert remaining args to floats - must be unique list        
        thresholds_list = [float(t) for t in options.thresholds]
        import sets
        thresholds_list = list(sets.Set(thresholds_list))
        thresholds_list.sort()
        
        li_objs={} # label image objects, dict of
        s = openimage(file_series_object.current()).data.shape # data array shape
        # Create label images
        for t in thresholds_list:
            # the last 4 chars are guaranteed to be .spt above
            mergefile="%s_t%d.flt"%(outfile[:-4], t)
            li_objs[t]=labelimage(shape = s, fileout = mergefile, spatial = corrfunc) 
        # Not sure why that was there (I think if glob was used)
        # files.sort()
        if options.dark!=None:
            print "Using dark (background)",options.dark,"with added offset",options.darkoffset
            darkimage= openimage(options.dark).data 
        else:
            darkimage=None
        if options.flood!=None:
            floodimage=openimage(options.flood).data
            cen0 = floodimage.shape[0]/6
            cen1 = floodimage.shape[0]/6
            middle = floodimage[cen0:-cen0, cen1:-cen1]
            nmid = middle.shape[0]*middle.shape[1]
            floodavg = Numeric.sum(
                Numeric.ravel(middle).astype(Numeric.Float32))/nmid
            print "Using flood",options.flood,"average value",floodavg
            if floodavg < 0.7 or floodavg > 1.3:
                print "Your flood image does not seem to be normalised!!!"
            
        else:
            floodimage=None
        start = time.time()
        print "File being treated in -> out, elapsed time"
        for filein in file_series_object:
            peaksearch( filein , corrfunc , thresholds_list , outfile , li_objs,
                     dark = darkimage, flood = floodimage , darkoffset = options.darkoffset)

        for t in thresholds_list:
            li_objs[t].finalise()
    # Catch --help exiting
    except SystemExit:
        sys.exit()
    except:
#      print "Usage: %s filename  outputfile first last spline threshold [threshold...]" % (sys.argv[0])
        parser.print_help()
        print "\n\nGot an error:"
        # Raise the exception to find out in more detail where we went wrong
        import traceback
        traceback.print_exc()
        sys.exit()

end = time.time()
t = end-reallystart
print "Total time = %f /s, per image = %f /s" % (t,t/len(file_series_object))

