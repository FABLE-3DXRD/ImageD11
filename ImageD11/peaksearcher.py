#!/usr/bin/python

# Separate this from the command line script.

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
global stop_now
stop_now = False

from math import sqrt
import sys , glob , os.path

from ImageD11 import blobcorrector
from ImageD11.labelimage import labelimage
import numpy

# Generic file format opener from fabio
from fabio.openimage import openimage
from fabio import file_series

# Global variables
OMEGA = 0
OMEGASTEP = 1.0
OMEGAOVERRIDE = False

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

def read_and_correct(filename,
                     dark = None,
                     flood = None,
                     darkoffset = 0):
    """
    Reads a file doing the dark and flood corrections
    """
    try:
        data_object = openimage(filename)
    except KeyboardInterrupt:
        raise
    except:
        sys.stdout.write(filename+" not found\n")
        return 0
    picture = data_object.data.astype(numpy.float32)
    if dark != None:
        # This is meant to be quicker
        picture = numpy.subtract( picture , dark, picture )
        data_object.data = picture
    if flood != None:
        picture = numpy.divide( picture, flood, picture )
        data_object.data = picture
    if darkoffset != None:
        picture = numpy.add( picture , darkoffset, picture )
        data_object.data = picture
        return data_object





def peaksearch( filename , data_object , 
                corrector , 
                thresholds , 
                outputfile , 
                labims ):
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

   
    
    picture = data_object.data.astype(numpy.float32)

    
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
        global OMEGA, OMEGASTEP, OMEGAOVERRIDE
        if not data_object.header.has_key("Omega") or OMEGAOVERRIDE:
            # Might have imagenumber or something??
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


def peaksearch_driver(options, args):
    """
    To be called with options from command line
    """
    ################## debugging still
    for a in args:
        print "arg:"+str(a)
    for o in options.__dict__.keys():
        print "option:",str(o),str(getattr(options,o))
    ###################


    if options.thresholds is None:
        raise ValueError("No thresholds supplied [-t 1234]")

    if len(args) == 0  and options.stem is None:
        raise ValueError("No files to process")

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
                digits = options.ndigits,
                padding = options.padding )
    # Output files:

    outfile =     options.outfile
    if outfile [-4:] != ".spt":
        outfile = outfile + ".spt"
        print "Your output file must end with .spt, changing to ",outfile

    # Omega overrides

    global OMEGA, OMEGASTEP, OMEGAOVERRIDE
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
        floodavg = numpy.mean(middle)
        print "Using flood",options.flood,"average value",floodavg
        if floodavg < 0.7 or floodavg > 1.3:
            print "Your flood image does not seem to be normalised!!!"
         
    else:
        floodimage=None
    start = time.time()
    print "File being treated in -> out, elapsed time"
    # If we want to do read-ahead threading we fill up a Queue object with data 
    # objects
    # THERE MUST BE ONLY ONE peaksearching thread for 3D merging to work
    # there could be several read_and_correct threads, but they'll have to get the order right,
    # for now only one
    if options.oneThread:
        for filein in file_series_object:
            data_object = read_and_correct( filein, darkimage, floodimage,
                                            options.darkoffset)
            peaksearch( filein, data_object , corrfunc , 
                         thresholds_list , outfile , li_objs )
        for t in thresholds_list:
            li_objs[t].finalise()
    else:
        print "Going to use threaded version!?!"
        try:
            import Queue, threading
            class read_all(threading.Thread):
                def __init__(self, q, file_series_obj, dark , flood, darkoffset):
                    self.q = q 
                    self.file_series_obj = file_series_obj
                    self.dark = dark
                    self.flood = flood
                    self.darkoffset = darkoffset
                    threading.Thread.__init__(self)
                def run(self):
                    global stop_now
                    try:
                        for filein in self.file_series_obj:
                            if stop_now:
                                break
                            data_object = read_and_correct(  filein, self.dark, self.flood,
                                                             self.darkoffset)
                            self.q.put((filein, data_object) , block = True)
                        self.q.put( (None, None) , block = True)
                    except KeyboardInterrupt:
                        print "Finishing from read_all"
                        sys.stdout.flush()
                        self.q.put( (None, None) , block = True)
                        stop_now = True
     
            class peaksearch_all(threading.Thread):
                def __init__(self, q, corrfunc, thresholds_list, outfile, li_objs ):
                    self.q = q
                    self.corrfunc = corrfunc
                    self.thresholds_list = thresholds_list
                    self.outfile = outfile
                    self.li_objs = li_objs
                    threading.Thread.__init__(self)

                def run(self):
                    global stop_now
                    try:
                        while 1:
                            if stop_now:
                                break
                            filein, data_object = self.q.get(block = True)
                            if data_object is None:
                                break
                            peaksearch( filein, data_object , self.corrfunc , 
                                        self.thresholds_list , self.outfile , self.li_objs )    
                        for t in self.thresholds_list:
                            self.li_objs[t].finalise()
                    except KeyboardInterrupt:
                        print "Finishing from peaksearcher"
                        sys.stdout.flush()
                        stop_now = True



            q = Queue.Queue(5) # maximum size
            reader = read_all(q, file_series_object, darkimage , floodimage, options.darkoffset)
            searcher = peaksearch_all(q, corrfunc, thresholds_list, outfile, li_objs )
            reader.start()
            searcher.start()
            reader.join()
            searcher.join()

        except ImportError:
            print "Probably no threading module present"
            raise
    


def get_options(parser):
        """ Add our options to a parser object """
        parser.add_option("-n", "--namestem", action="store",
            dest="stem", type="string", default="data",
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
        parser.add_option("--singleThread", action="store_true",
                          dest="oneThread", default=False, 
                          help="Do single threaded processing")
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
                          help="Is the image number to padded Y|N, e.g. "\
                    "should 1 be 0001 or just 1 in image name, default=Y")
        return parser

if __name__=="__main__":
    raise Exception("Please use the driver script peaksearch.py")

