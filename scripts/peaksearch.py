#! /bliss/users/blissadm/python/bliss_python/suse82/bin/python


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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  0211-1307  USA



"""
Script for peaksearching images from the command line

Uses the connectedpixels extension for finding blobs above a threshold
and the blobcorrector(+splines) for correcting them for spatial distortion

Defines one function (peaksearch) which might be reused
"""

import time, sys
# For benchmarking
reallystart = time.time()

from math import sqrt
import sys , glob , os.path
from ImageD11 import blobcorrector
from ImageD11 import opendata
from ImageD11 import connectedpixels
import Numeric
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
        print "%.2f/s"%(self.now-self.start)

class labelimage:
    def __init__(self,shape,fileout=sys.stdout,spatial=blobcorrector.perfect):
        self.shape=shape
        self.bl = Numeric.zeros(shape,Numeric.Int)
        self.lastbl = Numeric.zeros(shape,Numeric.Int)
        self.np = 0
        self.res = (0,)*10
        self.lastres = []
        self.lastnp=-1
        self.verbose=0
        self.nprop = 9 # number of accumulating properties
        self.corrector=spatial
        self.closed =None
        self.finalpeaks=[]
        self.outfile=open(fileout,"w")
        self.outfile.write("# xc yc omega npixels avg_intensity x_raw y_raw sigx sigy covxy\n")

    def peaksearch(self,data,threshold,dark=None,flood=None):
        """
        # Call the c extension to do the peaksearch, on entry:
        #
        # picture = 2D Numeric UInt16 array (of your data)
        # blobim  = 2D Numeric Int array of same shape as data
        # threshold = float - pixels above this number are put into objects
        # verbose = 0/1 - whether or not to print anything
        #
        # One return - np is the number of peaks found (number of
        # connnected objects) picture should be unchanged
        # blobim should be overwritten with integers giving the peak number for
        # that pixel
        #   or zero if the pixel was below the threshold
        #

        """
        self.bl = Numeric.zeros(self.shape,Numeric.Int)
        # print "labelimage.peaksearch",data[0,0],threshold
        optarg = {"verbose":self.verbose}
        if dark is not None: optarg["dark"]=dark
        if flood is not None: optarg["flood"]=flood
        self.np = connectedpixels.connectedpixels(data, self.bl, threshold,**optarg)
        return self.np
        

    def properties(self,data,omega=0.,dark=None,flood=None):
        """
        # Now get the properties of the identified objects, on entry:
        # picture = data, as before
        # blobim  = peak assignments
        # np      = the number of peaks to find
        # verbose = 0/1 - whether to print or not
        #
        #
        # This returns a series of 1D Numeric arrays:
        #  npi = Number of pixels in objects ...  for each blob pixel
        #          sum 1
        #  isum =  sum of intensity in blobs
        #          sum data[i][j]
        #  sumsq = sum of intensity squared
        #          sum data[i][j] * data[i][j]
        #  com0  = sum of intensity multiplied by fast index:
        #          sum i*data[i][j]
        #  com1  = sum of intensity multiplied by slow index:
        #          sum j*data[i][j]
        #  com00 = sum of intensity multiplied by fast index squared:
        #          sum i*i*data[i][j]
        #  com01 = sum of intensity multiplied by fast index and slow index:
        #          sum i*j*data[i][j]
        #  com00 = sum of intensity multiplied by slow index squared:
        #          sum j*j*data[i][j]
        #  com2  = sum of intensity multiplied by slowest index (normally rotation angle)
        #        = sum omega*isum
        #
        
        """
        if self.np>0:
            optarg = {"verbose":self.verbose}
            if dark is not None: optarg["dark"]=dark
            if flood is not None: optarg["flood"]=flood
            res = connectedpixels.blobproperties(data, self.bl, self.np,**optarg)
            try:
                self.res = res + (Numeric.ones(res[1].shape[0],Numeric.Float)*omega * res[1],)
            except:
#                print res,self.np
                raise
        else:
            self.res = (0,)*10
        return self.res

    def mergelast(self):
        """
        Merge the last two images searched
        """
        if self.lastnp < 0:
            # no previous to merge with, just save the current for next time
            self.lastbl = self.bl
            self.lastnp = self.np
            self.lastres = self.res
            return
        # print "Calling bloboverlaps",self.verbose
        ds = connectedpixels.bloboverlaps(self.bl, self.np, self.lastbl, self.lastnp, self.verbose)
        # ds is a disjoint set. 
        # print ds
        # already filled out self.res, self.lastres
        #
        # Move counts in self.lastres forwards to self.res
        # Once all checked move non-zeros in self.lastres to self.closed
        three_d_debug=False
        if three_d_debug:
            print len(ds),self.np,self.lastnp,ds
            print "np",self.np,"len(self.res[0])",len(self.res[0])
            print "lastnp",self.lastnp,"len(self.lastres[0])",len(self.lastres[0])
        for i in range(self.np+1,len(ds)):
#            print i
#            if ds[i]==i and i<self.np:
#                # this peak is on res and not overlapped, it will, eventually, move to lastres
#                continue
            if ds[i]!=i:
                if three_d_debug:
                    print "linking ds[i],i",ds[i],i,
                assert ds[i] < self.np, "Link on the same image"
                # propose to link these two spots...
                # print "l",i,ds[i],"  ",
                # bring the stuff forwards from last
                if three_d_debug:
                    peak1 = []
                    peak2 = []
                for j in range(self.nprop): # loop over properties
                    if three_d_debug: # debug
                        peak1.append(self.res[j][ds[i]])
                        peak2.append(self.lastres[j][i-self.np])
                    try:
                        self.res[j][ds[i]]=self.res[j][ds[i]]+self.lastres[j][i-self.np]
                    except:
                        print j,i,ds[i]
                        print ds
                        raise Exception("You should not be seeing this error!!")
                if three_d_debug:
                    print "Linking ds[i],[i-self.np]",i,ds[i],i-self.np
                    print self.to_string(peak1),
                    print self.to_string(peak2),                    
            else: # this peak is no longer active so write it out
                peak = []
                for j in range(self.nprop): # loop over properties
                    # self.res[property][peak]
                    try:
                        peak.append(self.lastres[j][i-self.np])
                    except:
                        print "indices",j,i,i-self.np,self.lastres
                        raise                    
                self.outputpeak(peak)
        # Now get the overlaps on the current frame
        #
        nbad =0
        for i in range(1,self.np):
                if ds[i]!=i :
                    # propose to link these two spots on the current active frame
                    nbad += 1
                    continue
                    # print "linking on current not implemented",i,ds[i]
                    peak1=[]; peak2=[]
                    
                    # write the correct number in the blob image
                    # ... expensive - could be done in C??
                    connectedpixels.update_blobs(self.bl, ds)
                    self.bl = Numeric.where(self.bl == i, ds[i], self.bl)
                    for j in range(self.nprop): # loop over properties
                        try:
                            self.res[j][ds[i]]=self.res[j][ds[i]]+self.res[j][i]
                            self.res[j][i]=0
                        except:
                            print j,i,ds[i]
                            print ds
                            raise Exception("You should not be seeing this error!!")
        if nbad != 0:
            print "unlinked pks",nbad,
        temp = self.lastbl
        self.lastbl = self.bl
        self.bl=temp
        self.lastnp = self.np
        self.lastres = self.res

        
    def outputpeak(self,peak):
        self.outfile.write(self.to_string(peak))

    def to_string(self,peak):
        """
        Convert the accumulators to more meaningful things
        self.outfile.write("# xc yc omega npixels avg_intensity x_raw y_raw sigx sigy covxy\n")
        """
        npi, isum, sumsq, com0, com1, com00, com01, com11, omisum = peak
        if npi < 1:
            return ""
        else:
            n    = npi
            # Average intensity
            avg  = isum/n
            if n==1:
                si   = sqrt(sumsq - n*avg*avg)
            else:
                si   = sqrt((sumsq - n*avg*avg)/(n-1.))
            # Standard dev on intensity
            c0   = com0/isum
            # Centre of mass in index 0
            c1   = com1/isum
            # Centre of mass in index 1
            # Covariances - try except to allow for zeros
            try:
                c00 = sqrt((com00/isum - c0*c0))
            except:
                c00 = 0. # this means a line of pixels and rounding errors
            try:
                c11 = sqrt((com11/isum - c1*c1))
            except:
                c11 = 0.
            try:
                c01 = (com01/isum - c0*c1)/c00/c11
            except:
                c01 = 0.
            # Spatial corrections,
            # c0c and c1c are the distortion corrected centre of mass :
            c0c, c1c = self.corrector.correct(c0, c1)
            om = omisum/isum
            if False:# abs(c0c-1229)<2 and abs(c1c-1486)<2:
                print "c0c,c1c",c0c,c1c
                print peak
                print "%f %f %f %d %f %f %f %f %f %f\n"%(c0c, c1c, om, n, avg, c0, c1, c00, c11, c01)
                raise Exception("who are you??")
            return "%f %f %f %d %f %f %f %f %f %f\n"%(c0c, c1c, om, n, avg, c0, c1, c00, c11, c01)


    def finalize(self):
        """
        Assume we finished the last image
        """
        # self.outfile.write("#final frame\n")
        for i in range(len(self.lastres[0])): # loop over peaks
            peak = []
            for j in range(len(self.lastres)): # loop over properties
                # self.res[property][peak]
                try:
                    peak.append(self.lastres[j][i])
                except:
                    print "indices",j,i,i-self.np,self.lastres
                    raise
            self.outputpeak(peak)
        self.outfile.close()

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
        print "%.2f/s"%(self.now-self.start)
        
def peaksearch(filename, outputfile, corrector, labims , thresholds,
               dark=None, flood=None):
    """
    Searches for peaks, arguments are:
    filename   : image filename to search for peaks in
    outputfile : ascii file to receive results
    corrector  : ImageD11.blobcorrector object to correct peak positions
    labim      : Numeric.Int array of same shape as data image to hold blob
                 assignments (dict[thresholds] of pairs)
    thresholds : List of threshold values to use (use more than one if you are
                 not sure what to use), any iterable will do, not just list
    """
    t = timer()
    f = open(outputfile,"aq") # Open the output file for appending
    # Assumes an edf file for now - TODO - should be independent
    data_object = opendata.opendata(filename)
    picture = data_object.data
    # picture is (hopefully) a 2D Numeric array of type UInt16
    # t.tick("read")
    #if dark != None:
        # Slows down timing but avoid overflows
        # print "Subtracting dark,",picture[0,0]
        # print picture.shape,picture.typecode(),dark.shape,dark.typecode()
        #picture.savespace(0) # avoid overflows - very slow and silly
        #picture = picture.astype(Numeric.Int)-dark  
        # print "Subtracted dark,",picture[0,0]
    # t.tick("dark")
    #if flood != None:
        #picture = picture/flood
    t.tick(filename+" io/cor") # Progress indicator
    #print "datatype for searching", picture.typecode()    #
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
        npks = 0
        #
        # Do the peaksearch
        np = labelim.peaksearch(picture,threshold,dark=dark,flood=flood)
        #
        try:
            omega = float(data_object.header["Omega"])
        except: # Might have imagenumber or something??
            omega = 0.
        props = labelim.properties(picture,omega)
        npi, isum, sumsq, com0, com1, com00, com01, com11, omisum = props
        # omisum is omega * isum
        # Now write results out for this threshold level
        f.write("\n#Threshold level %f\n" % (threshold))
        f.write( "# Number_of_pixels Average_counts    x   y     xc   yc      sig_x sig_y cov_xy\n")
        outstrfmt = "%d  %f    %f %f    %f %f    %f %f %f\n"
        for  i in range(len(npi)): # Loop over peaks
            if npi[i]>=1:    # Throw out one pixel peaks (div zero)
                npks = npks+1
                n    = npi[i]
                # Average intensity
                avg  = isum[i]/n
                if n>2:
                    si   = sqrt((sumsq[i] - n*avg*avg)/(n-1.))
                else:
                    si   = sqrt((sumsq[i] - n*avg*avg)/n)
                # Standard dev on intensity
                c0   = com0[i]/isum[i]
                # Centre of mass in index 0
                c1   = com1[i]/isum[i]
                # Centre of mass in index 1
                # Covariances - try except to allow for zeros
                try:
                    c00 = sqrt((com00[i]/isum[i] - c0*c0))
                except:
                    c00 = 0. # this means a line of pixels and rounding errors
                try:
                    c11 = sqrt((com11[i]/isum[i] - c1*c1))
                except:
                    c11 = 0.
                try:
                    c01 = (com01[i]/isum[i] - c0*c1)/c00/c11
                except:
                    c01 = 0.
                # Spatial corrections,
                # c0c and c1c are the distortion corrected centre of mass :
                try:
                    c0c, c1c = corrector.correct(c0, c1)
                except:
                    c0c, c1c = c0, c1
                f.write(outstrfmt % (n, avg, c0, c1, c0c, c1c,  c00, c11, c01))
        labelim.mergelast() # called here!!!
        print "T=%-5d n=%-5d;" % (int(threshold),npks),
        # Close the output file
    f.close()
    # Finish progress indicator for this file
    t.tock()
    sys.stdout.flush()
    return npks # Number of peaks found


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
            dest="format",default="edf", type="string",
            help="Image File format, eg edf or bruker" )
        parser.add_option("-f", "--first", action="store",
            dest="first", default=0, type="int",
            help="Number of first file to process, default=0")
        parser.add_option("-l", "--last", action="store",
            dest="last", type="int",
            help="Number of last file to process")
        parser.add_option("-o", "--outfile", action="store",
            dest="outfile",default="peaks.out", type="string",
            help="Output filename, default=peaks.out")
        parser.add_option("-d", "--darkfile", action="store",
            dest="dark", default=None,  type="string",
            help="Dark current filename, to be subtracted, default=None")
        parser.add_option("-D", "--darkfileoffset", action="store",
            dest="darkoffset", default=100, type="int",
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
        options , args = parser.parse_args()
        stem =        options.stem
        outfile =     options.outfile
        first =       options.first
        last =        options.last
        if options.perfect=="N":
            print "Using spatial from",options.spline
            corrfunc = blobcorrector.correctorclass(options.spline)
        else:
            print "Avoiding spatial correction"
            corrfunc = blobcorrector.perfect()

        # List comprehension - convert remaining args to floats - must be unique list        
        thresholds_list = [float(t) for t in options.thresholds]
        import sets
        thresholds_list = list(sets.Set(thresholds_list))
        thresholds_list.sort()
        
        # Generate list of files to process
        print "Input format is",options.format
        if options.format == "edf":
            files = ["%s%04d%s" % (stem,i,".edf") for i in range(first,last+1)]
            corrfunc.orientation = "edf"
        if options.format == "bruker":
            files = ["%s.%04d" % (stem,i) for i in range(first,last+1)]
            corrfunc.orientation = "bruker"
            if options.spline.find("spatial2k.spline")>0:
                if raw_input("Are you sure about the spline??")!="y":
                    sys.exit()
        # Make a blobimage the same size as the first image to process
        if len(files)==0:
            raise "No files found for stem %s" % (stem)
        li_objs={} # label image objects, dict of
        s = opendata.opendata(files[0]).data.shape # data array shape
        # Create label images
        for t in thresholds_list:
            mergefile="%s_merge_t%d"%(options.outfile,t)
            li_objs[t]=labelimage(shape=s,fileout=mergefile,spatial=corrfunc) 
        # Not sure why that was there (I think if glob was used)
        # files.sort()
        if options.dark!=None:
            print "Using dark (background)",options.dark,"with added offset",options.darkoffset
            darkimage= (opendata.opendata(options.dark).data.astype(Numeric.Int)-\
                       options.darkoffset).astype(Numeric.UInt16)
        else:
            darkimage=None
        if options.flood!=None:
            floodimage=opendata.opendata(options.flood).data.astype(Numeric.Float32)
            # Check flood normalisation
            m1 = floodimage.shape[0]/6
            m2 = floodimage.shape[1]/6
            middle = Numeric.ravel(floodimage[m1:-m1,m2:-m2])
            floodavg =  Numeric.sum(middle)/middle.shape[0]
            print "Using flood",options.flood,"average value",floodavg
            if floodavg < 0.7 or floodavg > 1.3:
                print "Your flood image does not seem to be normalised!!!"
        else:
            floodimage=None
        start = time.time()
        print "File being treated in -> out, elapsed time"
        for filein in files:
            peaksearch(filein, outfile, corrfunc, li_objs, \
                       thresholds_list,dark=darkimage,flood=floodimage)
        for t in thresholds_list:
            li_objs[t].finalize()
    except:
#      print "Usage: %s filename  outputfile first last spline threshold [threshold...]" % (sys.argv[0])
        print
        # Raise the exception to find out in more detail where we went wrong
        raise
        sys.exit()

end = time.time()
t = end-reallystart
print "Total time = %f /s, per image = %f /s" % (t,t/len(files))

