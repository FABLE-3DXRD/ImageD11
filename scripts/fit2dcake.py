#!/usr/bin/env python

from __future__ import print_function

"""
Fit2d caking python script

Uses a single input file to run fit2d from the command line in batch mode

Can generate parameter files from a .fit2d.def file

example usage:
  1) Run fit2d, getting good calibration parameters and caking parameters
  2) Run this script using the -c option to generate a parameter file
        eg: fit2dcake.py -cmypars.pars
  3) Have a look at the file mypars.pars in a text editor, fix it as you like
  4) Process your data using the parameter file you have generated
        eg: fit2dcake.py
"""
import logging, sys, os

import numpy as np


if sys.platform == 'win32':
    XVFB = "echo "
    FIT2D = """D:\\fit2d_12_077_i686_WXP.exe """
    REMOVE = "del "
else:
    if not os.path.exists("/tmp/.X1-lock"):
        XVFB = "Xvfb :942 -ac &"
    else:
        XVFB = "echo "
    FIT2D = "fit2d_12_081_i686_linux2.4.20"
    REMOVE = "rm -f "


def spr_to_median(filename, outputfilename):
    """
    Takes the median accross rows of an spr file.
    Intention is to do radial transform.
    Write out radially transformed as spr and take median.
    """
    data = []
    for line in open(filename,"r").readlines()[1:]:
        data.append([float(x) for x in  line.split()])
    data = np.array(data)
    # print data[512,:]
    data = np.sort(data, 1)
    # print data[512,:]
    x_axis = []
    for line in open(filename+"_y","r").readlines()[4:]:
        x_axis.append(float(line.split()[0]))
    # print len(x_axis)
    # Final result
    o = open(outputfilename, "w")
    # use central trimean (less noisy than median)
    d = data[:, data.shape[1]/4 : data.shape[1]*3/4]
    med = d.sum(axis=1) / d.shape[1]
    nsum = np.sum

    # VAR = E(X^2) - E(X)*E(X)
    EX2 = nsum(d*d, 1) / d.shape[1]
    EX = nsum(d, 1) / d.shape[1]
    VAR = EX2 - EX * EX
    try:
        err = np.sqrt(VAR)
    except ValueError:
        err = np.sqrt(np.where(VAR > 0. , VAR , 1))
    #print sum(data[:,2:-2],1).shape,
    #print sum(data[:,2:-2],0).shape
    for i in range(len(x_axis)):
        o.write("%f   %f    %f\n"%(x_axis[i], med[i], err[i]))
    logging.info("converted %s %s"%(filename, outputfilename))


def deg(x):
    """
    Convert string of radians to degrees helper function
    """
    from math import degrees
    return degrees(float(x))


class cakemacrogenerator:
    """
    Generates a macro to run fit2d to cake some images
    """

    # String to be used on the first time through the caking menu
    first_time_run_string = """GRAPHICAL COORDINATE
           1
 1.0E+00
 1.0E+00
           0
           0
           1
 2.0E00
 2.0E00
           1
 3.0E00
 3.0E00
"""

    # list of jobs to process after fit2d finishes
    do_when_finished = []

    # Options when reading in an image on the powder diffraction menu
    input_options_names = [
        "DARK CURRENT",
        "DC FILE",
        "FLAT-FIELD",
        "FF FILE",
        "FF SCALE",
        "FF MULTIPLIER",
        "SPATIAL DIS.",
        "SD FILE"
        ]

    input_options_values = {
        "DARK CURRENT" : "NO",
        "DC FILE"      : None,
        "FLAT-FIELD"   : "NO",
        "FF FILE"      : None,
        "FF SCALE"     : "NO",
        "FF MULTIPLIER": 1.0 ,
        "SPATIAL DIS." : "NO",
        "SD FILE"      : None,
        }

    # First caking menu parameters
    image_pars_names = [
        "X-PIXEL SIZE",
        "Y-PIXEL SIZE",
        "DISTANCE",
        "WAVELENGTH",
        "X-BEAM CENTRE",
        "Y-BEAM CENTRE",
        "TILT ROTATION",
        "ANGLE OF TILT"
        ]

    image_pars_values = {
        "X-PIXEL SIZE" : 46.77648 ,
        "Y-PIXEL SIZE" : 48.08150 ,
        "DISTANCE"     : 975.0000 ,
        "WAVELENGTH"   : 0.154981 ,
        "X-BEAM CENTRE": 2905.198 ,
        "Y-BEAM CENTRE": 1033.436 ,
        "TILT ROTATION": 154.3416 ,
        "ANGLE OF TILT": 0.899907 }


    # Second fit2d caking menu parameters
    integrate_pars_names = [
        "START AZIMUTH",
        "END AZIMUTH",
        "INNER RADIUS",
        "OUTER RADIUS",
        "SCAN TYPE",
        "1 DEGREE AZ",
        "AZIMUTH BINS",
        "RADIAL BINS",
        "CONSERVE INT.",
        "POLARISATION",
        "GEOMETRY COR.",
        "FACTOR",
        ]

    integrate_pars_values = {
        "START AZIMUTH" :  160.0   ,
        "END AZIMUTH"   : -160.0   ,
        "INNER RADIUS"  :  928.0   ,
        "OUTER RADIUS"  : 2825.0   ,
        "SCAN TYPE"     : "2-THETA",
        "1 DEGREE AZ"   : "NO"     ,
        "AZIMUTH BINS"  : 1        ,
        "RADIAL BINS"   : 2048     ,
        "CONSERVE INT." : "NO"     ,
        "POLARISATION"  : "NO"     ,
        "GEOMETRY COR." : "YES"    ,
        "FACTOR"        : 0.99
        }

    mask_pars_names = [
        "USE MASK",
        "MASK FILE",
        "DIM1_DATA",
        "DIM2_DATA"
        ]
    mask_pars_values = {
        "USE MASK"     : "NO",
        "MASK FILE"    : None,
        "DIM1_DATA"    : "2048",
        "DIM2_DATA"    : "2048"
        }


    # Parameters for writing out results
    input_extn = "edf"
    saving_format = "CHIPLOT"
    output_extn = "chi"

    def __init__(self):
        # Object holds a macro in a string
        self.macro = "I ACCEPT\n"
        self.macro = self.macro + "POWDER DIFFRACTION (2-D)\n"
        # Flag to see if you have been through the caking menu yet
        self.first_time_run = True
        # Flag to check parameters have been read in
        self.parsread = False

    def writepars(self, outparfile):
        """
        # Dump all dictionaries to a file as name value
        """
        pf = open(outparfile, "w")

        d = self.input_options_values
        keys = self.input_options_names
        pf.write("\n# Image correction parameters\n")
        for k in keys:
            pf.write("%s = %s\n"%(k, d[k]))

        d = self.integrate_pars_values
        keys = self.integrate_pars_names
        pf.write("\n# Integration parameters\n")
        for k in keys:
            pf.write("%s = %s\n"%(k, d[k]))

        d = self.image_pars_values
        keys = self.image_pars_names
        pf.write("\n# Experiment geometry parameters\n")
        for k in keys:
            pf.write("%s = %s\n"%(k, d[k]))

        d = self.mask_pars_values
        keys = self.mask_pars_names
        pf.write("\n# Mask and dimensions\n")
        for k in keys:
            pf.write("%s = %s\n"%(k, d[k]))

        pf.write("\n# I/O parameters\n")
        pf.write("input_extn = %s\n"%(self.input_extn))
        pf.write("saving_format = %s\n"%(self.saving_format))
        pf.write("output_extn = %s\n"%(self.output_extn))
        pf.close()

    def readpars(self, inparfile):
        """
        # Read in the parameters for image processing
        """
        pf = open(inparfile, "r").readlines()
        for line in pf:
            if line[0] == "#" or \
               len(line) < 3 or \
               line.find("=") < 0 :
                continue
            try:
                key, value = line.split(" = ")
                value = value.rstrip() # get rid of trailing \n
            except:
                print("$$$%s$$$"% (line))
                raise
            if value == "None":
                value = None
            self.setparameter(key, value)
        self.parsread = True
        self.checkpars()

    def checkpars(self):
        """
        Try to catch some (potentially fatal) errors
        """
        if int(self.mask_pars_values["DIM1_DATA"]) < \
           int(self.integrate_pars_values["RADIAL BINS"]) or \
           int(self.mask_pars_values["DIM2_DATA"]) < \
           int(self.integrate_pars_values["RADIAL BINS"]) or \
           int(self.mask_pars_values["DIM1_DATA"]) < \
           int(self.integrate_pars_values["AZIMUTH BINS"]) or \
           int(self.mask_pars_values["DIM2_DATA"]) < \
           int(self.integrate_pars_values["AZIMUTH BINS"]):
            print("Problem with array dimensions, check your input pars")
            print("RADIAL BINS", self.integrate_pars_values["RADIAL BINS"])
            print("AZIMUTH BINS", self.integrate_pars_values["AZIMUTH BINS"])
            print("DIM1_DATA", self.mask_pars_values["DIM1_DATA"])
            print("DIM2_DATA", self.mask_pars_values["DIM2_DATA"])
            sys.exit()




    def setparameter(self, key, value):
        """
        #        print "Trying to set",key,value
        """
        if key in list(self.integrate_pars_values.keys()):
            self.integrate_pars_values[key] = value
        elif key in list(self.image_pars_values.keys()):
            self.image_pars_values[key] = value
        elif key in list(self.input_options_values.keys()):
            self.input_options_values[key] = value
        elif key in list(self.mask_pars_values.keys()):
            self.mask_pars_values[key] = value
        elif key == "input_extn":
            self.input_extn = value
        elif key == "output_extn":
            self.output_extn = value
        elif key == "saving_format":
            self.saving_format = value
        else:
            print("Failure to set", key, value)

    # Dictionary to translate .fit2d.def file syntax into macro syntax

    fit2d_def_names_to_macro_names = {
        "DARK_CURRENT_CORRECTION" : "DARK CURRENT",
        "DARK_CURRENT_FILE"       : "DC FILE",
        "FF_CORRECTION"           : "FLAT-FIELD",
        "FF_FILE"                 : "FF FILE",
        "FF_SCALE"                : "FF SCALE",
        "FF_SCALER"               : "FF MULTIPLIER",
        "SD_CORRECTION"           : "SPATIAL DIS.",
        "SD_FILE"                 : "SD FILE",
        "X_PIXEL_SIZE"            : "X-PIXEL SIZE",  # meters
        "Y_PIXEL_SIZE"            : "Y-PIXEL SIZE",  #
        # "SCAN_TYPE"              : "SCAN TYPE",
             # Need to figure out how to interpret
        "CAKE_DEFAULT_1_DEGREE"   : "1 DEGREE AZ",
        "CAKE_START_AZIMUTH"      : "START AZIMUTH", # radians
        "CAKE_END_AZIMUTH"        : "END AZIMUTH",   # radians
        #"CAKE_INNER_LIMIT"        : "INNER RADIUS",  # ???
        #"CAKE_OUTER_LIMIT"        : "OUTER RADIUS",
        #"DETECTOR_ROTATION"       : "???",
        "X_BEAM_CENTRE"           : "X-BEAM CENTRE", # pixels, I think
        "Y_BEAM_CENTRE"           : "Y-BEAM CENTRE", #
        "WAVELENGTH"              : "WAVELENGTH",    # meters
        "SAMPLE_DISTANCE"         : "DISTANCE",      # meters
        "TILT_ROTATION"           : "TILT ROTATION", # radians
        "TILT_ANGLE"              : "ANGLE OF TILT",  # radians
        "DIM1_DATA"               : "DIM1_DATA",
        "DIM2_DATA"               : "DIM2_DATA"
        }



    # Convert units from .fit2d.def file to macro units

    fit2d_def_to_macro_converters = {
        "X_PIXEL_SIZE"          : lambda x : float(x)*1e6     , # to microns
        "Y_PIXEL_SIZE"          : lambda x : float(x)*1e6     , # to microns
        "CAKE_START_AZIMUTH"    : deg                         , # to degrees
        "CAKE_END_AZIMUTH"      : deg                         , #
        "WAVELENGTH"            : lambda x : float(x)*1e10    , # to angstrom
        "TILT_ROTATION"         : deg                         , #
        "TILT_ANGLE"            : deg                         , #
        "SAMPLE_DISTANCE"       : lambda x : float(x)*1000.     # to mm
        }


    def generate_pars_from_def_file(self, def_file):
        """
        Attempt to generate parameters from a .fit2d.def file
        So you have run fit2d, all worked and you cleanly exited to generate the def file
        This function picks up the parameters
        """
        def_file_obj = open(def_file, "r")
        while 1:
            try:
                key = def_file_obj.next().rstrip()
                value = def_file_obj.next().rstrip()
                if key in list(self.fit2d_def_names_to_macro_names.keys()):
                    mykey = self.fit2d_def_names_to_macro_names[key]
                    if key in list(self.fit2d_def_to_macro_converters.keys()):
                        # We have a conversion factor to apply
                        myvalue = self.fit2d_def_to_macro_converters[key](value)
                    elif value == "FALSE":
                        myvalue = "NO"
                    elif value == "TRUE":
                        myvalue = "YES"
                    else:
                        myvalue = value
                    self.setparameter(mykey, myvalue)
            except StopIteration:
                break

    def inputfile(self, filein):
        """
        Read in a file on the powder menu in fit2d
        """
        self.macro = self.macro + "INPUT\n" + filein + "\n"
        if self.input_extn.find("mccd") > -1:
            self.macro += "O.K.\n"
        if self.input_extn.find("cor")>-1:
            # hard wire for ID11, sorry
            self.macro += "X-PIXELS\n2048\nY-PIXELS\n2048\nDATA TYPE\nREAL\nSIGNED\nNO\nBYTE SWAP\nNO\n"
            self.macro += "STARTING BYTE\n8193\nO.K.\n"
        for name in self.input_options_names:
            value = self.input_options_values[name]
            if name.find("MASK") >= 0:
                continue
            if value is not None:
                self.macro = self.macro + name + "\n" + str(value) + "\n"
        self.macro = self.macro + "O.K.\n"

    def addmask(self):
        """
        Read in a mask file if requested in input file
        """
        #print "Trying to add mask"
        #print self.mask_pars_values["USE MASK"],
        # self.mask_pars_values["MASK FILE"]
        #print self.mask_pars_values["USE MASK"] == "YES",
        # self.mask_pars_values["MASK FILE"] is not None
        if self.mask_pars_values["USE MASK"] == "YES" and \
           self.mask_pars_values["MASK FILE"] is not None:
            self.macro = self.macro + "MASK\nLOAD MASK\n" +\
                        self.mask_pars_values["MASK FILE"] +\
                        "\nEXIT\n"


    def imagepars(self):
        """
        Set the image parameters after cake -> integrate
        """
        for name in self.image_pars_names:
            value = self.image_pars_values[name]
            if value is not None:
                self.macro = self.macro + name + "\n" + str(value) + "\n"
        self.macro = self.macro + "O.K.\n"

    def cakepars(self):
        """
        Set the caking parameters (second part of cake -> integrate)
        """
        for name in self.integrate_pars_names:
            value = self.integrate_pars_values[name]
            if value is not None:
                self.macro = self.macro + name + "\n" + str(value) + "\n"
        self.macro = self.macro + "O.K.\n"

    def ap(self, *args):
        """
        append to macro
        """
        for a in args:
            self.macro = self.macro + a

    def outputfile(self, filename):
        """
        Save the memory from fit2d (probably your results)
        TODO Make this format friendly - deal with other formats that have special options
        """
        if self.saving_format.find("SPREAD SHEET")>-1:
            self.do_when_finished.append(
                  [spr_to_median, filename, filename[:-3] + "med"]
                  )
            self.ap("EXIT\nEXIT\nIMAGE PROCESSING\nGEOMETRIC\nTRANSPOSE\n")
            self.ap("EXIT\nEXIT\nPOWDER DIFFRACTION\nOUTPUT\n")
            self.ap(self.saving_format)
            self.ap("\nYES\n", filename, "\n")
            self.ap("O.K.\n")
            # Now output a chi files to get the x axis and y axis
            self.macro=self.macro+"OUTPUT\nCHIPLOT\nFILE NAME\n"+filename+"_x\n"
            self.macro=self.macro+"OUTPUT ROWS\nYES\nO.K.\n"
            self.macro=self.macro+"OUTPUT\nCHIPLOT\nFILE NAME\n"+filename+"_y\n"
            self.macro=self.macro+"OUTPUT ROWS\nNO\nO.K.\n"

        elif self.saving_format.find("KLORA")>-1:
             self.macro=self.macro+"EXIT\nOUTPUT\n%s\n"%(self.saving_format)
             self.macro=self.macro+"NO\n"
             self.macro=self.macro+filename+"\n"
#            self.macro=self.macro+"O.K.\n"

        else:
            self.macro=self.macro+"EXIT\nOUTPUT\n%s\n"%(self.saving_format)
            self.macro=self.macro+"FILE NAME\n"
            self.macro=self.macro+filename+"\n"
            self.macro=self.macro+"O.K.\n"


    def cakeafile(self, filein, fileout):
        """
        Take one input file and convert it to one output file
        """
        if not self.parsread:
            raise Exception("Trying to process data with no parameters, give up!!")
        mac = self.macro # local copy
        self.macro=""
        # Read data in (assumes you are in powder diffraction menu)
        self.inputfile(filein)
        # Add mask
        self.addmask()
        # Select cake -> integrate
        self.macro=self.macro+"CAKE\n"
        if self.first_time_run:
            # Give annoying clicks
            self.macro=self.macro+self.first_time_run_string
            self.first_time_run=False
        self.macro=self.macro+"INTEGRATE\n"
        # fill out the first menu parameters (geometry)
        self.imagepars()
        # fill out the second menu parameters (transform)
        self.cakepars()
        # Save results
        self.outputfile(fileout)
        # Conventionally exchange and put previous back in
        self.macro=mac+self.macro+"EXCHANGE\n"

    def cakefileglob(self,stem_arg,replace=False,outstem=None):
        import glob
        filelist=glob.glob("%s????.%s"%(stem_arg,self.input_extn))
        if outstem is None:
            outstem = stem_arg
        outputfilelist = [s.replace(self.input_extn,self.output_extn) for s in filelist ]
        outputfilelist = [s.replace(stem_arg,outstem) for s in outputfilelist ]
        print(filelist)
        print(outputfilelist)
        if not replace:
            alreadydone = glob.glob(outstem+"*"+self.output_extn)
            # 1:1 mapping of input to output
        else:
            alreadydone = []
        for i in range(len(filelist)):
            if outputfilelist[i] not in alreadydone:
                print("%s %s" % (filelist[i], outputfilelist[i]))
                self.cakeafile(filelist[i],outputfilelist[i])




    def cakefileseries(self,stem_arg,first_arg,last_arg,outstem=None):
        """
        Cake a file sequence - eventually emulating fit2d run sequence command
        """
        for i in range(first_arg,last_arg+1):
            if self.input_extn.lstrip() == "":
                filein = "%s.%04d%s" % (stem_arg, i, self.input_extn)
            else:
                filein = "%s%04d.%s" % (stem_arg, i, self.input_extn)
            if outstem == None:
                outstem = stem_arg
            fileout = "%s%04d.%s" % (outstem, i, self.output_extn)
            self.cakeafile(filein,fileout)
            if i%100 == 0:
                print(i, end=' ')
                sys.stdout.flush()

    def run(self,show=True):
        """
        Run the generated macro
        """
        self.macro=self.macro+"EXIT\nEXIT FIT2d\nYES\n"
        import time, tempfile
        tmpfile = tempfile.NamedTemporaryFile(dir=".")
        tf = open( tmpfile.name + ".mac","w")
        tf.write(self.macro)
        tf.close()
        # Send the display to a local black hole
        if show != "SHOW":
            os.system(XVFB)
            time.sleep(1)
            try:
                displaywas=os.environ["DISPLAY"]
            except KeyError:
                displaywas=":0"
            try:
                os.environ["DISPLAY"]=os.environ["HOST"]+":942"
            except KeyError:
                os.environ["DISPLAY"]=":942"
        array=str(self.mask_pars_values["DIM1_DATA"])+"x"+str(self.mask_pars_values["DIM2_DATA"])
        cmd = FIT2D + " -dim%s -mac%s.mac "%(array,os.path.split(
            tmpfile.name)[-1])
        print("executing",cmd)
        os.system(cmd)
        os.system(REMOVE+" "+tmpfile.name+".mac")
        tmpfile.close()
        if show is None:
            os.environ["DISPLAY"]=displaywas
        print("Finalising")
        for job in self.do_when_finished:
            job[0](job[1],job[2])


if __name__=="__main__":

    # already imported sys

    from optparse import OptionParser
    parser=OptionParser(usage="%prog stem_name first last parfile")
    parser.add_option("-c", "--create-parameters", action="store",
          type="string", dest="cpars",
          help="Generate parameter file from DEFFILE file to save in CPARS")
    parser.add_option("-d", "--def-file", action="store", type="string",
                      dest="deffile",
          help="The fit2d .fit2d.def to read, defaults to $HOME/.fit2d.def")
    parser.add_option("-p", "--pars", action="store", type="string",
                      dest="parfile",
                      help="Input parameter file")
    parser.add_option("-s", "--show", action="store", type="string",
                      dest="show",# default="SHOW",
                      help="Show the flickering blue window")
    parser.add_option("-g", "--glob", action="store", type="string",
                      dest="glob",
                      help="Ignore numbers and glob for files, eg -g lab6_")
    parser.add_option("-b", "--debug", action="store", type="string",
                      dest="debug",
                      help="Debugging - write macro to file and stop")
    parser.add_option("-f", "--forceoverwrite", action="store", type="string",
                      dest="forceoverwrite",
                      help="Force overwrite when globbing")
    parser.add_option("-o", "--outstem", action="store", type="string",
                      dest="outstem",
                      help="Output file stem, if different from input")
    options , userargs = parser.parse_args()

    caker=cakemacrogenerator()

    if options.cpars is not  None:
        # already imported os
        if options.deffile == None:
            f = os.path.join( os.environ["HOME"] , ".fit2d.def" )
        else:
            f = options.deffile
        caker.generate_pars_from_def_file( f )
        caker.writepars(options.cpars)
        sys.exit()



    if options.glob is not None:
        caker.readpars(options.parfile)
        caker.cakefileglob(options.glob,options.forceoverwrite is not None,options.outstem)
        caker.run(show=options.show)
        sys.exit()

    try:
        stem = userargs[0]
        first = int(userargs[1])
        last = int(userargs[2])
        parfile = userargs[3]
    except IndexError:
        parser.print_usage()
        sys.exit()

    caker.readpars(parfile)

    caker.cakefileseries(stem,first,last,options.outstem)

    if options.debug is None:
        caker.run(show=options.show)
    else:
        open(options.debug,"w").write(caker.macro)

