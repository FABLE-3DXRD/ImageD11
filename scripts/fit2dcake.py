#!/usr/bin/env python

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

    # Options when reading in an image on the powder diffraction menu
    input_options_names=[
        "DARK CURRENT",
        "DC FILE",
        "FLAT-FIELD",
        "FF FILE",
        "FF SCALE",
        "FF MULTIPLIER",
        "SPATIAL DIS.",
        "SD FILE"
        ]
    input_options_values={
        "DARK CURRENT" : "NO",
        "DC FILE"      : None,
        "FLAT-FIELD"   : "NO",
        "FF FILE"      : None,
        "FF SCALE"     : "NO",
        "FF MULTIPLIER": 1.0 ,
        "SPATIAL DIS." : "NO",
        "SD FILE"      : None
        }


    # First caking menu parameters
    image_pars_names=[    
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


    # Second caking menu parameters
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
        "GEOMETRY COR."
        ]
    integrate_pars_values= {
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
        "GEOMETRY COR." : "YES"
        }


    # Parameters for writing out results
    input_extn = "edf"
    saving_format = "CHIPLOT"
    output_extn = "chi"

    def __init__(self):
        # Object holds a macro in a string
        self.macro="I ACCEPT\n"
        self.macro=self.macro+"POWDER DIFFRACTION (2-D)\n"
        # Flag to see if you have been through the caking menu yet
        self.first_time_run=True

    def writepars(self,parfile):
        # Dump all dictionaries to a file as name value
        pf=open(parfile,"w")
        
        d = self.input_options_values
        keys = self.input_options_names
        pf.write("\n# Image correction parameters\n")
        for k in keys:
            pf.write("%s = %s\n"%(k,d[k]))

        d = self.integrate_pars_values
        keys = self.integrate_pars_names
        pf.write("\n# Integration parameters\n")
        for k in keys:
            pf.write("%s = %s\n"%(k,d[k]))

        d = self.image_pars_values
        keys = self.image_pars_names 
        pf.write("\n# Experiment geometry parameters\n")
        for k in keys:
            pf.write("%s = %s\n"%(k,d[k]))

        pf.write("\n# I/O parameters\n")
        pf.write("input_extn = %s\n"%(self.input_extn))
        pf.write("saving_format = %s\n"%(self.saving_format))
        pf.write("output_extn = %s\n"%(self.output_extn))
        pf.close()

    def readpars(self,parfile):
        # Read in the parameters for image processing
        pf=open(parfile,"r").readlines()
        for line in pf:
            if line[0]=="#" or len(line)<3 or line.find("=")<0:
                continue
            try:
                key,value = line.split(" = ")
                value=value.rstrip() # get rid of trailing \n
            except:
                print "$$$%s$$$"%(line)
                raise
            if value == "None":
                value=None
            self.setparameter(key,value)


    def setparameter(self,key,value):
#        print "Trying to set",key,value
        if key in self.integrate_pars_values.keys():
            self.integrate_pars_values[key]=value  
        elif key in self.image_pars_values.keys():
            self.image_pars_values[key]=value
        elif key in self.input_options_values.keys():
            self.input_options_values[key]=value
        elif key == "input_extn":
            self.input_extn=value
        elif key == "output_extn":
            self.output_extn=value
        elif key == "saving_format":
            self.saving_format=value
        else:
            print "Failure to set",key,value

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
        # "SCAN_TYPE"              : "SCAN TYPE",     # Need to figure out how to interpret
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
        "TILT_ANGLE"              : "ANGLE OF TILT"  # radians
        }


    def deg(x):
        """
        Convert string of radians to degrees helper function
        """
        from math import degrees
        return degrees(float(x))

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

        
    def generate_pars_from_def_file(self,def_file):
        """
        Attempt to generate parameters from a .fit2d.def file
        So you have run fit2d, all worked and you cleanly exited to generate the def file
        This function picks up the parameters
        """
        f=open(def_file,"r")
        while 1:
            try:
                key=f.next().rstrip()
                value=f.next().rstrip()
                if key in self.fit2d_def_names_to_macro_names.keys():
                    mykey=self.fit2d_def_names_to_macro_names[key]
                    if key in self.fit2d_def_to_macro_converters.keys():
                        # We have a conversion factor to apply
                        myvalue =  self.fit2d_def_to_macro_converters[key](value)
                    elif value == "FALSE":
                        myvalue = "NO"
                    elif value == "TRUE":
                        myvalue = "YES"
                    else:
                        myvalue = value
                    self.setparameter(mykey,myvalue)
            except StopIteration:
                break
            
    def inputfile(self,filein):
        """
        Read in a file on the powder menu in fit2d
        """
        self.macro=self.macro+"INPUT\n"+filein+"\n"
        for name in self.input_options_names:
            value=self.input_options_values[name]
            if value!=None:
                self.macro=self.macro+name+"\n"+str(value)+"\n"
        self.macro=self.macro+"O.K.\n"
        
    def imagepars(self):
        """
        Set the image parameters after cake -> integrate
        """
        for name in self.image_pars_names:
            value=self.image_pars_values[name]
            if value!=None:
                self.macro=self.macro+name+"\n"+str(value)+"\n"
        self.macro=self.macro+"O.K.\n"

    def cakepars(self):
        """
        Set the caking parameters (second part of cake -> integrate)
        """
        for name in self.integrate_pars_names:
            value=self.integrate_pars_values[name]
            if value!=None:
                self.macro=self.macro+name+"\n"+str(value)+"\n"
        self.macro=self.macro+"O.K.\n"

    def outputfile(self,file):
        """
        Save the memory from fit2d (probably your results)
        TODO: Make this format friendly - deal with other formats
        """
        self.macro=self.macro+"EXIT\nOUTPUT\n%s\nFILE NAME\n"%(self.saving_format)
        self.macro=self.macro+file+"\n"
        self.macro=self.macro+"O.K.\n"

        
    def cakeafile(self,filein,fileout):
        """
        Take one input file and convert it to one output file
        """
        # Read data in (assumes you are in powder diffraction menu)
        self.inputfile(filein)
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
        # Conventionally exchange
        self.macro=self.macro+"EXCHANGE\n"

    def cakefileseries(self,stem,first,last):
        """
        Cake a file sequence - eventually emulating fit2d run sequence command
        """
        for i in range(first,last+1):
            filein ="%s%04d.%s"%(stem,i,self.input_extn)
            fileout="%s%04d.%s"%(stem,i,self.output_extn)
            self.cakeafile(filein,fileout)

    def run(self):
        """
        Run the generated macro
        """
        self.macro=self.macro+"EXIT\nEXIT FIT2d\nYES\n"
        open("fit2dcake.mac","w").write(self.macro)
        import os, time
        # Send the display to a local black hole
        os.system("Xvfb :1 &")
        time.sleep(1)
        displaywas=os.environ["DISPLAY"]
        os.environ["DISPLAY"]=os.environ["HOST"]+":1"
        os.system("fit2d_12_081_i686_linux2.4.20 -dim2048x2048 -macfit2dcake.mac")
        os.system("rm -f fit2dcake.mac")
        os.environ["DISPLAY"]=displaywas

if __name__=="__main__":

    import sys

    from optparse import OptionParser
    parser=OptionParser(usage="%prog stem_name first last parfile")
    parser.add_option("-c","--create-parameters",action="store",type="string",dest="cpars",
                      help="Generate parameter file from DEFFILE file to save in CPARS")
    parser.add_option("-d","--def-file",action="store",type="string",dest="deffile",
                      help="The fit2d .fit2d.def to read, defaults to $HOME/.fit2d.def")
    options , args = parser.parse_args()
    
    caker=cakemacrogenerator()

    if options.cpars != None:
        import os
        if options.deffile == None:
            f = os.path.join( os.environ["HOME"] , ".fit2d.def" )
        else:
            f=options.deffile
        caker.generate_pars_from_def_file( f )
        caker.writepars(options.cpars)
        sys.exit()
        
    try:
        stem = args[0]
        first = int(args[1])
        last = int(args[2])
        parfile = args[3]
    except:
        parser.print_usage()
        sys.exit()

    caker.readpars(parfile)

    caker.cakefileseries(sys.argv[1],int(sys.argv[2]),int(sys.argv[3]))

    caker.run()
        
    
