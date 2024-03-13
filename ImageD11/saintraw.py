
from __future__ import print_function

"""
Module for read saint raw files from bruker integration
"""
import numpy
from ImageD11 import columnfile


docsheader = """
...from the file sairefl._tp...
Copyright Bruker

ITEM                  FORMAT   DESCRIPTION
----                  ------   -----------

COMMENT                 '!'    An exclamation mark in the first column
                               indicates that the line contains a
                               comment.  (If the line is not a
                               comment, the record starts in the first
                               column, not the second; that is, the
                               first column is not "reserved" for a
                               possible exclamation mark)
"""

docs = """
IHKL                   3I4     HKL indices of the reflection      

  #IMNP                   3I4     PRESENT ONLY FOR MODULATED-STRUCTURE DATA
                               MNP indices of the reflection.  Unused
                               indices  (in 4- or 5-dimensional cases)
                               are written as zero

FI                     F8.x    Integrated intensity in photons, background-
                               subtracted, normalized to 1 min/deg, and
                               Lp corrected.  Number of places to right
                               of the decimal point is adjusted to
                               balance precision vs. overflow;  exponential
                               format is used for very large values.

SI                     F8.x    Estimated standard deviation of the
                               intensity.  Number of places to right
                               of the decimal point is adjusted to
                               balance precision vs. overflow;  exponential
                               format is used for very large values.

IBATNO                 I4      Batch number (for scaling) assigned to
                               the reflection

COSINES                6F8.5   Direction cosines of the incident and
                               diffracted beams, relative to the 
                               unit cell axes.  Order is XI,XD,YI,YD,
                               ZI,ZD where XI is the X-component of
                               the incident beam, XD is the X-component
                               of the diffracted beam, etc.  Used
                               for absorption correction.


MSTATUS                I3      1X,Z2   Status mask reserved for flagging
                               abnormal conditions.  There are
                               currently none defined.

XO                     F7.2    Observed X-pixel coordinate of the
                               intensity-weighted reflection centroid,
                               in reduced pixels (scale of 512x512)

YO                     F7.2    Observed Y-pixel coordinate of the
                               intensity-weighted reflection centroid,
                               in reduced pixels (scale of 512x512)

ZO                     F8.2    Observed frame number of the
                               intensity-weighted reflection centroid

XP                     F7.2    Predicted X-pixel of the reflection
                               in reduced pixels (scale of 512x512)

YP                     F7.2    Predicted Y-pixel of the reflection
                               in reduced pixels (scale of 512x512)

ZP                     F8.2    Predicted frame number of the reflection

CORLPAF                F6.3    Multiplicative correction for Lorentz
                               effect, polarization, air absorption,
                               and detector faceplate absorption,
                               already applied to integrated intensity,
                               FI

CORREL                 F5.2    The correlation coefficient between the
                               3-D profile observed for this reflection
                               and the corresponding model 3-D profile

ACCUMTIME              F7.2    Accumulated hours of exposure

SWING                  F7.2    Detector swing angle in degrees

ROT                    F7.2    Scan axis (omega or phi) setting in
                               degrees at which the reflection was
                               observed

IAXIS                  I2      Scan axis number (2=omega, 3=phi)

ISTL                   I5      SIN(THETA)/LAMBDA times 10,000 

PKSUM                  I9      Total raw peak counts in photons, with no
                               correction for background, normalization, 
                               or Lp

BGAVG                  I7      Normally, average BG per pixel in 
                               photons * 1000. In data sets where the 
                               background was reported to be very large 
                               ( > 1024 photons after 1 min/deg 
                               normalization), the program issues a warning 
                               message during integration and the scale of
                               1000 is omitted.

ROTREL                 F7.2    Rotation of the scan axis (omega or phi)
                               in degrees relative to the start of the
                               run

ICRYST                 I4      Crystal number (for scaling)

PKFRAC                 F6.3    Fraction of the profile volume nearest
                               this HKL.  1-PKFRAC is the fraction of
                               the intensity which had to be estimated
                               from the model profiles because of
                               overlap with neighboring spots

IKEY                   I11     The unique sort key for the group of
                               equivalents to which this HKL belongs.
                               IKEY = 1048576 * (H+511) + 1024 *
                               (K+511) + (L+511), where H,K and L are
                               the HKL indices transformed to the base
                               asymmetric unit.

IMULT                  I3      The point-group multiplicity of this HKL

CORLORENTZ             F6.3    Lorentz correction

XGEO                   F8.2    Spatially corrected X relative to beam
                               center, in pixels

YGEO                   F8.2    Spatially corrected Y relative to beam 
                               center, in pixels

CHI                    F8.3    Chi setting angle, degrees

OTHER_ANGLE            F8.3    The remaining setting angle, degrees.  This
                               will be phi if scans were in omega, or omega
                               if scans were in phi

ICOMPONENT             I4      Twin component number in SHELXTL HKLF 5
                               convention.  In a multi-component overlap,
                               ICOMPONENT is negated for all but the last
                               record of the overlap

"""

    
    

class saintraw(object):
    doc = docs
    titles  = []
    formats = {}
    helps   = {}
    def __init__(self, filename=None):
        """
        filename = filename to read in
        """
        self.parsedocs()
        if filename is not None:
            self.read(filename)
    
    def parsedocs(self):
        """
        Parse the saint documentation for the Bruker format
        """
        self.titles = []
        title = help = format = None
        for line in self.doc.split("\n"):
            if len(line.rstrip()) == 0:
                if title is not None:
                    self.formats[title] = format
                    self.helps[title]   = help
                    self.titles.append(title)
                    title = None
                    format = None
                continue
            if line[0] != " ":
                title, format = line.split()[0:2]
                help = " ".join(line.split()[2:])
            else:
                help = " ".join([help, line.lstrip()])

        alltitles = []
        slices = []
        funcs = []
        i = 0
        allformats = []
        for t in self.titles:
            f = self.formats[t]
            if f[0].isdigit():
                n = int(f[0])
                f = f[1:]
            else:
                n = 1
            if n > 1:
                for j in range(n):
                    alltitles.append( t + "_%d" % (j) )
                    allformats.append( f )
            else:
                alltitles.append( t )
                allformats.append( f )
            assert f[0] in ["I","F"]
            if f[0] == "I":
                for dummy in range(n):
                    funcs.append( int )
            if f[0] == "F":
                for dummy in range(n):
                    funcs.append( float )
            num = int(f[1:].split(".")[0])
            for dummy in range(n):
                slices.append( slice( i, i + num ) )
                i += num
        self.alltitles = alltitles
        self.allformats = allformats
        self.funcs = funcs
        self.slices = slices
        assert len(funcs) == len(slices)
        assert len(slices) == len(alltitles)
            
    def read(self, filename):
        """
        Read an ascii formatted saint reflection file
        """
        self.data = {}
        self.lines = open(filename,"r").readlines()
        for t in self.alltitles:
            self.data[t] = []
        zipped = list(zip(self.alltitles, self.slices, self.funcs))
        for line in self.lines:
            if line[0] == "!":
                # Comment line
                continue
            for t,s,f in zipped:
                # Parse this line
                try:
                    self.data[t].append( f( line[s] ) )
                except:
                    print(t,s,f)
                    raise

    def condition_filter(self, name, func):
        """
        Remove the peaks according to condition
        """
        assert len(self.lines) == len(self.data[name] )
        indices = numpy.compress( func( numpy.array( self.data[name]) ) , 
                                  list(range(len(self.lines))) )
        self.take( indices )                                

    def take(self, order):
        """
        Put the peaks in the order given in order (indices)
        """
        for t in list(self.data.keys()):
            self.data[t] = numpy.take( self.data[t],
                                       order) 
        self.lines = list( numpy.take( self.lines,
                                       order))
    
    def sort(self, name):
        """
        Sort according to a column in self.data
        """
        order = numpy.argsort( self.data[name] )
        self.take(order) 


    def write(self, filename):
        """
        Write an ascii formatted saint reflection file

        """
        outf = open(filename, "w")
        for line in self.lines:
            outf.write( line )
        # raise Exception("Not implemented writing yet!")

    def tocolumnfile(self):
        """
        Return a columnfile
        """
        cof = columnfile.newcolumnfile( self.alltitles )
        dlist = [ self.data[t] for t in self.alltitles ]
        cof.bigarray = numpy.array( dlist, float )
        cof.nrows = len( self.data[ self.alltitles[0] ] )
        cof.ncols = len( self.alltitles ) 
        cof.set_attributes()
        return cof

if __name__ == "__main__":

    import sys, time
    START = time.time()
    sra = saintraw()
    print("Making object", time.time() - START)
    
    START = time.time()
    sra.read(sys.argv[1])
    print("Reading", time.time() - START)
    
    print(len(sra.data['IHKL_0']))

    START = time.time()
    cra = sra.tocolumnfile()
    print(cra.bigarray.shape)
    print("Convert to colfile", time.time() - START)
