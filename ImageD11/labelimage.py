



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


from ImageD11 import blobcorrector, connectedpixels
# Names of property columns in array
from ImageD11.connectedpixels import s_1, s_I, s_I2, \
    s_fI, s_ffI, s_sI, s_ssI, s_sfI, s_oI, s_ooI, s_foI, s_soI, \
    bb_mn_f, bb_mn_s, bb_mx_f, bb_mx_s, bb_mn_o, bb_mx_o, \
    mx_I, mx_I_f, mx_I_s, mx_I_o

from math import sqrt

import sys

import numpy.oldnumeric as n


# These should match the definitions in 
# /sware/exp/saxs/doc/SaxsKeywords.pdf
def flip1(x,y): return x,y
def flip2(x,y): return -x,y   
def flip3(x,y): return x,-y
def flip4(x,y): return -x,-y
def flip5(x,y): return y,x
def flip6(x,y): return y,-x
def flip7(x,y): return -y,x
def flip8(x,y): return -y,-x





class labelimage:
    """
    For labelling spots in diffraction images
    """

    titles = "#  fc  sc  omega" 
    format = "  %.4f"*3
    titles += "  Number_of_pixels"
    format += "  %.0f"
    titles += "  avg_intensity  f_raw  s_raw  sigf  sigs  covsf"
    format += "  %.4f"*6
    titles += "  sigo  covso  covfo"
    format += "  %.4f"*3
    titles += "  sum_intensity"
    format += "  %.4f"
    titles += "  IMax_int  IMax_f  IMax_s  IMax_o"
    format += "  %.4f  %.0f  %.0f  %.4f"
    titles += "  Min_f  Max_f  Min_s  Max_s  Min_o  Max_o"
    format += "  %.0f"*4 + "  %.4f"*2
    titles += "  dety  detz"
    format += "  %.4f"*2
    titles += "\n"
    format += "\n"


    def __init__(self,
                 shape, 
                 fileout = sys.stdout,
                 spatial = blobcorrector.perfect(),
                 flipper = flip2 ):
        """
        Shape - image dimensions
        fileout - writeable stream for merged peaks
        spatial - correction of of peak positions
        """
        self.shape = shape
        self.bl = n.zeros(shape,n.Int)
        self.lastbl = n.zeros(shape,n.Int)
        self.np = 0
        self.lastres = None
        self.lastnp = "FIRST"
        self.verbose = 0
        self.corrector = spatial
        self.closed = None
        self.finalpeaks = []
        if hasattr(fileout,"write"):
            self.outfile = fileout
        else:
            self.outfile = open(fileout,"w")
        self.spot3d_ID = 0 # counter for printing
        self.outfile.write(self.titles)
        self.fs2yz = flipper


    def peaksearch(self,data,threshold,omega,dark=None,flood=None):
        """
        # Call the c extensions to do the peaksearch, on entry:
        #
        # data = 2D Numeric array (of your data)
        # threshold = float - pixels above this number are put into objects
        """
        corr_data = data
        if dark is not None   and     flood is None:
            corr_data = data - dark
        if dark is None       and     flood is not None: 
            corr_data = data / flood
        if dark is not None   and     flood is not None:
            corr_data = (data - dark)/flood
        self.np = connectedpixels.connectedpixels(corr_data, 
                                                  self.bl, 
                                                  threshold,
                                                  self.verbose)
        if self.np>0:
            self.res = connectedpixels.blobproperties(corr_data, 
                                                      self.bl, 
                                                      self.np,
                                                      omega=omega)
        else:
            self.res = None

    def mergelast(self):
        """
        Merge the last two images searches
        """
        if self.lastnp == "FIRST":
            # No previous image available, this was the first
            # Swap the blob images
            self.lastbl, self.bl = self.bl, self.lastbl
            self.lastnp = self.np
            self.lastres = self.res
            return
        ds = connectedpixels.bloboverlaps(self.lastbl,
                                          self.lastnp,
                                          self.lastres,
                                          self.bl,
                                          self.np,
                                          self.res,
                                          self.verbose)
        # lastres is now moved forward into res
        self.outputpeaks(self.lastres)
        self.lastnp = self.np   # This is array dim
        self.lastres = self.res # free old lastres I hope

    def outputpeaks(self, peaks):
        """
        Peaks are in Numeric arrays nowadays
        """
        for p in peaks:
            if p[s_1] < 0.1:
                # Merged with another
                continue
            avg_I = p[s_I]/p[s_1]
            # first moments
            fraw = p[s_fI]/p[s_I]
            sraw = p[s_sI]/p[s_I]
            oraw = p[s_oI]/p[s_I]
            # Spline correction
            fc, sc = self.corrector.correct(fraw, sraw)
            # second moments - the plus 1 is the zero width = 1 pixel
            ss = sqrt( p[s_ssI]/p[s_I] - sraw*sraw + 1 ) 
            ff = sqrt( p[s_ffI]/p[s_I] - fraw*fraw + 1 )
            oo = sqrt( p[s_ooI]/p[s_I] - oraw*oraw + 1 )
            sf = ( p[s_sfI]/p[s_I] - sraw*fraw )/ss/ff
            so = ( p[s_soI]/p[s_I] - sraw*oraw )/ss/oo
            fo = ( p[s_foI]/p[s_I] - fraw*oraw )/ff/oo
            dety, detz = self.fs2yz(fraw, sraw)
            self.outfile.write(self.format % (
                    fc, sc, oraw, 
                    p[s_1], avg_I,    
                    fraw, sraw,
                    ff, ss, sf, 
                    oo, so, fo, 
                    p[s_I],     
                    p[mx_I],p[mx_I_f],p[mx_I_s],p[mx_I_o], 
                    p[bb_mn_f],p[bb_mx_f],p[bb_mn_s],p[bb_mx_s],
                    p[bb_mn_o],p[bb_mx_o],
                    dety, detz ))
            self.spot3d_ID += 1
            
            
            
    def finalise(self):
        self.outputpeaks(self.lastres)



