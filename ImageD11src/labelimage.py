
from __future__ import print_function

"""
Class to wrap the cImageD11 c extensions for labelling
blobs in images.
"""



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


from ImageD11 import blobcorrector, cImageD11
# Names of property columns in array
from ImageD11.cImageD11 import s_1, s_I, s_I2,\
    s_fI, s_ffI, s_sI, s_ssI, s_sfI, s_oI, s_ooI, s_foI, s_soI, \
    bb_mn_f, bb_mn_s, bb_mx_f, bb_mx_s, bb_mn_o, bb_mx_o, \
    mx_I, mx_I_f, mx_I_s, mx_I_o, dety, detz, \
    avg_i, f_raw, s_raw, o_raw, f_cen, s_cen, \
    m_ss, m_ff, m_oo, m_sf, m_so, m_fo


from math import sqrt

import sys

import numpy as np


# These should match the definitions in
# /sware/exp/saxs/doc/SaxsKeywords.pdf
def flip1(x, y):
    """ fast, slow to dety, detz"""
    return  x,  y
def flip2(x, y):
    """ fast, slow to dety, detz"""
    return -x,  y
def flip3(x, y):
    """ fast, slow to dety, detz"""
    return  x, -y
def flip4(x, y):
    """ fast, slow to dety, detz"""
    return -x, -y
def flip5(x, y):
    """ fast, slow to dety, detz"""
    return  y,  x
def flip6(x, y):
    """ fast, slow to dety, detz"""
    return  y, -x
def flip7(x, y):
    """ fast, slow to dety, detz"""
    return -y,  x
def flip8(x, y):
    """ fast, slow to dety, detz"""
    return -y, -x





class labelimage:
    """
    For labelling spots in diffraction images
    """

    titles = "#  sc  fc  omega"
    format = "  %.4f"*3
    titles += "  Number_of_pixels"
    format += "  %.0f"
    titles += "  avg_intensity  s_raw  f_raw  sigs  sigf  covsf"
    format += "  %.4f"*6
    titles += "  sigo  covso  covfo"
    format += "  %.4f"*3
    titles += "  sum_intensity  sum_intensity^2"
    format += "  %.4f  %.4f"
    titles += "  IMax_int  IMax_s  IMax_f  IMax_o"
    format += "  %.4f  %.0f  %.0f  %.4f"
    titles += "  Min_s  Max_s  Min_f  Max_f  Min_o  Max_o"
    format += "  %.0f"*4 + "  %.4f"*2
    titles += "  dety  detz"
    format += "  %.4f"*2
    titles += "  onfirst  onlast  spot3d_id"
    format += "  %d  %d  %d"
    titles += "\n"
    format += "\n"


    def __init__(self,
                 shape,
                 fileout = sys.stdout,
                 spatial = blobcorrector.perfect(),
                 flipper = flip2,
                 sptfile = sys.stdout ):
        """
        Shape - image dimensions
        fileout - writeable stream for merged peaks
        spatial - correction of of peak positions
        """
        self.shape = shape  # Array shape
        if not hasattr(sptfile,"write"):
            self.sptfile = open(sptfile, "w")
        else:
            self.sptfile = sptfile # place for peaksearch to print - file object
        self.corrector = spatial  # applies spatial distortion
        self.fs2yz = flipper # generates y/z

        self.onfirst = 1    # Flag for first image in series
        self.onlast = 0     # Flag for last image in series
        self.blim = np.zeros(shape, np.int32)  # 'current' blob image
        self.npk = 0        #  Number of peaks on current
        self.res = None     #  properties of current

        self.threshold = None # cache for writing files

        self.lastbl = np.zeros(shape, np.int32)# 'previous' blob image
        self.lastres = None
        self.lastnp = "FIRST" # Flags initial state

        self.verbose = 0    # For debugging


        if hasattr(fileout,"write"):
            self.outfile = fileout
        else:
            self.outfile = open(fileout,"w")

        self.spot3d_id = 0 # counter for printing
        try:
            self.outfile.write(self.titles)
        except:
            print(type(self.outfile),self.outfile)
            raise




    def peaksearch(self, data, threshold, omega):
        """
        # Call the c extensions to do the peaksearch, on entry:
        #
        # data = 2D Numeric array (of your data)
        # threshold = float - pixels above this number are put into objects
        """
        d = data.astype(np.float32)
        self.labelpeaks( d, threshold )
        self.measurepeaks( d, omega )


    def labelpeaks(self, data, threshold):
        self.threshold = threshold
        self.npk = cImageD11.connectedpixels(data.astype(np.float32),
                                             self.blim,
                                             threshold,
                                             self.verbose)

    def measurepeaks(self, data, omega, blim=None):
        if blim is not None:
            self.npk = blim.max()
            self.blim = blim
        if self.npk > 0:
            self.res = cImageD11.blobproperties(data.astype(np.float32),
                                                self.blim,
                                                self.npk,
                                                omega=omega)
        else:
            # What to do?
            self.res = None
        
        

    def mergelast(self):
        """
        Merge the last two images searches
        """
        if self.lastnp == "FIRST":
            # No previous image available, this was the first
            # Swap the blob images
            self.lastbl, self.blim = self.blim, self.lastbl
            self.lastnp = self.npk
            self.lastres = self.res
            return
        if self.npk > 0 and self.lastnp > 0:
            # Thanks to Stine West for finding a bug here
            #
            self.npk = cImageD11.bloboverlaps(self.lastbl,
                                              self.lastnp,
                                              self.lastres,
                                              self.blim,
                                              self.npk,
                                              self.res,
                                              self.verbose)
        if self.lastnp > 0:
            # Fill out the moments of the "closed" peaks
            # print "calling blobmoments with",self.lastres
            cImageD11.blob_moments(self.lastres[:self.lastnp])
            # Write them to file
            self.outputpeaks(self.lastres[:self.lastnp])
        # lastres is now moved forward into res
        self.lastnp = self.npk   # This is array dim
        if self.npk > 0:
            self.lastres = self.res[:self.npk]  # free old lastres I hope
        else:
            self.lastres = None
        # Also swap the blob images
        self.lastbl, self.blim = self.blim, self.lastbl

    def output2dpeaks(self, file_obj):
        """
        Write something compatible with the old ImageD11 format
        which fabian is reading.
        This is called before mergelast, so we write self.npk/self.res
        """

        file_obj.write("# Threshold level %f\n"%( self.threshold))
        file_obj.write("# Number_of_pixels Average_counts    s   f     sc   fc      sig_s  sig_f  cov_sf  IMax_int\n")
        cImageD11.blob_moments(self.res)

        fs = "%d  "+ "%f  "*9 + "\n"
        for i in self.res[:self.npk]:
            if i[s_1] < 0.1:
                raise Exception("Empty peak on current frame")
            i[s_cen], i[f_cen] = self.corrector.correct(i[s_raw], i[f_raw])
            file_obj.write(fs % (i[s_1],  i[avg_i], i[s_raw], i[f_raw],
                                 i[s_cen], i[f_cen],
                                 i[m_ss], i[m_ff], i[m_sf], i[mx_I]))
        file_obj.write("\n")

    def outputpeaks(self, peaks):
        """
        Peaks are in Numeric arrays nowadays
        """
        for i in peaks:
            if i[s_1] < 0.1:
                # Merged with another
                continue

            # Spline correction
            i[s_cen], i[f_cen] = self.corrector.correct(i[s_raw], i[f_raw])
            i[dety], i[detz] = self.fs2yz(i[f_raw], i[s_raw])
            self.outfile.write(self.format % (
                    i[s_cen], i[f_cen], i[o_raw],
                    i[s_1], i[avg_i],
                    i[s_raw], i[f_raw],
                    i[m_ss], i[m_ff], i[m_sf],
                    i[m_oo], i[m_so], i[m_fo],
                    i[s_I],i[s_I2],
                    i[mx_I],i[mx_I_s],i[mx_I_f],i[mx_I_o],
                    i[bb_mn_s],i[bb_mx_s],i[bb_mn_f],i[bb_mx_f],
                    i[bb_mn_o],i[bb_mx_o],
                    i[dety], i[detz],
                    self.onfirst, self.onlast, self.spot3d_id ))
            self.spot3d_id += 1
        if self.onfirst > 0:
            self.onfirst = 0




    def finalise(self):
        """
        Write out the last frame
        """
        self.onlast = 1
        if self.lastres is not None:
            cImageD11.blob_moments(self.lastres)
            self.outputpeaks(self.lastres)
        #if hasattr(self.sptfile, "close"):
        #    self.sptfile.close()
        #     wonder what that does to stdout



