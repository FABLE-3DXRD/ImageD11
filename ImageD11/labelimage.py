



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
import numpy.oldnumeric as n

class labelimage:
    """
    For labelling spots in diffraction images
    """
    def __init__(self,
                 shape, 
                 fileout=sys.stdout,
                 spatial=blobcorrector.perfect):
        """
        Shape - image dimensions
        fileout - writeable stream for merged peaks
        spatial - correction of of peak positions
        """
        self.shape=shape
        self.bl = Numeric.zeros(shape,Numeric.Int)
        self.lastbl = Numeric.zeros(shape,Numeric.Int)
        self.np = 0
        self.lastres = None
        self.lastnp="FIRST"
        self.verbose=0
        self.nprop = 9 # number of accumulating properties
        self.corrector=spatial
        self.closed = None
        self.finalpeaks=[]
        self.outfile=open(fileout,"w")
        self.outfile.write(
            "# xc yc omega npixels avg_intensity x_raw y_raw sigx sigy covxy\n")

    def peaksearch(self,data,threshold,omega,dark=None,flood=None):
        """
        # Call the c extensions to do the peaksearch, on entry:
        #
        # data = 2D Numeric array (of your data)
        # threshold = float - pixels above this number are put into objects
        """
        corr_data = data
        if dark is not None and flood is None:
            corr_data = data - dark
        if dark is None flood is not None: 
            corr_data = data / flood
        if dark is not None and flood is not None:
            corr_data = (data - dark)/flood
        self.np = connectedpixels.connectedpixels(corr_data, 
                                                  self.bl, 
                                                  threshold,
                                                  omega,
                                                  self.verbose)

        if self.np>0:
            self.res = connectedpixels.blobproperties(corr_data, 
                                                 self.bl, 
                                                 self.np,
                                                 **optarg)
        else:
            self.res = None


    def mergelast(self):
        """
        Merge the last two images searches
        """
        if self.lastnp == "FIRST":
            # No previous image available, this was the first
            # Swap the blob images
            t = self.lastbl
            self.lastbl = self.bl
            self.bl = t
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
        for i in range(



