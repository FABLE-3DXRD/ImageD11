

 
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
Defines a class for correcting peak positions for spatial distortion
via a fit2d spline file
 
To think about doing - valid regions? What if someone uses a 1K spline
file for a 2K image etc?
"""
import logging
import Numeric as n
import bisplev, math

def readfit2dfloats(fp, nfl):
    """
    Interprets a 5E14.7 formatted fortran line
    """
    ret = []
    j = 0
    while j < nfl:
        i = 0
        thisline = fp.readline()
        # logging.debug("readfit2dfloats:"+thisline)
        while i < 5 * 14:
            # logging.debug(str(i)+ thisline[i:i+14])
            ret.append(float(thisline[i:i+14]) )
            j = j + 1
            i = i + 14
            if j == nfl: 
                break
        i = i + 1
    return ret

class correctorclass: #IGNORE:R0902
    """
    Applies a spatial distortion to a peak position using a fit2d splinefile
    """
    def __init__(self, argsplinefile, orientation="edf"):
        """
        Argument is the name of a fit2d spline file
        """
        self.splinefile = argsplinefile
        self.tolerance = 1e-5
        self.orientation = orientation
        self.pos_lut = None
        self.pixel_lut = None
        self.xmin = self.ymin = self.xmax = self.ymax = 0.0
        self.xsize = self.ysize = 1.0
        self.gridspacing = 0.0
        self.tck1 = None
        self.tck2 = None
        if self.splinefile is not None:
            self.readfit2dspline(self.splinefile)
        


    def correct(self, xin, yin):
        """
        Transform x,y in raw image coordinates into x,y of an
        idealised image. Returns a tuple (x,y), expects a
        pair of floats as arguments
        """
        if self.orientation == "edf":
            xcor = xin + bisplev.bisplev(yin, xin, self.tck2)
            ycor = yin + bisplev.bisplev(yin, xin, self.tck1)
        elif self.orientation == "bruker":
            # fit2d does a flip
            xpos = self.xmax - xin
            xcor = xin - bisplev.bisplev(yin, xpos, self.tck2)
            ycor = yin + bisplev.bisplev(yin, xpos, self.tck1)
        return xcor, ycor

    def make_pixel_lut(self, dims):
        """
        Generate an x and y image which maps the array indices into
        floating point array indices (to be corrected for pixel size later)

        returns 
        FIXME - check they are the right way around
                add some sort of known splinefile testcase
        """
        # Cache the value in case of multiple calls
        if self.pixel_lut is None:
            x_im = n.outerproduct(range(dims[0]), n.ones(dims[1]))
            y_im = n.outerproduct(n.ones(dims[1]), range(dims[0]))  
            self.pixel_lut = self.correct(x_im, y_im)
        return self.pixel_lut

    def make_pos_lut(self, dims):
        """
        Generate a look up table of pixel positions in microns
        # Cache the value in case of multiple calls
        returns ...
        """
        if self.pos_lut is None:
            if self.pixel_lut is None:
                self.make_pixel_lut(dims)                
            self.pos_lut = ( self.pixel_lut[0] * self.xsize, 
                             self.pixel_lut[1] * self.ysize )
        return self.pos_lut

    def distort(self, xin, yin):
        """
        Distort a pair of points xnew, ynew to find where they
        would be in a raw image

        Iterative algorithm...
        """
        yold = yin - bisplev.bisplev(yin, xin, self.tck1)
        xold = xin - bisplev.bisplev(yin, xin, self.tck2)
        # First guess, assumes distortion is constant
        yt = yin - bisplev.bisplev(yold, xold, self.tck1)
        xt = xin - bisplev.bisplev(yold, xold, self.tck2)
        # Second guess should be better
        error = math.sqrt((xt - xold) * (xt - xold) + (yt - yold) * (yt - yold))
        ntries = 0
        while error > self.tolerance:
            ntries = ntries + 1
            xold = xt
            yold = yt
            yt = yin - bisplev.bisplev(yold, xold, self.tck1)
            xt = xin - bisplev.bisplev(yold, xold, self.tck2)
            error = math.sqrt((xt - xold) * (xt - xold) + 
                              (yt - yold) * (yt - yold)   )
            # print error,xold,x,yold,y
            if ntries == 10:
                raise Exception("Error getting the inverse spline to converge")
        return xt, yt

    def test(self, xin, yin):
        """
        Checks that the correct and distort functions are indeed
        inversely related to each other
        """
        xtes, ytes = self.correct(xin, yin)
        xold, yold = self.distort(xtes, ytes)
        error = math.sqrt( (xin - xold) * (xin - xold) + 
                           (yin - yold) * (yin - yold))
        if error > self.tolerance:
            logging.error("Blobcorrector Test Failed!")
            raise Exception("Problem in correctorclass")






            # read the fit2d array into a tck tuple
    def readfit2dspline(self, name):
        """
        Reads a fit2d spline file into a scipy/fitpack tuple, tck
        A fairly long and dull routine...
        """
        fp = open(name, "r")
        # SPATIAL DISTORTION SPLINE INTERPOLATION COEFFICIENTS
        myline = fp.readline() 
        if myline[:7] != "SPATIAL":
            raise SyntaxError, name + \
                ": file does not seem to be a fit2d spline file"
        myline = fp.readline() # BLANK LINE
        myline = fp.readline() # VALID REGION
        myline = fp.readline() # the actual valid region, 
                               # assuming xmin,ymin,xmax,ymax
        logging.debug("xmin,ymin,xmax,ymax, read: "+myline)
        self.xmin, self.ymin, self.xmax, self.ymax = \
         [float(z) for z in myline.split()]
        myline = fp.readline() # BLANK
        myline = fp.readline() # GRID SPACING, X-PIXEL SIZE, Y-PIXEL SIZE
        myline = fp.readline()
        logging.debug("gridspace, xsize, ysize: "+myline)
        self.gridspacing, self.xsize, self.ysize = \
         [float(z) for z in  myline.split()]
        myline = fp.readline() # BLANK
        myline = fp.readline() # X-DISTORTION
        myline = fp.readline() # two integers nx1,ny1
        logging.debug("nx1, ny1 read: "+myline)
        nx1, ny1 = [int(z) for z in myline.split()]
        # Now follow fit2d formatted line 5E14.7
        tx1 = n.array(readfit2dfloats(fp, nx1), n.Float32)
        ty1 = n.array(readfit2dfloats(fp, ny1), n.Float32)
        c1 = n.array(readfit2dfloats(fp, (nx1 - 4) * (ny1 - 4)),  
                     n.Float32)
        myline = fp.readline() #BLANK
        myline = fp.readline() # Y-DISTORTION
        myline = fp.readline() # two integers nx2, ny2
        nx2 , ny2 = [int(z) for z in myline.split()]
        tx2 = n.array(readfit2dfloats(fp, nx2), n.Float32)
        ty2 = n.array(readfit2dfloats(fp, ny2), n.Float32)
        c2 = n.array(readfit2dfloats(fp, (nx2 - 4) * (ny2 - 4)), 
                     n.Float32)
        fp.close()
        # The 3 ,3 is the number of knots
        self.tck1 = (tx1, ty1, c1, 3, 3)
        self.tck2 = (tx2, ty2, c2, 3, 3)




class perfect(correctorclass):
    """
    To use on previously corrected when there is no splinefile
    Allows pixel size etc to be set
    """
    splinefile = "NO_CORRECTION_APPLIED"
    xsize = "UNKNOWN"
    ysize = "UNKNOWN"
    def __init__(self):
        correctorclass.__init__(self, None)
    def correct(self, xin, yin):
        """
        Do nothing - just return the same values
        """
        return xin, yin

#
#"""
#http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/OWENS/LECT5/node5.html
#
#Various interpolation schemes can be used. 
# A common one is bilinear interpolation, given by
#
#v(x,y) = c1x + c2y + c3xy + c4,
#
#where v(x,y) is the grey value at position (x,y).
#Thus we have four coefficients to solve for. We use the known grey values 
# of the 4 pixels
#surrounding the `come from' location to solve for the coefficients.
#
#We need to solve the equation
#
#v1   ( x1 y1 x1y1 1 ) c1
#v2 = ( x2 y2 x2y2 1 ) c2
#v3   ( x3 y3 x3y3 1 ) c3
#v4   ( x4 y4 x4y4 1 ) c4
#
#
#or, in short,
#[V] = [M][C],
#
#which implies
#[C] = [M]-1[V].
#
# This has to be done for every pixel location in the output image and 
# is thus a lot of computation!
# Alternatively one could simply use the integer pixel position closest 
# to the `come from location'.
# This is adequate for most cases.
#
#"""

#def unwarpimage(image, xpositions, ypositions):
#    """
#    xpositions/ypositions are floats giving pixel co-ords of the input image.
#
#    We need the positions of the pixels in the output image,
#     on the input image.
#
#    Hence, for now,
#    """
#    pass

