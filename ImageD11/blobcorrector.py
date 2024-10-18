from __future__ import print_function

 
# ImageD11 Software for beamline ID11
# Copyright (C) 2021  Jon Wright
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

2021 : added LUT method for eiger
"""
import logging, numpy, math
import numba
import fabio
from scipy.interpolate import bisplev

def readfit2dfloats(filep, nfl):
    """
    Interprets a 5E14.7 formatted fortran line
    filep = file object (has readline method)
    nfl   = the number of floats to read
    """
    ret = []
    j = 0
    while j < nfl:
        i = 0
        thisline = filep.readline()
        # logging.debug("readfit2dfloats:"+thisline)
        while i < 5 * 14:
            # logging.debug(str(i)+ thisline[i:i+14])
            ret.append(float(thisline[i:i+14]) )
            j = j + 1
            i = i + 14
            if j == nfl: 
                break
    return ret

class correctorclass: #IGNORE:R0902
    """
    Applies a spatial distortion to a peak position using a fit2d splinefile
    """
    def __init__(self, argsplinefile, orientation="edf"):
        """
        Argument is the name of a fit2d spline file
        """
        import warnings
        warnings.warn("For new data from ID11, better to use the dx,dy files instead of Fit2d spline", DeprecationWarning)
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
            self.dim = ( int(self.xmax-self.xmin), int(self.ymax-self.ymin))  # detector dimensions read from splinefile


    def correct(self, xin, yin):
        """
        Transform x,y in raw image coordinates into x,y of an
        idealised image. Returns a tuple (x,y), expects a
        pair of floats as arguments
        """
        if self.orientation == "edf":
            xcor = xin + bisplev(yin, xin, self.tck2)
            ycor = yin + bisplev(yin, xin, self.tck1)
        else: 
            # fit2d does a flip 
            raise Exception("Spline orientations must be edf, convert "
                            "your image to edf and remake the spline")
            # Unreachable code - we no longer accept this complexity
            # it means the spline file for ImageD11 bruker images
            # is not the same as for fit2d. 
            # xpos = self.xmax - xin
            # xcor = xin - bisplev(yin, xpos, self.tck2)
            # ycor = yin + bisplev(yin, xpos, self.tck1)
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
            x_im = numpy.outer(numpy.arange(dims[0]), numpy.ones(dims[1]))
            y_im = numpy.outer(numpy.ones(dims[0]), numpy.arange(dims[1]))
            # xcor is tck2
            x_im = numpy.add( x_im,
                              bisplev( numpy.arange(dims[1]),
                                         numpy.arange(dims[0]),
                                               self.tck2 ).T,
                              x_im)
            # ycor is tck1
            y_im = numpy.add( y_im,
                              bisplev( numpy.arange(dims[1]),
                                               numpy.arange(dims[0]),
                                               self.tck1 ).T,
                              y_im)
            self.pixel_lut = x_im, y_im
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


    def correct_px_lut(self, pks):
        """Transform x,y in raw image coordinates into x,y of an
        idealised image using a LUT. Faster than standard correct method for large peakfiles.
        pks : ImageD11 columnfile with raw coordinates s_raw, f_raw add new colulmns sc, fc with corrected coordinates"""
        
        assert self.dim is not None
        assert 's_raw' in pks.titles and 'f_raw' in pks.titles	

        # make pixel_lut  + substract xy grid coordinate (i,j) to keep only dx and  dy arrays.
        self.make_pixel_lut(self.dim)
        i,j = numpy.mgrid[ 0:self.dim[0], 0:self.dim[1] ]
        dx = self.pixel_lut[0] - i
        dy = self.pixel_lut[1] - j

        # get integer pixel index (si,fi) of each peak
        si = numpy.round(pks['s_raw']).astype(int)
        fi = numpy.round(pks['f_raw']).astype(int)
    
        # apply dx dy correction on s_raw / f_raw
        sc = (dx[ si, fi ] + pks.s_raw).astype(numpy.float32)
        fc = (dy[ si, fi ] + pks.f_raw).astype(numpy.float32)
        # add corrected arrays as new columns
        pks.addcolumn(sc, 'sc')
        pks.addcolumn(fc, 'fc')
    

    def distort(self, xin, yin):
        """
        Distort a pair of points xnew, ynew to find where they
        would be in a raw image

        Iterative algorithm...
        """
        yold = yin - bisplev(yin, xin, self.tck1)
        xold = xin - bisplev(yin, xin, self.tck2)
        # First guess, assumes distortion is constant
        ytmp = yin - bisplev(yold, xold, self.tck1)
        xtmp = xin - bisplev(yold, xold, self.tck2)
        # Second guess should be better
        error = math.sqrt((xtmp - xold) * (xtmp - xold) + 
                          (ytmp - yold) * (ytmp - yold)   )
        ntries = 0
        while error > self.tolerance:
            ntries = ntries + 1
            xold = xtmp
            yold = ytmp
            ytmp = yin - bisplev(yold, xold, self.tck1)
            xtmp = xin - bisplev(yold, xold, self.tck2)
            error = math.sqrt((xtmp - xold) * (xtmp - xold) + 
                              (ytmp - yold) * (ytmp - yold)   )
            # print error,xold,x,yold,y
            if ntries == 10:
                raise Exception("Error getting the inverse spline to converge")
        return xtmp, ytmp

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
        fin = open(name, "r")
        # SPATIAL DISTORTION SPLINE INTERPOLATION COEFFICIENTS
        myline = fin.readline() 
        if myline[:7] != "SPATIAL":
            raise SyntaxError(name + \
                ": file does not seem to be a fit2d spline file")
        fin.readline() # BLANK LINE
        fin.readline() # VALID REGION
        myline = fin.readline() # the actual valid region, 
                               # assuming xmin,ymin,xmax,ymax
        logging.debug("xmin,ymin,xmax,ymax, read: "+myline)
        self.xmin, self.ymin, self.xmax, self.ymax = \
         [float(z) for z in myline.split()]
        myline = fin.readline() # BLANK
        myline = fin.readline() # GRID SPACING, X-PIXEL SIZE, Y-PIXEL SIZE
        myline = fin.readline()
        logging.debug("gridspace, xsize, ysize: "+myline)
        self.gridspacing, self.xsize, self.ysize = \
         [float(z) for z in  myline.split()]
        fin.readline() # BLANK
        fin.readline() # X-DISTORTION
        myline = fin.readline() # two integers nx1,ny1
        logging.debug("nx1, ny1 read: "+myline)
        nx1, ny1 = [int(z) for z in myline.split()]
        # Now follow fit2d formatted line 5E14.7
        tx1 = numpy.array(readfit2dfloats(fin, nx1), numpy.float32)
        ty1 = numpy.array(readfit2dfloats(fin, ny1), numpy.float32)
        cf1 = numpy.array(readfit2dfloats(fin, (nx1 - 4) * (ny1 - 4)),  
                          numpy.float32)
        fin.readline() #BLANK
        fin.readline() # Y-DISTORTION
        myline = fin.readline() # two integers nx2, ny2
        nx2 , ny2 = [int(z) for z in myline.split()]
        tx2 = numpy.array(readfit2dfloats(fin, nx2), numpy.float32)
        ty2 = numpy.array(readfit2dfloats(fin, ny2), numpy.float32)
        cf2 = numpy.array(readfit2dfloats(fin, (nx2 - 4) * (ny2 - 4)), 
                     numpy.float32)
        fin.close()
        # The 3 ,3 is the number of knots
        self.tck1 = (tx1, ty1, cf1, 3, 3)
        self.tck2 = (tx2, ty2, cf2, 3, 3)




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
            x_im = numpy.outer(numpy.arange(dims[0]), numpy.ones(dims[1]))
            y_im = numpy.outer(numpy.ones(dims[0]), numpy.arange(dims[1]))
            self.pixel_lut = x_im, y_im
        return self.pixel_lut

@numba.njit(parallel=True)
def apply_lut_parallel(sr, fr, sc, fc, dx, dy):
    for i in numba.prange(len(sr)):
        si = int(numpy.round( sr[i] ))
        fi = int(numpy.round( fr[i] ))
        sc[i] = sr[i] + dy[ si, fi ]
        fc[i] = fr[i] + dx[ si, fi ]

class eiger_spatial(object):
    
    def __init__(self, 
                 dxfile="/data/id11/nanoscope/Eiger/spatial_20210415_JW/e2dx.edf",
                 dyfile="/data/id11/nanoscope/Eiger/spatial_20210415_JW/e2dy.edf",):
        self.dx = fabio.open(dxfile).data  # x == fast direction at ID11
        self.dy = fabio.open(dyfile).data  # y == slow direction
        assert self.dx.shape == self.dy.shape

    def __call__(self, pks, parallel=None):
        n = len(pks['s_raw'])
        if parallel is None:
            parallel = n > 1e6
        if parallel:
            sc = numpy.empty( n, float )
            fc = numpy.empty( n, float )
            apply_lut_parallel( pks['s_raw'], pks['f_raw'], sc, fc, self.dx, self.dy )
            pks['sc'] = sc
            pks['fc'] = fc
        else:
            si = numpy.round(pks['s_raw']).astype(int)
            fi = numpy.round(pks['f_raw']).astype(int)
            pks['fc'] = self.dx[ si, fi ] + pks['f_raw']
            pks['sc'] = self.dy[ si, fi ] + pks['s_raw']
        return pks
    
    def pixel_lut(self):
        """ returns (slow, fast) pixel postions of an image """
        s = self.dx.shape
        i, j = numpy.mgrid[ 0:s[0], 0:s[1] ]
        return self.dy + j, self.dx + i


def correct_cf_with_spline(cf, spline_file):
    """Creates a correctorclass from the spline file
       Corrects the columnfile with the spline file
       Returns the corrected columnfile"""
    corrector = correctorclass(spline_file)
    corrector.correct_px_lut(cf)
    return cf

def correct_cf_with_dxdyfiles(cf, dxfile, dyfile):
    """Corrects the columnfile with the dx/dy file
       Returns the corrected columnfile"""
    es = eiger_spatial( dxfile, dyfile )
    pkin = { 's_raw': cf['s_raw'], 'f_raw': cf['f_raw'] }
    pkout = es( pkin )
    cf.addcolumn( pkout['sc'], 'sc' )
    cf.addcolumn( pkout['fc'], 'fc' )
    return cf

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
