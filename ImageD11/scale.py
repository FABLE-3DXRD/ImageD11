
from __future__ import print_function

## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py




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
Class for scaling images with respect to each other due to moco
style instabilities.

Depends on there being a large background contribution to work with


Mathematically:
   y  =  a * x + b
   dy/da = x
   dy/db = 1
   Least squares problem with:
   matrix = ( dy/da*dy/da , dy/da*dy/db)
            ( dy/da*dy/db , dy/db*dy/db)
   rhs    = ( dy/da*y )
            ( dy/db*y )

Has the option to use only pixels in the image to scale to which
are above a threshold (eg the central circle on the bruker)
"""


import numpy, fabio

class scale:
    def __init__( self, im1, threshold = None):
        """
        Determines scale and offset values for images
        with respect to each other
        im1 = a * im2 + b
        returns a, b
        """
        lsqmat = numpy.zeros((2, 2), float)
        dyda   = numpy.ravel(im1).astype(float)
        self.threshold = threshold
        if threshold is None:
            self.indices = None
            self.notindices = None
        if threshold is not None:
            self.indices = numpy.compress(dyda > threshold,
                                            numpy.arange(dyda.shape[0]))
            self.notindices = numpy.compress(dyda <= threshold,
                                            numpy.arange(dyda.shape[0]))
            assert self.indices.shape[0] + self.notindices.shape[0] == \
                   dyda.shape[0], 'problem with threshold'
            dyda = numpy.take(dyda, self.indices)
        lsqmat[0, 0] = numpy.sum(dyda*dyda)
        lsqmat[1, 0] = lsqmat[0, 1] = numpy.sum(dyda)
        lsqmat[1, 1] = dyda.shape[0]
        self.dyda = dyda
        try:
            self.inverse = numpy.linalg.inv(lsqmat)
        except:
            print(lsqmat)
            raise
        

    def scaleimage(self, im2):
        """
        Return a copy of the image scaled to match the class
        """
        grad, off = self.scale(im2)
        new = im2/grad - off/grad
        new = numpy.where(new<0, 0, new)
        if self.notindices is None: 
            return new
        else:
            numpy.put(new, self.notindices, 0. )
            return new
    
    def scale(self, im2):
        """
        Fill out RHS and solve
        returns the scale to apply to the image stored in the class

        You probably want the scale to apply to the image you supply
        ...use scale image for that
        """
        if self.indices is None:
            rhs0 = numpy.sum(self.dyda * numpy.ravel(im2).astype(float))
            rhs1 = numpy.sum(numpy.ravel(im2).astype(float))
            ans = numpy.dot(self.inverse, [rhs0, rhs1])
            return ans[0], ans[1]
        else:
            usedata = numpy.take(numpy.ravel(im2) , self.indices)
            rhs0 = numpy.sum(self.dyda * usedata.astype(float))
            rhs1 = numpy.sum(usedata.astype(float))
            ans = numpy.dot(self.inverse, [rhs0, rhs1])
            return ans[0], ans[1]
            

def scaleseries( target, stem, first, last, 
                 thresh = None,
                 writeim = None ):
    """
    Scale a series of [bruker] images to the target
    TODO - make it work with fabio file series
    """
    # d0 = numpy.ravel(target.data.astype(float))
    scaler = scale(target.data, thresh)
    print("# Scaling with respect to:", sys.argv[1])
    if thresh is not None:
        print("# Using", scaler.indices.shape[0], "pixels above threshold")
    else:
        print("# Using all pixels")
    print("# Number Filename multiplier(t=" + str(thresh) + \
          ") offset multiplier(all) offset")
    if writeim is None:
        # we only look to see
        for i in range(first, last+1):
            name = "%s.%04d" % (stem, i)
            secondimage = fabio.open(name)
            a, b = scaler.scale(secondimage.data)
            print(i, name , a, b, end=' ')
    else: # we correct the image
        for i in range(first, last+1):
            name = "%s.%04d" % (stem, i)
            newname = "cor_%s.%04d" % (stem.split("/")[-1], i)
            secondimage = fabio.open(name)
            newdata = scaler.scaleimage(secondimage.data)
            # write out the file
            secondimage.data = newdata
            secondimage.write( newname )
            print(name, " -> ", newname)
            sys.stdout.flush()


if __name__ == "__main__":

    import sys
    FIRSTIMAGE = fabio.open(sys.argv[1])
    STEM = sys.argv[2]
    FIRST = int(sys.argv[3])
    LAST  = int(sys.argv[4])
    try:
        THRES = float(sys.argv[5])
    except:
        THRES = None
    try:
        WRIT = sys.argv[6]
    except:
        WRIT = None

    scaleseries( FIRSTIMAGE, STEM, FIRST, LAST, THRES, WRIT )
