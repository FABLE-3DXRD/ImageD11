


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


import Numeric, LinearAlgebra
from ImageD11 import opendata

class scale:
    def __init__( self, im1, threshold = None):
        """
        Determines scale and offset values for images
        with respect to each other
        im1 = a * im2 + b
        returns a, b
        """
        LsqMat = Numeric.zeros((2,2),Numeric.Float)
        dyda   = Numeric.ravel(im1).astype(Numeric.Float)
        self.threshold = threshold
        if threshold is None:
            self.indices = None
            self.notindices = None
        if threshold is not None:
            self.indices = Numeric.compress(dyda > threshold,
                                            Numeric.arange(dyda.shape[0]))
            self.notindices = Numeric.compress(dyda <= threshold,
                                            Numeric.arange(dyda.shape[0]))
            assert self.indices.shape[0] + self.notindices.shape[0] == \
                   dyda.shape[0], 'problem with threshold'
            dyda = Numeric.take(dyda,self.indices)
        LsqMat[0,0] = Numeric.sum(dyda*dyda)
        LsqMat[1,0] = LsqMat[0,1] = Numeric.sum(dyda)
        LsqMat[1,1] = dyda.shape[0]
        self.dyda = dyda
        try:
            self.Inverse = LinearAlgebra.inverse(LsqMat)
        except:
            print LsqMat
            raise
        

    def scaleimage(self,im2):
        """
        Return a copy of the image scaled to match the class
        """
        a,b = self.scale(im2)
        new = im2/a - b/a
        new = Numeric.where(new<0, 0, new)
        if self.notindices is None: 
            return new
        else:
            Numeric.put(new, self.notindices, 0. )
            return new
    
    def scale(self, im2):
        """
        Fill out RHS and solve
        returns the scale to apply to the image stored in the class

        You probably want the scale to apply to the image you supply
        ...use scale image for that
        """
        if self.indices is None:
            rhs0 = Numeric.sum(self.dyda * Numeric.ravel(im2).astype(Numeric.Float))
            rhs1 = Numeric.sum(Numeric.ravel(im2).astype(Numeric.Float))
            ans = Numeric.matrixmultiply(self.Inverse,[rhs0,rhs1])
            return ans[0],ans[1]
        else:
            usedata = Numeric.take(Numeric.ravel(im2),self.indices)
            rhs0 = Numeric.sum(self.dyda * usedata.astype(Numeric.Float))
            rhs1 = Numeric.sum(usedata.astype(Numeric.Float))
            ans = Numeric.matrixmultiply(self.Inverse,[rhs0,rhs1])
            return ans[0],ans[1]
            

def testscaleimage():
    im1 = Numeric.ones((10,10),Numeric.Float)
    im1[2:8,2:8] *= 2
    im1[4:6,4:6] *= 2
    im2 = im1 * 2. + 3.
    o = scale(im1)
    a , b = o.scale(im2)
    assert abs(a-2.)< 1e-6 and abs(b-3.) < 1e-6
    scaled = o.scaleimage(im2)
    diff = Numeric.ravel(scaled - im1)
    assert Numeric.sum(diff*diff) < 1e-6

if __name__=="__main__":
    import sys, time, glob
#    testscaleimage()
#    sys.exit()

    firstimage = opendata.opendata(sys.argv[1])
    stem = sys.argv[2]
    first = int(sys.argv[3])
    last  = int(sys.argv[4])
    try:
        threshold = float(sys.argv[5])
    except:
        threshold = None
    try:
        write = sys.argv[6]
    except:
        write = None
    d0 = Numeric.ravel(firstimage.data.astype(Numeric.Float))
    scaler = scale(firstimage.data,threshold)
    print "# Scaling with respect to:",sys.argv[1]
    if threshold is not None:
        print "# Using",scaler.indices.shape[0],"pixels above threshold"
    else:
        print "# Using all pixels"
    print "# Number Filename multiplier(t="+str(threshold)+") offset multiplier(all) offset"
    if write is None:
        # we only look to see
        for i in range(first,last+1):
            start = time.clock()
            name = "%s.%04d"%(stem,i)
            secondimage = opendata.opendata(name)
            a, b = scaler.scale(secondimage.data)
            print i, name , a, b,
    else: # we correct the image
        from ImageD11 import data
        for i in range(first,last+1):
            name = "%s.%04d"%(stem,i)
            newname = "cor_%s.%04d"%(stem.split("/")[-1],i)
            secondimage = opendata.opendata(name)
            newdata = scaler.scaleimage(secondimage.data)
            # write out the file
            dataobj = data.data(newdata, secondimage.header)
            opendata.writedata(newname,dataobj)
            print name," -> ",newname
            sys.stdout.flush()
