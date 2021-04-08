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
Class for determining the linearity correction for a frelon4m

Take a quadratic as a first guess

Mathematically:
   I  =  Io + a * Io * Io = Io ( 1 + a * Io )

   read series of images.

   try to scale them to one another.

   a is a constant for all images
   each image has it's own scale factor (exposure time)

   fit scale factor for each image
   fit a single 'a'
   fo

"""


import numpy as np
from ImageD11 import opendata

class linearity:
    def __init__( self , a=1e-6):
        """
        self.a = first guess at non linear term
        """
        self.a = 1e-5
        self.rawimages = []
        self.scaledimages = []
        self.scale_factors = []
        self.names = []
        self.average = None

    def anotherimage(self, ar, scale=1.0, name=None):
        """
        add another image to the class
        ar = data
        scale = guessed scale factor
        """
        t = ar.astype(float)
        self.rawimages.append(t)
        self.scaledimages.append(scale * t * ( 1.0 + t * self.a ))
        self.scale_factors.append(scale)
        if name is None:
            self.names.append(str(len(self.names)))
        else:
            self.names.append(name)



    def makeaverage(self):
        """
        try to get a better idea of the scale factors for each
        image
        """
        logging.info("Compute average")
        self.average = np.zeros(self.scaledimages[0].shape,
                                     float)
        for im in self.scaledimages:
            self.average = self.average + im
        self.average = self.average / len(self.scaledimages)


    def fittoaverage(self):
        """
        Try to fit the scale and a for each image to
        match the average

        chi2 = avg - calc
             = avg - scale * t * (1.0 + t * self.a )

        dc/ds = t * (1.0 + t * self.a ) = calc / scale
        dc/da = scale * t * t
        """
        logging.info("fittoaverage, self.a= "+str(self.a))
        logging.info("fittoaverage, [shift_s shift_a]")

        for i in range(len(self.rawimages)):
            t = self.rawimages[i]
            s = self.scale_factors[i]
            c = self.scaledimages[i]
            logging.debug("t.shape %s s %f %s c.shape"%(
                str(t.shape), s, str(c.shape)))
            dcds = c / s
            dcda = s * t * t
            dy = self.average - c
            logging.debug("dy.shape %s average.shape %s c.shape %s"%(
                str(dy.shape),str(self.average.shape), str(c.shape)))
            import sys
            sys.stdout.flush()
            shifts = lsq( dy , [dcds, dcda] )
            logging.info(self.names[i]+" "+str(shifts)+" "+str(s))





def lsq(diff, gradients):
    """
    difference
    gradients
    """
    nvar = len(gradients)
    lsqmat = np.zeros((nvar,nvar),float)
    rhs = np.zeros((nvar),float)
    for i in range(nvar):
        logging.debug(" lsq shapes: %d %d"%(gradients[i].ravel().shape[0],
                                            diff.ravel().shape[0]))
        try:
            rhs[i] = np.dot(gradients[i].ravel(), diff.ravel())
        except:
            print gradients[i].ravel().shape
            print diff.ravel().shape
            raise
        for j in range(i):
            lsqmat[i,j] = lsqmat[j,i] = \
                          np.dot(gradients[i].ravel(),
                                      gradients[j].ravel())
    inverse = np.linalg.inv(lsqmat)
    shifts = np.dot(inverse, rhs)
    return shifts


if __name__=="__main__":
    import sys, time, glob, logging
    log = logging.getLogger()
    log.setLevel(logging.INFO)
#    testscaleimage()
#    sys.exit()


    stem = sys.argv[1]
    first = int(sys.argv[2])
    last  = int(sys.argv[3])
    try:
        step = int(sys.argv[4])
    except:
        step = 1

    obj = linearity()

    for i in range( first , last + 1 , step ):
        name = opendata.makename(stem, i, ".edf")
        dataobj = opendata.opendata( name )
        try:
            logging.debug(str(dataobj.header["Integration"]))
            integration = float(dataobj.header["Integration"])
            scaling = 1.0 / integration
        except:
            raise
        logging.info("%s %f"%(name,integration))
        obj.anotherimage(dataobj.data, scaling, name)

    obj.makeaverage()
    obj.fittoaverage()


