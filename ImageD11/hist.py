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
Wrapper to ImageD11 histogram generators (c functions)

Eventually to be used to generating sensible image color maps
for plotedf

Also try to do a fast median filter based on pixel histograms
"""

from ImageD11 import _hist
import numpy.oldnumeric as Numeric
import logging

def hist(data,verbose=0):
    """
    Wrapper to ImageD11 histogram generators (c functions)
    returns histogram if data is a 2D Numeric.UInt16 array
    """
    if data.dtype.char == Numeric.UInt16:
        h = Numeric.zeros(pow(2,16)-1,Numeric.Int)
        _hist._hist(data,h,verbose)
        return h
    raise Exception("Sorry, not implemented")


def test_dvhist():
    import _hist
    gv= Numeric.zeros((10,3),Numeric.Float)
    h = Numeric.zeros((10,10,10),Numeric.Int)
    _hist._dvhist(gv,h,-1,1,-1,1,-1,1,1)

    logging.info("test_dvhist maximum.reduce(h)"% Numeric.maximum.reduce(h))


if __name__=="__main__":
    import sys
    from ImageD11 import opendata
    import time
    import logging
    start = time.time()
    d=opendata.openedf(sys.argv[1]).data
    logging.info("Opening time %f"% (time.time() - start))
    start=time.time()
    l1 = d.shape[0]/4
    l2 = d.shape[1]/4
    h1 = 3*d.shape[0]/4
    h2 = 3*d.shape[1]/4

    if len(sys.argv)>2:
        h2=hist(d[l1:h1,l2:h2],1)
    else:
        h2=hist(d[l1:h1,l2:h2])

    if len(sys.argv)>2:
        h1=hist(d,1)
    else:
        h1=hist(d)

    logging.info("Histogram time %f"% (time.time()-start))
    from matplotlib.pylab import plot,show
    plot(h1)
    plot(h2)
    show()
