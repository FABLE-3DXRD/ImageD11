


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

import Numeric as n
import bisplev, math

class correctorclass:
    """
    Applies a spatial distortion to a peak position using a fit2d splinefile
    """
    def __init__(self,splinefile):
        """
        Argument is the name of a fit2d spline file
        """
        self.splinefile=splinefile
        self.readfit2dspline(splinefile)
        self.tolerance=1e-5


    def correct(self,x,y):
        """
        Transform x,y in raw image coordinates into x,y of an
        idealised image. Returns a tuple (x,y), expects a
        pair of floats as arguments
        """
        ynew=bisplev.bisplev(y,x,self.tck1)+y
        xnew=bisplev.bisplev(y,x,self.tck2)+x
        return xnew, ynew

    def distort(self,xnew,ynew):
        """
        Distort a pair of points xnew, ynew to find where they
        would be in a raw image

        Iterative algorithm...
        """
        yold=ynew-bisplev.bisplev(ynew,xnew,self.tck1)
        xold=xnew-bisplev.bisplev(ynew,xnew,self.tck2)
        # First guess, assumes distortion is constant
        y=ynew-bisplev.bisplev(yold,xold,self.tck1)
        x=xnew-bisplev.bisplev(yold,xold,self.tck2)
        # Second guess should be better
        error=math.sqrt((x-xold)*(x-xold)+(y-yold)*(y-yold))
        ntries=0
        while error > self.tolerance:
            ntries=ntries+1
            xold=x
            yold=y
            y=ynew-bisplev.bisplev(yold,xold,self.tck1)
            x=xnew-bisplev.bisplev(yold,xold,self.tck2)
            error=math.sqrt((x-xold)*(x-xold)+(y-yold)*(y-yold))
            # print error,xold,x,yold,y
            if ntries == 10:
                raise "Error getting the inverse spline to converge"
        return x,y

    def test(self,x,y):
        """
        Checks that the correct and distort functions are indeed
        inversely related to each other
        """
        xnew,ynew = self.correct(x,y)
        xold,yold = self.distort(xnew,ynew)
        error = math.sqrt( (x-xold)*(x-xold) + (y-yold)*(y-yold))
        if error > self.tolerance:
            print "Failed!"
            raise "Problem in correctorclass"




    def readfit2dfloats(self,fp,n,debug=0):
        """
        Interprets a 5E14.7 formatted fortran line
        """
        vals=[]
        j=0
        while j < n:
            i=0
            line=fp.readline()
            while i < 5*14:
                if(debug):print line
                if(debug):print i,line[i:i+14]
                vals.append(float(line[i:i+14]) )
                j=j+1
                i=i+14
                if j == n: break
            i=i+1
        return vals


            # read the fit2d array into a tck tuple
    def readfit2dspline(self,name):
        """
        Reads a fit2d spline file into a scipy/fitpack tuple, tck
        A fairly long and dull routine...
        """
        kx=3
        ky=3
        fp=open(name,"r")
        line=fp.readline() # SPATIAL DISTORTION SPLINE INTERPOLATION COEFFICIENTS
        if line[:7] != "SPATIAL":
            raise SyntaxError, name+": file does not seem to be a fit2d spline file"
        line=fp.readline() # BLANK LINE
        line=fp.readline() # VALID REGION
        line=fp.readline() # the actual valid region, assume xmin,ymin,xmax,ymax
        vals=line.split()
        self.xmin=float(vals[0])
        self.ymin=float(vals[1])
        self.xmax=float(vals[3])
        self.ymax=float(vals[3])
        line=fp.readline() # BLANK
        line=fp.readline() # GRID SPACING, X-PIXEL SIZE, Y-PIXEL SIZE
        line=fp.readline()
        vals=line.split()
        self.gridspacing=float(vals[0])
        self.xsize=float(vals[1])
        self.ysize=float(vals[2])
        line=fp.readline() # BLANK
        line=fp.readline() # X-DISTORTION
        line=fp.readline() # two integers nx1,ny1
        vals=line.split()
        nx1=int(vals[0])
        ny1=int(vals[1])
        # Now follow fit2d formatted line 5E14.7
        tx1=n.array(self.readfit2dfloats(fp,nx1),n.Float32)
        ty1=n.array(self.readfit2dfloats(fp,ny1),n.Float32)
        c1 =n.array(self.readfit2dfloats(fp,(nx1-4)*(ny1-4)),n.Float32)
        line=fp.readline() #BLANK
        line=fp.readline() # Y-DISTORTION
        line=fp.readline() # two integers nx2, ny2
        vals=line.split()
        nx2=int(vals[0])
        ny2=int(vals[1])
        tx2=n.array(self.readfit2dfloats(fp,nx2),n.Float32)
        ty2=n.array(self.readfit2dfloats(fp,ny2),n.Float32)
        c2 =n.array(self.readfit2dfloats(fp,(nx2-4)*(ny2-4)),n.Float32)
        fp.close()
        self.tck1=(tx1,ty1,c1,kx,ky)
        self.tck2=(tx2,ty2,c2,kx,ky)

"""
http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/OWENS/LECT5/node5.html

Various interpolation schemes can be used. A common one is bilinear interpolation, given by

v(x,y) = c1x + c2y + c3xy + c4,

where v(x,y) is the grey value at position (x,y).
Thus we have four coefficients to solve for. We use the known grey values of the 4 pixels
surrounding the `come from' location to solve for the coefficients.

We need to solve the equation

v1   ( x1 y1 x1y1 1 ) c1
v2 = ( x2 y2 x2y2 1 ) c2
v3   ( x3 y3 x3y3 1 ) c3
v4   ( x4 y4 x4y4 1 ) c4


or, in short,
[V] = [M][C],

which implies
[C] = [M]-1[V].

This has to be done for every pixel location in the output image and is thus a lot of computation!
Alternatively one could simply use the integer pixel position closest to the `come from location'.
This is adequate for most cases.

"""

def unwarpimage(image, xpositions, ypositions):
    """
    xpositions/ypositions are floats giving pixel co-ords of the input image.

    We need the positions of the pixels in the output image, on the input image.

    Hence, for now,
    """
    pass

if __name__=="__main__":
    import sys,time
    from Numeric import *
    start=time.time()
    splinefile=sys.argv[1]
    infile=sys.argv[2]
    outfile=sys.argv[3]
    out=open(outfile,"w")
    npks=0
    tck1,tck2,xmin,xmax,ymin,ymax=readfit2dspline(splinefile)
    print "Time to read fit2d splinefile %f/s"%(time.time()-start)
    start=time.time()
    for line in open(infile,"r").readlines():
        try:
            vals = [float(v) for v in line.split()]
            if vals[0]>3:
                x=vals[2]
                y=vals[3]
                out.write(line[:-1])
                ynew=bisplev.bisplev(y,x,tck1)+y
                xnew=bisplev.bisplev(y,x,tck2)+x
                out.write(" ---> %f %f\n"%(xnew,ynew))
                npks=npks+1
        except ValueError:
            pass
        except:
            raise

    out.close()
    print "That took %f/s for %d peaks"%(time.time()-start,npks)
