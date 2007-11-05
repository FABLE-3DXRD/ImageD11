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

#
#
#


import numpy.oldnumeric as Numeric
from numpy.oldnumeric.linear_algebra import inverse
import math

import logging

def radians(x):
    return x*math.pi/180.

def degrees(x):
    return x*180./math.pi



def cross(a,b):
    """
    a x b has length |a||b|sin(theta)
    """
    return Numeric.array([ a[1]*b[2]-a[2]*b[1] ,a[2]*b[0]-b[2]*a[0], a[0]*b[1]-b[0]*a[1] ],Numeric.Float)

def unit(a):
    """
    Normalise vector a to unit length
    """
    try:
        return a/Numeric.sqrt(Numeric.dot(a,a))
    except:
        logging.error("cannot normalise to unit length a=%s"%(str(a)))
        raise




# Systematic absences

def P(h,k,l):
    return False

def A(h,k,l):
    return (k+l)%2 != 0

def B(h,k,l):
    return (h+l)%2 != 0

def C(h,k,l):
    return (h+k)%2 != 0

def I(h,k,l):
    return (h+k+l)%2 != 0

def F(h,k,l):
    return (h+k)%2!=0 or (h+l)%2!=0 or (k+l)%2!=0
                                                                                
def R(h,k,l):
    return (-h+k+l)%3 != 0


outif = {
    "P" : P ,
    "A" : I ,
    "B" : B ,
    "C" : C ,
    "I" : I ,
    "F" : F ,
    "R" : R}

def cellfromstring(s):
    items=s.split()
    # print items
    latt = [float(x) for x in items[0:6]]
    try:
        symm = items[6]
    except IndexError:
        symm = 'P'
    return unitcell(latt,symm)

class unitcell:
    # Unit cell stuff
    # Generate a list of peaks from a unit cell
    def __init__(self,lattice_parameters,symmetry="P",verbose=0):
        """
        Unit cell class
        supply a list (tuple etc) of a,b,c,alpha,beta,gamma
        optionally a symmetry, one of "P","A","B","C","I","F","R"
        """
        self.lattice_parameters=Numeric.array(lattice_parameters)
        if self.lattice_parameters.shape[0]!=6:
            raise Exception("You must supply 6 lattice parameters\n"+\
                            "      a,b,c,alpha,beta,gamma")
        self.symmetry=symmetry
        if self.symmetry not in ["P","A","B","C","I","F","R"]:
            raise Exception("Your symmetry "+self.symmetry+\
                            " was not recognised")
        # assigning a function here!
        self.absent=outif[self.symmetry]
        a = self.lattice_parameters[0]
        b = self.lattice_parameters[1]
        c = self.lattice_parameters[2]
        self.alpha=radians(self.lattice_parameters[3])
        ca= math.cos(radians(self.lattice_parameters[3]))
        cb= math.cos(radians(self.lattice_parameters[4]))
        cg= math.cos(radians(self.lattice_parameters[5]))
        if verbose==1: print "Unit cell",self.lattice_parameters
        self.g = Numeric.array( [[ a*a    ,  a*b*cg, a*c*cb ],
                                 [ a*b*cg ,  b*b   , b*c*ca ],
                                 [ a*c*cb ,  b*c*ca, c*c    ]],Numeric.Float)
        if verbose==1: print "Metric tensor\n",self.g
        try:
            self.gi = inverse(self.g)
        except:
            raise Exception("Unit cell was degenerate, could not determine"+\
                    "reciprocal metric tensor")
        if verbose==1: print "Reciprocal Metric tensor\n",self.gi
        self.as=Numeric.sqrt(self.gi[0,0])
        self.bs=Numeric.sqrt(self.gi[1,1])
        self.cs=Numeric.sqrt(self.gi[2,2])

        self.alphas=degrees(math.acos(self.gi[1,2]/self.bs/self.cs))
        self.betas =degrees(math.acos(self.gi[0,2]/self.as/self.cs))
        self.gammas=degrees(math.acos(self.gi[0,1]/self.as/self.bs))
        if verbose==1: print "Reciprocal cell"
        if verbose==1: print self.as,self.bs,self.cs,self.alphas,self.betas,self.gammas
        # Equation 3 from Busing and Levy
        self.B = Numeric.array ( [ [ self.as , self.bs*math.cos(radians(self.gammas)) , self.cs*math.cos(radians(self.betas)) ] ,
                           [       0 , self.bs*math.sin(radians(self.gammas)) , -self.cs*math.sin(radians(self.betas))*ca ],
                           [       0 ,                       0  ,      1./c ] ] , Numeric.Float)
        if verbose==1: print self.B
        if verbose==1: print Numeric.dot(Numeric.transpose(self.B),self.B)-self.gi # this should be zero
        self.hkls=None
        self.peaks=None
        self.limit=0

    def tostring(self):
        """
        Write out a line containing unit cell information
        """
        return "%f %f %f %f %f %f %s"%(self.lattice_parameters[0],
              self.lattice_parameters[1],
              self.lattice_parameters[2],
              self.lattice_parameters[3],
              self.lattice_parameters[4],
              self.lattice_parameters[5],
              self.symmetry)

    def anglehkls(self,h1,h2):
        """
        Compute the angle between reciprocal lattice vectors h1, h2
        """
        g1 = Numeric.dot(h1,Numeric.dot(self.gi,h1))
        g2 = Numeric.dot(h2,Numeric.dot(self.gi,h2))
        g12= Numeric.dot(h1,Numeric.dot(self.gi,h2))
        costheta = g12/math.sqrt(g1*g2)
        try:
            return degrees(math.acos(costheta)),costheta
        except:
            if abs(costheta-1) < 1e-6:
                return 0.,1.0
            if abs(costheta+1) < 1e-6:
                return 180.,-1.0
            print "Error in unit cell class determining angle"
            print "h1",h1,"h2",h2,"Costheta=",costheta
            raise

    def gethkls(self,dsmax):
        """
        Generate hkl list
        Argument dsmax is the d* limit (eg 1/d)
        Default of zero gives only the (000) reflection

        assumes [h|k|l] < 200
        """
        if dsmax == self.limit and self.peaks!=None:
            return self.peaks
        h=k=0
        l=1 # skip 0,0,0
        hs=ks=ls=1
        b=0
        peaks=[]
        while abs(h)<200: # H
            while abs(k)<200: # K
                while abs(l)<200: #L
                    ds=self.ds([h,k,l])
                    if ds < dsmax:
                        if not self.absent(h,k,l):
                            peaks.append([ds,(h,k,l)])
                        else:
                            pass
                        b=0
                    else:
                        if ls==1:
                            ls=-1
                            l=0
                        else:
                            ls=1
                            l=0
                            b=b+1
                            break
                    l=l+ls
                k=k+ks
                # l is always zero here
                if b>1:
                    if ks==1:
                        ks=-1
                        k=-1
                    else:
                        ks=1
                        k=0
                        b=b+1
                        break
            h=h+hs
            if b>3:
                if hs==1:
                    hs=-1
                    h=-1
                else:
                    hs=1
                    h=0
                    break

        peaks.sort()
        #for peak in peaks:
        #   print peak
        self.peaks=peaks
        self.limit=dsmax
        return peaks

    def ds(self,h):
        """ computes 1/d for this hkl = hgh """
        return math.sqrt(Numeric.dot(h,Numeric.dot(self.gi,h))) # 1/d or d*


    def makerings(self,limit,tol=0.001):
        """
        Makes a list of computed powder rings
        The tolerance is the difference in d* to decide
        if two peaks overlap
        """
        self.peaks=self.gethkls(limit) # [ ds, [hkl] ]
        self.ringds=[]   # a list of floats
        self.ringhkls={} # a dict of lists of integer hkl
        # Append first peak
        peak = self.peaks[0]
        self.ringds.append(peak[0])
        self.ringhkls[peak[0]] = [peak[1]]
        for peak in self.peaks[1:]:
            if abs(peak[0] - self.ringds[-1]) < tol:
                self.ringhkls[self.ringds[-1]].append(peak[1])
            else:
                self.ringds.append(peak[0])
                self.ringhkls[self.ringds[-1]]= [peak[1]]


    def orient(self,ring1,g1,ring2,g2,verbose=0):
        """
        Compute an orientation matrix using cell parameters and the indexing
        of two reflections

        Orientation matrix:
        A matrix such that h = A^1 g
        define 3 vectors t1,t2,t3
        t1 is parallel to first peak (unit vector along g1)
        t2 is in the plane of both   (unit vector along g1x(g1xg2))
        t3 is perpendicular to both  (unit vector along g1xg2)
        """
        costheta = Numeric.dot(g1,g2)/math.sqrt(Numeric.dot(g2,g2))/math.sqrt(Numeric.dot(g1,g1))
        if verbose==1: print "observed costheta",costheta
        best=5.
        for ha in self.ringhkls[self.ringds[ring1]]:
            for hb in self.ringhkls[self.ringds[ring2]]:
                ca=self.anglehkls(ha,hb)
                if verbose==1: print "%6.3f %.3f"% ca,ha,hb,best,ha,hb
                if abs(ca[1]-costheta) < best and ha!=hb:
                    h1=ha
                    h2=hb
                    # This was a serious bug!!! Best was assigned to ca[1]!!!
                    best=abs(ca[1]-costheta)
        if verbose==1:
            print "Assigning h1",h1,g1,self.ds(h1),math.sqrt(Numeric.dot(g1,g1)),self.ds(h1)-math.sqrt(Numeric.dot(g1,g1))
            print "Assigning h2",h2,g2,self.ds(h2),math.sqrt(Numeric.dot(g2,g2)),self.ds(h1)-math.sqrt(Numeric.dot(g1,g1))
            print "Cos angle calc",self.anglehkls(h1,h2),"obs",costheta
        try:
            h1c=Numeric.dot(self.B,h1)    
            h2c=Numeric.dot(self.B,h2)
            t1c=unit(h1c)
            t3c=unit(cross(h1c,h2c))
            t2c=unit(cross(h1c,t3c))
            t1g=unit(g1)
            t3g=unit(cross(g1,g2))
            t2g=unit(cross(g1,t3g))
            T_g = Numeric.transpose(Numeric.array([t1g,t2g,t3g]))  # Array are stored by rows and
            T_c = Numeric.transpose(Numeric.array([t1c,t2c,t3c]))  # these are columns
            U=Numeric.dot(T_g , inverse(T_c))
            UB=Numeric.dot(U,self.B)
            UBI=inverse(UB)
        except:
            logging.error("unitcell.orient h1 %s g1 %s h2 %s g2 %s self.B %s"%(
                     str(h1), str(g1), str(h2), str(g2), str(self.B)))
            import traceback
            traceback.print_exc()
        if verbose==1:
            print "UBI"
            print UBI
            print "Grain gi"
            print Numeric.dot(Numeric.transpose(UB),UB)
            print "Cell gi"
            print self.gi
            h=Numeric.dot(UBI,g1)
            print "(%9.3f, %9.3f, %9.3f)"%(h[0],h[1],h[2])
            h=Numeric.dot(UBI,g2)
            print "(%9.3f, %9.3f, %9.3f)"%(h[0],h[1],h[2])
        self.UBI=UBI
        self.UB=UB



if __name__=="__main__":
    import sys,time
    start=time.time()
    cell = unitcell([float(x) for x in sys.argv[1:7]],sys.argv[7])
    p=cell.gethkls(float(sys.argv[8]))
    n=len(p)
    print "Generated",n,"peaks in",time.time()-start,"/s"
    for k in p:
        try:
            print k[1],k[0],1/k[0],k[1]
        except:
            pass
