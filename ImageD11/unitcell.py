


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


from Numeric import *
from LinearAlgebra import inverse
import math

def radians(x):
    return x*math.pi/180.

def degrees(x):
   return x*180./math.pi



def cross(a,b):
   """
   a x b has length |a||b|sin(theta)
   """
   return array([ a[1]*b[2]-a[2]*b[1] ,a[2]*b[0]-b[2]*a[0], a[0]*b[1]-b[0]*a[1] ],Float)

def unit(a):
   """
   Normalise vector a to unit length
   """
   return a/sqrt(dot(a,a))
     



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

outif = {
    "P" : P ,
    "A" : I ,
    "B" : B ,
    "C" : C ,
    "I" : I ,
    "F" : F }

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
        optionally a symmetry, one of "P","A","B","C","I","F"
        """
        self.lattice_parameters=array(lattice_parameters)
        if self.lattice_parameters.shape[0]!=6:
            raise "You must supply 6 lattice parameters, a,b,c,alpha,beta,gamma"
        self.symmetry=symmetry
        if self.symmetry not in ["P","A","B","C","I","F"]:
            raise "Your symmetry "+self.symmetry+" was not recognised"
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
        self.g = array( [[ a*a    ,  a*b*cg, a*c*cb ],
                         [ a*b*cg ,  b*b   , b*c*ca ],
                         [ a*c*cb ,  b*c*ca, c*c    ]],Float)
        if verbose==1: print "Metric tensor\n",self.g
        try:
            self.gi = inverse(self.g)
        except:
            raise "Unit cell was degenerate, could not determine reciprocal metric tensor"
        if verbose==1: print "Reciprocal Metric tensor\n",self.gi
        self.as=sqrt(self.gi[0,0])
        self.bs=sqrt(self.gi[1,1])
        self.cs=sqrt(self.gi[2,2])
        
        self.alphas=degrees(math.acos(self.gi[1,2]/self.bs/self.cs))
        self.betas =degrees(math.acos(self.gi[0,2]/self.as/self.cs))
        self.gammas=degrees(math.acos(self.gi[0,1]/self.as/self.bs))
        if verbose==1: print "Reciprocal cell"
        if verbose==1: print self.as,self.bs,self.cs,self.alphas,self.betas,self.gammas
        # Equation 3 from Busing and Levy  
        self.B = array ( [ [ self.as , self.bs*cos(radians(self.gammas)) , self.cs*cos(radians(self.betas)) ] ,
                           [       0 , self.bs*sin(radians(self.gammas)) , -self.cs*sin(radians(self.betas))*ca ],
                           [       0 ,                       0  ,      1./c ] ] , Float)
        if verbose==1: print self.B
        if verbose==1: print matrixmultiply(transpose(self.B),self.B)-self.gi # this should be zero
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
        g1 = dot(h1,matrixmultiply(self.gi,h1))
        g2 = dot(h2,matrixmultiply(self.gi,h2))
        g12= dot(h1,matrixmultiply(self.gi,h2))
        costheta = g12/sqrt(g1*g2)
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
        """
        if dsmax == self.limit and self.peaks!=None:
           return self.peaks
        h=k=0
        l=1 # skip 0,0,0
        hs=ks=ls=1
        b=0
        peaks=[]
        while abs(h)<20: # H
            while abs(k)<20: # K
                while abs(l)<20: #L
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
        return math.sqrt(dot(h,matrixmultiply(self.gi,h))) # 1/d or d*


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

      costheta = dot(g1,g2)/sqrt(dot(g2,g2))/sqrt(dot(g1,g1))
      if verbose==1: print "observed costheta",costheta
      best=5.
      for ha in self.ringhkls[self.ringds[ring1]]:
          for hb in self.ringhkls[self.ringds[ring2]]:
             ca=self.anglehkls(ha,hb)
             if abs(ca[1]-costheta) < best and ha!=hb:
                   h1=ha
                   h2=hb
                   best=ca[1]
      if verbose==1:
         print "Assigning h1",h1,g1,self.ds(h1),sqrt(dot(g1,g1)),self.ds(h1)-sqrt(dot(g1,g1))
         print "Assigning h2",h2,g2,self.ds(h2),sqrt(dot(g2,g2)),self.ds(h1)-sqrt(dot(g1,g1))
         print "Cos angle calc",self.anglehkls(h1,h2),"obs",costheta
      h1c=matrixmultiply(self.B,h1)
      h2c=matrixmultiply(self.B,h2)
      t1c=unit(h1c)
      t3c=unit(cross(h1c,h2c))
      t2c=unit(cross(h1c,t3c))
      t1g=unit(g1)
      t3g=unit(cross(g1,g2))
      t2g=unit(cross(g1,t3g))
      T_g = transpose(array([t1g,t2g,t3g]))  # Array are stored by rows and 
      T_c = transpose(array([t1c,t2c,t3c]))  # these are columns
      U=matrixmultiply(T_g , inverse(T_c))
      UB=matrixmultiply(U,self.B)
      UBI=inverse(UB)
      if verbose==1:
         print "UBI"
         print UBI
         print "Grain gi"
         print matrixmultiply(transpose(UB),UB)
         print "Cell gi"
         print self.gi
         h=matrixmultiply(UBI,g1)
         print "(%9.3f, %9.3f, %9.3f)"%(h[0],h[1],h[2])
         h=matrixmultiply(UBI,g2)
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


        
        
        
