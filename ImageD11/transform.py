


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


from Numeric import *
import math as m

def degrees(x):
   """Convenience function"""
   return x*180.0/math.pi

def radians(x):
   """Convenience function"""
   return x*math.pi/180.0


def compute_tth_eta(peaks,yc,ys,ty,zc,zs,tz,dist):
   """
   Peaks is a 2 d array of x,y
   yc is the centre in y
   ys is the y pixel size
   ty is the tilt around y
   zc is the centre in z
   zs is the z pixel size
   tz is the tilt around z
   dist is the sample - detector distance
   """
   r1 = array( [ [ m.cos(tz) , m.sin(tz) , 0 ],
                 [ -m.sin(tz), m.cos(tz) , 0 ],
                 [         0 ,         0 , 1 ]],Float)
   r2 = array( [ [ m.cos(ty) , 0 , m.sin(ty) ],
                 [       0   , 1 , 0         ],
                 [-m.sin(ty) , 0 , m.cos(ty) ]],Float)
   r2r1=matrixmultiply(r1,r2)
   vec =array( [     zeros(peaks.shape[1]) , # place detector at zero, sample at -dist
                   (peaks[0,:]-yc)*ys      ,            # x in search
                   (peaks[1,:]-zc)*zs ]    , Float)     # y in search
   #print vec.shape
   rotvec=matrixmultiply(r2r1,vec)
   magrotvec=sqrt(sum(rotvec*rotvec,0))
   #print magrotvec.shape
   eta=degrees(arctan2(rotvec[2,:],rotvec[1,:])) # bugger this one for tilts
   # cosine rule a2 = b2+c2-2bccos(A)
   # a is distance from (0,0,0) to rotvec => magrotvec
   # b is distance from sample to detector
   # c is distance from sample to pixel
   a=magrotvec  # 1D
   b=dist            # 0D
   c=rotvec     # 3D
   #print c.shape
   c[0,:]=c[0,:]-dist# 3D
   c=sqrt(sum(c*c,0))# 1D
   #print a.shape,c.shape
   #   print a.shape,c.shape
   costwotheta = (b*b + c*c - a*a)/2/b/c
   #   print self.costwotheta.shape
   twothetarad=arccos(costwotheta)
   twotheta=degrees(twothetarad)
   return twotheta, eta



def compute_g_vectors(tth,eta,omega,wavelength, wedge = 0.0):
   """
   Generates spot positions in reciprocal space from twotheta, wavelength, omega and eta
   Assumes single axis vertical
   ... unless a wedge angle is specified
   """
   if wedge != 0.:
      print "Using a wedge angle of ",wedge
   tth=radians(tth)
   eta=radians(eta)
   om =radians(omega)
   # Compute k vector - the scattering in the laboratory
   c=cos(tth/2) # cos theta
   s=sin(tth/2) # sin theta
   ds=2*s/wavelength
   k=zeros((3,tth.shape[0]),Float)
   # x - along incident beam
   k[0,:] = -ds*s
   # y - towards door
   k[1,:] = -ds*c*sin(eta) # plus or minus sin or cos
   # z - towards roof
   k[2,:] =  ds*c*cos(eta)
   # G-vectors - rotate k onto the crystal axes
   g=zeros((3,ds.shape[0]),Float)
   t=zeros((3,ds.shape[0]),Float)
   # 
   # g =  R . W . k where:
   # R = ( cos(omega) , sin(omega), 0 )
   #     (-sin(omega) , cos(omega), 0 )
   #     (         0  ,         0 , 1 )
   #
   # W = ( cos(wedge) ,  0  ,  sin(wedge) )
   #     (         0  ,  1  ,          0  )
   #     (-sin(wedge) ,  0  ,  cos(wedge) )
   #
   if wedge != 0.0:
      c = cos(radians(wedge))
      s = sin(radians(wedge))
      t[0,:]= c * k[0,:]           + s * k[2,:]
      t[1,:]=            k[1,:]
      t[2,:]=-s * k[0,:]           + c * k[2,:]
      k=t
   g[0,:] = cos(om)*k[0,:]+sin(om)*k[1,:]
   g[1,:] =-sin(om)*k[0,:]+cos(om)*k[1,:]
   g[2,:] =                               k[2,:]
   return g


def uncompute_g_vectors(gv,wavelength, wedge=0.0):
   """
   Given g-vectors compute tth,eta,omega
   assert uncompute_g_vectors(compute_g_vector(tth,eta,omega))==tth,eta,omega
   """
   # k=array(g.shape,Float)
   # k[2,:]=g[2,:] # Z component in laboratory
   # Length gives two-theta
   #
   # Deal with wedge ... we had           g = W . R . k
   #                    now get     W-1 . g =     R . k
   # and proceed as before for R
   g=zeros(gv.shape,Float)
   c = cos(radians(wedge))
   s = sin(radians(wedge))
   g[0,:]= c * gv[0,:]           - s * gv[2,:]
   g[1,:]=            gv[1,:]
   g[2,:]= s * gv[0,:]           + c * gv[2,:]
   gv = g
   ds = sqrt(sum(gv*gv,0))
   s = ds*wavelength/2 # sin theta
   # Two theta from Bragg's law
   try:
      th = arcsin(s)
   except:
      print "Problem getting theta from sin(theta)"
      print "Maximum and minimum in sin(theta)=",maximum.reduce(s),minimum.reduce(s)
   c = cos(th) # cos theta
   # Two theta in degrees
   k=zeros(gv.shape,Float)
   k[2,:]=gv[2,:] # Last line above  == ds*c*cos(eta)
   k[0,:] = -ds*s # From above
   # Should check if this is true or not?
   bottom=ds*c
   try:
      coseta = gv[2,:]/bottom
   except: # Either a (0,0,0) peak or cos(theta) = 0 -> 90 degrees
      ind=compress(bottom<1e-6,range(ds.shape[0]))
      put(bottom,ind,1)
      coseta = gv[2,:]/bottom
      put(coseta,ind,1.)
   #print "max and min in coseta",maximum.reduce(coseta),minimum.reduce(coseta)
   coseta=clip(coseta,-1,1)
   tth = degrees(2*th)
   eta1 = degrees(arccos(coseta))
   eta2 = -eta1
   k[1,:] = -ds*c*sin(radians(eta1)) # Known also
   # Know k and g for each peak - need to solve for omega
   omega_crystal = arctan2(gv[0,:],gv[1,:])
   omega_laboratory = arctan2(k[0,:],k[1,:])
   omega1 = degrees(  omega_crystal - omega_laboratory )
   k[1,:] = -ds*c*sin(radians(eta2)) # Known also
   # Know k and g for each peak - need to solve for omega
   omega_crystal = arctan2(gv[0,:],gv[1,:])
   omega_laboratory = arctan2(k[0,:],k[1,:])
   omega2 = degrees(  omega_crystal - omega_laboratory )
   return tth,[eta1,eta2],[omega1,omega2]       

def uncompute_one_g_vector(gv,wavelength, wedge=0.0):
   """
   Given g-vectors compute tth,eta,omega
   assert uncompute_g_vectors(compute_g_vector(tth,eta,omega))==tth,eta,omega
   """
   t,e,o=uncompute_g_vectors(transpose(array([gv,gv])),wavelength,wedge=wedge)
   return t[0],[e[0][0],e[1][0]],[o[0][0],o[1][0]]       

   
   
   
