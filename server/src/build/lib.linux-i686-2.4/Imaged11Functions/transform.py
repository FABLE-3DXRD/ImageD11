


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

from math import pi
import time,math,unitcell,sys


def degrees(x):
    """Convenience function"""
    return x*180.0/pi

def radians(x):
    """Convenience function"""
    return x*pi/180.0


def compute_tth_eta(peaks,
                    yc,ys,ty,
                    zc,zs,tz,
                    dist,
                    detector_orientation=((1,0),(0,1)),
                    crystal_translation = None,
                    omega = None,         #       == phi at chi=90
                    axis_orientation1=0.0, # Wedge == theta on 4circ
                    axis_orientation2=0.0 #       == chi - 90
                    ):
    """
    Peaks is a 2 d array of x,y
    yc is the centre in y
    ys is the y pixel size
    ty is the tilt around y
    zc is the centre in z
    zs is the z pixel size
    tz is the tilt around z
    dist is the sample - detector distance
    detector_orientation is a matrix to apply to peaks arg to get ImageD11 convention
         (( 0, 1),( 1, 0)) for ( y, x)
         ((-1, 0),( 0, 1)) for (-x, y)
         (( 0,-1),(-1, 0)) for (-y,-x)
      etc...
    crystal_translation is the position of the grain giving rise to a diffraction spot
                 in x,y,z ImageD11 co-ordinates
                 x,y with respect to axis of rotation and or beam centre ??
                 z with respect to beam height, z centre
    omega data needed if crystal translations used
    """
    # Matrices for the tilt rotations
    r1 = array( [ [  cos(tz) , sin(tz) , 0 ],
                  [ -sin(tz) , cos(tz) , 0 ],
                  [    0     ,    0    , 1 ]],Float)
    r2 = array( [ [ cos(ty) , 0 , sin(ty) ],
                  [       0 , 1 ,   0     ],
                  [-sin(ty) , 0 , cos(ty) ]],Float)
    r2r1=matrixmultiply(r1,r2)
    # Peak positions in 3D space
    #  - apply detector orientation
    mypeaks = matrixmultiply(array(detector_orientation,Float),peaks)
    vec =array( [ zeros(mypeaks.shape[1])     , # place detector at zero, sample at -dist
                    (mypeaks[0,:]-yc)*ys      ,            # x in search
                    (mypeaks[1,:]-zc)*zs ]    , Float)     # y in search
    #print vec.shape
    # Position of diffraction spots in 3d space after detector tilts is:
    rotvec=matrixmultiply(r2r1,vec)
    # Scattering vectors
    if crystal_translation is None:
        magrotvec=sqrt(sum(rotvec*rotvec,0))
        #print magrotvec.shape

        eta=degrees(arctan2(rotvec[2,:],rotvec[1,:])) # bugger this one for tilts
        # cosine rule a2 = b2+c2-2bccos(A)
        # a is distance from (0,0,0) to rotvec => magrotvec
        # b is distance from sample to detector
        # c is distance from sample to pixel
        a=magrotvec  # 1D
        b=dist       # 0D
        c=rotvec     # 3D
        #print c.shape
        c[0,:]=c[0,:]+dist# 3D
        # print c
        c=sqrt(sum(c*c,0))# 1D
        #print a.shape,c.shape
        #   print a.shape,c.shape
        costwotheta = (b*b + c*c - a*a)/2/b/c
        #print costwotheta
        twothetarad=arccos(costwotheta)
        twotheta=degrees(twothetarad)
        return twotheta, eta
    else:
        # Compute positions of grains
        # expecting tx, ty, tz for each diffraction spot
        origin =  array([ -dist, 0, 0 ], Float)

        #
        # g =  R . W . k
        #  g - is g-vector w.r.t crystal
        #  k is scattering vector in lab
        #  so we want displacement in lab from displacement in sample
        #  shift =  W-1  R-1 crystal_translation
        #
        # R = ( cos(omega) , sin(omega), 0 )
        #     (-sin(omega) , cos(omega), 0 )
        #     (         0  ,         0 , 1 )
        #
        # W = ( cos(wedge) ,  0  ,  sin(wedge) )
        #     (         0  ,  1  ,          0  )
        #     (-sin(wedge) ,  0  ,  cos(wedge) )
        #
        # C = (         1  ,          0  ,       0     ) ??? Use eta0 instead
        #     (         0  ,   cos(chi)  ,  sin(chi)   )  ??? Use eta0 instead
        #     (         0  ,  -sin(chi)  ,  cos(chi)   )  ??? Use eta0 instead
        w=radians(axis_orientation1)
        WI = array( [ [ cos(w),         0, -sin(w)],
                      [      0,         1,       0],
                      [ sin(w),         0,  cos(w)] ] , Float)
        c=radians(axis_orientation2)
        CI = array( [ [      1,          0,       0],
                      [      0,     cos(c), -sin(c)],
                      [      0,     sin(c),  cos(c)] ] , Float)
        if omega.shape[0] != peaks.shape[1]:
            raise Exception("omega and peaks arrays must have same number of peaks")
        eta=zeros(omega.shape[0],Float)
        tth=zeros(omega.shape[0],Float)

        for i in range(omega.shape[0]):
            om = radians(omega[i])
            RI = array( [ [ cos(om), -sin(om), 0],
                          [ sin(om),  cos(om), 0],
                          [       0,        0, 1] ] , Float)
            rotate = matrixmultiply(WI,matrixmultiply(CI,RI))
            #print rotate.shape,origin.shape,crystal_translation.shape
            myorigin = origin + matrixmultiply(rotate , crystal_translation)
            scattering_vectors = s = rotvec[:,i] - myorigin
            #print i,s,
            eta[i]=degrees(arctan2(s[2],s[1]))
            mag_s=sqrt(sum(s*s,0))
            costth = s[0] / mag_s
            #print costth
            tth[i] = degrees(arccos(costth))

        return tth, eta

def compute_g_vectors(tth,eta,omega,wavelength, wedge = 0.0, chi=0.0):
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
    # C = (         1  ,         0  ,      0     )
    #     (         0  ,  cos(chi)  , sin(chi)   )
    #     (         0  , -sin(chi)  , cos(chi)   )
    #
    if wedge != 0.0:
        c = cos(radians(wedge))
        s = sin(radians(wedge))
        t[0,:]= c * k[0,:]           + s * k[2,:]
        t[1,:]=            k[1,:]
        t[2,:]=-s * k[0,:]           + c * k[2,:]
        k=t
    if chi != 0.0:
        c = cos(radians(chi))
        s = sin(radians(chi))
        t[0,:]= k[0,:]  
        t[1,:]=        c * k[1,:]    + s * k[2,:]
        t[2,:]=       -s * k[1,:]    + c * k[2,:]
        k=t
    g[0,:] = cos(om)*k[0,:]+sin(om)*k[1,:]
    g[1,:] =-sin(om)*k[0,:]+cos(om)*k[1,:]
    g[2,:] =                               k[2,:]
    return g


def uncompute_g_vectors(g,wavelength, wedge=0.0, chi=0.0):
    """
    Given g-vectors compute tth,eta,omega
    assert uncompute_g_vectors(compute_g_vector(tth,eta,omega))==tth,eta,omega
    """
    # print "new uncompute using wavelength",wavelength,"and wedge",wedge
    # k = W-1 R-1 g    ...or...      g = R W k
    #
    # |g| = 1 / d
    # sin(theta) = wavelength / 2 d = |g| wavelength / 2
    #
    ds = sqrt(sum(g*g,0))
    s = ds*wavelength/2.0 # sin theta
    #
    #     k0 = projection of scattering vector along x in lab
    #          must satisfy laue condition
    #        =  - sin(theta)/d
    #
    k0 = zeros(g.shape[0],Float)
    k0 = -ds*s
    #
    # this component: k = W-1 R-1 g
    # k = W-1 ( g0 cos(omega) - g1 sin(omega)  + 0  )
    #         ( g0 sin(omega) + g1 cos(omega)  + 0  )
    #         (       0       +      0         + g2 )
    #
    # k0 = cos(wedge)[g0 cos(omega) - g1 sin(omega)] - sin(wedge) g2
    #
    cw = cos(radians(wedge))
    sw = sin(radians(wedge))
    #
    # k0 = cw[g0 cos(omega) - g1 sin(omega)] - sw g2
    # (k0 + sw g2) / cw = g0 cos(omega) - g1 sin(omega)
    #
    lhs = ( k0 + sw * g[2,:] ) / cw
    #
    #      lhs = g0 cos(omega) - g1 sin(omega)
    #      lhs + g1 sin(omega) = g0 cos(omega)
    #      lhs^2 + 2 g1 lhs sin(omega) + g1^2 sin^2(omega) = g0^2 cos^2(omega)
    #      - g0^2 + lhs^2 + 2 g1 lhs sin(omega) + g1^2 sin^2(omega) + g0^2 sin^2(omega)  = 0
    #
    # Finally - solve the quadratic for sin(omega)
    # axx + bx + c = 0
    # x = {- b +/- sqrt(b*b-4ac) } / 2a
    #
    c = lhs * lhs - g[0,:]*g[0,:]
    b = 2. * g[1,:] * lhs
    a = g[1,:]*g[1,:] + g[0,:]*g[0,:]
    t = b*b - 4*a*c
    valid = where(t>0.0 , 1, 0)
    t = t * valid
    another_valid = where(a < 1e-12, 1, 0)
    a = where(another_valid == 1, 1, a)
    sinomega1 = (- b - sqrt(t)) / (a * 2.)   # a can clearly be zero
    sinomega2 = (- b + sqrt(t)) / (a * 2.)   # a can clearly be zero
    sinomega1 = sinomega1*(1-another_valid)
    sinomega2 = sinomega2*(1-another_valid)
    #
    # Now we need also cosomega to isolate the two choices of arcsin
    #
    #      lhs = g0 cos(omega) - g1 sin(omega)
    #      g0 cos(omega) - lhs =  g1 sin(omega)
    #      g0^2 cos^2(omega) - 2 lhs g0 cos(omega) + lhs^2 =  g1^2 sin^2(omega)
    #      (g0^2+g1^2) cos^2(omega) - 2 lhs g0 cos(omega) + lhs^2 - g1^2  = 0
    #
    c = lhs * lhs - g[1,:]*g[1,:]
    b = -2. * g[0,:] * lhs
    t = b*b - 4*a*c
    valid = where(t>0.0 , 1., 0.)
    t = t * valid
    another_valid = where(a < 1e-12, 1, 0)
    a = where(another_valid == 1, 1, a)
    cosomega1 = (- b - sqrt(t)) / (a * 2.)   # a can clearly be zero
    cosomega2 = (- b + sqrt(t)) / (a * 2.)   # a can clearly be zero
    cosomega1 = cosomega1*(1-another_valid)
    cosomega2 = cosomega2*(1-another_valid)

    #
    # Now we need to find out which solutions go together via sin^2 + cos^2 = 1
    # ... turns out it could be arctan2(sin1,cos1) or arctan2(sin1,cos2)
    #
    hypothesis_1 = abs(sinomega1*sinomega1+cosomega1*cosomega1 - 1.)
    hypothesis_2 = abs(sinomega1*sinomega1+cosomega2*cosomega2 - 1.)
    same_sign = where(abs(hypothesis_1 < hypothesis_2 ), 1, 0)
    omega1 = same_sign * arctan2(sinomega1,cosomega1) + (1.-same_sign) * arctan2(sinomega1,cosomega2)
    omega2 = same_sign * arctan2(sinomega2,cosomega2) + (1.-same_sign) * arctan2(sinomega2,cosomega1)
    #
    # Finally re-compute cosomega and sinomega due to the flipping possibility for getting eta
    cosomega1 = cos(omega1)
    sinomega1 = sin(omega1)
    cosomega2 = cos(omega2)
    sinomega2 = sin(omega2)
    #
    # Now we know R-1 and W-1 , so compute k = W-1 R-1 g
    #
    R1g1 = zeros(g.shape,Float)
    R1g2 = zeros(g.shape,Float)
    #
    R1g1[0,:] = cosomega1*g[0,:] - sinomega1*g[1,:]
    R1g1[1,:] = sinomega1*g[0,:] + cosomega1*g[1,:]
    R1g1[2,:] =                                     g[2,:]
    #
    R1g2[0,:] = cosomega2*g[0,:] - sinomega2*g[1,:]
    R1g2[1,:] = sinomega2*g[0,:] + cosomega2*g[1,:]
    R1g2[2,:] =                                     g[2,:]
    #
    #
    k_one =  zeros(g.shape,Float)
    k_two =  zeros(g.shape,Float)
    #
    # W-1 = ( cos(wedge) ,  0  , -sin(wedge) )
    #       (         0  ,  1  ,          0  )
    #       ( sin(wedge) ,  0  ,  cos(wedge) )
    #
    k_one[0,:] = cw * R1g1[0,:]                 - sw * R1g1[2,:]
    k_one[1,:] =                    R1g1[1,:]
    k_one[2,:] = sw * R1g1[0,:]                 + cw * R1g1[2,:]
    #
    k_two[0,:] = cw * R1g2[0,:]                 - sw * R1g2[2,:]
    k_two[1,:] =                    R1g2[1,:]
    k_two[2,:] = sw * R1g2[0,:]                 + cw * R1g2[2,:]
    #
    # k[1,:] = -ds*c*sin(eta)
    # ------    -------------   .... tan(eta) = -k1/k2
    # k[2,:] =  ds*c*cos(eta)
    #
    eta_one = arctan2(-k_one[1,:],k_one[2,:])
    eta_two = arctan2(-k_two[1,:],k_two[2,:])
    #
    #
    tth = degrees(arcsin(s)*2.) * valid
    eta1 = degrees(eta_one)*valid
    eta2 = degrees(eta_two)*valid
    omega1 = degrees(omega1)*valid
    omega2 = degrees(omega2)*valid
    return tth,[eta1,eta2],[omega1,omega2]

def old_uncompute_g_vectors(gv,wavelength, wedge=0.0):
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



def compute_lorentz_factors(args):
    """
    From Kabsch 1988 J. Appl. Cryst. 21 619

    Multiply the intensities by:
    Lorentz = | S.(u x S_0)| / |S|.|S_0|
    S = scattered vector
    S_0 = incident vector
    u = unit vector along rotation axis
    """
    pass

def compute_polarisation_factors(args):
    """
    From Kabsch 1988 J. Appl. Cryst. 21 619

    DIVIDE the intensities by:
    <sin2 psi> = ( 1 - 2p) [ 1 - (n.S/|S|^2) ] + p { 1 + [S.S_0/(|S||S_0|)^2]^2}

    p = degree of polarisation (sync = 1, tube = 0.5 , mono + tube in between)
        or "probability of finding electric field vector in plane having normal, n"
    S = scattered vector
    S_0 = incident vector
    n = normal to polarisation plane, typically perpendicular to S_0
    """
    pass

class transformer:
    """
    A class for transformation
    """
    def __init__(self, cell__a=2.9508, cell__b=2.9508, cell__c=4.6855, cell_alpha=90.0, cell_beta=90.0, cell_gamma=120.0,
                    cell_lattice="P", distance=303.7, fit_tolerance=0.05, omegasign=1, tilt_y=0.0, tilt_z=0.0,
                    wavelength=0.155, y_center=1005.768, y_size=48.08150, z_center=1002.832, z_size=46.77648):

        self.twotheta=None
        self.parameters={}
        self.wedge=0.0
        self.chi=0.0
        
        self.parameters['omegasign']=omegasign
        self.parameters['z-center']=z_center
        self.parameters['y-center']=y_center
        self.parameters['distance']=distance
        self.parameters['z-size']=z_size
        self.parameters['y-size']=y_size
        self.parameters['tilt-z']=tilt_z
        self.parameters['tilt-y']=tilt_y
        self.parameters['fit_tolerance']=fit_tolerance
        self.parameters['wavelength']=wavelength
        self.parameters['cell__a']=cell__a
        self.parameters['cell__b']=cell__b
        self.parameters['cell__c']=cell__c
        self.parameters['cell_alpha']=cell_alpha
        self.parameters['cell_beta']=cell_beta
        self.parameters['cell_gamma']=cell_gamma
        self.parameters['cell_lattice_[P,A,B,C,I,F]']=cell_lattice
        
        # it would make more sense to inherit the parameter object - will
        # have to think about this some more - how general is it?
        from ImageD11 import parameters
        self.parameterobj = parameters.parameters(cell__a=self.parameters['cell__a'],
                    cell__b=self.parameters['cell__b'], cell__c=self.parameters['cell__c'],
                    cell_alpha=self.parameters['cell_alpha'], cell_beta=self.parameters['cell_beta'],
                    cell_gamma=self.parameters['cell_gamma'], cell_lattice=self.parameters['cell_lattice_[P,A,B,C,I,F]'],
                    distance=self.parameters['distance'], fit_tolerance=self.parameters['fit_tolerance'],
                    omegasign=self.parameters['omegasign'], tilt_y=self.parameters['tilt-y'], tilt_z=self.parameters['tilt-z'],
                    wavelength=self.parameters['wavelength'], y_center=self.parameters['y-center'],
                    y_size=self.parameters['y-size'], z_center=self.parameters['z-center'], z_size=self.parameters['z-size'])


    def updatepars(self, cell__a, cell__b, cell__c, cell_alpha, cell_beta, cell_gamma,
                    cell_lattice, distance, fit_tolerance, omegasign, tilt_y, tilt_z,
                    wavelength, y_center, y_size, z_center, z_size):
        """
        Updates values of parameters, temporary used for 
        XML-RPC communication 
        """
        self.parameters['omegasign']=omegasign
        self.parameters['z-center']=z_center
        self.parameters['y-center']=y_center
        self.parameters['distance']=distance
        self.parameters['z-size']=z_size
        self.parameters['y-size']=y_size
        self.parameters['tilt-z']=tilt_z
        self.parameters['tilt-y']=tilt_y
        self.parameters['fit_tolerance']=fit_tolerance
        self.parameters['wavelength']=wavelength
        self.parameters['cell__a']=cell__a
        self.parameters['cell__b']=cell__b
        self.parameters['cell__c']=cell__c
        self.parameters['cell_alpha']=cell_alpha
        self.parameters['cell_beta']=cell_beta
        self.parameters['cell_gamma']=cell_gamma
        self.parameters['cell_lattice_[P,A,B,C,I,F]']=cell_lattice
        return "OK"

    def loadpars(self,filename=None):
        if filename is not None:
            self.parameterobj.loadparameters(filename)
        self.parameterobj.update_other(self)
        return "OK"

    def savepars(self,filename=None):
        self.parameterobj.update_yourself(self)
        if filename is not None:
            self.parameterobj.saveparameters(filename)
        return "OK"


    def loadfiltered(self,filename=None):
        #      if self.parent.finalpeaks!=None:
        #         from tkMessageBox import askyesno
        #         if not askyesno("Overwrite current peaks?","Loading new peaks
        #will overwrite the ones in memory now, are you sure?"):
        #            return
        f=open(filename,"r")
        #                   0123456789012
        line=f.readline()
        if line[0:12] !="# xc yc omega"[0:12]:
            print line
            print "Sorry! That does not seem to be a filter peaks file, output from the peaksearching menu option"
            return "Not OK"
        bigarray=[]
        for line in f.readlines():
            v=[float(z) for z in line.split()]
            bigarray.append(v)
        f.close()
        self.finalpeaks=transpose(array(bigarray))
        print self.finalpeaks.shape
        return "OK"
    
    def gety(self):
        """
        Sends x values to plot the x,y arrays
        """
        return self.finalpeaks[0,:].tolist()

    def getz(self):
        """
        Sends x values to plot the x,y arrays
        """
        return self.finalpeaks[1,:].tolist()
    
    def gettwotheta(self):
        """
        Sends x values to plot the x,y arrays
        """
        return self.twotheta.tolist()

    def geteta(self):
        """
        Sends x values to plot the x,y arrays
        """
        return self.eta.tolist()
    
    def gettheorytth(self):
        """
        Sends x values to plot the x,y arrays
        """
        return self.theorytth.tolist()
    

    def plotreta(self):
        self.peaks_xy = self.finalpeaks[0:2,:]
        self.x=self.peaks_xy[0,:]
        self.y=self.peaks_xy[1,:]
        self.omega=self.finalpeaks[2,:]
        print self.peaks_xy.shape
        try:
            self.twotheta, self.eta = compute_tth_eta( self.peaks_xy,
                  float(self.parameters['y-center']), # yc is the centre in y
                  float(self.parameters['y-size'])  , # ys is the y pixel size
                  float(self.parameters['tilt-y'])  , # ty is the tilt around y
                  float(self.parameters['z-center']), # zc is the centre in z
                  float(self.parameters['z-size'])  , # zs is the z pixel size
                  float(self.parameters['tilt-z'])  , # tz is the tilt around z
                  float(self.parameters['distance'])*1e3) # is the sample - detector distance
            self.ds = 2*sin(radians(self.twotheta)/2)/float(self.parameters['wavelength'])
        except:
            keys=self.parameters.keys()
            keys.sort()
            for key in keys:
                print self.parameters[key],type(self.parameters[key])
            raise
        return "OK"

    def addcellpeaks(self):
        #
        # Given unit cell, wavelength and distance, compute the radial positions
        # in microns of the unit cell peaks
        #
        a      =float(self.parameters['cell__a'])
        b      =float(self.parameters['cell__b'])
        c      =float(self.parameters['cell__c'])
        alpha  =float(self.parameters['cell_alpha'])
        beta   =float(self.parameters['cell_beta'])
        gamma  =float(self.parameters['cell_gamma'])
        lattice=self.parameters['cell_lattice_[P,A,B,C,I,F]']
        self.unitcell=unitcell.unitcell( [a,b,c,alpha,beta,gamma] ,  lattice)
        self.unitcell=self.unitcell
        if self.twotheta==None:
            self.twotheta, self.eta = compute_tth_eta( self.peaks_xy,
                  float(self.parameters['y-center']), # yc is the centre in y
                  float(self.parameters['y-size'])  , # ys is the y pixel size
                  float(self.parameters['tilt-y'])  , # ty is the tilt around y
                  float(self.parameters['z-center']), # zc is the centre in z
                  float(self.parameters['z-size'])  , # zs is the z pixel size
                  float(self.parameters['tilt-z'])  , # tz is the tilt around z
                  float(self.parameters['distance']*1e3)) # is the sample - detector distance
        # Find last peak in radius
        highest = maximum.reduce(self.twotheta)
        wavelength=float(self.parameters['wavelength'])
        ds = 2*sin(radians(highest)/2.)/wavelength
        self.dslimit=ds
        #print "highest peak",highest,"corresponding d*",ds
        self.theorypeaks=self.unitcell.gethkls(ds)
        self.unitcell.makerings(ds)
        self.theoryds=self.unitcell.ringds
        tths = [arcsin(wavelength*dstar/2)*2 for dstar in self.unitcell.ringds]
        self.theorytth=degrees(array(tths))
        return "OK"

    def gof(self,args):
        self.tolerance=float(self.parameters['fit_tolerance'])
        self.parameters['y-center']=args[0]
        self.parameters['z-center']=args[1]
        self.parameters['distance']=args[2]
#      self.parameters['wavelength']=args[3]
        self.parameters['tilt-y']=args[3]
        self.parameters['tilt-z']=args[4]
        self.twotheta, self.eta = compute_tth_eta( self.peaks_xy,
                                                         self.parameters['y-center'],
                                                         float(self.parameters['y-size']), # ys is the y pixel size
                                                         self.parameters['tilt-y'],
                                                         self.parameters['z-center'],
                                                         float(self.parameters['z-size'])  , # zs is the z pixel size
                                                         self.parameters['tilt-z'],
                                                         self.parameters['distance']*1e3
                                                         )
        w=float(self.parameters['wavelength']) # args[3] # what?
        gof=0.
        npeaks=0
        for i in range(len(self.tthc)):# (twotheta_rad_cell.shape[0]):
            self.tthc[i]=degrees(math.asin(self.fitds[i]*w/2)*2)
            diff=take(self.twotheta,self.indices[i]) - self.tthc[i]
#         print "peak",i,"diff",maximum.reduce(diff),minimum.reduce(diff)
            gof=gof+dot(diff,diff)
            npeaks=npeaks+len(diff)
        gof=gof/npeaks
        return gof*1e3


    def fit(self):
        import simplex
        # Assign observed peaks to rings
        w=float(self.parameters['wavelength'])
        self.indices=[]
        self.tthc=[]
        self.fitds=[]
        self.tolerance=float(self.parameters['fit_tolerance'])
        from itertools import izip
        tthmax = max(izip(self.twotheta, xrange(len(self.twotheta))))[0]
        print "Tolerance for assigning peaks to rings",self.tolerance,"max tth",tthmax
        for i in range(len(self.theoryds)):
            dsc=self.theoryds[i]
            tthcalc=math.asin(dsc*w/2)*360./math.pi # degrees
            if tthcalc>tthmax:
                break
#         print tthcalc
            logicals= logical_and( greater(self.twotheta, tthcalc-self.tolerance),
                                      less(self.twotheta, tthcalc+self.tolerance)  )

            if sum(logicals)>0:
#            print maximum.reduce(compress(logicals,self.twotheta)),minimum.reduce(compress(logicals,self.twotheta))
                self.tthc.append(tthcalc)
                self.fitds.append(dsc)
                ind=compress(logicals,range(self.twotheta.shape[0]))
                self.indices.append(ind)
#            print "Ring",i,tthcalc,maximum.reduce(take(self.twotheta,ind)),minimum.reduce(take(self.twotheta,ind))
#      if raw_input("OK?")[0] not in ["Y","y"]:
#         return
        guess=[ float(self.parameters['y-center']) ,
                float(self.parameters['z-center']) ,
                float(self.parameters['distance']) ,
  #              float(self.parameters['wavelength']) ,
                float(self.parameters['tilt-y']) ,
                float(self.parameters['tilt-z']) ]
        inc=[ .1 , .1 , .1 , radians(0.1) , radians(0.1) ]
        s=simplex.Simplex(self.gof,guess,inc)
        newguess,error,iter=s.minimize()
        self.parameters['y-center']=newguess[0]
        self.parameters['z-center']=newguess[1]
        self.parameters['distance']=newguess[2]
#      self.parameters['wavelength']=newguess[3]
        self.parameters['tilt-y']=newguess[3]
        self.parameters['tilt-z']=newguess[4]
        inc=[ .01 , .01 , .01 , 0.0001 , radians(0.01) , radians(0.01) ]
        guess=newguess
        s=simplex.Simplex(self.gof,guess,inc)
        newguess,error,iter=s.minimize()
        self.parameters['y-center']=newguess[0]
        self.parameters['z-center']=newguess[1]
        self.parameters['distance']=newguess[2]
#      self.parameters['wavelength']=newguess[3]
        self.parameters['tilt-y']=newguess[3]
        self.parameters['tilt-z']=newguess[4]
        return newguess
        
    def updatechiwedge(self, chi, wedge):
        """
        Allow the rotation axis to not be perpendicular to the beam
        """
        self.wedge = float(wedge)
        self.chi = float(chi)
        return "OK"
 
    def computegv(self):
        """
        Using self.twotheta, self.eta and omega angles, compute x,y,z of spot
        in reciprocal space
        """
        try:
            omegasign = float(self.parameters['omegasign'])
        except:
            omegasign = 1.
        self.gv = compute_g_vectors(self.twotheta,self.eta,self.omega*omegasign
                                              ,float(self.parameters['wavelength']),wedge=self.wedge,chi=self.chi)
        tthnew,etanew,omeganew=uncompute_g_vectors(self.gv,float(self.parameters['wavelength']),wedge=self.wedge)
      #  self.gv=self.gv
        print "Testing reverse transformations"
        for i in range(5):
            print "tth %8.3f eta %8.3f om %8.3f in "%(self.twotheta[i],self.eta[i],self.omega[i]),
            print "tth %8.3f [eta %8.3f om %8.3f or eta %8.3f om %8.3f]"%(tthnew[i],etanew[0][i],omeganew[0][i],etanew[1][i],omeganew[1][i])
        return "OK"      

    def savegv(self,filename=None):
        """
        Save g-vectors into a file
        Use crappy .ass format from previous for now (testing)
        """
        try:
            omegasign = float(self.parameters['omegasign'])
        except:
            omegasign = 1.
        f=open(filename,"w")
        f.write(self.unitcell.tostring())
        f.write("\n")
        f.write("# wavelength = %f\n"%( float(self.parameters['wavelength']) ) )
        f.write("# wedge = %f\n"%( float(self.wedge) ))
        f.write("# ds h k l\n")
        for peak in self.theorypeaks:
            f.write("%10.7f %4d %4d %4d\n"%(peak[0],peak[1][0],peak[1][1],peak[1][2]))
        order = argsort(self.twotheta)
        f.write("# xr yr zr xc yc ds phi omega\n")
        print maximum.reduce(self.omega),minimum.reduce(self.omega)
        for i in order:
            f.write("%f %f %f %f %f %f %f %f \n"%(self.gv[0,i],self.gv[1,i],self.gv[2,i],self.x[i],self.y[i],self.ds[i],self.eta[i],self.omega[i]*omegasign))
        f.close()
        return "OK"
    
    def write_graindex_gv(self,filename):
        
        try:
            omegasign = float(self.parameters['omegasign'])
        except:
            omegasign = 1.
        from Imaged11Functions import write_graindex_gv
        #self.parent.finalpeaks=transpose(array(bigarray))
        self.intensity = self.finalpeaks[3,:]*self.finalpeaks[4,:]
        print self.intensity.shape
        write_graindex_gv.write_graindex_gv(filename,
                                            self.gv,
                                            self.twotheta,
                                            self.eta,
                                            self.omega*omegasign,
                                            self.intensity,
                                            self.unitcell)
        return "OK"


if __name__=="__main__":
    # from indexing import mod_360
    def mod_360(theta, target):
        """
        Find multiple of 360 to add to theta to be closest to target
        """
        diff=theta-target
        while diff < -180:
            theta=theta+360
            diff=theta-target
        while diff > 180:
            theta=theta-360
            diff=theta-target
        return theta

    tth = array([   1,  2,  3,  4,  5,  6,  7,  8,  9, 10], Float)
    eta = array([  10, 40, 70,100,130,160,190,220,270,340], Float)
    om  = array([   0, 20, 40,100, 60,240,300, 20, 42, 99], Float)
#   tth = array([    4, 5, 7], Float)
#   eta = array([  100, 5, 190], Float)
#   om  = array([  100, 5, 300], Float)

    for wavelength in [ 0.1, 0.2, 0.3]:
        for wedge in [-10. , -5., 0., 5., 10.]:
            print "Wavelength",wavelength,"wedge",wedge
            print "tth, eta, omega   ...   tth, eta, omega   ... tth, eta, omega"
            gv = compute_g_vectors(tth,eta,om,wavelength,wedge)
            t,e,o = new_uncompute_g_vectors(gv,wavelength, wedge)
            for i in range(tth.shape[0]):
                print "%9.3f %9.3f %9.3f  "%(tth[i],eta[i],om[i]),
                print "%9.3f %9.3f %9.3f  "%(t[i],mod_360(e[0][i],eta[i]),mod_360(o[0][i],om[i])),
                print "%9.3f %9.3f %9.3f  "%(t[i],mod_360(e[1][i],eta[i]),mod_360(o[1][i],om[i])),
                # Choose best fitting
                e_eta1 = mod_360(e[0][i],eta[i]) - eta[i]
                e_om1  = mod_360(o[0][i], om[i]) -  om[i]
                score_1 = e_eta1*e_eta1 + e_om1*e_om1
                e_eta2 = mod_360(e[1][i],eta[i]) - eta[i]
                e_om2  = mod_360(o[1][i], om[i]) -  om[i]
                score_2 = e_eta2*e_eta2 + e_om2*e_om2
                print "%.5g %.5g"%(score_1,score_2)
