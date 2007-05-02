


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
Functions for transforming peaks
"""
import logging

from Numeric import *

try:
    # crazy debug 
    test = arccos(zeros(10,Float))
except:
    print dir()
    raise

from math import pi

def cross_product(a,b):
    """ returns axb for two len(3) vectors a,b"""
    assert len(a)==len(b)==3
    return array([a[1]*b[2]-a[2]*b[1],
                  a[2]*b[0]-a[0]*b[2],
                  a[0]*b[1]-a[1]*b[0]])

def degrees(x):
    """Convenience function"""
    return x*180.0/pi

def radians(x):
    """Convenience function"""
    return x*pi/180.0

def compute_tth_eta(peaks,y_center=0.,y_size=0.,tilt_y=0.,
                    z_center=0.,z_size=0.,tilt_z=0.,
                    tilt_x = 0.,
                    distance=0.,
                    # detector_orientation=((1,0),(0,1)),
                    o11 = 1.0 , o12 = 0.0, o21 = 0.0 , o22 = -1.0,
                    # crystal_translation = None,
                    t_x = 0.0, t_y = 0.0 , t_z = 0.0,
                    omega = None,         #       == phi at chi=90
                    wedge = 0.0, # Wedge == theta on 4circ
                    chi = 0.0, #       == chi - 90
                    **kwds): # last line is for laziness - 
                             # pass kwds you'd like to be ignored
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
    r1 = array( [ [  cos(tilt_z) ,-sin(tilt_z) , 0 ], # note this is r.h.
                  [  sin(tilt_z) , cos(tilt_z) , 0 ],
                  [    0         ,    0        , 1 ]],Float)
    r2 = array( [ [ cos(tilt_y) , 0 , sin(tilt_y) ],
                  [       0     , 1 ,   0     ],
                  [-sin(tilt_y) , 0 , cos(tilt_y) ]],Float)
    r3 = array( [ [  1 ,          0  ,       0     ],
                  [  0 ,  cos(tilt_x), sin(tilt_x) ],
                  [  0 , -sin(tilt_x), cos(tilt_x) ]],Float)
    r2r1=matrixmultiply(matrixmultiply(r1,r2),r3)
    # Peak positions in 3D space
    #  - apply detector orientation
    peaks_on_detector = array(peaks)
    peaks_on_detector[0,:] =  (peaks_on_detector[0,:]-z_center)*z_size
    peaks_on_detector[1,:] =  (peaks_on_detector[1,:]-y_center)*y_size
    # 
    detector_orientation = [[o11,o12],[o21,o22]]
    logging.debug("detector_orientation = "+str(detector_orientation))
    flipped = matrixmultiply(array(detector_orientation,Float),
                             peaks_on_detector)   
    # 
    vec =array( [ zeros(flipped.shape[1])     , # place detector at zero, 
                                                # sample at -dist
                    flipped[1,:]      ,        # x in search, frelon +z
                    flipped[0,:]]    , Float) # y in search, frelon -y 
    #print vec.shape
    # Position of diffraction spots in 3d space after detector tilts is:
    rotvec=matrixmultiply(r2r1,vec)
    # Scattering vectors
    if omega is None or ( t_x == 0. and t_y == 0 and t_z == 0):
        magrotvec=sqrt(sum(rotvec*rotvec,0))
        # tan(eta) = y/z
        eta=degrees(arctan2(rotvec[1,:],rotvec[2,:])) 
        # cosine rule a2 = b2+c2-2bccos(A)
        # a is distance from (0,0,0) to rotvec => magrotvec
        # b is distance from sample to detector
        # c is distance from sample to pixel
        a=magrotvec  # 1D
        b=distance   # 0D
        c=rotvec     # 3D
        #print c.shape
        c[0,:]=c[0,:]+distance# 3D
        # print c
        c=sqrt(sum(c*c,0))# 1D
        #print a.shape,c.shape
        #   print a.shape,c.shape
        costwotheta = (b*b + c*c - a*a)/2/b/c
        twothetarad=arccos(costwotheta)
        twotheta=degrees(twothetarad)

        return twotheta, eta
    else:
        # print "Using translations t_x %f t_y %f t_z %f"%(t_x,t_y,t_z)
        # Compute positions of grains
        # expecting tx, ty, tz for each diffraction spot
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
        w=radians(wedge)
        WI = array( [ [ cos(w),         0, -sin(w)],
                      [      0,         1,       0],
                      [ sin(w),         0,  cos(w)] ] , Float)
        c=radians(chi)
        CI = array( [ [      1,          0,       0],
                      [      0,     cos(c), -sin(c)],
                      [      0,     sin(c),  cos(c)] ] , Float)
        if omega.shape[0] != peaks.shape[1]:
            raise Exception(
                "omega and peaks arrays must have same number of peaks")
        eta=zeros(omega.shape[0],Float)
        tth=zeros(omega.shape[0],Float)
        t = zeros((3,omega.shape[0]),Float) # crystal translations
        # Rotations in reverse order compared to making g-vector
        # also reverse directions. this is trans at all zero to
        # current setting. gv is scattering vector to all zero
        om_r = radians(omega)
                # This is the real rotation (right handed, g back to k)
        t[0,:] = cos(om_r)*t_x - sin(om_r)*t_y
        t[1,:] = sin(om_r)*t_x + cos(om_r)*t_y
        t[2,:] =                                  t_z
        if wedge != 0.0:
            c = cos(radians(wedge))
            s = sin(radians(wedge))
            u = zeros(t.shape,Float)
            u[0,:]= c * t[0,:]           + -s * t[2,:]
            u[1,:]=            t[1,:]
            u[2,:]= s * t[0,:]           + c * t[2,:]
            t = u
        if chi != 0.0:
            c = cos(radians(chi))
            s = sin(radians(chi))
            u = zeros(t.shape,Float)
            u[0,:]= t[0,:]  
            u[1,:]=        c * t[1,:]    + -s * t[2,:]
            u[2,:]=        s * t[1,:]    + c * t[2,:]
            t = u
        myorigins = t.copy()
        myorigins[0,:] = myorigins[0,:] - distance
        # origin =  array([ -distance, 0, 0 ], Float)
        # scattering_vectors 
        s = rotvec - myorigins
        eta = degrees(arctan2(s[1],s[2]))
        mag_s = sqrt(sum(s*s,0))
        costth = s[0] / mag_s
        tth = degrees( arccos(costth) )
        return tth , eta
                      
    
def compute_tth_histo(finalpeaks,tth,no_bins=0,min_bin_ratio=1,
                    **kwds): # last line is for laziness - 
                             # pass kwds you'd like to be ignored
    tthsort = sort(tth)
    maxtth = tthsort[-1]
    logging.debug("maxtth=%f"%(maxtth))
    binsize = (maxtth+0.001)/no_bins
    tthbin = array(range(no_bins))*binsize # runs from 0->maxtth
    # vectorise this loop below from old Numeric manual
    #for t in tthsort:
    #    n= int(floor(t/binsize))
    #    histogram[n] = histogram[n] +1
    n = searchsorted(tthsort,tthbin) # position of bin in sorted
    n = concatenate([n,[len(tthsort)]])   # add on last position
    histogram = (n[1:] - n[:-1])*1.0        # this would otherwise be integer
    logging.debug("max(histogram) = %d"%(max(histogram)))
    histogram=histogram/max(histogram)
    
    # keeppeaks = [] 
    #for t in range(tth.shape[0]):
    #    if histogram[bins[t]]> min_bin_ratio:
    #        keeppeaks.append(t)
    #keeppeaks = tuple(keeppeaks)
    #finalpeaks = take(finalpeaks,keeppeaks,1)
    
    # Vectorised version
    bins = floor(tth/binsize).astype(Int) # bin for each two theta
    hpk = take(histogram, bins) # hist for each peak
    keepind = compress( hpk > float(min_bin_ratio) , range(tth.shape[0]) )
    print len(keepind),tth.shape[0]
    filteredpeaks = take(finalpeaks, keepind, 1) # beware of modifying array arg...

    return filteredpeaks
                             
def compute_g_vectors(tth, eta, omega, wavelength, wedge = 0.0, chi=0.0):
    """
    Generates spot positions in reciprocal space from 
      twotheta, wavelength, omega and eta
    Assumes single axis vertical
    ... unless a wedge angle is specified
    """
    #if wedge != 0.:
    #    print "Using a wedge angle of ",wedge
    tth=radians(tth)
    eta=radians(eta)

    om =radians(omega)
    # Compute k vector - the scattering in the laboratory
    c=cos(tth/2) # cos theta
    s=sin(tth/2) # sin theta
    ds=2*s/wavelength
    k=zeros((3,tth.shape[0]),Float)
    # x - along incident beam
    k[0,:] = -ds*s # this is negative x
    # y - towards door
    k[1,:] =  ds*c*sin(eta) # y direction
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
    # This is the reverse rotation (left handed, k back to g)
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
    #
    #
    #      - g0^2 + lhs^2 + 
    #        2 g1 lhs sin(omega) + 
    #          g1^2 sin^2(omega) + g0^2 sin^2(omega)  = 0
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
    omega1 = same_sign * arctan2(sinomega1,cosomega1) + \
        (1.-same_sign) * arctan2(sinomega1,cosomega2)
    omega2 = same_sign * arctan2(sinomega2,cosomega2) + \
        (1.-same_sign) * arctan2(sinomega2,cosomega1)
    #
    # Finally re-compute cosomega and sinomega due to the 
    # flipping possibility for getting eta
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
        print "Maximum and minimum in sin(theta)=",\
                 maximum.reduce(s),minimum.reduce(s)
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



def compute_lorentz_factors(tth, eta, omega, wavelength, wedge=0., chi=0.):
    """
    From Kabsch 1988 J. Appl. Cryst. 21 619

    Multiply the intensities by:
    Lorentz = | S.(u x So)| / |S|.|So|
    S = scattered vector
    So = incident vector
    u = unit vector along rotation axis
    """
    # So is along +x, the incident beam defines the co-ordinates in ImageD11
    # length is in reciprocal space units, 1/lambda
    So = [1./wavelength,0,0]
    #
    # u vector along rotation axis
    # starts as along z
    u = [0,0,1]
    # rotate by 
    # g =  R . W . k where:
    # R = ( cos(omega) , sin(omega), 0 )
    #     (-sin(omega) , cos(omega), 0 )
    #     (         0  ,         0 , 1 )
    #
    W =  [[ cos(wedge) ,  0  ,  sin(wedge) ],
          [         0  ,  1  ,          0  ],
          [-sin(wedge) ,  0  ,  cos(wedge) ]]
    #
    C =  [[         1  ,         0  ,      0     ],
          [         0  ,  cos(chi)  , sin(chi)   ],
          [         0  , -sin(chi)  , cos(chi)   ]]
    u = matrixmultiply(C,matrixmultiply(W,u))
    u_x_So = cross_product(u,So)
    # if DEBUG: print "axis orientation",u
    #
    # S = scattered vectors. Length 1/lambda.
    S = [ cos(radians(tth)/2.)*sin(radians(eta))/wavelength,
          cos(radians(tth)/2.)*cos(radians(eta))/wavelength,
          sin(radians(tth)/2.)/wavelength]
    S_dot_u_x_So = dot(S,u_x_So)
    mod_S = sqrt(dot(S,S))
    mod_So = sqrt(dot(So,So))
    try:
        lorentz = abs(S_dot_u_x_So)/mod_S/mod_So
    except:
        raise Exception("Please fix this div0 crap in lorentz")
    return lorentz


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

    In ImageD11 we normally expect to find:
    x axis along the beam
    z axis being up, and parallel to the normal n mentioned above
    """
    n = [0,0,1]
    

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
            print "tth, eta, omega   ...   " +\
                  "tth, eta, omega   ...   " +\
                  "tth, eta, omega"
            gv = compute_g_vectors(tth,eta,om,wavelength,wedge)
            t,e,o = uncompute_g_vectors(gv,wavelength, wedge)
            for i in range(tth.shape[0]):
                print "%9.3f %9.3f %9.3f  "%(tth[i],eta[i],om[i]),
                print "%9.3f %9.3f %9.3f  "%(t[i],
                             mod_360(e[0][i],eta[i]),
                             mod_360(o[0][i],om[i])),
                print "%9.3f %9.3f %9.3f  "%(t[i],
                             mod_360(e[1][i],eta[i]),
                             mod_360(o[1][i],om[i])),
                # Choose best fitting
                e_eta1 = mod_360(e[0][i],eta[i]) - eta[i]
                e_om1  = mod_360(o[0][i], om[i]) -  om[i]
                score_1 = e_eta1*e_eta1 + e_om1*e_om1
                e_eta2 = mod_360(e[1][i],eta[i]) - eta[i]
                e_om2  = mod_360(o[1][i], om[i]) -  om[i]
                score_2 = e_eta2*e_eta2 + e_om2*e_om2
                print "%.5g %.5g"%(score_1,score_2)
