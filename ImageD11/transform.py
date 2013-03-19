


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

import numpy as n
 
try:
    # crazy debug 
    test = n.arccos(n.zeros(10,n.float))
except:
    print dir()
    raise

from math import pi

def cross_product_2x2(a,b):
    """ returns axb for two len(3) vectors a,b"""
    assert len(a)==len(b)==3
    return n.array([a[1]*b[2]-a[2]*b[1],
                  a[2]*b[0]-a[0]*b[2],
                  a[0]*b[1]-a[1]*b[0]])

def degrees(x):
    """Convenience function"""
    return x*180.0/pi

def radians(x):
    """Convenience function"""
    return x*pi/180.0

def detector_rotation_matrix(tilt_x, tilt_y, tilt_z):
    """
    Return the tilt matrix to apply to peaks
    tilts in radians
    typically applied to peaks rotating around beam center
    """
    r1 = n.array( [ [  n.cos(tilt_z) ,-n.sin(tilt_z) , 0 ], # note this is r.h.
                    [  n.sin(tilt_z) , n.cos(tilt_z) , 0 ],
                    [    0         ,    0        , 1 ]],n.float)
    r2 = n.array( [ [ n.cos(tilt_y) , 0 , n.sin(tilt_y) ],
                    [       0     , 1 ,   0     ],
                    [-n.sin(tilt_y) , 0 , n.cos(tilt_y) ]],n.float)
    r3 = n.array( [ [  1 ,          0  ,       0     ],
                    [  0 ,  n.cos(tilt_x), -n.sin(tilt_x) ],
                    [  0 ,  n.sin(tilt_x), n.cos(tilt_x) ]],n.float)
    r2r1 = n.dot(n.dot(r3,r2),r1)
    return r2r1



def compute_xyz_lab(peaks,
                    y_center=0.,y_size=0.,tilt_y=0.,
                    z_center=0.,z_size=0.,tilt_z=0.,
                    tilt_x = 0.,
                    distance=0.,
                    # detector_orientation=((1,0),(0,1)),
                    o11 = 1.0 , o12 = 0.0, o21 = 0.0 , o22 = -1.0,
                    **kwds ):
    """
    Peaks is a 2 d array of x,y
    yc is the centre in y
    ys is the y pixel size
    ty is the tilt around y
    zc is the centre in z
    zs is the z pixel size
    tz is the tilt around z
    dist is the sample - detector distance
    detector_orientation is a matrix to apply to peaks arg to get 
    ImageD11 convention
         (( 0, 1),( 1, 0)) for ( y, x)
         ((-1, 0),( 0, 1)) for (-x, y)
         (( 0,-1),(-1, 0)) for (-y,-x)
      etc...
    """
    assert len(peaks)==2, "peaks must be a 2D array"
    # Matrix for the tilt rotations
    r2r1 = detector_rotation_matrix( tilt_x, tilt_y, tilt_z )
    # Peak positions in 3D space
    #  - apply detector orientation
    peaks_on_detector = n.array(peaks)
    peaks_on_detector[0,:] =  (peaks_on_detector[0,:]-z_center)*z_size
    peaks_on_detector[1,:] =  (peaks_on_detector[1,:]-y_center)*y_size
    # 
    detector_orientation = [[o11,o12],[o21,o22]]
    #logging.debug("detector_orientation = "+str(detector_orientation))
    flipped = n.dot(n.array(detector_orientation, n.float),
                             peaks_on_detector)   
    # 
    vec = n.array( [ n.zeros(flipped.shape[1])     , # place detector at zero, 
                                                     # sample at -dist
                     flipped[1,:]      ,             # x in search, frelon +z
                     flipped[0,:]]    , n.float)     # y in search, frelon -y 
    # Position of diffraction spots in 3d space after detector tilts about
    # the beam centre on the detector
    rotvec = n.dot(r2r1, vec)
    # Now add the distance (along x)
    rotvec[0,:] = rotvec[0,:] + distance
    return rotvec


def compute_tth_eta(peaks,
                    y_center=0.,y_size=0.,tilt_y=0.,
                    z_center=0.,z_size=0.,tilt_z=0.,
                    tilt_x = 0.,
                    distance=0.,
                    # detector_orientation=((1,0),(0,1)),
                    o11 = 1.0 , o12 = 0.0, o21 = 0.0 , o22 = -1.0,
                    t_x = 0.0, t_y = 0.0 , t_z = 0.0,
                    omega = None, #       == phi at chi=90
                    wedge = 0.0,  # Wedge == theta on 4circ
                    chi = 0.0,    #       == chi - 90
                    **kwds): # last line is for laziness - 
                             # pass kwds you'd like to be ignored
    """
    0/10 for style
    """
    peaks_xyz =  compute_xyz_lab( 
        peaks,
        y_center = y_center,y_size = y_size,tilt_y = tilt_y,
        z_center = z_center,z_size = z_size,tilt_z = tilt_z,
        tilt_x = tilt_x,
        distance = distance,
        # detector_orientation=((1,0),(0,1)),
        o11 = o11 , o12 = o12, o21 = o21 , o22 = o22 )
    
    tth, eta = compute_tth_eta_from_xyz(
        peaks_xyz ,
        t_x = t_x, t_y = t_y , t_z = t_z,
        omega = omega, 
        wedge = wedge, 
        chi = chi)   
        
    return tth, eta
    

def compute_tth_eta_from_xyz(peaks_xyz , omega , 
                             t_x = 0.0, t_y = 0.0 , t_z = 0.0,
                             #       == phi at chi=90
                             wedge = 0.0,  # Wedge == theta on 4circ
                             chi = 0.0,    #       == chi - 90
                             **kwds): # last line is for laziness - 
    """
    Peaks is a 3 d array of x,y,z peak co-ordinates
    crystal_translation is the position of the grain giving rise 
    to a diffraction spot
                 in x,y,z ImageD11 co-ordinates
                 x,y with respect to axis of rotation and or beam centre ??
                 z with respect to beam height, z centre
    omega data needed if crystal translations used
    """
    assert len(peaks_xyz) == 3
    # Scattering vectors
    if omega is None or ( t_x == 0. and t_y == 0 and t_z == 0):
        s1 = peaks_xyz
    else:
        # scattering_vectors 
        if len(omega) != len(peaks_xyz[0]):
            raise Exception(
                "omega and peaks arrays must have same number of peaks")
        s1 = peaks_xyz - compute_grain_origins( omega, wedge, chi,
                                                t_x, t_y, t_z)
    # CHANGED to HFP convention 4-9-2007
    eta = degrees(n.arctan2(-s1[1,:],s1[2,:]))
    s1_perp_x = n.sqrt(s1[1,:]*s1[1,:] + s1[2,:]*s1[2,:])
    tth = degrees( n.arctan2( s1_perp_x, s1[0,:]) )
    return tth , eta


def compute_xyz_from_tth_eta( tth , eta , omega , 
                              t_x = 0.0, t_y = 0.0 , t_z = 0.0,
                              #       == phi at chi=90
                              wedge = 0.0,  # Wedge == theta on 4circ
                              chi = 0.0,    #       == chi - 90
                              **kwds): # last line is for laziness - 
    """
    Given the tth, eta and omega, compute the xyz on the detector
    
    crystal_translation is the position of the grain giving rise 
    to a diffraction spot
                 in x,y,z ImageD11 co-ordinates
                 x,y with respect to axis of rotation and or beam centre ??
                 z with respect to beam height, z centre
                 
    omega data needed if crystal translations used
    """
    
    # xyz = unit vectors along the scattered vectors
    xyz = n.zeros( (3, tth.shape[0]), n.float)
    rtth = radians(tth)
    reta = radians(eta)
    xyz[0,:] =  n.cos( rtth )
    xyz[1,:] = -n.sin( rtth )*n.sin( reta )
    xyz[2,:] =  n.sin( rtth )*n.cos( reta )

    # Find some vectors in the fast, slow directions in the detector plane
    
    pks = n.array( [ ( 0, 0 ),
                     ( 1, 0 ),
                     ( 0, 1 ) ] )
    # ... not sure this is needed - puts the origin at the beam center
    # pks += n.array( [ float(kwds['z_center']) ,
    # float(kwds['y_center'])]  )

    # 3 peaks in the detector plane
    dxyzl = compute_xyz_lab(peaks, **kwds)

    # Detector plane normal
    norman = cross_product_2x2( dxyzl[1]-dxyzl[0],
                                dxyzl[2]-dxyzl[0] )

    # Equation of a plane gives:
    #   N dot ( P - Po ) == 0
    # ... where P is any point on the plane going through Po
    #     that is normal to P.
    # P is on our line, so:
    #     P = grain_origin + t * xyz
    # Also:
    #     N dot( grain_origin + t * xyz - Po ) == 0
    # ... solve for t ...
    
    grain_origins = compute_grain_origins(omega, **kwds )
    
    N_txyz = dot(norman, dxyzl[0] - grain_origins )
    N_xyz = dot( norman, xyz )

    # If N_xyz is zero anywhere the beam is || to plane and cannot
    # intersect.
    # expect a div0 and set to infinity
    N_xyz = n.where( N_xyz <= 0, 1e35, N_xyz )

    t = N_txyz / N_xyz
    
    # Now we have the teas we just wander along the lines
    xyz_det = t * xyz

    # And finally find those points in terms of our original vectors
    sc = dot( dxyzl[1] - dxyzl[0], xyz_det - dxyzl[0] )
    fc = dot( dxyzl[2] - dxyzl[0], xyz_det - dxyzl[0] )

    return fc, sc
    

def compute_grain_origins(omega, wedge = 0.0, chi = 0.0,
                          t_x = 0.0, t_y = 0.0, t_z = 0.0):
    """
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
    """
    w=radians(wedge)
    WI = n.array( [ [ n.cos(w),         0, -n.sin(w)],
                    [      0,           1,         0],
                    [ n.sin(w),         0,  n.cos(w)] ] , n.float)
    c=radians(chi)
    CI = n.array( [ [      1,            0,         0],
                    [      0,     n.cos(c), -n.sin(c)],
                    [      0,     n.sin(c),  n.cos(c)] ] , n.float)
    t   = n.zeros((3,omega.shape[0]),n.float) # crystal translations
    # Rotations in reverse order compared to making g-vector
    # also reverse directions. this is trans at all zero to
    # current setting. gv is scattering vector to all zero
    om_r = radians(omega)
    # This is the real rotation (right handed, g back to k)
    t[0,:] = n.cos(om_r)*t_x - n.sin(om_r)*t_y
    t[1,:] = n.sin(om_r)*t_x + n.cos(om_r)*t_y
    t[2,:] =                                  t_z
    if chi != 0.0:
        c = n.cos(radians(chi))
        s = n.sin(radians(chi))
        u = n.zeros(t.shape,n.float)
        u[0,:]= t[0,:]  
        u[1,:]=        c * t[1,:]    + -s * t[2,:]
        u[2,:]=        s * t[1,:]    +  c * t[2,:]
        t = u
    if wedge != 0.0:
        c = n.cos(radians(wedge))
        s = n.sin(radians(wedge))
        u = n.zeros(t.shape,n.float)
        u[0,:]= c * t[0,:]           + -s * t[2,:]
        u[1,:]=            t[1,:]
        u[2,:]= s * t[0,:]           +  c * t[2,:]
        t = u
    return t


                      
    
def compute_tth_histo(tth, no_bins = 100,
                    **kwds):
    """
    Compute a histogram of tth values
    
    Returns a normalised histogram (should make this a probability
    *and*
     For each datapoint, the number of other points in the same bin
    """
    tthsort = n.sort(tth)
    maxtth = tthsort[-1]
    mintth = tthsort[0]
    logging.debug("maxtth=%f , mintth=%f"%(maxtth, mintth))
    binsize = (maxtth-mintth)/(no_bins+1)
    tthbin = n.arange(mintth, maxtth+binsize, binsize)
    # print len(tthbin),tthbin[:10]
    nn = n.searchsorted(tthsort,tthbin) # position of bin in sorted
    nn = n.concatenate([nn,[len(tthsort)]])   # add on last position
    histogram = (nn[1:] - nn[:-1]).astype(n.float32) 
    # this would otherwise be integer
    logging.debug("max(histogram) = %d"%(max(histogram)))
    # Change from max
    # histogram = histogram/max(histogram)
    histogram = histogram/len(tth)
    # Vectorised version
    bins = n.floor((tth-mintth)/binsize).astype(n.int) # bin for each two theta
    #print "got bins",len(bins),len(tth),len(histogram)
    #print "bins",bins[:10]
    #print "tth",tth[:10]
    #print "histogram",histogram[:10]
    hpk = n.take(histogram, bins) # histogram value for each peak
    #print "hpk",hpk[:10]
    return tthbin, histogram, hpk
  
def compute_k_vectors(tth, eta, wvln):
    """
    generate k vectors - scattering vectors in laboratory frame
    """
    tth=radians(tth)
    eta=radians(eta)
    c=n.cos(tth/2) # cos theta
    s=n.sin(tth/2) # sin theta
    ds=2*s/wvln
    k=n.zeros((3,tth.shape[0]),n.float)
    # x - along incident beam
    k[0,:] = -ds*s # this is negative x
    # y - towards door
    k[1,:] = -ds*c*n.sin(eta) # CHANGED eta to HFP convention 4-9-2007
    # z - towards roof
    k[2,:] =  ds*c*n.cos(eta)
    return k

                             
def compute_g_vectors(tth, 
                      eta, 
                      omega, 
                      wvln, 
                      wedge = 0.0, 
                      chi   = 0.0):
    """
    Generates spot positions in reciprocal space from 
      twotheta, wavelength, omega and eta
    Assumes single axis vertical
    ... unless a wedge angle is specified
    """
    k = compute_k_vectors(tth, eta, wvln)
#    print k[:,0]
    return compute_g_from_k( k, omega, wedge, chi)

    
def compute_g_from_k( k, omega, wedge=0, chi=0):
    """
    Compute g-vectors with cached k-vectors
    """
    om =radians(omega)
    # G-vectors - rotate k onto the crystal axes
    g = n.zeros((3,k.shape[1]),n.float)
    t = n.zeros((3,k.shape[1]),n.float)
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
        c = n.cos(radians(wedge))
        s = n.sin(radians(wedge))
        t[0,:]= c * k[0,:]           + s * k[2,:]
        t[1,:]=            k[1,:]
        t[2,:]=-s * k[0,:]           + c * k[2,:]
        k=t.copy()
    if chi != 0.0:
        c = n.cos(radians(chi))
        s = n.sin(radians(chi))
        t[0,:]= k[0,:]  
        t[1,:]=        c * k[1,:]    + s * k[2,:]
        t[2,:]=       -s * k[1,:]    + c * k[2,:]
        k=t.copy()
    # This is the reverse rotation (left handed, k back to g)
    g[0,:] = n.cos(om)*k[0,:] + n.sin(om)*k[1,:]    
    g[1,:] =-n.sin(om)*k[0,:] + n.cos(om)*k[1,:]
    g[2,:] =                                     k[2,:]
    return g


def uncompute_g_vectors(g, wavelength, wedge=0.0, chi=0.0):
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
    ds = n.sqrt(n.sum(g*g,0))
    s = ds*wavelength/2.0 # sin theta
    #
    #     k0 = projection of scattering vector along x in lab
    #          must satisfy laue condition
    #        =  - sin(theta)/d
    #
    k0 = n.zeros(g.shape[0],n.float)
    k0 = -ds*s
    #
    # this component: k = W-1 R-1 g
    # k = W-1 ( g0 cos(omega) - g1 sin(omega)  + 0  )
    #         ( g0 sin(omega) + g1 cos(omega)  + 0  )
    #         (       0       +      0         + g2 )
    #
    # k0 = cos(wedge)[g0 cos(omega) - g1 sin(omega)] - sin(wedge) g2
    #
    cw = n.cos(radians(wedge))
    sw = n.sin(radians(wedge))
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
    valid = n.where(t>0.0 , 1, 0)
    t = t * valid
    another_valid = n.where(a < 1e-12, 1, 0)
    a = n.where(another_valid == 1, 1, a)
    sinomega1 = (- b - n.sqrt(t)) / (a * 2.)   # a can clearly be zero
    sinomega2 = (- b + n.sqrt(t)) / (a * 2.)   # a can clearly be zero
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
    valid = n.where(t>0.0 , 1., 0.)
    t = t * valid
    another_valid = n.where(a < 1e-12, 1, 0)
    a = n.where(another_valid == 1, 1, a)
    cosomega1 = (- b - n.sqrt(t)) / (a * 2.)   # a can clearly be zero
    cosomega2 = (- b + n.sqrt(t)) / (a * 2.)   # a can clearly be zero
    cosomega1 = cosomega1*(1-another_valid)
    cosomega2 = cosomega2*(1-another_valid)

    #
    # Now we need to find out which solutions go together via sin^2 + cos^2 = 1
    # ... turns out it could be arctan2(sin1,cos1) or arctan2(sin1,cos2)
    #
    hypothesis_1 = abs(sinomega1*sinomega1+cosomega1*cosomega1 - 1.)
    hypothesis_2 = abs(sinomega1*sinomega1+cosomega2*cosomega2 - 1.)
    same_sign = n.where(abs(hypothesis_1 < hypothesis_2 ), 1, 0)
    omega1 = same_sign * n.arctan2(sinomega1,cosomega1) + \
        (1.-same_sign) * n.arctan2(sinomega1,cosomega2)
    omega2 = same_sign * n.arctan2(sinomega2,cosomega2) + \
        (1.-same_sign) * n.arctan2(sinomega2,cosomega1)
    #
    # Finally re-compute cosomega and sinomega due to the 
    # flipping possibility for getting eta
    cosomega1 = n.cos(omega1)
    sinomega1 = n.sin(omega1)
    cosomega2 = n.cos(omega2)
    sinomega2 = n.sin(omega2)
    #
    # Now we know R-1 and W-1 , so compute k = W-1 R-1 g
    #
    R1g1 = n.zeros(g.shape,n.float)
    R1g2 = n.zeros(g.shape,n.float)
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
    k_one =  n.zeros(g.shape,n.float)
    k_two =  n.zeros(g.shape,n.float)
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
    eta_one = n.arctan2(-k_one[1,:],k_one[2,:])
    eta_two = n.arctan2(-k_two[1,:],k_two[2,:])
    #
    #
    tth = degrees(n.arcsin(s)*2.) * valid
    eta1 = degrees(eta_one)*valid
    eta2 = degrees(eta_two)*valid
    omega1 = degrees(omega1)*valid
    omega2 = degrees(omega2)*valid
    return tth,[eta1,eta2],[omega1,omega2]


def uncompute_one_g_vector(gv,wavelength, wedge=0.0):
    """
    Given g-vectors compute tth,eta,omega
    assert uncompute_g_vectors(compute_g_vector(tth,eta,omega))==tth,eta,omega
    """
    t,e,o = uncompute_g_vectors(
        n.transpose(
            n.array([gv,gv])),
        wavelength,
        wedge=wedge)

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
    u = n.matrixmultiply(C,n.matrixmultiply(W,u))
    u_x_So = cross_product_2x2(u,So)
    # if DEBUG: print "axis orientation",u
    #
    # S = scattered vectors. Length 1/lambda.
    S = n.array([ cos(radians(tth)/2.)*sin(radians(eta))/wavelength,
                cos(radians(tth)/2.)*cos(radians(eta))/wavelength,
                sin(radians(tth)/2.)/wavelength])
    try:    
        S_dot_u_x_So = n.dot(S,u_x_So)
    except:
        print S.shape, u_x_So.shape
    mod_S = n.sqrt(S*S)
    mod_So = n.sqrt(So*So)
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

    tth = n.array([   1,  2,  3,  4,  5,  6,  7,  8,  9, 10], n.float)
    eta = n.array([  10, 40, 70,100,130,160,190,220,270,340], n.float)
    om  = n.array([   0, 20, 40,100, 60,240,300, 20, 42, 99], n.float)

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
