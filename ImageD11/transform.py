
from __future__ import print_function

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
import numpy as np
from ImageD11 import gv_general, cImageD11
from numpy import radians, degrees
import fabio # for LUT

try:
    # crazy debug
    _ = np.arccos(np.zeros(10, float))
except:
    print(dir())
    raise

from math import pi


def cross_product_2x2(a, b):
    """ returns axb for two len(3) vectors a,b"""
    assert len(a) == len(b) == 3
    return np.array([a[1] * b[2] - a[2] * b[1],
                     a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]])




def detector_rotation_matrix(tilt_x, tilt_y, tilt_z):
    """
    Return the tilt matrix to apply to peaks
    tilts are in radians
    typically applied to peaks rotating around beam center
    """
    r1 = np.array([[np.cos(tilt_z), -np.sin(tilt_z), 0],  # note this is r.h.
                   [np.sin(tilt_z), np.cos(tilt_z), 0],
                   [0,    0, 1]], float)
    r2 = np.array([[np.cos(tilt_y), 0, np.sin(tilt_y)],
                   [0, 1,   0],
                   [-np.sin(tilt_y), 0, np.cos(tilt_y)]], float)
    r3 = np.array([[1,          0,       0],
                   [0,  np.cos(tilt_x), -np.sin(tilt_x)],
                   [0,  np.sin(tilt_x), np.cos(tilt_x)]], float)
    r2r1 = np.dot(np.dot(r3, r2), r1)
    return r2r1


def compute_xyz_lab(peaks,
                    y_center=0., y_size=0., tilt_y=0.,
                    z_center=0., z_size=0., tilt_z=0.,
                    tilt_x=0.,
                    distance=0.,
                    # detector_orientation=((1,0),(0,1)),
                    o11=1.0, o12=0.0, o21=0.0, o22=-1.0,
                    **kwds):
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
      
    kwds are not used (but lets you pass in a dict with other things in it)
    """
    assert len(peaks) == 2, "peaks must be a 2D array"
    # Matrix for the tilt rotations
    r2r1 = detector_rotation_matrix(tilt_x, tilt_y, tilt_z)
    # Peak positions in 3D space
    #  - apply detector orientation
    peaks_on_detector = np.array(peaks)
    peaks_on_detector[0, :] = (peaks_on_detector[0, :] - z_center) * z_size
    peaks_on_detector[1, :] = (peaks_on_detector[1, :] - y_center) * y_size
    #
    detector_orientation = [[o11, o12], [o21, o22]]
    # logging.debug("detector_orientation = "+str(detector_orientation))
    flipped = np.dot(np.array(detector_orientation, float),
                     peaks_on_detector)
    #
    vec = np.array([np.zeros(flipped.shape[1]),  # place detector at zero,
                    # sample at -dist
                    flipped[1, :],             # x in search, frelon +z
                    flipped[0, :]], float)     # y in search, frelon -y
    # Position of diffraction spots in 3d space after detector tilts about
    # the beam centre on the detector
    rotvec = np.dot(r2r1, vec)
    # Now add the distance (along x)
    rotvec[0, :] = rotvec[0, :] + distance
    return rotvec


def compute_tth_eta(peaks,
                    y_center=0., y_size=0., tilt_y=0.,
                    z_center=0., z_size=0., tilt_z=0.,
                    tilt_x=0.,
                    distance=0.,
                    # detector_orientation=((1,0),(0,1)),
                    o11=1.0, o12=0.0, o21=0.0, o22=-1.0,
                    t_x=0.0, t_y=0.0, t_z=0.0,
                    omega=None,  # == phi at chi=90
                    wedge=0.0,  # Wedge == theta on 4circ
                    chi=0.0,  # == chi - 90
                    **kwds):  # spare args are ignored
    """
    Finds x,y,z co-ordinates of peaks in the laboratory frame
    Computes tth/eta from these (in degrees)
    
    kwds are not used (left for convenience if you have a parameter dict)
    """
    peaks_xyz = compute_xyz_lab(
        peaks,
        y_center=y_center, y_size=y_size, tilt_y=tilt_y,
        z_center=z_center, z_size=z_size, tilt_z=tilt_z,
        tilt_x=tilt_x,
        distance=distance,
        # detector_orientation=((1,0),(0,1)),
        o11=o11, o12=o12, o21=o21, o22=o22)

    tth, eta = compute_tth_eta_from_xyz(
        peaks_xyz,
        t_x=t_x, t_y=t_y, t_z=t_z,
        omega=omega,
        wedge=wedge,
        chi=chi)

    return tth, eta


def compute_tth_eta_from_xyz(peaks_xyz, omega,
                             t_x=0.0, t_y=0.0, t_z=0.0,
                             #       == phi at chi=90
                             wedge=0.0,  # Wedge == theta on 4circ
                             chi=0.0,  # == chi - 90
                             **kwds):  # last line is for laziness -
    """
    Peaks is a 3 d array of x,y,z peak co-ordinates
    crystal_translation is the position of the grain giving rise to a diffraction spot
    in x,y,z ImageD11 co-ordinates
         x,y is with respect to the axis of rotation (usually also beam centre).
         z with respect to beam height, z centre
    omega data are needed if crystal translations are used
    
    computed via the arctan recipe.
    
    returns tth/eta in degrees
    """
    assert len(peaks_xyz) == 3
    # Scattering vectors
    if omega is None or (t_x == 0. and t_y == 0 and t_z == 0):
        s1 = peaks_xyz
    else:
        # scattering_vectors
        if len(omega) != len(peaks_xyz[0]):
            raise Exception(
                "omega and peaks arrays must have same number of peaks")
        s1 = peaks_xyz - compute_grain_origins(omega, wedge, chi,
                                               t_x, t_y, t_z)
    # CHANGED to HFP convention 4-9-2007
    eta = np.degrees(np.arctan2(-s1[1, :], s1[2, :]))
    s1_perp_x = np.sqrt(s1[1, :] * s1[1, :] + s1[2, :] * s1[2, :])
    tth = np.degrees(np.arctan2(s1_perp_x, s1[0, :]))
    return tth, eta


def compute_sinsqth_from_xyz(xyz):
    """ Computes sin(theta)**2
    x,y,z = co-ordinates of the pixel in cartesian space
    
    if you need grain translations then use func(peaks_xyz - compute_grain_origins(...) )

    seems to be competitive with arctan2 (docs/sintheta_squared_geometry.ipynb)
    
    returns sin(theta)**2 
    """
    # R = hypotenuse of component normal to incident beam (defines x). e.g. y*y+z*z
    R = xyz[1]*xyz[1] + xyz[2]*xyz[2]
    # Q = hypotenuse along the scattered beam, e.g. x*x+y*y+z*z
    Q = xyz[0]*xyz[0] + R
    # if Q == 0 then this is undefined
    sinsqth = 0.5*R/( Q + xyz[0]*np.sqrt(Q) ) 
    return sinsqth


def sinth2_sqrt_deriv(xyz):
    """ sin(theta)**2 from xyz, and derivatives w.r.t x,y,z
    """
    x,y,z = xyz
    R = z*z + y*y
    Q = R + x*x
    SQ = np.sqrt(Q)
    R2 = R/2
    # at x==y==0 this is undefined. ac
    rQ_xSQ = 1/(Q + x*SQ)
    sinth2 = R2*rQ_xSQ
    # some simplification and collecting terms from expressions above to get:
    sr = sinth2*rQ_xSQ
    p = (x / SQ + 2)*sr  # p should be in the range 3sr -> 2sr for x/x to 0/sqrt(R)
    t = (rQ_xSQ - p)     #
    sinth2_dx =  -(SQ*sr+x*p)
    sinth2_dy =   y*t
    sinth2_dz =   z*t
    return sinth2, sinth2_dx, sinth2_dy, sinth2_dz


def compute_xyz_from_tth_eta(tth, eta, omega,
                             t_x=0.0, t_y=0.0, t_z=0.0,
                             #       == phi at chi=90
                             wedge=0.0,  # Wedge == theta on 4circ
                             chi=0.0,  # == chi - 90
                             **kwds):  # last line is for laziness -
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
    xyz = np.zeros((3, tth.shape[0]), float)
    rtth = np.radians(tth)
    reta = np.radians(eta)
    xyz[0, :] =  np.cos(rtth)
    #  eta = np.degrees(np.arctan2(-s1[1, :], s1[2, :]))
    xyz[1, :] = -np.sin(rtth) * np.sin(reta)
    xyz[2, :] =  np.sin(rtth) * np.cos(reta)

    # Find vectors in the fast, slow directions in the detector plane
    pks = np.array([(1, 0),
                    (0, 1),
                    (0, 0) ], float).T
    dxyzl = compute_xyz_lab(pks, **kwds)
    # == [xpos, ypos, zpos] shape (3,n)
    #
    # This was based on the recipe from Thomas in Acta Cryst ...
    #  ... Modern Equations of ...

    ds = dxyzl[:,0] - dxyzl[:,2]  # 1,0 in plane is (1,0)-(0,0)
    df = dxyzl[:,1] - dxyzl[:,2]  # 0,1 in plane
    dO = dxyzl[:,2]               # origin pixel

    # Cross products to get the detector normal
    # Thomas uses an inverse matrix, but then divides out the determinant anyway
    det_norm = np.cross( ds, df )

    # Scattered rays on detector normal
    norm = np.dot( det_norm, xyz )
    # Check for divide by zero
    msk = (norm == 0)
    needmask = False
    if msk.sum()>0:
        norm += msk
        needmask = True

    # Intersect ray on detector plane
    sc = np.dot( np.cross( df, dO ), xyz ) / norm
    fc = np.dot( np.cross( dO, ds ), xyz ) / norm

    if (t_x != 0) or (t_y != 0) or (t_z != 0):
        go = compute_grain_origins(omega,
                                   wedge=wedge, chi=chi,
                                   t_x=t_x, t_y=t_y, t_z=t_z)
        # project these onto the detector face to give shifts
        sct = (  xyz * np.cross( df, go.T ).T ).sum(axis=0) / norm
        fct = (  xyz * np.cross( go.T, ds ).T ).sum(axis=0) / norm
        sc -= sct
        fc -= fct

    if needmask:
        fc = np.where( msk, 0, fc )
        sc = np.where( msk, 0, sc )

    return fc, sc


def compute_grain_origins(omega, wedge=0.0, chi=0.0,
                          t_x=0.0, t_y=0.0, t_z=0.0):
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
    w = np.radians(wedge)
    WI = np.array([[np.cos(w),         0, -np.sin(w)],
                   [0,           1,         0],
                   [np.sin(w),         0,  np.cos(w)]], float)
    c = np.radians(chi)
    CI = np.array([[1,            0,         0],
                   [0,     np.cos(c), -np.sin(c)],
                   [0,     np.sin(c),  np.cos(c)]], float)
    t = np.zeros((3, omega.shape[0]), float)  # crystal translations
    # Rotations in reverse order compared to making g-vector
    # also reverse directions. this is trans at all zero to
    # current setting. gv is scattering vector to all zero
    om_r = np.radians(omega)
    # This is the real rotation (right handed, g back to k)
    t[0, :] = np.cos(om_r) * t_x - np.sin(om_r) * t_y
    t[1, :] = np.sin(om_r) * t_x + np.cos(om_r) * t_y
    t[2, :] = t_z
    if chi != 0.0:
        c = np.cos(np.radians(chi))
        s = np.sin(np.radians(chi))
        u = np.zeros(t.shape, float)
        u[0, :] = t[0, :]
        u[1, :] = c * t[1, :] + -s * t[2, :]
        u[2, :] = s * t[1, :] + c * t[2, :]
        t = u
    if wedge != 0.0:
        c = np.cos(np.radians(wedge))
        s = np.sin(np.radians(wedge))
        u = np.zeros(t.shape, float)
        u[0, :] = c * t[0, :] + -s * t[2, :]
        u[1, :] = t[1, :]
        u[2, :] = s * t[0, :] + c * t[2, :]
        t = u
    return t


def compute_tth_histo(tth, no_bins=100, weight = False, weights = None,
                      **kwds):
    """
    Compute a histogram of tth values
    Uses numpy's histogram rather that doing it by hand as above
    New feature: weight by something (peak intensity for instance), send true for weight and weights values

    Returns a normalised histogram (should make this a probability
    *and*
     For each datapoint, the corresponding histogram weight
 
    Updated and modernized 2021-02-11 S. Merkel
    """
    maxtth = tth.max()
    mintth = tth.min()
    logging.debug("Histogram: maxtth=%f , mintth=%f, bins=%d" % (maxtth, mintth, no_bins))
    if (weight):
        logging.debug("Weighted histogram")
        histogram,binedges = np.histogram(tth, bins=no_bins, weights=weights, density=True)
    else:
        logging.debug("Un-weighted histogram")
        histogram,binedges = np.histogram(tth, bins=no_bins, density=True)
    tthbin = 0.5 *(binedges[:-1] + binedges[1:])
    histogram = histogram/histogram.sum()
    # histogram value for each peak
    # len(hpk) = number of peaks
    # Tried to use numpy's digitize but failed. Edges are treated differently between np.histogram and np.digitize (both are inclusive in np.histogram)
    # Tried many combinations and gave up
    binsize = (maxtth - mintth) / (no_bins-1)
    bins = np.floor((tth - mintth) / binsize).astype(np.int32)
    hpk = np.take(histogram, bins)
    return tthbin, histogram, hpk


def compute_k_vectors(tth, eta, wvln):
    """
    generate k vectors - scattering vectors in laboratory frame
    """
    tth = np.radians(tth)
    eta = np.radians(eta)
    c = np.cos(tth / 2)  # cos theta
    s = np.sin(tth / 2)  # sin theta
    ds = 2 * s / wvln
    k = np.zeros((3, tth.shape[0]), float)
    # x - along incident beam
    k[0, :] = -ds * s  # this is negative x
    # y - towards door
    k[1, :] = -ds * c * np.sin(eta)  # CHANGED eta to HFP convention 4-9-2007
    # z - towards roof
    k[2, :] = ds * c * np.cos(eta)
    return k


def compute_g_vectors(tth,
                      eta,
                      omega,
                      wvln,
                      wedge=0.0,
                      chi=0.0):
    """
    Generates spot positions in reciprocal space from
      twotheta, wavelength, omega and eta
    Assumes single axis vertical
    ... unless a wedge angle is specified
    """
    k = compute_k_vectors(tth, eta, wvln)
    return compute_g_from_k(k, omega, wedge, chi)


def compute_g_from_k(k, omega, wedge=0, chi=0):
    """
    Compute g-vectors with cached k-vectors
    """
    om = np.radians(omega)
    # G-vectors - rotate k onto the crystal axes
    g = np.zeros((3, k.shape[1]), float)
    t = np.zeros((3, k.shape[1]), float)
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
        c = np.cos(np.radians(wedge))
        s = np.sin(np.radians(wedge))
        t[0, :] = c * k[0, :] + s * k[2, :]
        t[1, :] = k[1, :]
        t[2, :] = -s * k[0, :] + c * k[2, :]
        k = t.copy()
    if chi != 0.0:
        c = np.cos(np.radians(chi))
        s = np.sin(np.radians(chi))
        t[0, :] = k[0, :]
        t[1, :] = c * k[1, :] + s * k[2, :]
        t[2, :] = -s * k[1, :] + c * k[2, :]
        k = t.copy()
    # This is the reverse rotation (left handed, k back to g)
    g[0, :] = np.cos(om) * k[0, :] + np.sin(om) * k[1, :]
    g[1, :] = -np.sin(om) * k[0, :] + np.cos(om) * k[1, :]
    g[2, :] = k[2, :]
    return g


def uncompute_g_vectors(g, wavelength, wedge=0.0, chi=0.0):
    """
    Given g-vectors compute tth,eta,omega
    assert uncompute_g_vectors(compute_g_vector(tth,eta,omega))==tth,eta,omega
    """
    if wedge == chi == 0:
        post = None
    else:
        post = gv_general.wedgechi( wedge=wedge, chi=chi )
    omega1, omega2, valid = gv_general.g_to_k(
        g, wavelength,axis=[0,0,-1], pre=None, post=post )
    # we know g, omega. Compute k as ... ?
    if post is None:
        pre = None
    else:
        pre = gv_general.chiwedge( wedge=wedge, chi=chi ).T
    k_one = gv_general.k_to_g( g, omega1, axis=[0,0,1],
                               pre = pre, post=None)
    k_two = gv_general.k_to_g( g, omega2, axis=[0,0,1],
                               pre = pre, post=None)
    #
    # k[1,:] = -ds*c*sin(eta)
    # ------    -------------   .... tan(eta) = -k1/k2
    # k[2,:] =  ds*c*cos(eta)
    #
    eta_one = np.arctan2(-k_one[1, :], k_one[2, :])
    eta_two = np.arctan2(-k_two[1, :], k_two[2, :])
    #
    #
    ds = np.sqrt(np.sum(g * g, 0))
    s = ds * wavelength / 2.0  # sin theta
    tth = np.degrees(np.arcsin(s) * 2.) * valid
    eta1 = np.degrees(eta_one) * valid
    eta2 = np.degrees(eta_two) * valid
    omega1 = omega1 * valid
    omega2 = omega2 * valid
    return tth, [eta1, eta2], [omega1, omega2]


def uncompute_one_g_vector(gv, wavelength, wedge=0.0):
    """
    Given g-vectors compute tth,eta,omega
    assert uncompute_g_vectors(compute_g_vector(tth,eta,omega))==tth,eta,omega
    """
    t, e, o = uncompute_g_vectors(
        np.transpose(
            np.array([gv, gv])),
        wavelength,
        wedge=wedge)

    return t[0], [e[0][0], e[1][0]], [o[0][0], o[1][0]]


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
    So = [1. / wavelength, 0, 0]
    #
    # u vector along rotation axis
    # starts as along z
    u = [0, 0, 1]
    # rotate by
    # g =  R . W . k where:
    # R = ( cos(omega) , sin(omega), 0 )
    #     (-sin(omega) , cos(omega), 0 )
    #     (         0  ,         0 , 1 )
    #
    W = [[np.cos(wedge),  0,  np.sin(wedge)],
         [0,  1,          0],
         [-np.sin(wedge),  0,  np.cos(wedge)]]
    #
    C = [[1,         0,      0],
         [0,  np.cos(chi), np.sin(chi)],
         [0, -np.sin(chi), np.cos(chi)]]
    u = np.dot(C, np.dot(W, u))
    u_x_So = cross_product_2x2(u, So)
    # if DEBUG: print "axis orientation",u
    #
    # S = scattered vectors. Length 1/lambda.
    S = np.array([np.cos(np.radians(tth) / 2.) * np.sin(np.radians(eta)) / wavelength,
                  np.cos(np.radians(tth) / 2.) * np.cos(np.radians(eta)) / wavelength,
                  np.sin(np.radians(tth) / 2.) / wavelength])
    try:
        S_dot_u_x_So = np.dot(S, u_x_So)
    except:
        print(S.shape, u_x_So.shape)
    mod_S = np.sqrt(S * S)
    mod_So = np.sqrt(So * So)
    try:
        lorentz = abs(S_dot_u_x_So) / mod_S / mod_So
    except:
        raise Exception("Please fix this div0 crap in lorentz")
    return lorentz


def compute_polarisation_factors(args):
    """
    From Kabsch 1988 J. Appl. Cryst. 21 619

    DIVIDE the intensities by:
    <sin2 psi> = (1 - 2p) [ 1 - (n.S/|S|^2) ] + p { 1 + [S.S_0/(|S||S_0|)^2]^2}

    p = degree of polarisation (sync = 1, tube = 0.5 , mono + tube in between)
        or "probability of finding electric field vector in plane having
            normal, n"
    S = scattered vector
    S_0 = incident vector
    n = normal to polarisation plane, typically perpendicular to S_0

    In ImageD11 we normally expect to find:
    x axis along the beam
    z axis being up, and parallel to the normal n mentioned above
    """
    n = [0, 0, 1]

class Ctransform(object):
    pnames = ( "y_center", "z_center", "y_size", "z_size",
               "distance", "wavelength","omegasign",
               "tilt_x","tilt_y","tilt_z",
               "o11", "o12", "o21", "o22",
               "wedge", "chi" )
    def __init__(self, pars ):
        
        """ To Do ...:
                origin = xyz(0,0), ds = xyz(1,0), df = xyz(0,1)
                xyz(s,f) = origin + ds*s + df*f
        """
        self.pars={}
        for p in self.pnames:
            self.pars[p] = pars[p] # copy
        self.reset()
        
    def reset(self):
        p = self.pars
        self.distance_vec = np.array( (p['distance'],0.,0.))
        self.dmat = detector_rotation_matrix( p['tilt_x'], p['tilt_y'], p['tilt_z'])
        self.fmat = np.array( [[ 1,   0,   0],
                               [ 0, p['o22'], p['o21']],
                               [ 0, p['o12'], p['o11']]] )
        self.rmat = np.dot(self.dmat, self.fmat).ravel()
        self.cen = np.array( ( p["z_center"], p["y_center"], p["z_size"], p["y_size"] ))
                                
    def sf2xyz(self, sc, fc, out=None):
        """
        computes the xlab, ylab, zlab co-ordinates of spots on the detector
        input: sc, fc = corrected spot positions, shape(n)
        output: xlylzl array shape (n,3)
        """
        assert len(sc) == len(fc)
        if out is None:
            out = np.empty( (len(sc),3), float)
        # Not used! t = np.array( (tx,ty,tz) )
        cImageD11.compute_xlylzl( sc, fc, self.cen, self.rmat, self.distance_vec, out)
        return out
    
    def xyz2gv(self, xyz, omega, tx=0, ty=0, tz=0, out=None ):
        """
        computes g vectors from xlylzl co-ordinates on the detector
        input:  xlylzl array shape (n,3)
        output: gve shape (n,3)
        """
        assert len(omega) == len(xyz)
        if out is None:
            out = np.empty( (len(xyz),3), float)
        cImageD11.compute_gv( xyz,
                             omega, 
                             self.pars['omegasign'], 
                             self.pars['wavelength'], 
                             self.pars['wedge'], 
                             self.pars['chi'], 
                             np.array((tx,ty,tz)), 
                             out)
        return out
    
    def sf2gv( self, sc, fc, omega, tx=0, ty=0, tz=0, out=None ):
        """
        computes g vectors in the sample frame
        input sc, fc, omega = spot coordinates, shape(n)
        tx, ty, tz = scalar grain centre of mass
        returns: gve shape(n,3)
        """
        xyz = self.sf2xyz( sc, fc )
        return self.xyz2gv( xyz, omega, tx, ty, tz, out )

    def xyz2geometry( self, xyz, omega, tx=0, ty=0, tz=0, out=None):
        """
        Takes x,y,z laborator spot positions.
        Computes everything else needed for updateGeometry in columnfile
        """
        assert len(omega) == len(xyz)
        assert xyz.shape[1] == 3
        if out is None:
            out = np.empty( (len(xyz), 6), float )
        # (xlylzl,omega,omegasign,wvln,wedge,chi,t,out,ng)
        # compute_geometry(xlylzl,omega,omegasign,wvln,wedge,chi,t,out,ng)
        cImageD11.compute_geometry( xyz,
                                   omega,
                                   self.pars['omegasign'],
                                   self.pars['wavelength'],
                                   self.pars['wedge'],
                                   self.pars['chi'],
                                   np.array((tx, ty, tz)),
                                   out )
        return out
        
class PixelLUT( object ):

    """ A look up table for a 2D image to store pixel-by-pixel values
    """
    
    # parameters that can be used to create this LUT
    pnames = ( "y_center", "z_center", "y_size", "z_size",
               "distance", "wavelength", "omegasign",
               "tilt_x","tilt_y","tilt_z",
               "o11", "o12", "o21", "o22",
               "wedge", "chi", "dxfile", "dyfile", "spline", "shape" )
                 
    def __init__( self, pars ):
        """
        pars is a dictionary containing the calibration parameters
        """
        self.pars = {}
        for p in self.pnames:
            if p in pars:
                self.pars[p] = pars[p] # make a copy
        if 'dxfile' in pars:
            # slow/fast coordinates on image at pixel centers
            self.df = fabio.open( pars['dxfile'] ).data
            self.ds = fabio.open( pars['dyfile'] ).data
            self.shape = s = self.ds.shape # get shape from file
            self.pars['shape'] = s
            slow, fast = np.mgrid[ 0:s[0], 0:s[1] ]
            self.sc = slow + self.ds
            self.fc = fast + self.df
        elif 'spline' in pars: # need to test this...
            from ImageD11 import blobcorrector
            b = blobcorrector.correctorclass( self.pars['spline'] )
            s = int(b.ymax - b.ymin), int(b.xmax - b.xmin)
            if 'shape' in self.pars:      # override. Probabl
                s = self.pars['shape'] 
            self.shape = s
            self.fc, self.sc = b.make_pixel_lut( s ) 
            slow, fast = np.mgrid[ 0:s[0], 0:s[1] ]
            self.df = self.fc - fast
            self.ds = self.sc - slow
        else:
            s = self.shape = self.pars.get("shape")
            self.sc, self.fc = np.mgrid[0:s[0], 0:s[1]]
            self.df = None 
            self.ds = None
    
        self.xyz = compute_xyz_lab( (self.sc.ravel(), self.fc.ravel()), **self.pars )
        self.sinthsq = compute_sinsqth_from_xyz( self.xyz )
        self.tth, self.eta = compute_tth_eta_from_xyz(self.xyz, None, **self.pars)
        # scattering angles:
        self.tth, self.eta = compute_tth_eta( (self.sc.ravel(), self.fc.ravel()), **self.pars )
        # scattering vectors:
        self.k = compute_k_vectors( self.tth, self.eta, self.pars.get('wavelength') )
        self.sinthsq.shape = s
        self.tth.shape = s
        self.eta.shape = s
        self.k.shape = (3, s[0], s[1])
        self.xyz.shape = (3, s[0], s[1])
    
    def spatial(self, sraw, fraw):
        """ applies a spatial distortion to sraw, fraw (for peak centroids) """
        if self.df is None:
            return sraw, fraw
        else:
            si = np.round(sraw.astype(int)).clip( 0, self.shape[1] - 1 )
            fi = np.round(fraw.astype(int)).clip( 0, self.shape[1] - 1 )
            sc = sraw + self.ds[ si, fi ]
            fc = fraw + self.df[ si, fi ]
            return sc, fc
        
    def __repr__(self):
        """ print yourself in a way we can use for eval """
        sp = "\n".join( [ "%s : %s,"%(repr(p), repr(self.pars[p])) for p in self.pnames
                         if p in self.pars ] )
        return "PixelLUT( { %s } )"%(sp)    
    
    
    
if __name__ == "__main__":
    # from indexing import mod_360
    def mod_360(theta, target):
        """
        Find multiple of 360 to add to theta to be closest to target
        """
        diff = theta - target
        while diff < -180:
            theta = theta + 360
            diff = theta - target
        while diff > 180:
            theta = theta - 360
            diff = theta - target
        return theta

    tth = np.array([1,  2,  3,  4,  5,  6,  7,  8,  9, 10], float)
    eta = np.array([10, 40, 70, 100, 130, 160, 190, 220, 270, 340], float)
    om = np.array([0, 20, 40, 100, 60, 240, 300, 20, 42, 99], float)

    for wavelength in [0.1, 0.2, 0.3]:
        for wedge in [-10., -5., 0., 5., 10.]:
            print("Wavelength", wavelength, "wedge", wedge)
            print("tth, eta, omega   ...   " +\
                  "tth, eta, omega   ...   " +\
                  "tth, eta, omega")
            gv = compute_g_vectors(tth, eta, om, wavelength, wedge)
            t, e, o = uncompute_g_vectors(gv, wavelength, wedge)
            for i in range(tth.shape[0]):
                print("%9.3f %9.3f %9.3f  " % (tth[i], eta[i], om[i]), end=' ')
                print("%9.3f %9.3f %9.3f  " % (t[i],
                                               mod_360(e[0][i], eta[i]),
                                               mod_360(o[0][i], om[i])), end=' ')
                print("%9.3f %9.3f %9.3f  " % (t[i],
                                               mod_360(e[1][i], eta[i]),
                                               mod_360(o[1][i], om[i])), end=' ')
                # Choose best fitting
                e_eta1 = mod_360(e[0][i], eta[i]) - eta[i]
                e_om1 = mod_360(o[0][i], om[i]) - om[i]
                score_1 = e_eta1 * e_eta1 + e_om1 * e_om1
                e_eta2 = mod_360(e[1][i], eta[i]) - eta[i]
                e_om2 = mod_360(o[1][i], om[i]) - om[i]
                score_2 = e_eta2 * e_eta2 + e_om2 * e_om2
                print("%.5g %.5g" % (score_1, score_2))
