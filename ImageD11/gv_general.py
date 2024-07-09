
from __future__ import print_function


## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py




# ImageD11 Software for beamline ID11
# Copyright (C) 2005-2007  Jon Wright
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
Find a generalised way to compute g-vectors backwards and forwards
"""
import math, logging
from numpy.linalg import det, inv
import numpy as np

# Use python -v myscript.py args
# print("gv_general from ",__file__)

def angmod(a):
    return np.arctan2( np.sin(a), np.cos(a) )

class rotation_axis:
    """
    Represents an axis - so far includes:
        axis/angle
        matrix (direction cosines)

    potentially can add quaternions and euler etc later?
    """
    def __init__(self, direction, angle = 0.0):
        """
        direction is the rotation direction
        angle is the angle of rotation (degrees)
        """
        self.direction = np.asarray(direction)
        assert self.direction.shape == (3,) , \
            "direction.shape != 3, is it a vector??"
        mag = np.dot(self.direction, self.direction)
        if abs(mag - 1.0) > 1e-5 :
            self.direction = self.direction / mag
            logging.warning("Non-normalised direction vector "+str(direction))
        self.angle = angle
        self.matrix = self.to_matrix()

    def rotate_vectors(self, vectors, angles = None):
        """
        Given a list of vectors, rotate them to new ones
        Angle is either self.angle, or 1 per vector, in degrees
        
        http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/rotation.html
        p' = a . p a + (p - a . p a) cos q + a * p sin q
           = p cos q + a . p a (1-cos q) + a * p sin q 
           point p, axis a, angle q 
           
           
        http://mathworld.wolfram.com/RotationFormula.html
        r' = r cos(t) + n(n.r)(1-cos(t)) + rxn sin(t)
        """
        p = np.array(vectors, float)
        assert p.shape[0] == 3
        if angles is None:
            self.to_matrix()
            return np.dot(self.matrix, p)
        q = np.radians(angles)
        assert q.shape[0] == p.shape[1]
        rp = p * np.cos(q)
        assert rp.shape == p.shape
        a_dot_p = np.dot( self.direction, p)
        apa = np.outer(self.direction, a_dot_p)
        rp += apa*(1-np.cos(q))
        rp += np.sin(q)*np.transpose( np.cross( 
                self.direction, np.transpose(p)) )
        return rp
    
    def rotate_vectors_inverse(self, vectors, angles = None):
        """
        Same as rotate vectors, but opposing sense
        """
        p = np.array(vectors, float)
        assert p.shape[0] == 3
        if angles is None:
            self.to_matrix()
            return np.dot(self.inversematrix, p)
        # Just reverse the angles here
        else:
            return self.rotate_vectors(vectors, -angles)

    def to_matrix(self):
        """
        Returns a 3x3 rotation matrix
        R = I3cos(t) + ww^T (1-cos(t)) - W sin(t)
        t = angle
        W = vector with hat operator = [ 0  -wz  wy 
                                         wz  0  -wx
                                        -wy  wx  0  ]
        """ 
        # FIXME - check for caching
        dx, dy, dz = self.direction
        e = np.array([self.direction]*3)
        w = np.transpose(np.array( [ [ 0, -dz, dy ] , 
                                   [dz, 0, -dx] , 
                                   [-dy, dx, 0] ], float))
        st = math.sin( math.radians( self.angle ))
        ct = math.cos( math.radians( self.angle ))
        self.matrix = np.identity(3, float)*ct - st * w  + \
                      (1 - ct)*e*np.transpose(e)
        self.inversematrix = inv(self.matrix)
        return self.matrix

def axis_from_matrix( m ):
    """
    Factory method to create a rotation_axis object from a 
    direction cosine matrix
    http://en.wikipedia.org/wiki/Rotation_representation_%28mathematics%29
    """
    # should check for m being a pure rotation
    d = det(m)
    assert abs(d - 1.0)<1e-6, "pure rotation? det(m) = %f"%(d)
    arg = (m[0,0]+m[1,1]+m[2,2]-1.)/2.0
    if arg == 1:
        # identity matrix - oh bugger - direction is undefined
        angle_rad = 0.0
        direc = np.array([0,0,1],float)
    else:
        angle_rad = math.acos(arg)  
        direc = np.array([m[2,1] - m[1,2],
                         m[0,2] - m[2,0],
                         m[1,0] - m[0,1] ], float )
        direc = direc / np.sqrt(np.dot(direc, direc))
    o = rotation_axis( direc , math.degrees( angle_rad ) )
    if not (abs(o.matrix - m) < 1e-5).all():
        print("o.matrix\n",o.matrix)
        print("m\n",m)
        raise Exception("error in axis_from_matrix")
    return o

rotate_identity = rotation_axis( np.array([0,0,1],float) , 0.0 )


def k_to_g( k , # scattering vectors in the laboratory
            angles , # eg samome in degrees
            axis = None,
            # eg z axis 
            pre = None , 
            post = None ):
    """
    Computes g = pre . rot(axis, angle) . post . k
    Typically in ImageD11 post = [wedge][chi] and pre = identity
       since g = Omega.Chi.Wedge.k
       or    k = Wedge-1.Chi-1.Omega-1.g
       that is to say omega is directly under the sample and 
       wedge is at the bottom of the stack
    """
    if axis is None:
        raxis = rotation_axis( [0,0,1] )
    else:
        raxis = rotation_axis( axis )
    assert k.shape[1] == angles.shape[0] , \
        "Number of vectors and scan axis must be same"+\
        str(k.shape)+str(angles.shape)
    assert k.shape[0] == 3 , "k must be a [:,3] array"
    if post is not None:
        pk =  np.dot(post, k )
        rpk = raxis.rotate_vectors(pk , angles)
    else:
        rpk = raxis.rotate_vectors(k , angles)
    if pre is not None:
        return np.dot(pre, rpk )
    else:
        return rpk

    
    
def g_to_k( g,  # g-vectors [3,:]
            wavelength, # Needed - depends on curve of Ewald sphere!
            axis = np.array([0,0,1], float),
            pre = None,
            post = None ):
    """
    Get the k and angles given the g in 
         g = pre . rot(axis, angle) . post . k
         g = omega . chi . wedge . k
    ...where k will satisfy the Laue equations
    
    The recipe is:
        Find the components of g along and perpendicular to the rotation axis
            co-ords are a0 = axis, 
                        a1 = axis x g, 
                        a2 = a1 x axis = ( axis x g ) x axis
        Rotated vector will be a0 + a1 sin(angle) + a2 cos(angle)
        Laue condition says that [incident beam].[k] = k.sin(theta)
                                                     = 2 sin^2(theta) / lambda 
                                                     = sin(t)/d = |g|.sin(t)
                                                     = |g|*|g|*lambda / 2
        Apply any post rotations to the incident beam 
        Apply any pre-rotations to the g vector
        |g| = [a0 + a1 sin(angle) + a2 cos(angle) ] . [incident beam]
        => solve for angle

        http://www.ping.be/~ping1339/gonio.htm
        a sin(u) + b cos(u) = c
        let: tan(u') = -b/a
        and: A = a / cos(u')
        Finally you'll find:
             A.sin(u-u') = c
             A.cos(pi/2 - u - u') = c
    """
    assert g.shape[0] == 3
    npeaks = g.shape[1]
    # First deal with the pre and post rotations
    if pre is not None:
        # rg = pre.rotate_vectors_inverse(g)
        rg = np.dot( pre, g )
    else:
        rg = g
    assert rg.shape == g.shape
    # beam = np.zeros(rg.shape, float)
    # beam[0,:] = -1./wavelength
    # beam[1,:] = 0.
    # beam[2,:] = 0.
    beam = np.array( [ -1.0/wavelength, 0, 0] , float )
    if post is not None:
        # rb = post.rotate_vectors(beam)
        rb = np.dot( post.T, beam )
    else:
        rb = beam
    assert rb.shape == beam.shape
    # Find the components of g with respect to our rotation axis
    # a1 = perpendicular to both axis and g
    a1 = np.transpose(np.cross(axis, rg.T))
    # a2 perpendicular to axis, along g
    a2 = np.transpose(np.cross(a1.T, axis))
    # projection of g along axis
    a0 = rg - a2
    assert a0.shape == a1.shape == a2.shape == g.shape
    # Dot product with incident beam
    rbda0 = np.sum(rb[:,np.newaxis] * a0, 0)
    rbda1 = np.sum(rb[:,np.newaxis] * a1, 0)
    rbda2 = np.sum(rb[:,np.newaxis] * a2, 0)
    assert rbda0.shape == rbda1.shape == rbda2.shape == (npeaks,)
    modg = np.sqrt(np.sum(g * g, 0))
    kdotbeam = -modg*modg/2.
    # print kdotbeam,"uyou"
    # k.b = rbda0 + rbda1.sin(t) + rbda2.cos(t)
    # a = rbda1
    # b = rbda2
    # c = kdotbeam - rbda0
    # From wikipedia: 
    # http://en.wikipedia.org/wiki/List_of_trigonometric_identities#Linear_combinations
    # a.sin(x) + b.cos(x) = sqrt(a*a+b*b) sin(x+p)
    # with p = atan(b/a)
    phi = np.arctan2(rbda2,rbda1)
    den = np.sqrt(rbda1*rbda1+rbda2*rbda2)
    msk = (den <= 0)
    quot = (kdotbeam - rbda0)/(den + msk)
    valid = (~msk) & ( quot >= -1) & ( quot <= 1)
    quot = np.where( valid, quot, 0 ) 
    x_plus_p = np.arcsin(quot)
    sol1 = x_plus_p + phi
    sol2 = math.pi - x_plus_p + phi
    # k 
    return np.degrees(angmod(sol1)), np.degrees(angmod(sol2)), valid


def wedgemat(w):
    cw = np.cos(np.radians(w))
    sw = np.sin(np.radians(w))
    W = np.array([[ cw ,  0  , sw ],
                  [ 0  ,  1  , 0  ],
                  [-sw ,  0  , cw ]])
    return W

def chimat(c):
    cc = np.cos(np.radians(c))
    sc = np.sin(np.radians(c))
    C = np.array([[ 1  ,   0  ,  0   ],
                  [ 0  ,  cc  , sc   ],
                  [ 0  , -sc  , cc   ]])
    return C

    
def wedgechi(wedge=0., chi=0.):
    """
    in transform:
    g = omega . chi . wedge .k
    """
    return np.dot( wedgemat(wedge), chimat(chi) )

def chiwedge(chi=0., wedge=0.):
    """
    in transform:
    g = omega . chi . wedge .k
    """
    return np.dot( chimat(chi), wedgemat(wedge) )


