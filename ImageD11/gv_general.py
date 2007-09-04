


# ImageD11 Software for beamline ID11
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
Find a generalised way to compute g-vectors backwards and forwards
"""
import math, logging, LinearAlgebra
import Numeric as n
from ImageD11.transform import degrees

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
        self.direction = n.array(direction)
        assert self.direction.shape == (3,) , "direction.shape != 3, is it a vector??"
        mag = n.dot(self.direction, self.direction)
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
        p = n.array(vectors, n.Float)
        assert p.shape[0] == 3
        if angles is None:
            self.to_matrix()
            return n.matrixmultiply(self.matrix, p)
        q = n.array(angles)*math.pi/180.
        assert q.shape[0] == p.shape[1]
        rp = p * n.cos(q)
        assert rp.shape == p.shape
        a_dot_p = n.matrixmultiply( self.direction, p)
        apa = n.outerproduct(self.direction, a_dot_p)
        rp += apa*(1-n.cos(q))
        rp += n.sin(q)*n.transpose( n.cross_product( self.direction, n.transpose(p)) )
        return rp
    
    def rotate_vectors_inverse(self, vectors, angles = None):
        """
        Same as rotate vectors, but opposing sense
        """
        p = n.array(vectors, n.Float)
        assert p.shape[0] == 3
        if angles is None:
            self.to_matrix()
            return n.matrixmultiply(self.inversematrix, p)
        # Just reverse the angles here
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
        e = n.array([self.direction]*3)
        w = n.transpose(n.array( [ [ 0, -dz, dy ] , 
                                   [dz, 0, -dx] , 
                                   [-dy, dx, 0] ], n.Float))
        st = math.sin( math.radians( self.angle ))
        ct = math.cos( math.radians( self.angle ))
        self.matrix = n.identity(3, n.Float)*ct - st * w  + \
                      (1 - ct)*e*n.transpose(e)
	self.inversematrix = LinearAlgebra.inverse(self.matrix)
        return self.matrix

def axis_from_matrix( m ):
    """
    Factory method to create a rotation_axis object from a direction cosine matrix
    http://en.wikipedia.org/wiki/Rotation_representation_%28mathematics%29
    """
    # should check for m being a pure rotation
    d = LinearAlgebra.determinant(m)
    assert abs(d - 1.0)<1e-6, "pure rotation? det(m) = %f"%(d)
    arg = (m[0,0]+m[1,1]+m[2,2]-1.)/2.0
    if arg == 1:
        # identity matrix - oh bugger - direction is undefined
        angle_rad = 0.0
        dir = n.array([0,0,1],n.Float)
    else:
        angle_rad = math.acos(arg)  
        dir = n.array([m[2,1] - m[1,2],
                       m[0,2] - m[2,0],
                       m[1,0] - m[0,1] ], n.Float )
    o = rotation_axis( dir , math.degrees( angle_rad ) )
    assert o.matrix == m
    return o

rotate_identity = rotation_axis( n.array([0,0,1],n.Float) , 0.0 )


def k_to_g( k , # scattering vectors in the laboratory
            angles , # eg samome in degrees
            axis = rotation_axis( n.array([0,0,1],n.Float) , 0.0 ) , # eg z axis 
            pre = rotate_identity , 
            post = rotate_identity ):
    """
    Computes g = pre . rot(axis, angle) . post . k
    Typically in ImageD11 pre = [wedge][chi] and post = identity
    """
    assert k.shape[1] == angles.shape[0] , "Number of vectors and scan axis must be same"+str(k.shape)+str(angles.shape)
    assert k.shape[0] == 3 , "k must be a [:,3] array"
    #print "angles=",angles[0]
    #print "k.shape\n",k.shape,k[:,0]
    pk =  post.rotate_vectors_inverse( k )
    #print "pk.shape\n" ,pk.shape,pk[:,0]
    rpk = axis.rotate_vectors_inverse(pk , angles)
    #print "rpk.shape\n",rpk.shape,rpk[:,0]
    g =   pre.rotate_vectors_inverse( rpk )
    #print "g.shape\n",g.shape,g[:,0]
    return g
    
    
def g_to_k( g,  # g-vectors [3,:]
            wavelength, # Needed - depends on curve of Ewald sphere!
            axis = n.array([0,0,1], n.Float),
            pre = None,
            post = None ):
    """
    Get the k and angles given the g in 
         g = pre . rot(axis, angle) . post . k
    ...where k will satisfy the Laue equations
    
    The recipe is:
        Find the components of g along and perpendicular to the rotation axis
            co-ords are a0 = axis, a1 = axis x g, a2 = a1 x axis = ( axis x g ) x axis
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
        rg = pre.rotate_vectors_inverse(g)
    else:
        rg = g
    assert rg.shape == g.shape
    beam = n.zeros(rg.shape, n.Float)
    beam[0,:] = -1.
    beam[1,:] = 0.
    beam[2,:] = 0.
    if post is not None:
        rb = post.rotate_vectors(beam)
    else:
        rb = beam
    assert rb.shape == g.shape
    # Find the components of g with respect to our rotation axis
    # a1 = perpendicular to both axis and g
    a1 = n.transpose(n.cross_product(axis, g, axis1=0, axis2=0))
    # a2 perpendicular to axis, along g
    a2 = n.transpose(n.cross_product(a1, axis, axis1=0, axis2=0))
    # projection of g along axis
    a0 = g - a2
    assert a0.shape == a1.shape == a2.shape == g.shape
    # Dot product with incident beam
    rbda0 = n.sum(rb * a0, 0) 
    rbda1 = n.sum(rb * a1, 0) 
    rbda2 = n.sum(rb * a2, 0)
    assert rbda0.shape == rbda1.shape == rbda2.shape == (g.shape[1],)
    modg = n.sqrt(n.sum(g * g, 0))
    sin_theta = modg*wavelength/2.
    kdotbeam = -sin_theta*modg
    # k.b = rbda0 + rbda1.sin(t) + rbda2.cos(t)
    a = rbda1
    b = rbda2
    c = kdotbeam - rbda0
    # From wikipedia: 
    # http://en.wikipedia.org/wiki/List_of_trigonometric_identities#Linear_combinations
    # a.sin(x) + b.cos(x) = sqrt(a*a+b*b) sin(x+p)
    # with p = atan(b/a)
    p = n.arctan2(b,a)
    k = n.sqrt(a*a+b*b)
    quot = c/k
    valid = n.logical_and( n.greater_equal( quot, -1) ,
                           n.less_equal(    quot,  1) )
    quot = n.where( valid, quot, 0 ) 
    x_plus_p = n.arcsin(quot)
    sol1 = x_plus_p + p
    sol2 = math.pi - x_plus_p + p

    return degrees(sol1), degrees(sol2), valid
    
    
