


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
import math
import Numeric as n

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
        self.angle = angle
        self.matrix = self.to_matrix()

    def rotate_vectors(self, vectors, angles = None):
        """
        Given a list of vectors, rotate them to new ones
        Angle is either self.angle, or 1 per vector, in degrees
        
        http://mathworld.wolfram.com/RotationFormula.html
        r' = r cos(t) + n(n.r)(1-cos(t)) + rxn sin(t)
        """
        v = n.array(vectors,n.Float)
        assert v.shape[1] == 3
        if angles is None:
            self.to_matrix()
            return n.transpose(n.matrixmultiply(self.matrix, n.transpose(v)))
        a = n.array(angles)*math.pi/180.
        assert a.shape[0] == v.shape[0]
        rp = n.cos(a)*n.transpose(v)
        t = n.matrixmultiply( v, self.direction)
        rp += t * (1-n.cos(a)) * n.transpose(v) 
        rp += n.sin(a)*n.transpose( n.cross_product(self.direction, v) )
        return n.transpose(rp)
        
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
        w = n.array( [ [ 0, -dz, dy ] , [dz, 0, -dx] , [-dy, dx, 0] ], n.Float)
        st = math.sin( math.radians( self.angle ))
        ct = math.cos( math.radians( self.angle ))
        self.matrix = n.identity(3, n.Float)*ct - st * w  + \
                      (1 - ct)*e*n.transpose(e)
        return self.matrix
        
        
def axis_from_matrix( m ):
    """
    Factory method to create a rotation_axis object from a direction cosine matrix
    http://en.wikipedia.org/wiki/Rotation_representation_%28mathematics%29
    """
    angle_rad = math.acos((m[0,0]+m[1,1]+m[2,2]-1.)/2.0)   
    ts = 2.0 * math.sin(angle_rad)
    dir = n.array([ m[1,2] - m[2,1],
                    m[2,0] - m[0,2],
                    m[0,1] - m[1,0] ], n.Float )
    o = rotation_axis( dir , math.degrees( angle_rad ) )
    assert o.matrix == m
    return o

def k_to_g( k , # scattering vectors in the laboratory
            angles , # eg samome in degrees
            axis = n.array([0,0,1],n.Float) , # eg z axis 
            pre = n.identity(3) , 
            post = n.identity(3) ):
    """
    Computes g = pre . rot(axis, angle) . post . k
    Typically in ImageD11 pre = [wedge][chi] and post = identity
    """
    assert len(k) == len(angle) , "Number of vectors and scan axis must be same"
    pk = matrixmultiply( post , k )
    r = rotation_axis(axis, 0.)
    rpk = j
    g = matrixmultiply( pre , rpk )
