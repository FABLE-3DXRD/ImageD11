
from __future__ import print_function




# ImageD11 Software for beamline ID11
# Copyright (C) 2008  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  0211-1307  USA

"""
If you google for "fit plane points" you can find similar things. This 
uses an SVD to fit the plane to the points. The story goes:

    center of mass is in the plane
    vectors from the c.o.m to the points are in a space
    svd finds principal vectors for that space
    null vector (or smallest) is perpendicular to plane holding points

Then plane computations go via n.(r-c)=0

Enjoy!
"""


# Fit a least squares plane to some data
# Data are likely to be the edge of a box from an integration
# ... but could be anything (eg: g-vectors etc)

import numpy as np

class plane:

    def __init__(self, normal, point):
        """
        Normal - vector normal to the surface
        point  - a point on the surface

        Equation of plane is normal.(x - point) = 0
        """
        assert len(point) == len(normal) == 3, "3D space please"
        assert abs (np.dot(normal , normal) - 1) < 1e-6, \
                "Normalise normal please"
        self.point = point
        self.normal = normal

    def calc(self, points, index=2):
        """
        points - [x, y, z]
        index  - compute x or y or z
        """
        assert index in [0,1,2],"3D space please"
        assert isinstance(points, np.ndarray),"Please use numpy array for points"
        assert points.shape[1] == 3
        assert abs(self.normal[index]) > 1e-6, "Sorry, your plane is parallel"
        points[:, index ] = 0.
        # Compute vectors from self.point to points
        vecs = points - self.point
        # Projections along normal
        proj = np.dot(self.normal, vecs.T)
        # Now determine points[index, :] such that proj == 0 everywhere
        points[:, index ] = -proj/self.normal[index]
        return points



def fitplane(x, y, z):
    """
    x, y and z are the co-ordinates of the points to fit to
    """
    assert len(x) == len(y) == len(z)
    cx = np.mean( x ) 
    cy = np.mean( y ) 
    cz = np.mean( z ) 
    skinny = np.array( [np.asarray(x) - cx ,
                        np.asarray(y) - cy ,
                        np.asarray(z) - cz ] )     
    u, s, vh = np.linalg.svd( skinny , full_matrices = 0 )
    # print "Through",cx,cy,cz,"Normal is u",u[:,-1]
    return plane( u[:,-1].copy(), np.array([cx, cy, cz]))


def test():
    p = fitplane( [ 1 , 0 , 0 ]    ,
                  [ 0 , 1 , 0 ]    ,
                  [ 0 , 0 , 1 ]    )
    for i in range(3):
        assert abs(abs(p.normal[i]) - 1/np.sqrt(3))< 1e-6
    assert np.allclose( p.normal[0], p.normal[1])
    assert np.allclose( p.normal[0], p.normal[2])

    assert (p.calc( np.array([[ -1, 0, 0]] ) ) == np.array( [-1, 0, 2] )).all()
    assert (p.calc( np.array([[  2, 0, 0]] ) ) == np.array( [ 2, 0,-1] )).all()
    assert (p.calc( np.array([[  0, 2, 0]] ) ) == np.array( [ 0, 2,-1] )).all()
    assert (p.calc( np.array([[  0,-1, 0]] ) ) == np.array( [ 0,-1, 2] )).all()
    
    print("Tests look OK")

if __name__ == "__main__":
    test()

