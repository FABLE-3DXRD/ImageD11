
from __future__ import print_function

# ImageD11 Software for beamline ID11
# Copyright (C) 2020  Jon Wright
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

import numpy as np, math

# helpers 
def e6_to_symm(e):
    """ Follows the eps convention from xfab 
    eps = [e11, e12, e13, e22, e23, e33]
    """
    e11, e12, e13, e22, e23, e33 = e
    return np.array(((e11,e12,e13),
                     (e12,e22,e23),
                     (e13,e23,e33)))

def symm_to_e6(m):
    """ Follows the eps convention from xfab 
    eps = [e11, e12, e13, e22, e23, e33]
    """
    return np.array( ( m[0,0], m[0,1], m[0,2],
                               m[1,1], m[1,2],
                                       m[2,2] ) )

def cell_to_B( cell ):
    from ImageD11.unitcell import unitcell
    B = unitcell( cell ).B
    return B

class DeformationGradientTensor( object ):
    
    def __init__(self, ubi, ub0):
        """
        see docs/DeformationGradientTensor.ipynb
        F = dot( ubi.T, ub0.T )
        F = ui.bi.b0.u0
        """
        if hasattr(ubi, 'ubi'): # you passed in and ImageD11.grain
            ubi = ubi.ubi
        assert ubi.shape == (3,3)
        if hasattr( ub0, 'ub'): 
            # you passed in and ImageD11.grain
            ub0 = ub0.ub
        assert ub0.shape == (3,3)
        self.F = np.dot( ubi.T, ub0.T )
        self._svd = None
        self._vrs = None

    @property
    def SVD(self):
        """ Returns the singular value decomposition of F """
        if self._svd is None:
            self._svd = np.linalg.svd(self.F)
        return self._svd

    @property
    def VRS(self):
        """ Returns the Polar decomposition of F=V.R=R.S 
        with R as a rotation matrix and V, S as symmetric
        """
        if self._vrs is None:
            w,sing,vh = self.SVD
            R = np.dot( w, vh )
            S = np.dot( vh.T, np.dot( np.diag(sing),  vh ) )
            V = np.dot( w   , np.dot( np.diag(sing), w.T ) )
            self._vrs = V,R,S
        return self._vrs

    @property
    def U(self):
        """ Returns the Busing and Levy U matrix relating ubi to ub0 """
        v,r,s = self.VRS
        return r
    
    def finite_strain_ref(self, m=0.5):
        """
        Returns the finite strain in the reference co-ordinate system
        if the polar decomposition gives : 
        F = V.R = R.S = ui.bi.b0.u0
        Ft.F removes ui final orientation effect
        Em = (S^m - I )/2m
        """
        m2 = int(round(m*2))
        assert np.allclose( m2 * 0.5, m )
        if m2 == 0:
            u,s,vt = self.SVD # == F
            logs = np.diag( np.log( s ) )
            logFFT = np.dot( np.dot( vt.T, logs ), vt )
            Em = logFFT # empirically, no division by 2 for matching in small strain limit
        elif (m2 % 2) == 0: # No SVD in this path
            m = int(round(m))
            Cm = np.linalg.matrix_power( np.dot( self.F.T, self.F ), m )
            Em = (Cm - np.eye(3))/m2
        else:
            V,R,S = self.VRS # == F
            Um = np.linalg.matrix_power(S, m2)
            Em = (Um - np.eye(3))/m2
        return Em
    
    def finite_strain_lab(self, m=0.5):
        """
        Returns the finite strain in the lab co-ordinate system
        if the polar decomposition gives : 
        F = V.R = R.S = ui.bi.b0.u0
        F.Ft removes u0 initial orientation effect
        Em = (V^m - I )/2m
        """
        m2 = int(round(m*2))
        assert np.allclose( m2 * 0.5, m )
        if m2 == 0:
            u,s,vt = self.SVD # == F
            logs = np.diag( np.log( s ) )
            logFTF = np.dot( np.dot( u, logs ), u.T )
            em = logFTF # *0.5
        elif (m2 % 2) == 0: # No SVD in this path
            m = int(round(m))
            Bm = np.linalg.matrix_power( np.dot( self.F, self.F.T ), m )
            em = (Bm - np.eye(3))/m2
        else:
            V,R,S = self.VRS # == F
            Vm = np.linalg.matrix_power(V, m2)
            em = (Vm - np.eye(3))/m2
        return em
    
            


