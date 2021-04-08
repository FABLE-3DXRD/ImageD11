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
Tools for processing UBI matrices together with data

Uses indexing objects to do refinement

Mostly unimplemented dummies - see also refinegrains

"""

import numpy as np

class ubitool:
    def __init__(self,ubilist=None,ubifile=None,obsdata=None):
        self.ubilist=ubilist
        self.ubifile=ubifile
        self.obsdata=obsdata
    def get_unit_cell(self,ubi=None):
        """
        Convert UBI representation to give unit cell
        """
        pass
    def get_orientation(self,ubi=None):
        """
        Give orientation matrix independant of unit cell
        """
        pass
    def validate_ubi(self,ubi=None):
        """
        Find the number of peaks versus hkl_tol when refining
        This should plateau when the tolerance is correct and the
        data are good enough.
        """
        pass
    def validate_ubi_collection(self):
        """
        Check for duplicate orientations
        """
        pass
    def validate_peak_assignements(self):
        """
        Make sure each hkl is only assigned to one peak.
        Offer to merge if there are duplicates
        """
        pass
    def refine_translations(self,ubi=None):
        """
        Compute an offset in x/y/z for the origin of
        the grain with respect to the centre of rotation
        """
        pass
    def read_ubi_file(self,filename):
        """
        Get ubi matrices from a file
        """
        i=0; u = np.zeros((3,3),float)
        for line in open(filename,"r").readlines():
            uij = [float(x) for x in line.split()]
            if len(uij)==3:
                u[i,0]=uij[0] ; u[i,1]=uij[1] ; u[i,2]=uij[2]
                i=i+1
            else:
                self.ubilist.append(u)
                u = np.zeros((3,3),float)
