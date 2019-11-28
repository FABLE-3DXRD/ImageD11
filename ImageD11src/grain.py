## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py


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

import numpy as np, math
import ImageD11.indexing
import xfab.tools

# helpers : put these into xfab.tools at some point?
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



class grain:
    def __init__(self,ubi,translation=None, **kwds):
        if translation is None:
            # If translation has not been read from ubi file make it 
            # be None to avoid confusion with a grain which is known
            # to be at [0,0,0]
            self.translation = None
        else:
            self.translation = np.array(translation,np.float)
        self.set_ubi(ubi)

    def set_ubi(self, ubi):
        """ Update the orientation and clear cached values """
        self.ubi = np.array(ubi,np.float)
        self.clear_cache()        

    def clear_cache(self):
        # We will cache a bunch of things and access them
        # via properties. Usually we return copies of these
        # so we avoid people overwriting cached values by accident.
        # e.g. cell = grain.unitcell
        #      cell[4:] = 90 # degree angle for orthorhombic
        #      oups : grain._unitcell was corrupted!
        self._ub = None
        self._u = None
        self._rod = None
        self._mt = None
        self._rmt = None
        self._unitcell = None

    @property 
    def ub(self):
        """ The UB matrix from Busing and Levy
        columns are the reciprocal space lattice vectors
        """
        if self._ub is None:
            self._ub = np.linalg.inv(self.ubi)
        return self._ub.copy()

    @property
    def u(self):
        """ The orientation matrix (U) from Busing and Levy """
        if self._u is None:
            self._u = xfab.tools.ubi_to_u(self.ubi)
        return self._u.copy()

    @property
    def Rod(self):
        """ A Rodriguez vector. 
        Length proportional to angle, direction is axis"""
        if self._rod is None:
            self._rod = xfab.tools.u_to_rod( self.u )
        return self._rod.copy()

    @property
    def mt(self):
        """Metric tensor """
        if self._mt is None:
            self._mt = np.dot( self.ubi, self.ubi.T )
        return self._mt.copy()

    @property
    def rmt(self):
        """ Reciprocal metric tensor """
        if self._rmt is None:
            self._rmt = np.linalg.inv( self.mt )
        return self._rmt.copy()
        
    @property
    def unitcell(self):
        """ a,b,c,alpha,beta,gamma """
        if self._unitcell is None:
            self._unitcell = xfab.tools.ubi_to_cell( self.ubi )
        return self._unitcell.copy()

    def eps_grain(self, dzero_cell):
        """ dzero_cell = [a,b,c,alpha,beta,gamma]
        Returns eps = [e11, e12, e13, e22, e23, e33]
        ... in the grain system
        """
        B_xfab = xfab.tools.form_b_mat( self.unitcell )
        eps = xfab.tools.b_to_epsilon( B_xfab, dzero_cell )
        return eps       

    def eps_sample_matrix(self, dzero_cell):
        """ dzero_cell = [a,b,c,alpha,beta,gamma]
        Returns eps as a symmetric matrix
        ... in the sample system (z usually up axis, x along the beam at omega=0)
        """
        Eg = self.eps_grain_matrix( dzero_cell )
        U = self.u
        Es = np.dot( U,  np.dot( Eg, U.T ) )
        return Es

    def eps_grain_matrix(self, dzero_cell):
        """ dzero_cell = [a,b,c,alpha,beta,gamma]
        Returns eps as a symmetric matrix
        ... in the grain system 
        """
        eps = self.eps_grain( dzero_cell )
        return e6_to_symm( eps )

    def eps_sample(self, dzero_cell):
        """ dzero_cell = [a,b,c,alpha,beta,gamma]
        Returns eps as [e11, e12, e13, e22, e23, e33]
        ... in the sample system (z usually up axis, x along the beam at omega=0)
        """
        Es = self.eps_sample_matrix( dzero_cell )
        return symm_to_e6( Es )

        


    
def write_grain_file(filename, list_of_grains):
    f = open(filename, "w")
    for g in list_of_grains:
        t = g.translation
        f.write("#translation: %g %g %g\n"%(t[0],t[1],t[2]))
        if hasattr(g,"name"):
            f.write("#name %s\n"%(g.name.rstrip()))
        if hasattr(g,"intensity_info"):
            f.write("#intensity_info %s\n"%(g.intensity_info.rstrip()))
        if hasattr(g,"npks"):
            f.write("#npks %d\n"%(int(g.npks)))
        if hasattr(g,"nuniq"):
            f.write("#nuniq %d\n"%(int(g.nuniq)))
        if hasattr(g,"Rod"):
            try:
                f.write("#Rod %f %f %f\n"%tuple([float(r) for r in g.Rod]))
            except:
                f.write("#Rod %s"%(g.Rod))
        f.write("#UBI:\n")
        u = g.ubi
        # More than float32 precision
        f.write("%.9g %.9g %.9g\n"  %(u[0,0],u[0,1],u[0,2]))
        f.write("%.9g %.9g %.9g\n"  %(u[1,0],u[1,1],u[1,2]))
        f.write("%.9g %.9g %.9g\n\n"%(u[2,0],u[2,1],u[2,2]))
    f.close()

def read_grain_file(filename):
    """read ubifile and return a list of ubi arrays """
    f = open(filename, "r")
    grainsread = []
    u = []
    t = None
    p = {}
    for line in f:
        if line.find("#translation:")==0:
            t = [ float(x) for x in line.split()[1:]]
            continue
        if line[0] == "#" and line.find("UBI")<0:
            k,v=line[1:].split(" ",1)
            p[k]=v
            continue
        if line[0] == "#" and line.find("intensity_info")>-1:
            p["intensity_info"] = line.split("intensity_info")[1].rstrip()
        if line.find("#")==0: continue
        vals = [ float(x) for x in line.split() ]
        if len(vals) == 3:
            u = u + [vals]
        if len(u)==3:
            grainsread.append( grain(u, t) )
            for k in ["name","npks","nuniq","intensity_info"]: # Rod - is recomputed when needed
                if k in p:
                    setattr(grainsread[-1], k, p[k])
            p={}
            u = []
            t = None
    f.close()
    return grainsread
