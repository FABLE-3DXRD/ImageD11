
from __future__ import print_function
# ImageD11_v0.4 Software for beamline ID11
# Copyright (C) 2015  Jon Wright
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


import numpy as np
from ImageD11 import cImageD11, grain, transform, unitcell

import math, time, sys, logging

def unit(a): return a/np.sqrt(np.dot(a,a))

def get_tth( ds, wvln):
    return 2 * np.degrees( np.arcsin( ds * wvln / 2 ) )


class indexer:
    """
    A class for searching for orientation matrices

    Base the errors on detector pixels and omega steps
    """
    def __init__(self,
                 transformpars,
                 colfile,
                 ):
        """
        """
        self.transformpars = transformpars
        self.loadcolfile(colfile)        
        self.setcell()
        self.reset()

    def setcell(self):
        cp = [self.transformpars.get("cell_%s"%(s)) for s in
              "_a _b _c alpha beta gamma".split()]
        self.unitcell = unitcell.unitcell(cp, 
                        self.transformpars.get("cell_lattice_[P,A,B,C,I,F,R]"))
        

    def loadcolfile(self, colfile):
        self.cf = colfile
        self.updatecolfile()
        for label in "gx gy gz omega tth eta xl yl zl".split():
            assert label in self.cf.titles, label

    def updatecolfile(self):
       if "xl" not in self.cf.titles:
          if "sc" in self.cf.titles:
             pks = self.cf.sc, self.cf.fc
          elif "xc" in self.cf.titles: 
             pks = self.cf.xc, self.cf.yc
          else:
             raise Exception("peaks file misses xc or sc")
          xl,yl,zl = transform.compute_xyz_lab( pks,
                                                **self.transformpars.parameters)
          self.cf.addcolumn(xl,"xl")
          self.cf.addcolumn(yl,"yl")
          self.cf.addcolumn(zl,"zl")
       peaks_xyz = np.array((self.cf.xl,self.cf.yl,self.cf.zl))
       om = self.cf.omega
       sign = self.transformpars.get("omegasign")
       tth, eta = transform.compute_tth_eta_from_xyz(
          peaks_xyz, om*sign,
          **self.transformpars.parameters)
       self.cf.addcolumn( tth, "tth", )
       self.cf.addcolumn( eta, "eta", )
       gx, gy, gz = transform.compute_g_vectors(
           tth, eta, om*sign,
           wvln  = self.transformpars.get("wavelength"),
           wedge = self.transformpars.get("wedge"),
           chi   = self.transformpars.get("chi") )
       self.cf.addcolumn(gx, "gx")
       self.cf.addcolumn(gy, "gy")
       self.cf.addcolumn(gz, "gz")           
       self.cf.addcolumn( np.sqrt( gx * gx + 
                                   gy * gy + 
                                   gz * gz ),
                          "modg")
       
    
    def reset(self):
        """ when pars change or colfile changes etc """
        if "ring" not in self.cf.titles:
            self.cf.addcolumn( -np.ones(self.cf.nrows), "ring" )
            self.cf.addcolumn( 42*np.ones(self.cf.nrows), "ringerr" )
        if "labels" not in self.cf.titles:
            self.cf.addcolumn( -np.ones(self.cf.nrows), "labels" )

       

    def tthcalc(self, hkls = None):
        if hkls is None:
            dslimit = np.maximum.reduce(self.cf.modg)
            wvln = self.transformpars.get( "wavelength" )
            ttol = self.transformpars.get( "fit_tolerance" ) # tth units
            # at limit this ttol make which ds range?
            tthlim = 2 * np.degrees( np.arcsin( dslimit * wvln / 2 ) )
            dstol  = 2 * np.sin( np.radians( (tthlim + ttol) / 2 ) ) / wvln - dslimit
            self.unitcell.makerings(dslimit+dstol/2, tol = dstol)
        else:
            self.unitcell.peaks = [ (self.unitcell.ds(h),h) for h in hkls ]
            self.unitcell.makerings( dslimit, tol=dstol, newpeaks=False )
        self.unitcell.ringtth = [ 2 * np.degrees( np.arcsin( ds * wvln / 2 ) )
                                  for ds
                                  in self.unitcell.ringds ]

    def assigntorings(self):
        """
        Assign the g-vectors to hkl rings
        """
        dsr = self.unitcell.ringds
        self.cf.ring[:]=-1
        self.cf.ringerr[:]=42.
        tol = self.transformpars.get( "fit_tolerance" )
        for i,tth in enumerate(self.unitcell.ringtth):
           diff = self.cf.tth - tth
           adiff = np.abs(diff)
           mask = (adiff < tol) & (adiff < self.cf.ringerr)
           if mask.sum()>0:
              self.cf.ring[mask] = i
              self.cf.ringerr[mask] = diff[mask]
        # Report on assignments
        print("Ring     (  h,  k,  l) Mult  total indexed to_index  ")
        # try reverse order instead
        dsr = self.unitcell.ringds
        for j in range(len(dsr))[::-1]:
            ind = self.cf.ring == j
            n_indexed  = (self.cf.labels[ind] >  -0.5).sum()
            n_to_index = (self.cf.labels[ind] <  -0.5).sum()
            h=self.unitcell.ringhkls[dsr[j]][0]
            print("Ring %-3d (%3d,%3d,%3d)  %3d  %5d  %5d  %5d"%(\
                j,h[0],h[1],h[2],len(self.unitcell.ringhkls[dsr[j]]),
                     ind.sum(),n_indexed,n_to_index))


    def pairs(self, hkl1, hkl2, cos_tol = 0.02, hkl_tol = 0.05):
        """
        We only look for reflection pairs matching a single hkl pairing
        """
        import time
        start = time.time()
        w = self.transformpars.get( "wavelength" )
        tth1 = get_tth(self.unitcell.ds(hkl1), w) 
        tth2 = get_tth(self.unitcell.ds(hkl2), w)
        tthtol = self.transformpars.get( "fit_tolerance" )
        allinds = np.arange( self.cf.nrows )
        ind1 = allinds[ abs(self.cf.tth - tth1) < tthtol ]
        ind2 = allinds[ abs(self.cf.tth - tth2) < tthtol ]
        angle, cosangle = self.unitcell.anglehkls( hkl1, hkl2 )
        g = np.array( (self.cf.gx, self.cf.gy, self.cf.gz), float )
        n = g/self.cf.modg
        gvf = g.T.copy()
        n1 = n[:,ind1]
        n2 = n[:,ind2]
        pairs = []
        j = np.arange( n2.shape[1] )
        # from unitcell.orient
        h1c=unit(np.dot( self.unitcell.B, hkl1 ))
        h2c=unit(np.dot( self.unitcell.B, hkl2 ))
        t1c=unit(h1c)
        t3c=unit(np.cross(h1c,h2c))
        t2c=unit(np.cross(h1c,t3c))
        T_c =  np.array( [t1c, t2c, t3c] )
        T_g = np.zeros((3,3))
        for i,this_n1 in enumerate(n1.T):
           cosa = np.dot( this_n1, n2 )
           goodones = j[abs(cosa - cosangle) < cos_tol]
           # t1 is along g1
           # t2 is plane of both: g1x(g1xg2)
           # t3 is perpendicular to both
           if i%100==0:
              print(i,time.time()-start,len(n1.T))
           for k in goodones:
               this_n2 = n2[:,k]
               T_g[0] = this_n1
               T_g[2] = unit( np.cross( this_n1, this_n2 ) )
               T_g[1] = unit( np.cross( this_n1, T_g[2]) )
               U = np.dot( T_g.T, T_c)
               ubi =  np.dot( U, self.unitcell.B) 
#               npks = cImageD11.score(ubi,gvf,0.1)
               pairs.append( (ind1[i], ind2[k], U, ubi ) )
#               print npks, ubi
        self.pairs=pairs
        print(time.time()-start,"for",len(pairs),n1.shape, n2.shape)
        return pairs



 


if __name__=="__main__":
   from ImageD11.columnfile import columnfile
   from ImageD11.parameters import read_par_file
   import sys
   p = read_par_file(sys.argv[1])
   c = columnfile( sys.argv[2] )
   i = indexer( p, c )
   i.tthcalc()
   i.assigntorings()
   hkl1 = [int(h) for h in sys.argv[3].split(",")]
   hkl2 = [int(h) for h in sys.argv[4].split(",")]
   i.pairs(hkl1, hkl2)
