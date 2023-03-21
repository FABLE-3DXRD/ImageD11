
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
from ImageD11 import grain, transform, cImageD11, indexing, unitcell, refinegrains
import scipy.optimize
from scipy.spatial.transform import Rotation
import math, time, sys, logging

def unit(a):
    """ Normalise vector """
    return a/np.sqrt(np.dot(a,a))

def get_tth( ds, wvln):
    """ Convert 1/d to Bragg angle """
    return 2 * np.degrees( np.arcsin( ds * wvln / 2 ) )


class indexer:
    """
    A class for searching for orientation matrices
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
        """ Sets the unit cell parameters for indexing """
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
        self.peaks_xyzT = np.array((self.cf.xl,self.cf.yl,self.cf.zl)).T.copy()
        om = self.cf.omega
        osign = self.transformpars.get("omegasign")
        tth, eta = transform.compute_tth_eta_from_xyz(
            self.peaks_xyzT.T, om*osign,
            **self.transformpars.parameters)
        self.cf.addcolumn( tth, "tth", )
        self.cf.addcolumn( eta, "eta", )
        self.cf.addcolumn( refinegrains.lf(tth,eta), 'lf')
        gx, gy, gz = transform.compute_g_vectors(
            tth, eta, om*osign,
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
        self.mok = np.zeros( self.cf.nrows, bool )
        

    def tthcalc(self, hkls = None):
        """
        Computes the twotheta for a unit cell given a list of 
        hkls or from the unit cell generating a list of hkls.
        """
        if hkls is None:
            dslimit = self.cf.modg.max()
            wvln = self.transformpars.get( "wavelength" )
            ttol = self.transformpars.get( "fit_tolerance" ) # tth units
            # at limit this ttol make which ds range?
            tthlim = 2 * np.degrees( np.arcsin( dslimit * wvln / 2 ) )
            dstol  = 2 * np.sin( np.radians( (tthlim + ttol) / 2 )
                                 ) / wvln - dslimit
            self.unitcell.makerings(dslimit+dstol/2, tol = dstol)
        else:
            self.unitcell.peaks = [ (self.unitcell.ds(h),h) for h in hkls ]
            self.unitcell.makerings( dslimit, tol=dstol ) # , newpeaks=False )
        self.unitcell.ringtth = [ 2 * np.degrees( np.arcsin( ds * wvln / 2 ) )
                                  for ds
                                  in self.unitcell.ringds ]

    def assigntorings(self):
        """
        Assign the g-vectors to hkl rings that are in self.unitcell
        """
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
                
                
    def printringassign(self):
        # Report on assignments
        print("Ring     (  h,  k,  l) Mult  total indexed to_index  ")
        # try reverse order instead
        dsr = self.unitcell.ringds
        for j in range(len(dsr))[::-1]:
            ind = self.cf.ring == j
            n_indexed  = (self.cf.labels[ind] >  -0.5).sum()
            n_to_index = (self.cf.labels[ind] <  -0.5).sum()
            h=self.unitcell.ringhkls[dsr[j]][0]
            print("Ring %-3d (%3d,%3d,%3d)  %3d  %5d  %5d  %5d  %.4f"%(\
                j,h[0],h[1],h[2],len(self.unitcell.ringhkls[dsr[j]]),
                     ind.sum(),n_indexed,n_to_index,
                     self.unitcell.ringtth[j]  ))
        print("Total peaks",self.cf.nrows,"assigned",(self.cf.ring>=0).sum())

        
    def rings_2_use(self, rings = None, multimin = 12 ):
        """Filter rings as having low multiplicity for indexing searches
        
        Give rings = [list of rings to use]
        Or multimin = all rings with low multiplicity
        """
        ring = self.cf.ring.astype(int)
        if rings is None:
            mok = np.zeros( (self.cf.nrows,), bool )
            for i, ds in enumerate( self.unitcell.ringds ):
                mult = len(self.unitcell.ringhkls[ds])
                if mult <= multimin:
                    mok[ ring == i ] = True
        else:
            mok = np.zeros( (self.cf.nrows,), bool )
            for r in rings:
                mok[ ring == r ] = True
        self.mok = mok        

    def search1d(self, gvec_id, 
                 hkl = None,
                 angstart = -180, 
                 angend = 180,
                 nang = 3600,
                 tol = 0.1,
                ):
        """
        gvec_id = an integer for a row of self.colfile 
        hkl  = hkl indices to assign to the peak (None 
            means guess from self.cf.ring)
        This is inspired from Bernier's fibre texture method (citation:
        https://github.com/HEXRD/hexrd/blob/master/hexrd/findorientations.py
        )        
        """
        gv = np.array( (self.cf.gx[gvec_id], 
                        self.cf.gy[gvec_id], 
                        self.cf.gz[gvec_id]), float)
        print(gvec_id)
        if hkl is None:
            ring = self.cf.ring[gvec_id]
            ds = self.unitcell.ringds[int(ring)]
            hkls = self.unitcell.ringhkls[ ds ]
            hkl = np.array( hkls[0], int )
            print("Choosing",hkl,"from hkls", hkls)
        assert hkl.shape == (3,)
        g0 = np.dot( self.unitcell.B, hkl ) # non rotated g
        # normalised vectors
        n0 = g0/np.linalg.norm( g0 )
        nobs = gv/np.linalg.norm( gv )
        cosa = np.dot( n0, nobs )
        ang  = np.arccos(cosa)
        # if the vectors are already parallel ?
        if ang < np.radians(0.001):
            u0 = np.eye(3)
        else:
            vec  = np.cross( n0, nobs )
            sina = np.linalg.norm( vec )
            u0 = Rotation.from_rotvec( ang * vec / sina ).as_matrix()
        ub0 = np.dot( u0, self.unitcell.B )
        ubi0 = np.linalg.inv( ub0 )
        # Now we want to rotate around nobs
        allgve = np.array( (self.cf.gx,self.cf.gy,self.cf.gz) ).T.copy()
        scores = []
        ubis = []
        gc = np.dot( ub0, hkl )
        angs = np.radians( np.linspace( angstart, angend, nang) )
        ubis = [ np.dot(ubi0, 
                Rotation.from_rotvec( nobs * a ).as_matrix())
                for a in angs ]
        scores = [ cImageD11.score( ubi, allgve, tol )
                  for ubi in ubis ]
        return scores, ubis
    
    def choose( self, gnum, tol):
        isig = self.cf.npixels * self.cf.avg_intensity * self.cf.lf
        idp = np.argmax( self.mok * (self.cf.labels<0) * isig )
        angstart = 0
        angend = 360
        nang = 3600
        s,ubis=self.search1d(idp, tol=tol,
                             angstart = angstart, angend=angend, nang=nang)
        matfit = ubis[np.argmax(s)].copy()
        tfit = np.zeros(3)
        tol = 0.05
        inds, hkls = self.assign( matfit, tfit, tol )
        matfit, tfit = self.refine( matfit, tfit, inds, hkls, tol )
        inds, hkls = self.assign( matfit, tfit, tol )
        print( inds.shape, tfit ) 
        print( indexing.ubitocellpars( matfit ) )
        self.cf.labels[ inds ] = gnum
        self.grains[ gnum ] = grain.grain( matfit, tfit )
    
        
    def pairs(self, hkl1, hkl2, cos_tol = 0.02, hkl_tol = 0.1):
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
        print("Angle, cosangle",angle,cosangle,hkl1,hkl2)
        assert angle > 1 and angle < 179, "hkls are parallel"
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
           #if i%100==0:
           #   print i,time.time()-start,len(n1.T)
           for k in goodones:
               this_n2 = n2[:,k]
               T_g[0] = this_n1
               T_g[2] = unit( np.cross( this_n1, this_n2 ) )
               T_g[1] = unit( np.cross( this_n1, T_g[2]) )
               U = np.dot( T_g.T, T_c)
               ub =  np.dot( U, self.unitcell.B)
               ubi = np.linalg.inv( ub )
               ubio = ubi.copy()
               npks = cImageD11.score(ubi,gvf,hkl_tol)
               pairs.append( (ind1[i], ind2[k], U, ubi ) )
               print(npks, end=' ') 
               
               ubi, trans = self.refine( ubi, np.zeros(3,float), tol=hkl_tol )
               inds, hkls = self.assign( ubi, trans, hkl_tol )
               ubi, trans = self.refine( ubi, trans, inds = inds, hkls= hkls, tol=hkl_tol )
               print(npks, ubi)
               print("cell: ",6*"%.6f "%(  indexing.ubitocellpars(ubi) ))
               print("position: ",trans)
               print()
        self.pairscache=pairs
        print(time.time()-start,"for",len(pairs),n1.shape, n2.shape)
        return pairs


    def assign(self, ubi, translation, tol):
        gv = np.zeros(self.peaks_xyzT.shape, float )
        cImageD11.compute_gv( self.peaks_xyzT,
                    self.cf.omega,
                    self.transformpars.get('omegasign'),
                    self.transformpars.get('wavelength'),
                    self.transformpars.get('wedge'),  
                    self.transformpars.get('chi'),
                    translation,
                    gv)
        hkl = np.dot( ubi, gv.T )
        hkli = np.floor( hkl + 0.5 )
        e = (hkl - hkli)
        e = (e*e).sum(axis=0)
        inds = np.compress( e < tol*tol, np.arange(len(hkl[0])) )
        return inds , hkli[:,inds]


    NC = 0
    def gof( self, p, *args ):
        self.NC += 1
        try:
            hkls, inds, peaks_xyzT, gobs, omega = args
        except:
            print(args, len(args))
            raise
        p.shape = 4,3
        ub = p[:3]
        t = p[3]
        gcalc = np.dot( ub, hkls )
        cImageD11.compute_gv( peaks_xyzT,
                              omega,
                              self.transformpars.get('omegasign'),
                              self.transformpars.get('wavelength'),
                              self.transformpars.get('wedge'),  
                              self.transformpars.get('chi'),
                              t,
                              gobs)
        #print gobs
        #print gcalc
        #print (gobs-gcalc).ravel()
        #1/0
        e = (gcalc - gobs.T).ravel()
        p.shape = 12,
 #       print p-0.1,(e*e).sum()
        return e#(e*e).sum()
        
    def refine(self, ubi, translation=(0,0,0) ,
               inds = None, hkls = None,
               tol = 0.05):
        """
        Fit ubi and translation
        ihkls = array of peak_index, h, k, l
        tol = hkl tolerance to use in fitting
        """
        import time
        self.NC =0
        start = time.time()
        if inds is None:
            inds , hkls = self.assign(  ubi, translation, tol )
        ub = np.linalg.inv( ubi )
        x0 = np.array( list( ub.ravel() ) + list(translation ) )
        fun = self.gof
        args = (hkls, inds, self.peaks_xyzT[inds],
                np.zeros(self.peaks_xyzT[inds].shape ),
                self.cf.omega[inds])
        def Dfun( x0, *args ):
            epsilon = np.ones(12)*1e-6
            epsilon[-3:] = 1.
            return deriv( x0, fun, epsilon, *args)
        res, ier = scipy.optimize.leastsq( fun, x0, args, )#
                                            # Dfun, col_deriv=True)
        ub = np.reshape(res[:9], (3,3))
        t = res[-3:]
        ubi = np.linalg.inv( ub )
        return ubi, t
 
def deriv( xk, f, epsilon, *args):
    f0 = f(*((xk,) + args))
    grad = np.zeros((len(xk),len(f0)), float)
    ei = np.zeros((len(xk),), float)
    for k in range(len(xk)):
        ei[k] = 1.0
        d = epsilon * ei
        grad[k] = (f(*((xk + d,) + args)) - f0) / d[k]
        ei[k] = 0.0
    return grad

if __name__=="__main__":
    from ImageD11.columnfile import columnfile
    from ImageD11.parameters import read_par_file
    import ImageD11.grain
    import sys
    p = read_par_file(sys.argv[1])
    c = columnfile( sys.argv[2] )
    i = indexer( p, c )
    if sys.argv[3][:3] == "fit":
        #                                            0                   1                  2           3    4             5
        # \test\simul_1000_grains>python ..\..\ImageD11\indexer.py Al1000\Al1000.par Al1000\Al1000.flt fit allgrid.map allgridfitscipy.map
        gl = ImageD11.grain.read_grain_file( sys.argv[4] )
        inds = np.arange( len(gl), dtype=int )
        allhkls = np.array( (c.h, c.k, c.l) )
        for k,g in enumerate(gl):
            if len(sys.argv[3]) == len("fit"):
                for j, tol in enumerate([0.05,0.02,0.01,0.0075]):
                    inds, hkls = i.assign( g.ubi, g.translation, tol )
                    ubi, t = i.refine(   g.ubi, 
                                         translation=g.translation,
                                         inds = inds, hkls = hkls,
                                         tol = tol)  
                    g.translation = t
                    g.set_ubi( ubi )
            else:
                inds = np.compress( c.labels == k , inds )
                hkls = np.compress( c.labels == k , allhkls )
                for j in range(3):
                    ubi, t = i.refine(   g.ubi, 
                                         translation=g.translation,
                                         inds = inds, hkls = hkls,
                                         tol = tol)  
                    g.translation = t
                    g.set_ubi( ubi )
            print(k, len(inds),6*"%.8f "%(indexing.ubitocellpars(ubi)))
            print("\t",t)
        ImageD11.grain.write_grain_file( sys.argv[5], gl )
    else:
        i.updatecolfile()
        i.tthcalc()
        print("Calling assign")
        i.assigntorings()
        hkl1 = [int(h) for h in sys.argv[3].split(",")]
        hkl2 = [int(h) for h in sys.argv[4].split(",")]
        print("Calling pairs")
        i.pairs(hkl1, hkl2)


   # e.g. python indexer.py ../test/demo/eu3.pars ../test/demo/eu.flt 2,0,0 0,0,4
