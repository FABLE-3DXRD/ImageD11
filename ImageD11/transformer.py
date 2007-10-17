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

import numpy.oldnumeric as n
import logging
import math
from ImageD11 import transform, unitcell, columnfile
from ImageD11.parameters import par, parameters

PARAMETERS = [
     par( "omegasign", 1.0, 
          helpstring = "Sign of the rotation about z (normally +1 for right handed)",
          vary=False, 
          can_vary=False),
     par('z_center',  1024.0,
          helpstring = "Beam centre in vertical, pixels",
          vary=True, 
          can_vary=True,
          stepsize = 1.0),
     par('y_center',  1024.0,
          helpstring = "Beam centre in horizontal, pixels",
          vary=True, 
          can_vary=True,
          stepsize = 1.0),
     par('distance',   50000.0,
          helpstring = "sample detector distance, same units as pixel size",
          vary=True, 
          can_vary=True,
          stepsize = 1.0),
     par('z_size',   48.08150,
          helpstring = "pixel size in vertical, same units distance",
          vary=False, 
          can_vary=False), # this could actually vary - a bit crazy?
     par('y_size',   46.77648,
          helpstring = "pixel size in horizontal, same units as distance",
          vary=False, 
          can_vary=False), # this could actually vary - a bit crazy?
     par('tilt_z',    0.0,
          helpstring = "detector tilt, right handed around z",
          vary=True, 
          can_vary=True,
          stepsize = transform.radians(0.1) ),
     par('tilt_y',    0.0,
          helpstring = "detector tilt, right handed around y",
          vary=True, 
          can_vary=True,
          stepsize = transform.radians(0.1) ),
     par('tilt_x',    0.0,
          helpstring = "detector tilt, right handed around x",
          vary=False, 
          can_vary=True,
          stepsize = transform.radians(0.1) ),
     par('fit_tolerance', 0.05,
          helpstring = "tolerance to decide which peaks to use",
          vary=False, 
          can_vary=False),
     par('wavelength', 0.155,
          helpstring = "wavelength, normally angstrom, same as units unit cell ",
          vary=False, 
          can_vary=True,# but you'll be lucky!
          stepsize = 0.0001 ),
     par('wedge',0.0,
          helpstring = "wedge, rotation around y under omega",
          vary=False, 
          can_vary=True,
          stepsize = transform.radians(0.1) ),
     par('chi', 0.0,
          helpstring = "wedge, rotation around x under omega",
          vary=False, 
          can_vary=True,
          stepsize = transform.radians(0.1) ),
     par('cell__a' , 4.1569,
          helpstring = "unit cell par, same units as wavelength",
          vary=False, 
          can_vary=True,
          stepsize = 0.01 ),
     par('cell__b' , 4.1569,
          helpstring = "unit cell par, same units as wavelength",
          vary=False, 
          can_vary=True,
          stepsize = 0.01 ),
     par('cell__c' , 4.1569,
          helpstring = "unit cell par, same units as wavelength",
          vary=False, 
          can_vary=True,
          stepsize = 0.01 ),
     par('cell_alpha' , 90.0,
          helpstring = "unit cell par, degrees",
          vary=False, 
          can_vary=True,
          stepsize = 0.01 ),
     par('cell_beta' , 90.0,
          helpstring = "unit cell par, degrees",
          vary=False, 
          can_vary=True,
          stepsize = 0.01 ),
     par('cell_gamma' , 90.0,
          helpstring = "unit cell par, degrees",
          vary=False, 
          can_vary=True,
          stepsize = 0.01 ),
     par('cell_lattice_[P,A,B,C,I,F,R]', "P",
          helpstring = "lattice centering type. Try P if you are not sure",
          vary=False, 
          can_vary=False),
     par('o11' , 1, 
          helpstring = "detector flip element +1 for frelon & quantix",
          vary=False, 
          can_vary=False ),
     par('o12' , 0, 
          helpstring = "detector flip element 0 for frelon & quantix",
          vary=False, 
          can_vary=False),
     par('o21' , 0, 
          helpstring = "detector flip element 0 for frelon & quantix",
          vary=False, 
          can_vary=False),
     par('o22' , -1, 
          helpstring = "detector flip element -1 for frelon & +1 for quantix",
          vary=False, 
          can_vary=False),
     par('t_x' , 0,
          helpstring = "crystal translation, units as distance/pixels",
          vary=False, 
          can_vary=True,
          stepsize = 1.),
     par('t_y' , 0,
          helpstring = "crystal translation, units as distance/pixels",
          vary=False, 
          can_vary=True,
          stepsize = 1.),
     par('t_z' , 0,
          helpstring = "crystal translation, units as distance/pixels",
          vary=False, 
          can_vary=True,
          stepsize = 1.),
     ]



class transformer:
    """
    Handles the algorithmic, fitting and state information for 
    fitting parameters to give experimental calibrations
    """
    def __init__(self):
        """
        Nothing passed in ?
        It should really be inheriting from the parameter object.
        """
        self.twotheta = None
        self.finalpeaks = None
        self.gv = None
        self.unitcell = None
        # this sets defaults according to class dict.
        self.parameterobj = parameters()
        for p in PARAMETERS:
            self.parameterobj.addpar(p)
        self.pars = self.parameterobj.get_parameters()
        self.colfile = None
   
    def updateparameters(self):
        self.pars = self.parameterobj.get_parameters()

    def get_variable_list(self):
        return self.parameterobj.get_variable_list()

    def getvars(self):
        """ decide what is refinable """
        return self.parameterobj.varylist
    
    def setvars(self, varlist):
        """ set the things to refine """
        self.parameterobj.varylist = varlist

    def loadfiltered(self, filename):
        """
        Read in 3D peaks from peaksearch
        """
        self.colfile = columnfile.columnfile(filename)
        self.finalpeaks = self.colfile.bigarray
        if   ((self.colfile.titles[0:3] == ["sc","fc","omega"]) or
              (self.colfile.titles[0:3] == ["xc","yc","omega"])):
            logging.warning("titles are %s"%(str(self.colfile.titles[0:3])))
        self.peaks_xy = self.finalpeaks[0:2,:]
        self.x = self.finalpeaks[0,:]
        self.y = self.finalpeaks[1,:]
        self.omega = self.finalpeaks[2,:]
        if hasattr(self.colfile, "spot3d_id"):
            self.spot3d_id = self.colfile.spot3d_id
        else:
            self.spot3d_id = range(len(self.x))

    def loadfileparameters(self,filename):
        """ Read in beam center etc from file """
        self.parameterobj.loadparameters(filename)

    def saveparameters(self,filename):
        """ Save beam center etc to file """
        self.parameterobj.saveparameters(filename)

    def applyargs(self,args):
        """ for use with simplex/gof function, alter parameters """
        self.parameterobj.set_variable_values(args)

    def compute_tth_eta(self):
        """ Compute the twotheta and eta for peaks previous read in """
        self.twotheta, self.eta = transform.compute_tth_eta( self.peaks_xy,
                                                             omega=self.omega,
                 **self.parameterobj.get_parameters() ) 
        #self.ds =
        
    def compute_tth_histo(self):
        """ Compute the histogram over twotheta for peaks previous read in
        and filter peaks out falling in bins with less than min_bin_ratio """
        self.finalpeaks = transform.compute_tth_histo(self.finalpeaks,
                                                      self.twotheta , 
                                  **self.parameterobj.get_parameters())
        self.peaks_xy = self.finalpeaks[0:2,:]
        self.x = self.finalpeaks[0,:]
        self.y = self.finalpeaks[1,:]
        self.omega = self.finalpeaks[2,:]
        self.compute_tth_eta()

    def gof(self,args):
        """ Compute how good is the fit of obs/calc peak positions in tth """
        self.applyargs(args)
        # Here, pars is a dictionary of name/value pairs to pass to compute_tth_eta
        self.compute_tth_eta()
        w=self.wavelength
        gof=0.
        npeaks=0
        for i in range(len(self.tthc)):# (twotheta_rad_cell.shape[0]):
            self.tthc[i]=transform.degrees(math.asin(self.fitds[i]*w/2)*2)
            diff=n.take(self.twotheta,self.indices[i]) - self.tthc[i]
#         print "peak",i,"diff",maximum.reduce(diff),minimum.reduce(diff)
            gof=gof+n.sum(diff*diff)
            npeaks=npeaks+len(diff)
        gof=gof/npeaks
        return gof*1e3

    def fit(self,tthmax=180):
        """ Apply simplex to improve fit of obs/calc tth """
        import simplex
        if self.theoryds==None:
            self.addcellpeaks()
        # Assign observed peaks to rings
        self.wavelength = None
        self.indices=[]  # which peaks used
        self.tthc=[]     # computed two theta values
        self.fitds=[]    # hmm?
        self.fit_tolerance=1.
        pars = self.parameterobj.get_parameters()
        w = float(pars['wavelength'])
        self.wavelength = w
        self.fit_tolerance = float(pars['fit_tolerance'])
        print "Tolerance for assigning peaks to rings",self.fit_tolerance,"max tth",tthmax
        for i in range(len(self.theoryds)):
            dsc=self.theoryds[i]
            tthcalc=math.asin(dsc*w/2)*360./math.pi # degrees
            if tthcalc>tthmax:
                break
#         print tthcalc
            logicals= n.logical_and( n.greater(self.twotheta, tthcalc-self.fit_tolerance),
                                     n.less(self.twotheta, tthcalc+self.fit_tolerance)  )

            if sum(logicals)>0:
#            print maximum.reduce(compress(logicals,self.twotheta)),minimum.reduce(compress(logicals,self.twotheta))
                self.tthc.append(tthcalc)
                self.fitds.append(dsc)
                ind=n.compress(logicals,range(self.twotheta.shape[0]))
                self.indices.append(ind)
#            print "Ring",i,tthcalc,maximum.reduce(take(self.twotheta,ind)),minimum.reduce(take(self.twotheta,ind))
#      if raw_input("OK?")[0] not in ["Y","y"]:
#         return
        guess = self.parameterobj.get_variable_values()
        inc = self.parameterobj.get_variable_stepsizes()
        s=simplex.Simplex(self.gof,guess,inc)
        newguess,error,niter=s.minimize()
        inc=[v/10 for v in inc]
        guess=newguess
        s=simplex.Simplex(self.gof,guess,inc)
        newguess,error,niter=s.minimize()
        self.parameterobj.set_variable_values(newguess)
        self.gof(newguess) 
        print newguess


    def addcellpeaks(self,limit=None):
        """
        Adds unit cell predicted peaks for fitting against
        Optional arg limit gives highest angle to use

        Depends on parameters:
            'cell__a','cell__b','cell__c', 'wavelength' # in angstrom
            'cell_alpha','cell_beta','cell_gamma' # in degrees
            'cell_lattice_[P,A,B,C,I,F,R]'
        """
        #
        # Given unit cell, wavelength and distance, compute the radial positions
        # in microns of the unit cell peaks
        #
        pars = self.parameterobj.get_parameters()
        cell = [ pars[name] for name in ['cell__a','cell__b','cell__c',
                                        'cell_alpha','cell_beta','cell_gamma']]
        lattice=pars['cell_lattice_[P,A,B,C,I,F,R]']
        self.unitcell=unitcell.unitcell( cell ,  lattice)
        if self.twotheta==None:
            self.compute_tth_eta()
        # Find last peak in radius
        if limit is None:
            highest = n.maximum.reduce(self.twotheta)
        else:
            highest = limit
        self.wavelength=pars['wavelength']
        ds = 2*n.sin(transform.radians(highest)/2.)/self.wavelength
        self.dslimit = ds
        #print "highest peak",highest,"corresponding d*",ds
        self.theorypeaks = self.unitcell.gethkls(ds)
        self.unitcell.makerings(ds)
        self.theoryds = self.unitcell.ringds
        tths = [n.arcsin(self.wavelength*dstar/2)*2 
                 for dstar in self.unitcell.ringds]
        self.theorytth = transform.degrees(n.array(tths))

    def computegv(self):
        """
        Using self.twotheta, self.eta and omega angles, compute x,y,z of spot
        in reciprocal space
        """
        pars = self.parameterobj.get_parameters()
        self.wavelength = 1.0 # dumb
        self.omegasign = 1.0
        self.wedge=0.
        self.chi=0.
        for a in ['wavelength','omegasign','wedge','chi']:
            if pars.has_key(a):
                setattr(self,a,pars[a])
        
        if self.twotheta is None:
            self.compute_tth_eta()
        self.gv = transform.compute_g_vectors(self.twotheta,
                self.eta, 
                self.omega*self.omegasign, # eek
                self.wavelength,
                wedge=self.wedge,
                chi=self.chi)
        tthnew,etanew,omeganew=transform.uncompute_g_vectors(self.gv,
                self.wavelength,wedge=self.wedge)
        #
        #print "Testing reverse transformations"
        #for i in range(5):
        #    print "tth %8.3f eta %8.3f om %8.3f in "%(self.twotheta[i],self.eta[i],self.omega[i]),
        #    print "tth %8.3f [eta %8.3f om %8.3f or eta %8.3f om %8.3f]"%(tthnew[i],etanew[0][i],omeganew[0][i],etanew[1][i],omeganew[1][i])

    def savegv(self,filename):
        """
        Save g-vectors into a file
        Use crappy .ass format from previous for now (testing)
        """
        #        self.parameterobj.update_other(self)
        if self.gv is None:
            self.computegv()
        if self.unitcell is None:
            self.addcellpeaks(0.)
        f=open(filename,"w")
        f.write(self.unitcell.tostring())
        f.write("\n")
        f.write("# wavelength = %f\n"%( float(self.wavelength) ) )
        f.write("# wedge = %f\n"%( float(self.wedge) ))
        f.write("# ds h k l\n")
        for peak in self.theorypeaks:
            f.write("%10.7f %4d %4d %4d\n"%(peak[0],peak[1][0],peak[1][1],peak[1][2]))
        order = n.argsort(self.twotheta)
        f.write("# xr yr zr xc yc ds phi omega spot3d_id\n")
        print n.maximum.reduce(self.omega),n.minimum.reduce(self.omega)
        ds = 2*n.sin(transform.radians(self.twotheta/2))/self.wavelength
        for i in order:
            f.write("%f %f %f %f %f %f %f %f %d \n"%(self.gv[0,i],self.gv[1,i],self.gv[2,i],
                self.x[i],self.y[i],ds[i],self.eta[i],self.omega[i]*self.omegasign, 
                self.spot3d_id[i]))
        f.close()


    def write_graindex_gv(self,filename):
        from ImageD11 import write_graindex_gv
        if self.gv is None:
            self.computegv()
        self.intensity = self.finalpeaks[3,:]*self.finalpeaks[4,:]
        write_graindex_gv.write_graindex_gv(filename,
                                            self.gv,
                                            self.twotheta,
                                            self.eta,
                                            self.omega*self.omegasign,
                                            self.intensity,
                                            self.unitcell)
        
