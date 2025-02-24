
from __future__ import print_function

# ImageD11_v1.1 Software for beamline ID11
# Copyright (C) 2005 - 2008  Jon Wright
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

import numpy
import logging
import math
import sys, os
from ImageD11 import transform, unitcell, columnfile
from ImageD11.parameters import par, parameters

PARAMETERS = [
     par("omegasign", 1.0,
          helpstring="Sign of the rotation about z " + \
              "(normally +1 for right handed)",
          vary=False,
          can_vary=False),
     par('z_center', 1024.0,
          helpstring="Beam centre in vertical, pixels",
          vary=True,
          can_vary=True,
          stepsize=1.0),
     par('y_center', 1024.0,
          helpstring="Beam centre in horizontal, pixels",
          vary=True,
          can_vary=True,
          stepsize=1.0),
     par('distance', 50000.0,
          helpstring="sample detector distance, same units as pixel size",
          vary=True,
          can_vary=True,
          stepsize=100.0),
     par('z_size', 48.08150,
          helpstring="pixel size in vertical, same units distance",
          vary=False,
          can_vary=True,
          stepsize=0.1), # this could actually vary - a bit crazy?
     par('y_size', 46.77648,
          helpstring="pixel size in horizontal, same units as distance",
          vary=False,
          can_vary=True,
         stepsize=0.1), # this could actually vary - a bit crazy?
     par('tilt_z', 0.0,
          helpstring="detector tilt, right handed around z",
          vary=True,
          can_vary=True,
          stepsize=transform.radians(0.1)),
     par('tilt_y', 0.0,
          helpstring="detector tilt, right handed around y",
          vary=True,
          can_vary=True,
          stepsize=transform.radians(0.1)),
     par('tilt_x', 0.0,
          helpstring="detector tilt, right handed around x",
          vary=False,
          can_vary=True,
          stepsize=transform.radians(0.1)),
     par('fit_tolerance', 0.05,
          helpstring="tolerance to decide which peaks to use",
          vary=False,
          can_vary=False),
     par('wavelength', 0.155,
          helpstring="wavelength, normally angstrom, " + \
             "same as units unit cell ",
          vary=False,
          can_vary=True, # but you'll be lucky!
          stepsize=0.00001),
     par('wedge', 0.0,
          helpstring="wedge, rotation around y under omega",
          vary=False,
          can_vary=True,
          stepsize=transform.radians(0.1)),
     par('chi', 0.0,
          helpstring="wedge, rotation around x under omega",
          vary=False,
          can_vary=True,
          stepsize=transform.radians(0.1)),
     par('cell__a' , 4.1569,
          helpstring="unit cell par, same units as wavelength",
          vary=False,
          can_vary=True,
          stepsize=0.01),
     par('cell__b' , 4.1569,
          helpstring="unit cell par, same units as wavelength",
          vary=False,
          can_vary=True,
          stepsize=0.01),
     par('cell__c' , 4.1569,
          helpstring="unit cell par, same units as wavelength",
          vary=False,
          can_vary=True,
          stepsize=0.01),
     par('cell_alpha' , 90.0,
          helpstring="unit cell par, degrees",
          vary=False,
          can_vary=True,
          stepsize=0.01),
     par('cell_beta' , 90.0,
          helpstring="unit cell par, degrees",
          vary=False,
          can_vary=True,
          stepsize=0.01),
     par('cell_gamma' , 90.0,
          helpstring="unit cell par, degrees",
          vary=False,
          can_vary=True,
          stepsize=0.01),
     par('cell_lattice_[P,A,B,C,I,F,R]', "P",
          helpstring="lattice centering type. Try P if you are not sure",
          vary=False,
          can_vary=False),
     par('o11' , 1,
          helpstring="detector flip element +1 for frelon & quantix",
          vary=False,
          can_vary=False),
     par('o12' , 0,
          helpstring="detector flip element 0 for frelon & quantix",
          vary=False,
          can_vary=False),
     par('o21' , 0,
          helpstring="detector flip element 0 for frelon & quantix",
          vary=False,
          can_vary=False),
     par('o22' , -1,
          helpstring="detector flip element -1 for frelon & +1 for quantix",
          vary=False,
          can_vary=False),
     par('t_x' , 0,
          helpstring="crystal translation, units as distance/pixels",
          vary=False,
          can_vary=True,
          stepsize=1.),
     par('t_y' , 0,
          helpstring="crystal translation, units as distance/pixels",
          vary=False,
          can_vary=True,
          stepsize=1.),
     par('t_z' , 0,
          helpstring="crystal translation, units as distance/pixels",
          vary=False,
          can_vary=True,
          stepsize=1.),
     par('no_bins', 10000,
         helpstring="Number of bins to use in histogram based filters",
         vary=False,
         can_vary=False),
     par('min_bin_prob', 1e-5,
         helpstring="Number of bins to use in histogram based filters",
         vary=False,
         can_vary=False),
     par('weight_hist_intensities', False,
         helpstring="If True or 1, weight histograms by peak intensities. If False or 0, histogram by number of peaks.",
         vary=False,
         can_vary=False),
     ]



class transformer:
    """
    Handles the algorithmic, fitting and state information for 
    fitting parameters to give experimental calibrations
    """
    def __init__(self, parfile=None, fltfile=None):
        """
        Nothing is passed in
        ...will need to loadfileparameters and also peaks
        """
        self.unitcell = None
        # this sets defaults according to class dict.
        self.parameterobj = parameters()
        for p in PARAMETERS:
            self.parameterobj.addpar(p)
        # Interesting - is this an alias for the dict?
        self.pars = self.parameterobj.get_parameters()
        self.colfile = None
        self.xname = 'sc'
        self.yname = 'fc'
        self.omeganame = 'omega'
        self.theoryds = None
        if parfile is not None:
            self.loadfileparameters(parfile)
        if fltfile is not None:
            self.loadfiltered(fltfile)

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
    
    def setfiltered(self, colfile):
        """
        Set self.colfile as colfile
        """
        self.colfile = colfile
        if ( ("sc" in self.colfile.titles) and
             ("fc" in self.colfile.titles) and
             ("omega" in self.colfile.titles)):
            self.setxyomcols("sc", "fc", "omega")
        if (self.colfile.titles[0:3] == ["sc", "fc", "omega"]):
            self.setxyomcols("sc", "fc", "omega")
        if (self.colfile.titles[0:3] == ["xc", "yc", "omega"]):
            self.setxyomcols("xc", "yc", "omega")
        if "spot3d_id" not in self.colfile.titles:
            self.colfile.addcolumn(list(range(self.colfile.nrows)),
                                    "spot3d_id")
    
    def loadfiltered(self, filename):
        """
        Read in 3D peaks from peaksearch
        """
        colfile = columnfile.columnfile(filename)
        self.setfiltered(colfile)

    def setxyomcols(self, xname, yname, omeganame):
        self.xname = xname
        self.yname = yname
        self.omeganame = omeganame
        #logging.warning("titles are %s  %s  %s" % (self.xname,
        #                                         self.yname,
        #                                         self.omeganame))

    def getcols(self):
        return self.colfile.titles

    def loadfileparameters(self, filename, phase_name=None):
        """ Read in beam center etc from file """
        self.parameterobj.loadparameters(filename, phase_name=phase_name)

    def saveparameters(self, filename):
        """ Save beam center etc to file """
        self.parameterobj.saveparameters(filename)

    def applyargs(self, args):
        """ for use with simplex/gof function, alter parameters """
        self.parameterobj.set_variable_values(args)

    def getcolumn(self, name):
        """Return the data"""
        return self.colfile.getcolumn(name)

    def addcolumn(self, col, name):
        """Return the data"""
        return self.colfile.addcolumn(col, name)

    def compute_tth_eta(self):
        """ Compute the twotheta and eta for peaks previous read in """
        if None in [self.xname, self.yname]:
            raise Exception("No peaks loaded")
        peaks = [self.getcolumn(self.xname),
                 self.getcolumn(self.yname)]
        peaks_xyz = transform.compute_xyz_lab(peaks,
                                                **self.parameterobj.get_parameters())
        # Store these in the columnfile
        self.addcolumn(peaks_xyz[0], "xl")
        self.addcolumn(peaks_xyz[1], "yl")
        self.addcolumn(peaks_xyz[2], "zl")
        # Get the Omega name?
        omega = self.getcolumn(self.omeganame)
        tth, eta = transform.compute_tth_eta_from_xyz(
            peaks_xyz,
            omega,
            **self.parameterobj.get_parameters())
        self.addcolumn(tth , "tth")
        self.addcolumn(eta , "eta")
        return tth, eta

    def compute_histo(self, colname):
        """ Compute the histogram over twotheta for peaks previous read in
        Filtering is moved to a separate function
        
        colname is most-often "tth"
        
        other parameters are set in the parameter object
        no_bins = number of bins
        weight_hist_intensities: True or False
            False: histogram by number of measured peaks
            True: weight by peak intensities 
        """
        if colname not in self.colfile.titles:
            raise Exception("Cannot find column " + colname)
        weight = self.parameterobj.get("weight_hist_intensities")
        if (weight):
            bins, hist, hpk = transform.compute_tth_histo(self.getcolumn(colname), weight = True, weights = self.getcolumn("sum_intensity"), **self.parameterobj.get_parameters())
        else:
            bins, hist, hpk = transform.compute_tth_histo(self.getcolumn(colname),
                                             **self.parameterobj.get_parameters())
        self.addcolumn(hpk, colname + "_hist_prob")
        return bins, hist

    def compute_tth_histo(self):
        """ Give hardwire access to tth """
        if "tth" not in self.colfile.titles:
            self.compute_tth_eta()
        return self.compute_histo("tth")

    def filter_min(self, col, minval):
        """
        and filter peaks out falling in bins with less than min_prob
        """
        if "tth_hist_prob" not in self.colfile.titles:
            self.compute_tth_histo()
        mask = self.colfile.getcolumn("tth_hist_prob") > minval
        logging.info("Number of peaks before filtering = %d" % (
            self.colfile.nrows))
        self.colfile.filter(mask)
        logging.info("Number of peaks after filtering = %d" % (
            self.colfile.nrows))

    def tth_entropy(self):
        """
        Compute the entropy of the two theta values
        ... this may depend on the number of tth bins (perhaps?)
        ... to be used for cell parameter free line straightening
        S = \sum -p log p
        """
        if "tth_hist_prob" not in self.colfile.titles:
            self.compute_tth_histo()
        hpk = self.getcolumn("tth_hist_prob")
        lp = numpy.log(hpk)
        entropy = numpy.sum(-hpk * lp)
        return entropy

    def gof(self, args):
        """ Compute how good is the fit of obs/calc peak positions in tth """
        self.applyargs(args)
        # 
        if self.update_fitds:
            cell = unitcell.unitcell_from_parameters( self.parameterobj )
        # Here, pars is a dictionary of name/value pairs to pass to compute_tth_eta
        tth, eta = self.compute_tth_eta()
        w = self.parameterobj.get("wavelength")
        gof = 0.
        npeaks = 0
        for i in range(len(self.tthc)):# (twotheta_rad_cell.shape[0]):
            if self.update_fitds:
                b4 = self.fitds[i]
                self.fitds[i] = cell.ds( self.fithkls[i] )
            self.tthc[i] = transform.degrees(math.asin(self.fitds[i] * w / 2) * 2)
            diff = numpy.take(tth, self.indices[i]) - self.tthc[i]
#         print "peak",i,"diff",maximum.reduce(diff),minimum.reduce(diff)
            gof = gof + numpy.sum(diff * diff)
            npeaks = npeaks + len(diff)
        gof = gof / npeaks
        return gof * 1e6

    def fit(self, tthmin=0, tthmax=180):
        """ Apply simplex to improve fit of obs/calc tth """
        tthmin = float(tthmin)
        tthmax = float(tthmax)
        from . import simplex
        if self.theoryds == None:
            self.addcellpeaks()
        # Assign observed peaks to rings
        self.wavelength = None
        self.indices = []  # which peaks used
        self.tthc = []     # computed two theta values
        self.fitds = []    # hmm?
        self.fithkls = []
        self.fit_tolerance = 1.
        pars = self.parameterobj.get_parameters()
        w = float(pars['wavelength'])
        self.wavelength = w
        self.fit_tolerance = float(pars['fit_tolerance'])
        print("Tolerance for assigning peaks to rings", \
            self.fit_tolerance, ", min tth", tthmin, ", max tth", tthmax)
        tth, eta = self.compute_tth_eta()
        for i in range(len(self.theoryds)):
            dsc = self.theoryds[i]
            tthcalc = math.asin(dsc * w / 2) * 360. / math.pi # degrees
            if tthcalc > tthmax:
                break
            elif tthcalc < tthmin:
                continue
            logicals = numpy.logical_and(numpy.greater(tth,
                                               tthcalc - self.fit_tolerance),
                                     numpy.less(tth ,
                                            tthcalc + self.fit_tolerance))

            if sum(logicals) > 0:
                self.tthc.append(tthcalc)
                self.fitds.append(dsc)
                self.fithkls.append( self.unitcell.ringhkls[dsc][0] )
                ind = numpy.compress(logicals, list(range(len(tth))))
                self.indices.append(ind)
        self.update_fitds = False
        for p in self.parameterobj.varylist:
            if p.startswith('cell'):
                self.update_fitds = True
        guess = self.parameterobj.get_variable_values()
        inc = self.parameterobj.get_variable_stepsizes()
        if len(guess) == 0:
            # There is nothing to fit.
            logging.warning("You try to fit with no variables!?")
            return None
        s = simplex.Simplex(self.gof, guess, inc)
        newguess, error, niter = s.minimize()
        inc = [v / 10 for v in inc]
        guess = newguess
        s = simplex.Simplex(self.gof, guess, inc)
        newguess, error, niter = s.minimize()
        self.parameterobj.set_variable_values(newguess)
        self.wavelength = self.parameterobj.get("wavelength")
        self.gof(newguess)
        print(newguess)
        if self.update_fitds:
            self.addcellpeaks()



    def addcellpeaks(self, limit=None):
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
        cell = [ pars[name] for name in ['cell__a', 'cell__b', 'cell__c',
                                        'cell_alpha', 'cell_beta', 'cell_gamma']]
        lattice = pars['cell_lattice_[P,A,B,C,I,F,R]']
        if "tth" not in self.colfile.titles:
            self.compute_tth_eta()
        # Find last peak in radius
        if limit is None:
            highest = numpy.maximum.reduce(self.getcolumn("tth"))
        else:
            highest = limit
        w = pars['wavelength']
        ds = 2 * numpy.sin(transform.radians(highest) / 2.) / w
        self.dslimit = ds
        self.unitcell = unitcell.unitcell(cell, lattice)
        # If the space group is provided use xfab for generate unique hkls
        if 'cell_sg' in pars:
            self.theorypeaks = self.unitcell.gethkls_xfab(ds, pars['cell_sg'])
            tths = []
            self.theoryds = []
            for i in range(len(self.theorypeaks)):
                tths.append(2 * numpy.arcsin(w * self.theorypeaks[i][0] / 2.))
                self.theoryds.append( self.theorypeaks[i][0] )
        else:
            # HO:  I have removed this part as it seems redundant ringds also calls gethkls 
            # JPW: It was not redundant. theorypeaks is not defined anywhere else and you 
            #      can't write a g-vector file without it.
            self.theorypeaks = self.unitcell.gethkls(ds)
            self.unitcell.makerings(ds)
            self.theoryds = self.unitcell.ringds
            tths = [numpy.arcsin(w * dstar / 2) * 2
                    for dstar in self.unitcell.ringds]
        self.theorytth = transform.degrees(numpy.array(tths))


    def computegv(self):
        """
        Using tth, eta and omega angles, compute x,y,z of spot
        in reciprocal space
        """
        if "tth" not in self.colfile.titles:
            self.compute_tth_eta()
        pars = self.parameterobj.get_parameters()
        if "omegasign" in pars:
            om_sgn = pars["omegasign"]
        else:
            om_sgn = 1.0
        gv = transform.compute_g_vectors(
            self.getcolumn("tth"),
            self.getcolumn("eta"),
            self.getcolumn("omega") * om_sgn,
            self.parameterobj.get('wavelength'),
            wedge=self.parameterobj.get('wedge'),
            chi=self.parameterobj.get('chi'))

        self.addcolumn(gv[0], "gx")
        self.addcolumn(gv[1], "gy")
        self.addcolumn(gv[2], "gz")

    def getaxis(self):
        """
        Compute the rotation axis
        This handles omegasign in a more elegant way
        """
        # unit vector along z
        v = numpy.array([0, 0, 1], float)
        p = self.parameterobj.get("omegasign")
        return v * p

    def savegv(self, filename):
        """
        Save g-vectors into a file
        Use crappy .ascii format from previous for now (testing)
        """
        #        self.parameterobj.update_other(self)
        self.colfile.updateGeometry( self.parameterobj )
        if self.unitcell is None:
            self.addcellpeaks()
        f = open(filename, "w")
        f.write(self.unitcell.tostring())
        f.write("\n")
        pars = self.parameterobj.get_parameters()
        f.write("# wavelength = %f\n" % (float(pars['wavelength'])))
        f.write("# wedge = %f\n" % (float(pars['wedge'])))
        # Handle the axis direction somehow
        f.write("# axis %f %f %f\n" % tuple(self.getaxis()))
        # Put a copy of all the parameters in the gve file
        pkl = list(pars.keys())
        pkl.sort()
        for k in pkl:
            f.write("# %s = %s \n" % (k, pars[k]))
        f.write("# ds h k l\n")
        for peak in self.theorypeaks:
            f.write("%10.7f %4d %4d %4d\n" % (peak[0],
                                              peak[1][0],
                                              peak[1][1],
                                              peak[1][2]))
        tth = self.getcolumn("tth")
        ome = self.getcolumn("omega")
        eta = self.getcolumn("eta")
        gx = self.getcolumn("gx")
        gy = self.getcolumn("gy")
        gz = self.getcolumn("gz")
        x = self.getcolumn(self.xname)
        y = self.getcolumn(self.yname)
        spot3d_id = self.getcolumn("spot3d_id")
        xl = self.getcolumn("xl")
        yl = self.getcolumn("yl")
        zl = self.getcolumn("zl")
        order = numpy.argsort(tth)
        f.write("#  gx  gy  gz  xc  yc  ds  eta  omega  spot3d_id  xl  yl  zl\n")
        print(numpy.maximum.reduce(ome), numpy.minimum.reduce(ome))
        ds = 2 * numpy.sin(transform.radians(tth / 2)) / pars["wavelength"]
        fmt = "%f "*8 + "%d " + "%f "*3 + "\n"
        for i in order:
            f.write(fmt % (gx[i], gy[i], gz[i],
                           x[i], y[i],
                           ds[i], eta[i], ome[i],
                           spot3d_id[i],
                           xl[i], yl[i], zl[i]))
        f.close()

    def write_colfile(self, filename):
        """
        Save out the column file (all info stored we hope)
        """
        self.colfile.parameters = self.parameterobj
        self.colfile.writefile(filename)

    def write_graindex_gv(self, filename):
        from ImageD11 import write_graindex_gv
        if ("gx" not in self.colfile.titles):
            self.computegv()
        gv = [ self.getcolumn("gx"),
               self.getcolumn("gy"),
               self.getcolumn("gz") ]

        if ("sum_intensity" in self.colfile.titles):
            ints = self.getcolumn("sum_intensity")
        elif ("avg_intensity" in self.colfile.titles) and \
             ("Number_of_pixels" in self.colfile.titles):
            ints = self.getcolumn("sum_intensity") * \
                   self.getcolumn("Number_of_pixels")
        elif ("avg_intensity" in self.colfile.titles) and \
             ("npixels" in self.colfile.titles):
            ints = self.getcolumn("sum_intensity") * \
                   self.getcolumn("npixels")
        else:
            ints = numpy.zeros(self.colfile.nrows)

        pars = self.parameterobj.get_parameters()
        if "omegasign" in pars:
            om_sgn = pars["omegasign"]
        else:
            om_sgn = 1.0

        write_graindex_gv.write_graindex_gv(filename,
                                            numpy.array(gv),
                                            self.getcolumn("tth"),
                                            self.getcolumn("eta"),
                                            self.getcolumn("omega") * om_sgn,
                                            ints,
                                            self.unitcell)


    def write_pyFAI(self, filename, tthmin=0, tthmax=180):
        """
        Write file for Jerome Kieffer's pyFAI fitting routine to use and run the refinment...
        """
        try:
            import pyFAI
            import pyFAI.peakPicker
        except ImportError:
            logging.error("pyFAI module not available ... ")
            return
        controlpoints = pyFAI.peakPicker.ControlPoints()
        tthmin = float(tthmin)
        tthmax = float(tthmax)
        self.addcellpeaks()
        # Assign observed peaks to rings
        self.wavelength = None
        self.indices = []  # which peaks used
        self.tthc = []     # computed two theta values
        pars = self.parameterobj.get_parameters()
        w = float(pars['wavelength'])
        tol = float(pars['fit_tolerance'])
        tth, eta = self.compute_tth_eta()
        # Loop over calc peak positions
        z = self.getcolumn('s_raw')
        y = self.getcolumn('f_raw')
        for i, dsc in enumerate(self.theoryds):
            tthcalc = math.asin(dsc * w / 2) * 360. / math.pi # degrees
            if tthcalc > tthmax:
                break
            elif tthcalc < tthmin:
                continue
            ind = abs(tth - tthcalc) <= tol
            if len(ind) > 0:
                controlpoints.append_2theta_deg(list(zip(z[ind], y[ind])), tthcalc)
        controlpoints.save(filename)


        # There is no spline, hence return
        return

        #if len(controlpoints) == 0:
        #    logging.error("THe number of control point found in null !!! skipping optimization ")
        #    return
        #
        #if "spline" in pars:
        #    print("with spline")
        #    #TODO: where is the spline file stored ???
        #    geoRef = pyFAI.geometryRefinement.GeometryRefinement(controlpoints.getList(), dist=0.1, splineFile=self.splineFile)
        #else:
        #    geoRef = pyFAI.geometryRefinement.GeometryRefinement(controlpoints.getList(), dist=0.1, pixel1=pars["z_size"] * 1e-6, pixel2=pars["y_size"] * 1e-6)
        #geoRef.wavelength = w * 1e-10
        #previous = sys.maxint
        #while previous > geoRef.chi2():
        #    previous = geoRef.chi2()
        #    geoRef.refine2(1000000)
        #    print geoRef
        #geoRef.save(os.path.splitext(filename)[0] + ".poni")

    def save_tth_his(self, filename, bins, hist):
        """
        Saves a 2 theta histogram into a file
        bins and histogram should have been calculated previously
        """
        #        self.parameterobj.update_other(self)

        f = open(filename, "w")
        f.write("# Peaks: %s\n" % self.colfile.filename)
        f.write("# N. bins: %d\n" % len(bins))
        f.write("# 2tth intensity\n")
        for i in range(0,len(hist)):
            f.write("%10.7f  %.7g \n" % (bins[i], hist[i]))
        f.close()
