
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


import numpy as np
#try:
#    from Tkinter import *
#except:
#    # python3?
#    from tkinter import *

from .listdialog import listdialog, columnchooser
from . import twodplot


import logging

class guitransformer:

    def __init__(self,parent,quiet="No"):
        """
        Parent is a hook to features of the parent gui
        """
        self.quiet=quiet
        self.parent=parent
        self.menuitems = ( "Transformation", 0,
            [ ( "Load filtered peaks", 0, self.loadfiltered),
              ( "Plot y/z", 5, self.plotyz     ),
              ( "Load parameters", 1, self.loadfileparameters),
              ( "Edit parameters", 0, self.editparameters),
              ( "Plot tth/eta", 0, self.plotreta ),
              ( "Add unit cell peaks",0, self.addcellpeaks),
              ( "Fit",0, self.fit),
              ( "Save parameters", 0, self.saveparameters),
#              ( "Plot selected columns", 0, self.plotcols),
              ( "Plot tth histogram", 0, self.plothisto ),
              ( "Export tth histogram", 0, self.savehisto ), 
              ( "Filter peaks based on tth histogram", 0, self.filterhisto ),
              ( "Compute g-vectors", 0, self.computegv),
              ( "Save g-vectors", 0, self.savegv),
              ( "Save new colfile", 0, self.savecolfile),
              ( "Write graindex finalpeaks.log",0, self.write_graindex_gv),
              ( "Write pyFAI data file", 0, self.write_pyFAI)
              ] )

    def loadfiltered(self):
        filename=self.parent.opener.show(title=
                 "File containing filtered peaks",
                  filetypes=[("filtered peaks", "*.flt"),
                             ("All Files ", "*")])
        self.parent.guicommander.execute("transformer",
                "loadfiltered",filename)

    def loadfileparameters(self):
        filename=self.parent.opener.show(title=
                "File containing detector parameters",
                filetypes=[("parameter files", "*.par"),
                           ("parameter files", "*.pars"),
                           ("parameter files", "*.prm"),
                           ("All Files ", "*")])
        self.parent.guicommander.execute("transformer",
                "loadfileparameters",filename)

    def saveparameters(self,filename=None):
        if filename==None:
            filename=self.parent.saver.show(title=
                         "File to save detector parameters")
        self.parent.guicommander.execute("transformer",
                "saveparameters",filename)


    def editparameters(self):
        """
        Gets a copy of the parameter object
        Allows user to edit parameters
        """
        self.parent.guicommander.execute("transformer","updateparameters")
        pars = self.parent.guicommander.getdata("transformer","pars")
        vars = self.parent.guicommander.execute("transformer","getvars")
        possvars = self.parent.guicommander.execute("transformer",
                                                    "get_variable_list")
        logging.debug("possible variables "+str(possvars))
        # wtf?
        logic = {}
        for v in possvars:
            if v in vars:
                logic[v]=1
            else:
                logic[v]=0
        logging.debug("transformer pars: %s"% (str(pars)))
        d = listdialog(self.parent,items=pars,title="Detector parameters",
                       logic=logic)
        self.parent.guicommander.execute("transformer",
                                         "parameterobj.set_parameters",
                                         d.result)
        # wtf d.fv
        vars = []
        print("d.fv",d.fv)
        for v in possvars:
            logging.debug(str(v)+" "+str(d.fv[v]))
            if d.fv[v]==1:
                vars.append(v)
        logging.debug("vars: "+str(vars))
        self.parent.guicommander.execute("transformer",
                "parameterobj.set_varylist",vars)

    def plotyz(self):
        """
        Plots the x,y arrays being used
        """
        xname = self.parent.guicommander.getdata("transformer","xname")
        yname = self.parent.guicommander.getdata("transformer","yname")
        x = self.parent.guicommander.execute("transformer",
                                             "getcolumn", xname )
        y = self.parent.guicommander.execute("transformer",
                                             "getcolumn", yname )
        self.parent.twodplotter.hideall()
        self.parent.twodplotter.adddata(
            ( "Filtered peaks",
              twodplot.data(
                  x, y,
                  { "xlabel" : xname,
                    "ylabel" : yname,
                    "title"  : "Peak positions in array",
                    'plotopts' : {'color':'g',
                                  'marker':'.',
                                  'markersize': 1,
                                  'linestyle' : 'none',
                                  'alpha':0.8}
                }
              )))

    def chooseyz(self):
        """
        choose the columns to use for x / y on detector
        """
        pass

    def plotcols(self):
        names = self.parent.guicommander.execute("transformer","getcols")
        print(names)
        d = columnchooser(self.parent, names)
        print(d.result)


    def fit(self):
        tthmin = self.parent.twodplotter.a.get_xlim()[0]
        tthmax = self.parent.twodplotter.a.get_xlim()[1]
        self.parent.guicommander.execute("transformer","fit",tthmin,tthmax)
        self.plotreta()

    def plotreta(self):
        self.parent.guicommander.execute("transformer","compute_tth_eta")
        tth = self.parent.guicommander.execute("transformer","getcolumn","tth")
        eta = self.parent.guicommander.execute("transformer","getcolumn","eta")
        self.parent.twodplotter.adddata(
              ( "2Theta/Eta",
                 twodplot.data(
                    tth,
                    eta,
                    {"xlabel":"TwoTheta / degrees",
                     "ylabel":"Azimuth / degrees",
                     "title" :"Peak positions",
                     'plotopts' : {'color':'g',
                                   'marker':'.',
                                   'markersize': 1,
                                   'linestyle' : 'none',
                                   'alpha':0.8}
}
                     )))

    def plothisto(self, nbins = None):
        if nbins is None:
            nbins = self.parent.guicommander.execute("transformer",
                                                     "parameterobj.get",
                                                     "no_bins")
            doweight = self.parent.guicommander.execute("transformer",
                                                     "parameterobj.get",
                                                     "weight_hist_intensities")
            d = listdialog( self.parent,
                            items={"no_bins": nbins, "weight_hist_intensities": doweight},
                            title="Histogram - no of bins")

            nbins = int(d.result['no_bins'])
            doweight = int(d.result['weight_hist_intensities'])

            self.parent.guicommander.execute("transformer",
                                             "parameterobj.set_parameters",
                                             d.result)
        else:
            self.parent.guicommander.execute("transformer",
                                             "parameterobj.set",
                                             "no_bins", nbins)

        bins, hist = self.parent.guicommander.execute("transformer",
                                                      "compute_tth_histo")
        self.parent.twodplotter.adddata(
              ( "2Theta/Eta",
                 twodplot.data(
                    bins,
                    hist,
                    {"xlabel":"TwoTheta / degrees",
                     "ylabel":"No in bin",
                     "title" :"TwoTheta histogram",

                     }
                     )))

    def savehisto(self, nbins = None):
        if nbins is None:
            nbins = self.parent.guicommander.execute("transformer",
                                                     "parameterobj.get",
                                                     "no_bins")
            doweight = self.parent.guicommander.execute("transformer",
                                                     "parameterobj.get",
                                                     "weight_hist_intensities")
            d = listdialog( self.parent,
                            items={"no_bins": nbins, "weight_hist_intensities": doweight},
                            title="Histogram - no of bins")

            nbins = int(d.result['no_bins'])
            doweight = int(d.result['weight_hist_intensities'])

            self.parent.guicommander.execute("transformer",
                                             "parameterobj.set_parameters",
                                             d.result)
        else:
            self.parent.guicommander.execute("transformer",
                                             "parameterobj.set",
                                             "no_bins", nbins)

        bins, hist = self.parent.guicommander.execute("transformer",
                                                      "compute_tth_histo")
        filename=self.parent.saver.show(title="File to save histogram")
        self.parent.guicommander.execute("transformer","save_tth_his",filename,bins,hist)


    def filterhisto(self):
        """ Call plot histo, then filter on it """
        nbins = self.parent.guicommander.execute("transformer",
                                                 "parameterobj.get",
                                                 "no_bins")

        min_bin_prob = self.parent.guicommander.execute("transformer",
                                                        "parameterobj.get",
                                                        "min_bin_prob")
        
        doweight = self.parent.guicommander.execute("transformer",
                                                     "parameterobj.get",
                                                     "weight_hist_intensities")

        d=listdialog(self.parent,items={
            "no_bins": nbins,
            "min_bin_prob": min_bin_prob,  "weight_hist_intensities": doweight},
                     title="Histogram filter")


        self.parent.guicommander.execute("transformer",
                                         "parameterobj.set_parameters",
                                         d.result)

        min_bin_prob = self.parent.guicommander.execute("transformer",
                                                        "parameterobj.get",
                                                        "min_bin_prob")

        self.plothisto(nbins)
        self.parent.guicommander.execute("transformer",
                                         "filter_min",
                                         "tth_hist_prob",
                                         min_bin_prob)


    def addcellpeaks(self):
        self.parent.guicommander.execute("transformer","addcellpeaks")
        tth=self.parent.guicommander.getdata("transformer","theorytth")
        self.parent.twodplotter.adddata(
              ( "HKL peaks",
                 twodplot.data(
                         tth,
                         np.zeros(tth.shape[0]),
                     {'plotopts' : {'color':'r',
                                    'marker':'|',
                                    'markersize': 50,
                                    'linestyle' : 'none',
                                    'alpha':1.0}
                          }
                 )))

    def computegv(self):
        self.parent.guicommander.execute("transformer","computegv")

    def savegv(self):
        filename=self.parent.saver.show(title="File to save gvectors")
        self.parent.guicommander.execute("transformer","savegv",filename)

    def savecolfile(self):
        filename=self.parent.saver.show(title="File to save newcolfile")
        self.parent.guicommander.execute("transformer",
                                         "write_colfile",
                                         filename)


    def write_graindex_gv(self):
        filename=self.parent.saver.show(title="File for graindex, try finalpeaks.log")
        self.parent.guicommander.execute("transformer","write_graindex_gv",filename)



    def write_pyFAI(self):
        tthmin = self.parent.twodplotter.a.get_xlim()[0]
        tthmax = self.parent.twodplotter.a.get_xlim()[1]
        filename=self.parent.saver.show(title="File for pyFAI, try data.py")
        self.parent.guicommander.execute("transformer","write_pyFAI",filename,
                                         tthmin,tthmax)

