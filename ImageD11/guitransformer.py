

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

from Numeric import *
from Tkinter import *

from listdialog import listdialog

from ImageD11 import twodplot

class guitransformer:

    def __init__(self,parent,quiet="No"):
        """
        Parent is a hook to features of the parent gui
        """
        self.quiet=quiet
        self.parent=parent
        self.nbins = 10000
        self.min_bin_ratio = 0.10
        self.menuitems = ( "Transformation", 0,
                           [ ( "Load filtered peaks", 0, self.loadfiltered),
                             ( "Plot y/z", 5, self.plotyz     ),
                             ( "Load parameters", 1, self.loadfileparameters),
                             ( "Edit parameters", 0, self.editparameters),
                             ( "Plot tth/eta", 0, self.plotreta ),
                             ( "Add unit cell peaks",0, self.addcellpeaks),
                             ( "Fit",0, self.fit),
                             ( "Save parameters", 0, self.saveparameters),
                             ( "Plot tth histogram", 0, self.plothisto ),
                             ( "Filter peaks based on tth histogram", 0, self.filterhisto ),

#                             ( "Set axis orientation", 0, self.setaxisorientation),
                             ( "Compute g-vectors", 0, self.computegv),
                             ( "Save g-vectors", 0, self.savegv),
                             ( "Write graindex finalpeaks.log",0, self.write_graindex_gv)
                           ] )

    def loadfiltered(self):
        filename=self.parent.opener.show(title= \
                "File containing filtered peaks")
        self.parent.guicommander.execute("transformer",
                "loadfiltered",filename)

    def loadfileparameters(self):
        filename=self.parent.opener.show(title= \
                "File containing detector parameters")
        self.parent.guicommander.execute("transformer",
                "loadfileparameters",filename)

    def saveparameters(self,filename=None):
        if filename==None:
            filename=self.parent.saver.show(title="File to save detector parameters")
        self.parent.guicommander.execute("transformer",
                "saveparameters",filename)


    def editparameters(self):
        """
        Gets a copy of the parameter object
        Allows user to edit parameters
        """
        self.parent.guicommander.execute("transformer","updateparameters")
        pars = self.parent.guicommander.getdata("transformer","pars")
        import logging
        logging.debug("transformer pars: %s"% (str(pars)))
        d=listdialog(self.parent,items=pars,title="Detector parameters")
        self.parent.guicommander.execute("transformer",
                "parameterobj.set_parameters",d.result)
        
    def plotyz(self):
        """
        Plots the x,y arrays being used
        """
        d=self.parent.guicommander.getdata("transformer","finalpeaks")
        self.parent.twodplotter.hideall()
        self.parent.twodplotter.adddata(
              ( "Filtered peaks",
                 twodplot.data(
                    d[0,:],
                    d[1,:],
                    { "xlabel" : "y", 
                      "ylabel" : "z", 
                      "title"  : "Peak positions on detector"} ) ) )


    def fit(self):
        tthmax = self.parent.twodplotter.a.get_xlim()[1]
        self.parent.guicommander.execute("transformer","fit",tthmax)
        # maybe update plot?

    def plotreta(self):
        self.parent.guicommander.execute("transformer","compute_tth_eta")
        tth = self.parent.guicommander.getdata("transformer","twotheta")
        eta = self.parent.guicommander.getdata("transformer","eta")
#       self.parent.twodplotter.hideall()
        self.parent.twodplotter.adddata(
              ( "2Theta/Eta",
                 twodplot.data(
                    tth,
                    eta,
                    {"xlabel":"TwoTheta / degrees",
                     "ylabel":"Azimuth / degrees",
                     "title" :"Peak positions"}
                     )))

    def plothisto(self):
        d=listdialog(self.parent,items={"no_bins": self.nbins},title="Histogram - no of bins")
        self.nbins = d.result['no_bins']
        self.parent.guicommander.execute("transformer","parameterobj.set_parameters",d.result)
        tth = self.parent.guicommander.getdata("transformer","twotheta")
        self.parent.twodplotter.adddata(
              ( "2Theta/Eta",
                 twodplot.data(
                    None,
                    tth,
                    #self.histogram,
                    {"xlabel":"TwoTheta / degrees",
                     "ylabel":"No in bin",
                     "title" :"TwoTheta histogram",
                     "plottype" :"hist"}
                     )))

    def filterhisto(self):
        d=listdialog(self.parent,items={"no_bins": self.nbins, "min_bin_ratio": self.min_bin_ratio},title="Histogram filter")
        self.parent.guicommander.execute("transformer","parameterobj.set_parameters",d.result)
        self.parent.guicommander.execute("transformer","compute_tth_histo")


    def addcellpeaks(self):
        self.parent.guicommander.execute("transformer","addcellpeaks")
        tth=self.parent.guicommander.getdata("transformer","theorytth")
        self.parent.twodplotter.adddata(
              ( "HKL peaks",
                 twodplot.data(
                         tth,
                         zeros(tth.shape[0]),
                         {'color':'r',
                          'pointtype':'|'
                          }
                 )))

    def computegv(self):
        self.parent.guicommander.execute("transformer","computegv")

    def savegv(self):       
        filename=self.parent.saver.show(title="File to save gvectors")
        self.parent.guicommander.execute("transformer","savegv",filename)

    def write_graindex_gv(self):
        filename=self.parent.saver.show(title="File for graindex, try finalpeaks.log")
        self.parent.guicommander.execute("transformer","write_graindex_gv",filename)

