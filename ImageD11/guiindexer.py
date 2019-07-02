

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


"""
Tkinter wrapper for ImageD11 indexing object

All communication should be via parent guicommander object

Owner of the plot3d window
"""

from .listdialog import listdialog
from . import twodplot
import threading
try:
    import Tkinter as Tk
except:
    import tkinter as Tk


class run_idx:
    def __init__(self, parent, idxer):
        self.top = Tk.Toplevel(parent)
        self.label = Tk.Label(self.top, text="Found so far %8d"%(0))
        self.label.pack()
        self.btn = Tk.Button(self.top, text="Stop", command=self.stop)
        self.btn.pack()
        self.idxer = idxer
        self.thr = threading.Thread( target=self.idxer.score_all_pairs )
        self.thr.start()
        self.top.after(1000, self.update)

    def stop(self):
        self.label.configure(text="Got a stop, so far %8d"%(len(self.idxer.ubis)))
        self.btn.configure(text='Stopping')
        self.idxer.stop=True
        self.thr.join()
        self.top.destroy()

    def update(self):
        self.label.configure( text="Tested %8d of %8d\nFound so far %8d"%(
            self.idxer.tried,
            self.idxer.npairs,
            len(self.idxer.ubis)) )
        if self.idxer.tried == self.idxer.npairs:
            self.top.destroy()
        else:
            self.top.after(1000, self.update)


class guiindexer:

    def __init__(self,parent):
        """
        Parent (arg) is a hook to features of the parent gui
        sets up indexing menuitems
        """
        self.parent=parent
# peaks are in      self.parent.finalpeaks
        self.menuitems = ( "Indexing", 0,
                           [ ( "Load g-vectors", 0, self.loadgv),
                             ( "Plot x/y/z", 5, self.plotxyz),
                             ( "Load parameters", 1, self.loadfileparameters),
                             ( "Edit parameters", 0, self.editparameters),
                             ( "Assign peaks to powder rings", 0, self.assignpeaks),
                             ( "Make Friedel pair file", 5, self.makefriedel),
                             ( "Generate trial orientations",0, self.find),
                             ( "Score trial orientations",0, self.scorethem),
                             ( "Auto-find",2, self.autofind),
                             ( "Histogram fit quality",0, self.histogram_drlv_fit),
                             ( "Save parameters", 0, self.saveparameters),
                             ( "Save UBI matrices", 5, self.saveubis),
                             ( "Write out indexed peaks",0,self.saveindexing),
                             ( "Reset indexer",0,self.reset),
                           ] )
        self.plot3d=None

    def autofind(self):
        self.parent.guicommander.commandscript += "myindexer.score_all_pairs()\n"
        run_idx( self.parent, self.parent.guicommander.objects['indexer'])


    def loadgv(self):
        """ see indexing.readgvfile """
        filename=self.parent.opener.show(
            title="File containing g-vectors",
            filetypes=[ ("Gvector files", "*.gve"),
                        ("Gvector files", "*.gv") ] )
        self.parent.guicommander.execute("indexer","readgvfile",filename)

    def saveubis(self):
        """ see indexing.saveubis """
        filename=self.parent.saver.show(title="File to save UBIS")
        self.parent.guicommander.execute("indexer","saveubis",filename)

    def makefriedel(self):
        """ see indexing.friedelpairs """
        filename=self.parent.saver.show(title="File to save Friedelpairs")
        self.parent.guicommander.execute("indexer","friedelpairs",filename)

    def scorethem(self):
        """ see indexing.scorethem """
        self.parent.guicommander.execute("indexer","scorethem")

    def histogram_drlv_fit(self):
        """
        Calls indexer.histogram_drlv_fit
        Plots indexer.bins versus indexer.histogram
        """
        self.parent.guicommander.execute("indexer","histogram_drlv_fit")
        x=self.parent.guicommander.getdata("indexer","bins")
        y=self.parent.guicommander.getdata("indexer","histogram")
        self.parent.twodplotter.plotitems={} # clears plot
        import matplotlib.cm
        for yline in range(y.shape[0]):
            color = matplotlib.cm.jet( yline*1.0/y.shape[0] )
            print("yline, color",yline,color)
            self.parent.twodplotter.plotitems["drlv histogram"+str(yline)]=twodplot.data(
                     x[1:],y[yline,:],
                   {"xlabel" : "drlv",
                    "ylabel" : "freq",
                    "title"  : "drlv histogram",
                    'plotopts' : { "linestyle" : "-",
                                   "marker" : "o",
                                   "markersize" : 1,
                                   "alpha" : 0.8,
                                   "color" : color,
                                   }
                    }
                   ) # data
        self.parent.twodplotter.replot()


    def assignpeaks(self):
        """ see indexing.assigntorings """
        self.parent.guicommander.execute("indexer","assigntorings")

    def loadfileparameters(self):
        """ see indexing.loadpars and parameters.loadpars """
        filename=self.parent.opener.show(
            title="File containing indexing parameters",
            filetypes = [ ("Parameter files", "*.prm"),
                          ("Parameter files", "*.par") ] )
        self.parent.guicommander.execute("indexer","loadpars",filename)

    def saveparameters(self):
        """ see indexing.savepars and parameters.savepars """
        filename=self.parent.saver.show(title="File to save indexing parameters")
        self.parent.guicommander.execute("indexer","savepars",filename)


    def editparameters(self):
        """
        Has the indexing object update its parameter object
               eg : savepars(None)
        gets a copy of the parameter object
        Allows user to edit parameters
        Has the indexing object update itself from the parameter object
               eg : loadpars(None)
        """
        # First make the indexer update its parameters object
        self.parent.guicommander.execute("indexer","updateparameters") # no filename arg
        # Now borrow a copy to read them and edit
        pars = self.parent.guicommander.getdata("indexer","pars")
        d=listdialog(self.parent,items=pars,title="Indexing parameters")
        self.parent.guicommander.execute("indexer","parameterobj.set_parameters",d.result) 
        self.parent.guicommander.execute("indexer","loadpars") # and use them 
        



    def plotxyz(self):
        """
        Gets gv from indexing object
        Plots the x,y,z (gv) array in a 3D opengl window
        """
        import logging
        try:
            from . import plot3d
        except:
            import traceback
            traceback.print_last()
            logging.warning("You might have a PyOpenGl problem??")
            return
        gv = self.parent.guicommander.getdata("indexer","gv")
        if gv is not None:
            if self.plot3d==None:
                self.plot3d = plot3d.plot3d(self.parent,gv)
                self.plot3d.go()
                logging.debug("self.plot3d " + str(self.plot3d))
            else:
                self.plot3d.changedata(gv)

    def find(self):
        """ see indexing.find """
        self.parent.guicommander.execute("indexer","find")

    def saveindexing(self):
        """ see indexing.saveindexing """
        filename=self.parent.saver.show(title="File to save indexing output")
        self.parent.guicommander.execute("indexer","saveindexing",filename)


    def reset(self):
        """ see indexing.reset """
        self.parent.guicommander.execute("indexer","reset")
        
