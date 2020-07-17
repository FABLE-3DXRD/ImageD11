
# Get Strain/Stress from ImageD11 UBI/map files
# Copyright (C) 2015 Younes ELHACHI
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


class guisolver:

    def __init__(self,parent):
        """
        Parent (arg) is a hook to features of the parent gui
        sets up eps_sig_solver menuitems
        """
        self.parent=parent

        self.menuitems = ( "Strain/Stress", 0,
                           [ ( "Load ubis", 0, self.loadubis),
                             ( "Load parameters", 1, self.loadfileparameters),
                             ( "Edit parameters", 0, self.editparameters),
                             ( "Compute and save strain and stress", 0, self.compute_save_epsig),
                             ( "Save parameters", 0, self.saveparameters)
                           ] )


    def loadubis(self):
        """ see eps_sig_solver.loadmap """
        filename=self.parent.opener.show(
            title="File containing ubi matrices",
            filetypes=[ ("Grain map files", "*.map"),
                        ("UBI files", "*.ubi") ] )
        self.parent.guicommander.execute("solver","loadmap",filename)

    def compute_save_epsig(self):
        """ see eps_sig_solver.compute_eps_sig """
        filename=self.parent.saver.show(title="File to save strain and stress")
        self.parent.guicommander.execute("solver","compute_write_eps_sig",filename)

    def loadfileparameters(self):
        """ see eps_sig_solver.loadpars and parameters.loadpars """
        filename=self.parent.opener.show(
            title="File containing unit cell and elastic constants",
            filetypes = [ ("Parameter files", "*.prm"),
                          ("Parameter files", "*.par") ] )
        self.parent.guicommander.execute("solver","loadpars",filename)

    def saveparameters(self):
        """ see eps_sig_solver.savepars and parameters.savepars """
        filename=self.parent.saver.show(title="File to save unit cell and elastic constants")
        self.parent.guicommander.execute("solver","savepars",filename)


    def editparameters(self):
        """
        Has the eps_sig_solver object update its parameter object
               eg : savepars(None)
        gets a copy of the parameter object
        Allows user to edit parameters
        Has the eps_sig_solver object update itself from the parameter object
               eg : loadpars(None)
        """
        # First make the eps_sig_solver update its parameters object
        self.parent.guicommander.execute("solver","updateparameters") # no filename arg
        # Now borrow a copy to read them and edit
        pars = self.parent.guicommander.getdata("solver","pars")
        d=listdialog(self.parent,items=pars,title="unit cell and elastic constants")
        self.parent.guicommander.execute("solver","parameterobj.set_parameters",d.result) 
        self.parent.guicommander.execute("solver","loadpars") # and use them 
        
