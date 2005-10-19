

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

import time,math,sys

import factors


class guifactors:
                           
   def __init__(self,parent):
      """
      Parent is a hook to features of the parent gui
      """
      self.parent=parent
# peaks are in      self.parent.finalpeaks
      self.menuitems = ( "Factors", 0,
                         [ ( "Load chi files", 0, self.loadchis),
                           ( "Plot a chi", 0, self.plotchi     ),
			   ( "Generate SVD", 0, self.svd ),
			   ( "Select number of factors", 0, self.setnfactors),
			   ( "Generate differences", 9, self.dataminuscalc)
                         ] ) 
      self.factors=factors.factors()
	
   def loadchis(self):
      filename=self.parent.opener.show(title="Name of a chi file")
      # globs by default
      self.factors.loadchis(filename)

   def plotchi(self):
      # which one?
      self.parent.twodplotter.adddata(
      	( "Chi file", twodplot.data(self.factors.x,self.factors.y[i],
	{ "xlabel":"probably 2theta",
	  "ylabel":"probably intensity",
	  "title":"chifile %d"%(i) } ) )

   def svd(self):
      self.factors.svd()
      # Give scree plot

   def setnfactors(self):
      self.factors.setnfactors(i)

   def dataminuscalc(self):
      self.factors.generate()
      
