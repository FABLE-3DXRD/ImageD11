

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

import time,math,unitcell,sys

import indexing

from listdialog import listdialog


class guiindexer:
                           
   def __init__(self,parent):
      """
      Parent is a hook to features of the parent gui
      """
      self.parent=parent
# peaks are in      self.parent.finalpeaks
      self.menuitems = ( "Indexing", 0,
                         [ ( "Load g-vectors", 0, self.loadgv),
                           ( "Plot x/y/z", 5, self.plotxyz     ),
                           ( "Load parameters", 1, self.loadfileparameters),
                           ( "Edit parameters", 0, self.editparameters),
                           ( "Assign peaks to powder rings", 0, self.assignpeaks),
                           ( "Make Friedel pair file", 5, self.makefriedel),
                           ( "Generate trial orientations",0, self.find),
                           ( "Score trial orientations",0, self.scorethem),
                           ( "Histogram fit quality",0, self.histogram_drlv_fit),
                           ( "Save parameters", 0, self.saveparameters),
                           ( "Save UBI matrices", 5, self.saveubis),
                           ( "Write out indexed peaks",0,self.saveindexing )
                         ] ) 
      self.parameters={}
      self.indexer=indexing.indexer(unitcell=self.parent.unitcell,gv=self.parent.gv)
      self.loadparameters()
      self.plot3d=None

   def loadgv(self):
#      if self.parent.gv!=None:
#         from tkMessageBox import askyesno
#         if not askyesno("Overwrite current g-vectors?","Loading new will 
#overwrite the ones in memory now, are you sure?"):
#            return
      filename=self.parent.opener.show(title="File containing g-vectors")
      self.indexer.readgvfile(filename)

   def saveubis(self):
      filename=self.parent.saver.show(title="File to save UBIS")
      self.indexer.saveubis(filename)

   def makefriedel(self):
      filename=self.parent.saver.show(title="File to save Friedels")
      self.indexer.friedelpairs(filename)



   def scorethem(self):
      self.indexer.scorethem()

   def histogram_drlv_fit(self):
      self.indexer.histogram_drlv_fit()
      x=self.indexer.bins
      y=self.indexer.histogram
      self.parent.twodplotter.plotitems={}
      from ImageD11 import twodplot
      for yline in range(y.shape[0]):
         self.parent.twodplotter.plotitems["drlv histogram"+str(yline)]=twodplot.data(
                  x[1:],y[yline,:],
                {"xlabel" : "drlv",
                 "ylabel" : "freq",
                 "title"  : "drlv histogram",
                 "pointtype" : "-"
                 }
                ) # data
      self.parent.twodplotter.replot()

      
   def assignpeaks(self):
      self.indexer.assigntorings()

   def loadfileparameters(self):
      filename=self.parent.opener.show(title="File containing indexing parameters")
      self.loadparameters(filename)

   def loadparameters(self,filename=None):   
      if filename==None:
         self.parameters['cosine_tol']=self.indexer.cosine_tol
         self.parameters['hkl_tol']=self.indexer.hkl_tol
         self.parameters['ring_1']=self.indexer.ring_1
         self.parameters['ring_2']=self.indexer.ring_2
         self.parameters['minpks']=self.indexer.minpks
         self.parameters['uniqueness']=self.indexer.uniqueness
         self.parameters['ds_tol']=self.indexer.ds_tol
         self.parameters['wavelength']=self.indexer.wavelength
      else:
         try:
            lines = open(filename,"r").readlines()
            for line in lines:
               name,value=line.split()
               self.parameters[name]=value
            self.indexer.cosine_tol=float(self.parameters['cosine_tol'])
            self.indexer.hkl_tol= float(self.parameters['hkl_tol'])
            self.indexer.ring_1=int(self.parameters['ring_1'])
            self.indexer.ring_2=int(self.parameters['ring_2'])
            self.indexer.minpks=int(self.parameters['minpks'])
            self.indexer.uniqueness=float(self.parameters['uniqueness'])
            self.indexer.ds_tol=float(self.parameters['ds_tol'])
            self.indexer.wavelength=float(self.parameters['wavelength'])
         except IOError:
            from tkMessageBox import showinfo
            showinfo("Sorry","Could not open/interpret your file")
               

   def saveparameters(self,filename=None):
      if filename==None:
         filename=self.parent.saver.show(title="File to save indexing parameters")
      f=open(filename,"w")
      keys=self.parameters.keys()
      keys.sort()
      for key in keys:
         try:
            f.write("%s %f\n"%(key,self.parameters[key]))
         except:
            f.write("%s %s\n"%(key,self.parameters[key]))
      f.close()

   def saveparameters(self,filename=None):
      if filename==None:
         filename=self.parent.saver.show(title="File to save indexing parameters")
      f=open(filename,"w")
      keys=self.parameters.keys()
      keys.sort()
      for key in keys:
         try:
            f.write("%s %f\n"%(key,self.parameters[key]))
         except:
            f.write("%s %s\n"%(key,self.parameters[key]))
      f.close()




   def editparameters(self):
      self.parameters['cosine_tol']=self.indexer.cosine_tol
      self.parameters['hkl_tol']=self.indexer.hkl_tol
      self.parameters['ring_1']=self.indexer.ring_1
      self.parameters['ring_2']=self.indexer.ring_2
      self.parameters['minpks']=self.indexer.minpks
      self.parameters['uniqueness']=self.indexer.uniqueness
      self.parameters['ds_tol']=self.indexer.ds_tol
      self.parameters['wavelength']=self.indexer.wavelength
      self.parameters['eta_range']=self.indexer.eta_range
      d=listdialog(self.parent,items=self.parameters,title="Indexing parameters")
      self.parameters=d.result
      self.indexer.cosine_tol=float(self.parameters['cosine_tol'])
      self.indexer.hkl_tol= float(self.parameters['hkl_tol'])
      self.indexer.ring_1=int(self.parameters['ring_1'])
      self.indexer.ring_2=int(self.parameters['ring_2'])
      self.indexer.minpks=int(self.parameters['minpks'])
      self.indexer.ds_tol=float(self.parameters['ds_tol'])
      self.indexer.uniqueness=float(self.parameters['uniqueness'])
      self.indexer.eta_range=float(self.parameters['eta_range'])
      try:
         self.indexer.wavelength=float(self.parameters['wavelength'])
      except:
         self.indexer.wavelength=None
      
      
   def plotxyz(self):
      """
      Plots the x,y arrays being used
      """
      import plot3d
      if self.indexer.gv!=None:
         if self.plot3d==None:
            self.plot3d = plot3d.plot3d(self.parent,self.indexer.gv)
            self.plot3d.go()
            print self.plot3d
         else:
            self.plot3d.changedata(self.indexer.gv)


   def find(self):
      self.indexer.find()


   def saveindexing(self,filename=None):
      if filename==None:
         filename=self.parent.saver.show(title="File to save indexing parameters")
      self.indexer.saveindexing(filename)
