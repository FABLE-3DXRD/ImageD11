

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

import transform

from listdialog import listdialog



class guitransformer:
                           
   def __init__(self,parent,quiet="No"):
      """
      Parent is a hook to features of the parent gui
      """
      self.quiet=quiet
      self.parent=parent
# peaks are in      self.parent.finalpeaks
      self.menuitems = ( "Transformation", 0,
                         [ ( "Load filtered peaks", 0, self.loadfiltered),
                           ( "Plot y/z", 5, self.plotyz     ),
                           ( "Load parameters", 1, self.loadfileparameters),
                           ( "Edit parameters", 0, self.editparameters),
                           ( "Plot tth/eta", 0, self.plotreta ),
                           ( "Add unit cell peaks",0, self.addcellpeaks),
                           ( "Fit",0, self.fit),
                           ( "Save parameters", 0, self.saveparameters),
                           ( "Set axis orientation", 0, self.setaxisorientation),
                           ( "Compute g-vectors", 0, self.computegv),
                           ( "Save g-vectors", 0, self.savegv),
                           ( "Write graindex finalpeaks.log",0, self.write_graindex_gv)
                         ] ) 
      self.twotheta=None
      self.parameters={}
      self.wedge=0.0
      self.loadparameters()

   def loadfiltered(self,filename=None):
      #      if self.parent.finalpeaks!=None:
      #         from tkMessageBox import askyesno
      #         if not askyesno("Overwrite current peaks?","Loading new peaks 
      #will overwrite the ones in memory now, are you sure?"):
      #            return
      if filename==None: 
         filename=self.parent.opener.show(title="File containing filtered peaks")
      f=open(filename,"r")
      #                   0123456789012
      line=f.readline()
      if line[0:12] !="# xc yc omega"[0:12]:
         print line
         from tkMessageBox import showinfo
         showinfo("Sorry","That does not seem to be a filter peaks file, output from the peaksearching menu option")
         return
      bigarray=[]
      for line in f.readlines():
         v=[float(z) for z in line.split()]
         bigarray.append(v)
      f.close()
      self.parent.finalpeaks=transpose(array(bigarray))
      print self.parent.finalpeaks.shape

   def loadfileparameters(self,filename=None):
      if filename==None:
         filename=self.parent.opener.show(title="File containing detector parameters")
      self.loadparameters(filename)

   def loadparameters(self,filename=None):   
      if filename==None:
         self.parameters['z-center']=1002.832
         self.parameters['y-center']=1005.768
         self.parameters['distance']=303.7
         self.parameters['z-size']=46.77648
         self.parameters['y-size']=48.08150
         self.parameters['tilt-z']=0.0
         self.parameters['tilt-y']=0.0
         self.parameters['fit_tolerance']=0.05
         self.parameters['wavelength']=0.155
         self.parameters['cell__a']=2.9508
         self.parameters['cell__b']=2.9508
         self.parameters['cell__c']=4.6855
         self.parameters['cell_alpha']=90.0
         self.parameters['cell_beta']=90.0
         self.parameters['cell_gamma']=120.0
         self.parameters['cell_lattice_[P,A,B,C,I,F]']="P"
      else:
         try:
            lines = open(filename,"r").readlines()
            for line in lines:
               name,value=line.split()
               self.parameters[name]=value
         except IOError:
            from tkMessageBox import showinfo
            showinfo("Could not open your file")
               

   def saveparameters(self,filename=None):
      if filename==None:
         filename=self.parent.saver.show(title="File to save detector parameters")
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
      d=listdialog(self.parent,items=self.parameters,title="Detector parameters")
      self.parameters=d.result
      
   def plotyz(self):
      """
      Plots the x,y arrays being used
      """
      import twodplot
      self.parent.twodplotter.hideall() 
      self.parent.twodplotter.adddata( 
            ( "Filtered peaks",
               twodplot.data(
                  self.parent.finalpeaks[0,:],
                  self.parent.finalpeaks[1,:],
                  {"xlabel" : "y", "ylabel" : "z", "title" : "Peak positions on detector"} ) ) )
           

   def gof(self,args):
      self.tolerance=float(self.parameters['fit_tolerance'])
      self.parameters['y-center']=args[0]
      self.parameters['z-center']=args[1]
      self.parameters['distance']=args[2]
#      self.parameters['wavelength']=args[3]
      self.parameters['tilt-y']=args[3]
      self.parameters['tilt-z']=args[4]
      self.twotheta, self.eta = transform.compute_tth_eta( self.peaks_xy,
                                                           self.parameters['y-center'],
                                                           float(self.parameters['y-size']), # ys is the y pixel size
                                                           self.parameters['tilt-y'],
                                                           self.parameters['z-center'],
                                                           float(self.parameters['z-size'])  , # zs is the z pixel size
                                                           self.parameters['tilt-z'],
                                                           self.parameters['distance']*1e3
                                                           )
      w=float(self.parameters['wavelength']) # args[3] # what?
      gof=0.
      npeaks=0
      for i in range(len(self.tthc)):# (twotheta_rad_cell.shape[0]):
         self.tthc[i]=transform.degrees(math.asin(self.fitds[i]*w/2)*2)
         diff=take(self.twotheta,self.indices[i]) - self.tthc[i]
#         print "peak",i,"diff",maximum.reduce(diff),minimum.reduce(diff)
         gof=gof+dot(diff,diff)
         npeaks=npeaks+len(diff)
      gof=gof/npeaks
      return gof*1e3

   def fit(self):
      import simplex
      if self.theoryds==None:
         from tkMessageBox import showinfo
         showinfo("Try again","You need to have added the unitcell peaks already")
         return
      # Assign observed peaks to rings
      w=float(self.parameters['wavelength'])
      self.indices=[]
      self.tthc=[]
      self.fitds=[]
      self.tolerance=float(self.parameters['fit_tolerance'])
      tthmax = self.parent.twodplotter.a.get_xlim()[1]
      print "Tolerance for assigning peaks to rings",self.tolerance,"max tth",tthmax
      for i in range(len(self.theoryds)):
         dsc=self.theoryds[i]
         tthcalc=math.asin(dsc*w/2)*360./math.pi # degrees
         if tthcalc>tthmax:
            break
#         print tthcalc
         logicals= logical_and( greater(self.twotheta, tthcalc-self.tolerance),
                                   less(self.twotheta, tthcalc+self.tolerance)  ) 
         
         if sum(logicals)>0:
#            print maximum.reduce(compress(logicals,self.twotheta)),minimum.reduce(compress(logicals,self.twotheta))
            self.tthc.append(tthcalc)
            self.fitds.append(dsc)
            ind=compress(logicals,range(self.twotheta.shape[0]))
            self.indices.append(ind)
#            print "Ring",i,tthcalc,maximum.reduce(take(self.twotheta,ind)),minimum.reduce(take(self.twotheta,ind))
#      if raw_input("OK?")[0] not in ["Y","y"]:
#         return
      guess=[ float(self.parameters['y-center']) , 
              float(self.parameters['z-center']) ,
              float(self.parameters['distance']) , 
#              float(self.parameters['wavelength']) ,
              float(self.parameters['tilt-y']) ,
              float(self.parameters['tilt-z']) ]
      inc=[ .1 , .1 , .1 , transform.radians(0.1) , transform.radians(0.1) ]
      s=simplex.Simplex(self.gof,guess,inc)
      newguess,error,iter=s.minimize()
      self.parameters['y-center']=newguess[0]
      self.parameters['z-center']=newguess[1]
      self.parameters['distance']=newguess[2]
#      self.parameters['wavelength']=newguess[3]
      self.parameters['tilt-y']=newguess[3]
      self.parameters['tilt-z']=newguess[4]
      inc=[ .01 , .01 , .01 , 0.0001 , transform.radians(0.01) , transform.radians(0.01) ]
      guess=newguess
      s=simplex.Simplex(self.gof,guess,inc)
      newguess,error,iter=s.minimize()
      self.parameters['y-center']=newguess[0]
      self.parameters['z-center']=newguess[1]
      self.parameters['distance']=newguess[2]
#      self.parameters['wavelength']=newguess[3]
      self.parameters['tilt-y']=newguess[3]
      self.parameters['tilt-z']=newguess[4]      
      print newguess

   def plotreta(self):
      self.peaks_xy = self.parent.finalpeaks[0:2,:]
      self.x=self.peaks_xy[0,:]
      self.y=self.peaks_xy[1,:]
      self.omega=self.parent.finalpeaks[2,:]
      print self.peaks_xy.shape
      try:
         self.twotheta, self.eta = transform.compute_tth_eta( self.peaks_xy,
               float(self.parameters['y-center']), # yc is the centre in y
               float(self.parameters['y-size'])  , # ys is the y pixel size
               float(self.parameters['tilt-y'])  , # ty is the tilt around y
               float(self.parameters['z-center']), # zc is the centre in z
               float(self.parameters['z-size'])  , # zs is the z pixel size
               float(self.parameters['tilt-z'])  , # tz is the tilt around z
               float(self.parameters['distance'])*1e3) # is the sample - detector distance
         self.ds = 2*sin(transform.radians(self.twotheta)/2)/float(self.parameters['wavelength'])
      except:
         keys=self.parameters.keys()
         keys.sort()
         for key in keys:
            print self.parameters[key],type(self.parameters[key])
         raise
      import twodplot
#      self.parent.twodplotter.hideall()
      self.parent.twodplotter.adddata(
            ( "2Theta/Eta",
               twodplot.data(
                  self.twotheta,
                  self.eta,
                  {"xlabel":"TwoTheta / degrees",
                   "ylabel":"Azimuth / degrees",
                   "title" :"Peak positions"}
                   )))


   def addcellpeaks(self):
      #
      # Given unit cell, wavelength and distance, compute the radial positions
      # in microns of the unit cell peaks
      #
      a      =float(self.parameters['cell__a'])
      b      =float(self.parameters['cell__b'])
      c      =float(self.parameters['cell__c'])
      alpha  =float(self.parameters['cell_alpha'])
      beta   =float(self.parameters['cell_beta'])
      gamma  =float(self.parameters['cell_gamma'])
      lattice=self.parameters['cell_lattice_[P,A,B,C,I,F]']
      self.unitcell=unitcell.unitcell( [a,b,c,alpha,beta,gamma] ,  lattice)
      self.parent.unitcell=self.unitcell
      if self.twotheta==None:
         self.twotheta, self.eta = transform.compute_tth_eta( self.peaks_xy,
               float(self.parameters['y-center']), # yc is the centre in y
               float(self.parameters['y-size'])  , # ys is the y pixel size
               float(self.parameters['tilt-y'])  , # ty is the tilt around y
               float(self.parameters['z-center']), # zc is the centre in z
               float(self.parameters['z-size'])  , # zs is the z pixel size
               float(self.parameters['tilt-z'])  , # tz is the tilt around z
               float(self.parameters['distance']*1e3)) # is the sample - detector distance
      # Find last peak in radius
      highest = maximum.reduce(self.twotheta)
      wavelength=float(self.parameters['wavelength'])
      ds = 2*sin(transform.radians(highest)/2.)/wavelength
      self.dslimit=ds
      #print "highest peak",highest,"corresponding d*",ds
      self.theorypeaks=self.unitcell.gethkls(ds)
      self.unitcell.makerings(ds)
      self.theoryds=self.unitcell.ringds
      tths = [arcsin(wavelength*dstar/2)*2 for dstar in self.unitcell.ringds]
      self.theorytth=transform.degrees(array(tths))
      import twodplot
      self.parent.twodplotter.adddata(
            ( "HKL peaks",
               twodplot.data(
                       self.theorytth,
                       zeros(self.theorytth.shape[0]),
                       {'color':'r',
                        'pointtype':'+'}
               )))
      
   def setaxisorientation(self):
      """
      Allow the rotation axis to not be perpendicular to the beam
      """
      p = { "wedge" : self.wedge }
      d=listdialog(self.parent,items=p,title="Axis Orientation")
      self.wedge = float(d.result["wedge"])
      print "set self.wedge to",self.wedge


            
   def computegv(self):
      """
      Using self.twotheta, self.eta and omega angles, compute x,y,z of spot
      in reciprocal space
      """
      self.gv = transform.compute_g_vectors(self.twotheta,self.eta,self.omega,float(self.parameters['wavelength']),wedge=self.wedge)
      tthnew,etanew,omeganew=transform.uncompute_g_vectors(self.gv,float(self.parameters['wavelength']),wedge=self.wedge)
      self.parent.gv=self.gv
      print "Testing reverse transformations"
      for i in range(5):
         print "%8.3f %8.3f %8.3f "%(self.twotheta[i],self.eta[i],self.omega[i]),
         print "%8.3f %8.3f %8.3f %8.3f %8.3f"%(tthnew[i],etanew[0][i],omeganew[0][i],etanew[1][i],omeganew[1][i])
      return

   def savegv(self,filename=None):
      """
      Save g-vectors into a file
      Use crappy .ass format from previous for now (testing)
      """
      if filename==None:
         filename=self.parent.saver.show(title="File to save gvectors")
      f=open(filename,"w")
      f.write(self.unitcell.tostring())
      f.write("\n")
      f.write("# wavelength = %f\n"%( float(self.parameters['wavelength']) ) )
      f.write("# wedge = %f\n"%( float(self.wedge) ))
      f.write("# ds h k l\n")
      for peak in self.theorypeaks:
           f.write("%10.7f %4d %4d %4d\n"%(peak[0],peak[1][0],peak[1][1],peak[1][2]))
      order = argsort(self.twotheta)
      f.write("# xr yr zr xc yc ds phi omega\n")
      print maximum.reduce(self.omega),minimum.reduce(self.omega)
      for i in order:
           f.write("%f %f %f %f %f %f %f %f \n"%(self.gv[0,i],self.gv[1,i],self.gv[2,i],self.x[i],self.y[i],self.ds[i],self.eta[i],self.omega[i]))
      f.close()


   def write_graindex_gv(self):
      filename=self.parent.saver.show(title="File for graindex, try finalpeaks.log")
      from ImageD11 import write_graindex_gv
      #self.parent.finalpeaks=transpose(array(bigarray))
      self.intensity = self.parent.finalpeaks[3,:]*self.parent.finalpeaks[4,:]
      print self.intensity.shape
      write_graindex_gv.write_graindex_gv(filename,
                                          self.gv,
                                          self.twotheta,
                                          self.eta,
                                          self.omega,
                                          self.intensity,
                                          self.unitcell)
      
