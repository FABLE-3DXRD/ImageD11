from ImageD11.guimaker import GuiMaker
from Tkinter import *


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




if __name__=="__main__":
   def help():
      from tkMessageBox import showinfo
      showinfo("Help","Sorry, no help for you!\nPlease try harder")
   ImageD11credits = """Thanks to:
                       
                       All of the Fable team, which includes at least:
                          Andy Goetz, Gavin Vaughan,
                          Soren Schmidt, Henning Poulsen, Larry Margulies
                       ...and others who should remind me to mention them

                       John Hunter for the matplotlib plotting
                       
                       Anyone who tests me and gives useful feedback
                          
                       Jon Wright, for writing me!
                      """
   from ImageD11 import license
   license = license.license


   def credits():
      from tkMessageBox import showinfo
      showinfo("Credits",ImageD11credits)

   def showlicense():
      from ScrolledText import ScrolledText
      win = ScrolledText(Toplevel(),width=100)
      win.insert(END,license)
#      print dir(win)
      win.pack(expand=1,fill=BOTH)
      win.focus_set()
#      showinfo("License",license)

   import tkFileDialog,os

   class TestGuiMaker(GuiMaker):
      def start(self):
          from tkMessageBox import showinfo


          startmessage = """
ImageD11 version 0.4, Copyright (C) 2005 Jon Wright
ImageD11 comes with ABSOLUTELY NO WARRANTY; for details select help, license.
This is free software, and you are welcome to redistribute it under certain conditions

Please send useful feedback to wright@esrf.fr
"""

          startmessage += """
Stuff to do:

   Implement the image peaksearching in the gui (maybe display opengl images??)
   Separate the peak merging/reading from the gui for standalone scripting
   Same for the transformations - once parameters are known gui is not needed
   Tidy the mulitple threshold consequences
   Implement those filters based in intensity etc

   Sort out a file format which tracks all information properly?
"""
          showinfo("Welcome to ImageD11_v0.4",startmessage)
#             from tkMessageBox import showinfo
#             showinfo("Sorry","You had to say yes then if you want to use the program")
#             sys.exit()

          from ImageD11 import guipeaksearch 
          self.peaksearcher = guipeaksearch.guipeaksearcher(self)
          self.finalpeaks=None
          from ImageD11 import guitransformer
          self.unitcell=None
          self.gv=None
          self.transformer  = guitransformer.guitransformer(self)
          from ImageD11 import guiindexer
          self.indexer = guiindexer.guiindexer(self)
          import sys
          self.menuBar = [ ( "File", 0,
                              [ ( "Print" , 0, self.printplot ),
                                ( "Exit", 1, sys.exit) ] ) ,
                           self.peaksearcher.menuitems,
                           self.transformer.menuitems,
                           self.indexer.menuitems,
                           ( "Plotting", 0, 
                              [ ( "Autoscale", 0, self.autoscaleplot),
                                ( "Clear plot",0, self.clearplot),
                              ] ) ,
                           
                           ( "Help", 0, 
                              [ ( "Help Me!", 0, help) ,
                                ( "Credits" , 0, credits) ,
                                ( "License" , 0, showlicense)
                                ] ) ]

      def printplot(self): self.twodplotter.printplot()
  
      def autoscaleplot(self):
         self.twodplotter.autoscale()
      def clearplot(self):
         self.twodplotter.plotitems={}
         self.twodplotter.replot()

      def makeWidgets(self):
         # Get size of TopLevels window
         self.opener=tkFileDialog.Open()
         self.saver=tkFileDialog.SaveAs()
         from ImageD11 import twodplot
         self.twodplotter=twodplot.twodplot(self)
         self.twodplotter.pack(side=RIGHT, expand=1, fill=BOTH)
        
   root = Tk()
   root.wm_title("ImageD11")
   TestGuiMaker()
   root.mainloop()
