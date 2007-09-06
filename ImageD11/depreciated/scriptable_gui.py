#!/bliss/users/blissadm/python/bliss_python/suse82/bin/python

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
Graphical user interface for ImageD11

Uses Tkinter (comes with python, normally)
Also depends eventually on matplotlib (publication quality 2d plotting)
and on OpenGl for plotting things in 3D (reciprocal space)

Hopefully the gui and the gruntwork can be completely separated, eventually.
"""

raise Exception("""depreciated script 
Use the normal ImageD11_gui.py script and click help -> history to get a script""")



from ImageD11.guimaker import GuiMaker
# GuiMaker is for building up the windows etc

from Tkinter import *




# This does not do anything unless you call it as a program:

if 1:


    # Help message - TODO - proper online help
    def help():
        from tkMessageBox import showinfo
        showinfo("Help","Sorry, no help for you!\nPlease try harder")


    # Credits message
    ImageD11credits = """Thanks to:

                        All of the Fable team, which includes at least:
                           Andy Goetz, Gavin Vaughan,
                           Soren Schmidt, Henning Poulsen, Larry Margulies
                        ...and others who should remind me to mention them

                        John Hunter for the matplotlib plotting

                        Anyone who tests me and gives useful feedback

                        Jon Wright, for writing me!
                       """

    def credits():
        from tkMessageBox import showinfo
        showinfo("Credits",ImageD11credits)


    # GPL is stored in ImageD11/license.py as a string to be
    # displayed in the GUI if the user asks to see it

    from ImageD11 import license
    license = license.license

    def showlicense():
        from ScrolledText import ScrolledText
        win = ScrolledText(Toplevel(),width=100)
        win.insert(END,license)
        win.pack(expand=1,fill=BOTH)
        win.focus_set()


    import tkFileDialog,os


    # Inherits from the GuiMaker and uses functions defined above
    class TestGuiMaker(GuiMaker):
        def start(self):
            """
            Override the GuiMaker start
            These are things to do when the gui starts
            eg: show a message about the license and list of things to do
            Then build the actual gui
            """
            from tkMessageBox import showinfo
            import ImageD11
            startmessage = """
  ImageD11 version %s, Copyright (C) 2005 Jon Wright
  ImageD11 comes with ABSOLUTELY NO WARRANTY; for details select help, license.
  This is free software, and you are welcome to redistribute it under certain conditions

  Please send useful feedback to wright@esrf.fr
  """%(ImageD11.__version__)

            startmessage += """
  Stuff to do:

     Implement the image peaksearching in the gui (maybe display opengl images??)
     Separate the peak merging/reading from the gui for standalone scripting
     Same for the transformations - once parameters are known gui is not needed
     Tidy the mulitple threshold consequences
     Implement those filters based in intensity etc

     Sort peaks in output file by integer hkl
     Allow merging algorith to be more flexible in tolerance, perhaps decide
     overlap as a multiple of FWHM observed.

     Sort out a file format which tracks all information properly?
  """
#          showinfo("Welcome to ImageD11_v0.4",startmessage)


            # For the peaksearch menu
            from ImageD11 import guipeaksearch
            self.peaksearcher = guipeaksearch.guipeaksearcher(self,quiet="yes")

            # self.finalpeaks is what the peaksearchmenu is meant to generate
            # and what the transformer should transform
            self.finalpeaks=None

            # For the transformation menu
            from ImageD11 import guitransformer

            # Unit cell is for generating theoretical peak positions
            self.unitcell=None

            # self.gv holds the g-vectors, after they are generated
            self.gv=None
            self.transformer  = guitransformer.guitransformer(self,quiet="yes")

            # For the indexing - supposed to generate orientations from the
            # unitcell and g-vectors
            from ImageD11 import guiindexer
            self.indexer = guiindexer.guiindexer(self)

            # sys is for sys.exit
            import sys

            # Configure the menubar (lists of Tuples of (name, underline_char, command_or_submenu) )
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

        # The twodplot object should be taking care of it's own menu
        # Stop doing it here - TODO
        def printplot(self):
            """
            Print the 2D plot (probably to a file?)
            """
            self.twodplotter.printplot()

        def autoscaleplot(self):
            """
            Autoscale the plot
            """
            self.twodplotter.autoscale()

        def clearplot(self):
            """
            Clear out the twodplot
            """
            self.twodplotter.plotitems={}
            self.twodplotter.replot()

        def makeWidgets(self):
            """
            Draw the gui and initialise some hidden things
            """
            # TODO Get size of TopLevels window and position it in a sensible way
            #
            # Opening and saving file widgets, normally hidden, they remember where
            # you were for working directories
            self.opener=tkFileDialog.Open()
            self.saver=tkFileDialog.SaveAs()
            #
            # Draw the twodplot
            from ImageD11 import twodplot
            self.twodplotter=twodplot.twodplot(self,quiet="yes")
            self.twodplotter.pack(side=RIGHT, expand=1, fill=BOTH)


if __name__=="__main__":
    # Start up Tkinter
    root = Tk()
    root.wm_title("ImageD11")
    # Instantiate an object of the class TestGuiMaker
    o=TestGuiMaker()
    # Thats it!
    import sys
    o.peaksearcher.readpeaks(sys.argv[1])
    o.peaksearcher.harvestpeaks()
    o.peaksearcher.mergepeaks()
    o.peaksearcher.filter()
    o.peaksearcher.savepeaks(sys.argv[1][:-3]+"pks")
    #o.transformer.parameterobj.loadparameters(sys.argv[2])
    o.clearplot()
    o.transformer.plotreta()
    o.transformer.addcellpeaks()
    o.transformer.computegv()
    o.transformer.savegv(sys.argv[1][:-3]+"gve")
#   root.mainloop()
    sys.exit()
