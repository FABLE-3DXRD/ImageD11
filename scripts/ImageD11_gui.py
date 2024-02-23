#!/usr/bin/env python
from __future__ import print_function


# ImageD11_v1.0 Software for beamline ID11
# Copyright (C) 2005-2007  Jon Wright
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

try:
    from Tkinter import *
    import tkFileDialog as filedialog
    from tkMessageBox import showinfo
    from ScrolledText import ScrolledText
except:
    from tkinter import *
    import tkinter.filedialog as filedialog
    from tkinter.messagebox import showinfo
    from tkinter.scrolledtext import ScrolledText

import logging
import sys
import os

# GuiMaker is for building up the windows etc

from ImageD11.tkGui.guimaker import GuiMaker
from ImageD11.tkGui import twodplot,  guipeaksearch, guitransformer, guiindexer
# guisolver
from ImageD11 import __version__, guicommand
from ImageD11.license import license

# This does not do anything unless you call it as a program:

if __name__ == "__main__":
    # get the output!
    # Set up the logging stuff
    console = logging.StreamHandler(sys.stdout)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(levelname)-8s : %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    console.setLevel(logging.DEBUG)
    root = logging.getLogger('')
    root.addHandler(console)
    root.setLevel(logging.DEBUG)  # sh

    # Help message - TODO - proper online help
    def help():
        hlp = """Try pointing your web browser at:
   http://fable.sourceforge.net/index.php/ImageD11
   See also the html documentation for programmers, somewhere like:
   file:///c:/python24/Lib/site-packages/ImageD11/doc/ImageD11.html
"""
        showinfo("Help", "Please try harder\n"+hlp)

    # Credits message
    ImageD11credits = """Thanks to:

                        All of the Fable team, which includes at least:
                           Andy Gotz, Gavin Vaughan, Henning O. Sorensen,
                           Soren Schmidt, Henning Poulsen, Larry Margulies
                           Erik Knudsen,
                        ...and others who should remind me to mention them

                        Tine Knudsen who bravely helped commission the
                        introduction of a second rotation axis (wedge).

                        Benoit Mallard for his assistance with some extreme
                        programming to debug the transformation module.

                        Younes ElHachi for adding the eps_sig calculations.

                        John Hunter for the matplotlib plotting.

                        All of the pyopengl, Numeric, numpy and python teams

                        Anyone who tests me and gives useful feedback

                        Jon Wright, for writing me!
                       """

    def credits():
        showinfo("Credits", ImageD11credits)

    # GPL is stored in ImageD11/license.py as a string to be
    # displayed in the GUI if the user asks to see it


    def showlicense():
        win = ScrolledText(Toplevel(), width=100)
        win.insert(END, license)
        win.pack(expand=1, fill=BOTH)
        win.focus_set()

    # Inherits from the GuiMaker and uses functions defined above
    class TestGuiMaker(GuiMaker):
        guicommand.RETURN_NUMERICS = True
        guicommander = guicommand.guicommand()

        def start(self):
            """
            Override the GuiMaker start
            These are things to do when the gui starts
            eg: show a message about the license and list of things to do
            Then build the actual gui
            """
            startmessage = """
  ImageD11 version %s, Copyright (C) 2005-2017 Jon Wright
  ImageD11 comes with ABSOLUTELY NO WARRANTY; for details select help,
  license. This is free software, and you are welcome to redistribute it
  under certain conditions

  Please send useful feedback to wright@esrf.fr
  """ % (__version__)

            startmessage += """
  You are using version %s

  There have been lots of changes recently!

  I would also be happily surprised if it is currently working.
  """ % (__version__)
            showinfo("Welcome to ImageD11 " + __version__,
                     startmessage)

            # For the peaksearch menu

            self.peaksearcher = guipeaksearch.guipeaksearcher(self)


            self.transformer = guitransformer.guitransformer(self)

            # For the indexing - supposed to generate orientations from the
            # unitcell and g-vectors
            self.indexer = guiindexer.guiindexer(self)

#            self.solver = guisolver.guisolver(self)

            # Configure the menubar (lists of Tuples of (name,
            # underline_char, command_or_submenu) )
            self.menuBar = [("File", 0,
                             [("Print", 0, self.printplot),
                              ("Exit", 1, sys.exit)]),
                            self.peaksearcher.menuitems,
                            self.transformer.menuitems,
                            self.indexer.menuitems,
#                            self.solver.menuitems,
                            ("Plotting", 0,
                             [("Autoscale", 0, self.autoscaleplot),
                              ("Clear plot", 0, self.clearplot),
                              ]),
                            ("Help", 0,
                             [("Help Me!", 0, help),
                              ("History", 1, self.history),
                              ("Credits", 0, credits),
                              ("License", 0, showlicense)
                              ])]

        # The twodplot object should be taking care of it's own menu
        # Stop doing it here - TODO
        def history(self):
            win = ScrolledText(Toplevel(), width=100)
            history = self.guicommander.gethistory()
            win.insert(END, history)
            win.pack(expand=1, fill=BOTH)
            win.focus_set()

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
            self.twodplotter.clear()

        def makeWidgets(self):
            """
            Draw the gui and initialise some hidden things
            """
            # TODO Get size of TopLevels window and position it in a
            # sensible way
            #
            # Opening and saving file widgets, normally hidden, they
            # remember where you were for working directories
            self.opener = filedialog.Open()
            self.saver = filedialog.SaveAs()
            #
            # Draw the twodplot
            self.twodplotter = twodplot.twodplot(self)
            self.twodplotter.pack(side=RIGHT, expand=1, fill=BOTH)

    # Start up Tkinter
    root = Tk()
    root.wm_title("ImageD11")
    # Instantiate an object of the class TestGuiMaker
    TestGuiMaker()
    # Thats it!
    root.mainloop()
