

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

try:
    from tkMessageBox import showinfo, askyesno
except:
    # python3?
    from tkinter.messagebox import showinfo, askyesno

from . import twodplot


class guipeaksearcher:
    """
    Tkinter wrapper for ImageD11 guipeakmerge object

    Note that peaksearching is done via the command line

    All communication should be via parent guicommander object
    """
    def __init__(self,parent,quiet="No"):
        """
        Parent is a hook to features of the parent gui
        sets up peakmerging menu items + message for searching
        """
        self.parent=parent
        self.quiet=quiet
        #print "I am in quiet mode",quiet
        self.menuitems = ( "PeakSearching", 4,
                           [ ( "Search raw images" , 0, self.searchraw ),
                             ( "Read pks file", 0, self.readpeaks ),
                             ( "Harvest peaks", 0, self.harvestpeaks),
                             ( "Merge peaks",  0, self.mergepeaks ),
                             ( "Filter peaks", 0, self.filter     ),
                             ( "Save good peaks",2, self.savepeaks) ] )

    def searchraw(self):
        """ Explains to user about the command line script """
        showinfo("Sorry","""Not implemented for gui, please use the jonpeaksearch script on amber/crunch for now""")

    def readpeaks(self,filename=None):
        """
        Runs peakmerger.readpeaks
        gets names for first and last image for plot title
        gets omega angles and image numbers for plot
        plots image number versus omega angle to see if you did a scan
        """
        if filename is None:
            filename = self.parent.opener.show(title="Peaksearch results file")
        import os
        if os.path.isfile(filename):
            # apply this to peakmerge
            self.parent.guicommander.execute("peakmerger","readpeaks",filename)
        else:
            showinfo("Sorry, bad filename %s"%(filename))
        imagenumbers = self.parent.guicommander.getdata("peakmerger","imagenumbers")
        omegas = self.parent.guicommander.getdata("peakmerger","omegas")
        images = self.parent.guicommander.getdata("peakmerger","images")
        first = images[0].name
        last = images[-1].name
        self.parent.twodplotter.adddata(
              ("number versus omega",
                 twodplot.data(
                  imagenumbers,omegas,
                  {"xlabel" : "imagenumber",
                   "ylabel" : "Omega",
                   "title"  : first+"..."+last }
                  ) # data
                 )# tuple to plotter
              ) # call

    def harvestpeaks(self):
        """
        gets range from 2d plot (image numbers and omegas)
        calls peakmerger.harvestpeaks(image_number_range,omega_range)
        """
        # Now we need to select the range of peaks to use
        if self.quiet=="No":
            ans = askyesno("Have you selected a sensible range of images?","""
   Use the mouse to select the range of image numbers and omega angles that
   you want to use from the graph on the screen.

   If all is OK, then say yes now and we will try to harvest the peaks.
   Otherwise, say no, select the right range and come back "harvestpeaks" again
   """       )
            # FIXME? ans is ignored anyway
        # We now have the ranges in imagenumber and omega from
        # the graph
        numlim=self.parent.twodplotter.a.get_xlim()
        omlim=self.parent.twodplotter.a.get_ylim()
        self.parent.guicommander.execute("peakmerger","harvestpeaks",numlim=numlim,omlim=omlim)
        if self.quiet=="No":
            npks = len(self.parent.guicommander.getdata("peakmerger","allpeaks") )
            showinfo("Harvested peaks","You have a total of %d peaks,"%(npks)+
               "no peaks have been merged")

    def mergepeaks(self):
        """ calls peakmerger.mergepeaks and reports number of peaks to user """
        self.parent.guicommander.execute("peakmerger","mergepeaks")
        nmerged = len(self.parent.guicommander.getdata("peakmerger","merged"))
        if self.quiet=="No":
            showinfo("Finished merging peaks","You have a total of "+str(nmerged)+" after merging")
        return

    def filter(self):
        """ calls peakmerger.filter (does very little)
        plots x and y final peak positions

        TODO implement filters!!!
        """
        self.parent.guicommander.execute("peakmerger","filter")
        peaks = self.parent.guicommander.getdata("peakmerger","finalpeaks")
        if self.quiet=="No":
            self.parent.twodplotter.hideall() # get rid of image number versus omega plot
            self.parent.twodplotter.adddata(
                  ( "Filtered peaks",
                     twodplot.data(
                        peaks[0,:],
                        peaks[1,:],
                        {"xlabel" : "x",
                         "ylabel" : "y",
                         "title"  : "Peak positions on detector"}
                         ) ) )
        return "nothing"
        # Need to filter based on x,y
        # also based on intensity
        # also based on shape

    def savepeaks(self,filename=None):
        """ see peakmerger.savepeaks """
        if filename==None:
            filename=self.parent.saver.show(title="Filtered peak positions")
        self.parent.guicommander.execute("peakmerger","savepeaks",filename)
