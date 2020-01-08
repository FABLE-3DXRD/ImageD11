============
Program List
============

This is a list of the different things which can be run from a command line.

Image Related
=============
peaksearch.py
  Peaksearches the images, see the detailed help.

recoveromega.py
  Attempts to repair the problem of having no header information about sample orientation. Largely superceeded by command line options to peaksearch.py nowadays.

fit2dcake.py
  Drives the fit2d program to do radial integration. See the detailed help and also the official fit2d site.

edfheader.py
  Reads the header information from ID11/ESRF edf files.

id11_summarize.py
  Runs through a set of images trying to figure out what happened during an experiment.

bgmaker.py
  Runs through a series of images finding the minimum.

edf2bruker.py
  Converts ID11/ESRF edf files into bruker format for processing with Bruker tools (eg smart, saint etc).

powderimagetopeaks.py
  Draws a grid onto a powder image to cut the rings up into spots. The resulting image can then be peaksearched and used for calibration.

Graphical applications 
=======================

ImageD11_gui.py
  Your one stop solution for peak transformation, calibration and indexing.

plotedf.py
  OpenGl based image display. Should be obsoleted. See also imageviewer and fabian.

rubber.py
  Tk based Image display. Should be obsoleted. See also imageviewer and fabian.

Grain refinement
================

fitgrain.py
  Fits a grain (eg ubi) and diffractometer parameters.

filtergrain.py
  Selects peaks which belong to a grain for fitgrain.

filterout.py
  Tries to remove peaks belonging to a grain from a filtered peaks file so you can then index the rest.

makemap.py
  Refines the position and orientation of grains with diffractometer pars fixed to make a 3D centre of mass map.

plotgrainhist.py
  Shows how well a set of grains fit a set of peaks in an flt file, via a list of parameters.

ubi2cellpars.py
  Computes the unit cell parameters for the UBI matrices in a ubi file

pars_2_sweeper.py
  Converts the ImageD11 grain refinement parameters into something closert to format which is needed by  Soren Schmidt's grainsweeper program.

Miscellaneaous
==============

index_unknown.py
  A new indexing program aimed at small numbers of unknown single crystals.

ImageD11_2_shelx.py
  Rather untested, unlikely to work yet. Gives integrated intensities for peaks in shelx format. You are better off with fabric.

ImageD11Server.py
  Obsolete - an xmlrpc server for the Java Gui's which now use jepp instead.

The Java/Eclipse/RCP based graphical interfaces to ImageD11:

These optional packages offer a richer user experience. Currently 
all features which are contained in ImageD11 are thought to be
accessible without using these graphical interfaces. The hope is 
that these interfaces make the program easier to use for a beginner,
and still offer the full range of features and bugs for the experts.
If you are running jobs via a condor batch processing system, you
probably still want the command line versions. Beware that the was
an early version of the Java based gui which was obsoleted after 
a change from xmlrpc based communication to using JNI hosted
python interpreter.
