Post-processing
===============

The 3D merge is now done during the peaksearch. This section is 
obsolete 
After peaksearching all of your diffraction images you will have 
reduced your data 
analysis problem from a multi-gigabyte scale to something which is 
hopefully less only a 
few megabytes. The output from the peaksearch.py script should 
contain all the 
information you could possibly need, but currently ImageD11 only 
provides tools to use a 
small part of that information. Specifically it processes peak 
positions in terms of 
xc,yc and Omega (with the possible addition of wedge, manually, 
later).

Currently this processing can only be carried out in the gui, but 
it is a very high priority to make this possible via a command line 
script (FIXME!!!). Start up the graphical interface via the script 
"ImageD11_gui.py", either by typing that command at a unix terminal 
or double clicking it under windows (windows users are likely to 
find it installed somewhere like 
c:\python24\scripts\ImageD11_gui.py). 
Eventually a message giving a list of things to do should appear on 
the screen. You are warmly invited to work on any of those tasks! 
Click the "OK" button to get rid of it and the graphical interface 
should appear after a few seconds. If it is still not there after a 
few minutes then either you have a very slow computer or you have 
some installation or configuration problems. Try to start the 
program again from a command line, READ the error message and act 
accordingly. 
This page describes the use of the "PeakSearching" menu. The first 
item (search raw images) does not do anything, see 
ImageD11:peaksearching for carrying out a peaksearch, or 
alternatively implement a graphical interface for carrying out that 
task. The second item (read .pks file) is the item you want to 
begin with. It will read the output of a peaksearch run and then 
display a plot of image number versus omega. This plot allows you 
to interpret whether or not a single continuous rotation was 
carried out and to select (via a mouse zoom) the range of images 
you wish to use. If your plot is garbage or you get an error 
message then probably you do not have good Omega values in your 
peaksearch output. Use the script recoveromega.py to fix that. 
Having selected the range of images you then need to ask the 
program to Harvest peaks, which means to go and read in all of the 
peaks from the output of the peaksearch for the range of images you 
want to use. It will tell you how many peaks it found. 
Once those peaks are read from the file they probably should be 
merged via the merge peaks menu option. This command tells the 
program to compute weighted averages of peaks on adjacent frames 
and is a poor mans version of a 3 dimensional peaksearch. Whether 
or not peaks are merged currently depends on the lines in file 
guipeaksearch.py::

  class peak:
     def __init__(self,line,omega,threshold,num,tolerance=2):
        self.TOLERANCE = tolerance # Pixel separation for combining peaks
        ...
FIXME Clearly this should be configurable by the average user from the 
interface and the whole processing business should not contain any 
reference to Tk! 
The next menu option is filter peaks which currently does not really do 
anything useful except to plot the peak positions on the screen. This will 
give you an idea how many spots you have observed, and perhaps indicate 
potential problems (eg: peak positions due to the direct beam). 
You are now ready to save good peaks into an output file. This output file 
contains all the information from the peaksearch merged together for peaks 
on adjacent images.::

  # xc yc omega npixels avg_intensity x_raw y_raw sigx sigy covxy
  1105.310387 1934.819179 -38.750000 102.000000 2022.323529 1119.384182 1953.811254 1.815357 3.546957 -0.082630
  372.770569 1469.674753 -38.750000 21.000000 1326.047619 377.557906 1468.206162 1 .761871 1.075446 -0.469539

  You can filter this file as you like. For example, using the awk program 
(standard for most unix systems, also available for windows [find a link]). 
For example, to remove peaks containing less the 10 pixels, try to following:

  awk '($4>10){print}' < originalpeaks.flt > filteredpeaks.flt
(You might have to adjust the first line of the file to be the same as 
before) 
The filtered peaks are now ready to me read into the transform menu, see 
ImageD11:Calibration.