=============
Peaksearching
=============

Peaksearching in two dimensional images can be carried out in a relatively automated way using the peaksearch.py script. There is a command line interface and the underlying code can also be used via the fable gui.

Giving the command "peaksearch.ph -h" should print the following help message::

 options:
 -h, --help            show this help message and exit
 -n STEM, --namestem=STEM
                       Name of the files up the digits part, eg mydata in
                       mydata0000.edf
 -F FORMAT, --format=FORMAT
                       Image File format, eg edf or bruker
 -S OMEGASTEP, --step=OMEGASTEP
                       Step size in Omega when you have no header info
 -T OMEGA, --start=OMEGA
                       Start position in Omega when you have no header info
 -f FIRST, --first=FIRST
                       Number of first file to process, default=0
 -l LAST, --last=LAST  Number of last file to process
 -o OUTFILE, --outfile=OUTFILE
                       Output filename, gets extension .spt
 -d DARK, --darkfile=DARK
                       Dark current filename, to be subtracted, default=None
 -D DARKOFFSET, --darkfileoffset=DARKOFFSET
                       Constant to subtract from dark to avoid overflows,
                       default=100
 -s SPLINE, --splinefile=SPLINE
                       Spline file for spatial distortion,
                       default=/data/opid11/inhouse/Frelon2K/spatial2k.spline
 -p PERFECT, --perfect_images=PERFECT
                       Ignore spline Y|N, default=N
 -O FLOOD, --flood=FLOOD
                       Flood file, default=None
 -t THRESHOLDS, --threshold=THRESHOLDS
                       Threshold level, you can have several

					   
Usually data images are collected as a numerical series while the sample rotates and are named using the convention stem0000.edf, stem0001.edf etc. The first three options "-n", "-f" and "-l" specify the series of images in terms of stem name start and end number. Currently the script will automatically assume ".edf" or ".edf.gz" extensions and corresponding format. If you want for example ".tif.gz" say "-F .tif.gz". 
It should accept the following formats::

 ".edf",".edf.gz",".edf.bz2",
 ".tif",".tif.gz",".tif.bz2",
 ".mccd",".mccd.gz",".mccd.bz2"
 
Others are present in spirit, but not perhaps via a user interface. The "-d" option is optional, for a dark current image that can be subtracted. If you use the Frelon2K CCD camera at ESRF the please allow for the offset of 1000 ADU if you do not decide to subtract a dark current. The "-D" option alters the offset to try to avoid overflow when subtraction gives a negative number (it is though to be marginally faster to remain with unsigned short data than going over to int).
The "-s" option specifies a fit2d spline file characterising the spatial distortion of the camera. The default is probably OK for recent data from ID11, but more generally you will want to change this option. Be sure to use the pixel size data from this file when computing the detector calibration parameters. In the case where no spatial corrections are needed then please be sure to specify the "-p Y" option to indicate there is no need for spatial distortion correction (eg a nice image plate). Corrections are applied to peaks positions, NOT to the raw images. 
Threshold levels are used to determine which pixels contribute to diffraction spots. A segmentation of the image will be made where all pixels above the threshold are the collected into connected objects. Try to choose levels which are significantly above the background level. You may use more than one level if you like, eg "-t 1050 -t 5000" for levels 1050 ADU and 5000 ADU. The output file contains all of the information contained in the header as well as the obtained peak positions. You can specify a more logical name for this file if you do not like "peaks.out". 
Some examples::

  peaksearch.py -n plat_a -f 50 -l 90 -o peaks_50_90.out -s spatial.spline -t 5000

This will peaksearch images::

  plat_a0050.edf
  plat_a0051.edf
  ...
  plat_a0090.edf
  
A spatial distortion computed from the file spatial.spline (in the current directory) will be applied. A threshold level of 5000 is used to decide which pixels contribute to peaks. 
Example output::

  # File plat_a0050.edf
  # Processed on Wed Nov 30 17:09:34 2005
  # Spatial correction from /data/opid11/inhouse/Frelon2K/spatial2k.spline
  # SPLINE X-PIXEL-SIZE 46.776480
  # SPLINE Y-PIXEL-SIZE 48.081500
  # HeaderID = EH:000001:000000:000000
  # y1 = -3.88889e-05
  # Instron strain = 0.002
  [... lots of header information snipped ...]
  # Omega = -39

  #Threshold level 1200.000000
  # Number_of_pixels Average_counts    x   y     xc   yc      sig_x sig_y cov_xy
  119  2371.050420    377.755365 1467.049632    372.961306 1468.529932    3.122391 2.406393 -0.363226
  44  1600.363636    424.880254 1857.962651    420.424252 1852.106077    1.570551 2.092762 -0.350194
  41  1495.170732    1118.681952 1952.197416    1104.657513 1933.297804    1.19973 6 2.575913 -0.063616

  # File plat_a0050.edf
  [ ...etc... ]
  
  
The numbers in the output file have the following meanings:

**Number_of_pixels** is the number of pixels in the connected object.

**Average_counts** is the average counts in each pixel in that object. If the camera has a constant offset (in this case about 1000) then you need to subtract that number before computing intensities. In this case intensity is approximately number_of_pixels*(Average_counts-1000.0)

**x** and **y** are the position of the centre of mass of the peak in the uncorrected image (first moments of the intensity distributions in the blobs is the mean). **xc** and **yc** are the x and y positions corrected for spatial distortion.

**sig_x** and **sig_y** are the second moments of the intensity distribution in the blob. Something like the width, but off by some factor of twopi perhaps (FIXME).

**cov_xy** is the covariance (the third one of the second moments). It ranges between -1 and 1 with a value of 0 for a circular peak and +1 and -1 refering to elliptical shapes rotated by 90 degrees from each other (FIXME a picture would help). The definition of x and y are in terms of the fast and slow array indices as the image comes into memory, so it depends on how the image was stored in the file and the routine which read it in as to what you might finally get. This should later be better defined by the use of a better image file format, like ImageCIF, which defines such things. For most experiments at ID11 the position in the rotation scan should be specified by a value "Omega" in the image headers. If this is not the case there is a script called recoveromega.py which will read a text file containing lines with omega values followed by filenames to recover the appropriate information.

How does it work?
=================

For now, we are talking about the most up to date SVN version (0.8.1), a new release 
will be made soon. The algorithm used is based on a disjoint set, which is described in "Introduction to Algorithms" by Cormen, Leiserson., Rivest and Stein. An image is scanned row by row and each pixel is compared to a threshold value. If the pixel is above the threshold, it will be labelled as a peak. To determine the labels the pixel is compared to the previous pixel, and the pixels on the previous row. If one of these pixels is already labelled, the current pixel takes the same label. When there is a disagreement about labelling, the two labels are made the same using the "disjoint set".

| oooooooooooo    o = pixel below threshold
| oo1ooooo2ooo    1,2 = labels
| o111ooo222oo    X pixel where label 1 and label 2 must be made equivalent
| oo1111X

Within the sourcecode the routine "connectedpixels" is a compiled extension in C which takes a data image and computes an image of peak assignment labels. It optionally applies the dark and flood corrections when thresholding. The connectivity in the 2D image is therefore::

  0  0  0
  0  X  0
  0  0  0
  
To compute the properties of each connected object the routine "blobproperties" is then called with the data image and label image (this could be more efficient). It scans through the image forming properties on a pixel by pixel basis. In order to merge peaks on adjacent frames, the label images are compared. When exactly overlapping pixels are both labelled, then these labels should be made equivalent. The code which does it is called "bloboverlaps".

Files which are produced
========================

| Historically:
| 
| name [default: peaks.out] : 2D peaks, header info
| name_merge_t500 [default: peaks.out_merge_t500] : 3D peaks for threshold at 500
| 
| In the Brave New World (tm):
| 
| name.spt [add ".spt" to users request] : 2D peaks, header info
| foreach threshold:
| 
| name_t500.flt [add "_t500.flt" to users request] : 3D peaks for threshold at 500
| "users request" to default to stem_first_last in place of peaks.out.
| ouch - forgot the first and last ... c'est la vie
  
Using peaksearch for programmers
--------------------------------

Well, it used to be simpler than it is now. For the full pleasures of 3D  peaksearching in a program designed for fast 2D, one should read the ImageD11.labelimage code. The quick way to get an idea what thresholding will do is this::

  data = fabio.openimage("my_lovely_data.edf")
  blobs = numpy.zeros(data.shape, numpy.int)
  threshold = 1500.0
  npks = connectedpixels.connectedpixels(data, blobs, threshold, verbose=0)

On exit from connectedpixels the npks is the number of peaks found and blobs is an array of integer peak assignments (the argument is supposed to be modified during the function call). You can then have a look at blobs to see where the peaks are. In practice the peaksearch code will be doing a dark 
and flood correction too, which can add to the confusion.

For the 3D version it could be something like this::

  with h5py.File('mydata.h5','r') as h:
      frames = h['/entry/scan/detector/data'][()]
      omegas = h['/entry/scan/omega'][()]
  nframes, rows, cols = frames.shape
  lio = labelimage.labelimage( (rows, cols), peaksfile )
  for i in np.argsort(omegas):
      lio.peaksearch(frames[i], threshold, omegas[i])
      lio.mergelast()
  lio.finalise()




Various dataset input formats - python class input format
---------------------------------------------------------

Image input was handled by fabio, the code attempts to build a fabio file series from the command line.
At some point in 2018 the ESRF control system upgraded to using hdf5 files and the peaksearch.py code
has never quite adapted. A workaround is to have a python class defining the data for peaksearching,
see sandbox/hdfscan.py in the source code.

Example cases of input image data include:

#. edf files with omega angles in headers. The image number order is the order to peasearch::
   data0000.edf ... Omega=0
   data0001.edf ... Omega=0.5
   ...
   data0360.edf ... Omega=180.

#. edf files with no angles in headers (ftomo). The image number order may be wrong due to interlacing (--interlaced and --iflip options). You have to supply the omega angles in command line arguments.
  
#. bliss files (2018). Interlaced frames in two distant folders::
   Omega = scanfolder/interlaced_1_1/data.h5::/measurement/rot_master/mean_pos:rot_mean
   scanfolder/interlaced_1_1/Frelon/interlaced_1_1_Frelon0000.edf
   scanfolder/interlaced_1_1/Frelon/interlaced_1_1_Frelon0359.edf
   Omega = scanfolder/interlaced_1_2/data.h5::/measurement/rot_master/mean_pos:rot_mean
   scanfolder/interlaced_1_2/Frelon/interlaced_1_2_Frelon0000.edf
   scanfolder/interlaced_1_2/Frelon/interlaced_1_2_Frelon0359.edf 

#. hdf5 files (2020). Frames in order.

#. hdf5 files (2020). Interlaced rewind.

#. hdf5 files (2020). Interlaced zigzag.
   
#. hdf5 files (2020). Interlaced forwards (e.g. mod 360).


   

