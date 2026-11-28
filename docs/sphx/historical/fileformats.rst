============
File Formats
============
Unfortunately ImageD11 has evolved with some very poor file formats. These should be replaced as soon as possible with something more powerful, but in the meantime, here is some documentation about what is actually there right now. Mostly the files are ascii following a simplistic model of:

 # name = value

 or

 # name value

 or

 # name1  name2    name3
 value1   value2   value3

Image Formats
=============
FIXME - it uses fabio.open() 

Input/Output : 2D Peak search results
=====================================

The peaksearch.py script saves results in an ascii file which is read in again in the graphical interface (peaksearch menu). For each image processed, the file contains the header information on lines beginning with '#' a series of threshold levels and peaks. The file is typically named peaks.out.

Example with most of the header information snipped::

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
  44

The Omega header entry is currently the only thing which is relied on. If your raw images did not contain Omega header entries then a helper script recoveromegas.py is available to help to insert these values. 
The numerical values for each peak found in the image which are recorded are:

:Number_of_pixels: the number of pixels in the blob (units: pure integer)
:Average_counts: the average pixel intensity in the blob (units : ADU)
:x: the centre of mass peak position in the raw image in the fast/slow (CHECK) array direction (units: pixels)
:y: the centre of mass peak position in the raw image in the fast/slow (CHECK) array direction (units: pixels)
:xc: the x value corrected for spatial distortion, using the spline file mentioned in the header (units: pixels)
:yc: the x value corrected for spatial distortion, using the spline file mentioned in the header (units: pixels)
:sig_x: the second moment of the intensity distribution in the x direction (eg: Gaussian sigma, related to the peak width) (units: pixels)
:sig_y: the second moment of the intensity distribution in the y direction (eg: Gaussian sigma, related to the peak width) (units: pixels)
:cov_xy: the x/y covariance of the intensity distribution. Ranges from -1 to 1 with zero being a circular peak and extreme values representing ellipses. (units: pure float)
 
| Written by : scripts/peaksearch.py scripts/recoveromega.py 
| Read by : ImageD11/guipeaksearch.py 

Input/Output : Merged peaks
===========================

Versions > 1.0

Many columns, they are:

:sc: distortion corrected centre of mass, slow index direction
:fc: distortion corrected centre of mass, fast index direction
:omega: omega (rotation sequence) centre of mass
:Number_of_pixels: ... in the 3D connected object
:avg_intensity: average intensity of the pixels in the object
:s_raw: raw position (no spatial), slow direction
:f_raw: raw position (no spatial), fast direction
:sigs: second moment, slow direction
:sigf: second moment, fast direction
:covsf: covariance for fast/slow
:sigo: second moment, omega direction (out of plane of image)
:covso: covariance for slow/omega
:covfo: covariance for fast/omega
:sum_intensity: total intensity for pixels above threshold
:sum_intensity^2: summed intensity squared [not useful?]
:IMax_int: maximum pixel
:IMax_s: array index slow direction for maximum pixel
:IMax_f: array index fast direction for maximum pixel
:IMax_o: array index omega direction for maximum pixel
:Min_s: minimum pixel position in slow direction
:Max_s: maximum pixel position in slow direction
:Min_f: minimum pixel position in fast direction
:Max_f: maximum pixel position in fast direction
:Min_o: minimum pixel position in omega direction
:Max_o: maximum pixel position in omega direction
:dety: unlikely to be correct flipped direction
:detz: unlikely to be correct flipped direction
:onfirst: blob is present on first image
:onlast: blob is present on last image
:spot3d_id: line number - 1 (titles) in the file on creation by peaksearch.py

The peaksearching menu in the gui reads in the previous peaksearch output and merges peaks which have adjacent Omega values (and centre of mass within some tolerance) to give a new file which contains only the peaks. Typically named peaks.flt::

 # xc yc omega npixels avg_intensity x_raw y_raw sigx sigy covxy
 1020.355425 1458.043233 -37.802212 170.000000 2696.852941 1024.230029 1458.756529 2.060744 2.430557 -0.124729
 1060.621433 1456.483538 -37.750000 4.000000 1227.250000 1064.493991 1457.499083 0.499964 0.499999 -0.003485
 665.542804 728.847852 -37.750000 6.000000 1243.166667 672.161282 732.332082 1.060509 0.743424 0.767337

Column labels are as before but with the addition of an omega column giving the centre of mass of the blob in omega too. 

| Written by : ImageD11/guipeaksearch.py 
| Read by : ImageD11/guitransformer.py ImageD11/refinegrains.py

gvectorfile
============

Scattering vectors computed by the transformation module, normally via the guitransformer module... The first line of the file contains the unit cell parameters and lattice centering (one of P,A,B,C,I,F). The wavelength and wedge angle are expected on the next lines (needed to compute ideal two theta, omega and azimuth angles for computed orientations). A list of computed d* values expected h,k,l values for the unit cell follow, I don't think the program uses them anymore. Note that d* is 1/d-spacing in Angstrom. The actual scattering vectors follow::

 78.712000 78.712000 78.712000 90.000000 90.000000 90.000000 I
 # wavelength = 0.939500
 # wedge = -0.051189
 # axis 0.000000 0.000000 1.000000
 [... the remaining part of the parameter information is snipped ...]
 # ds h k l
  0.0179669   -1   -1    0
  0.0179669   -1    0   -1
  0.0179669   -1    0    1
  0.0179669   -1    1    0
  0.0179669    0   -1   -1
  0.0179669    0   -1    1
  0.0179669    0    1   -1
  0.0179669    0    1    1
  0.0179669    1   -1    0
  0.0179669    1    0   -1
  0.0179669    1    0    1
  0.0179669    1    1    0
  0.0254091   -2    0    0
  0.0254091    0   -2    0
  0.0254091    0    0   -2
  [ lots of hkl's snipped ]
 #  gx  gy  gz  xc  yc  ds  eta  omega  spot3d_id  xl  yl  zl
 -0.006049 -0.006821 -0.011411 1553.000000 1506.534500 0.014605 141.377789 41.000000 35325 211440.500505 -1811.094736 -2266.916293
 -0.012846 0.008203 -0.009504 1564.869800 1510.194700 0.017962 121.948198 122.022400 105544 211439.925254 -3027.977786 -1888.286052
 0.013156 0.009757 0.007648 1503.726500 1543.564400 0.018077 -64.969351 54.000300 46271 211436.535556 3254.006490 1519.487824
 0.000474 0.018005 -0.001647 1500.521000 1525.574600 0.018086 -95.223069 1.993700 3308 211438.759075 3578.198430 -327.094270
 -0.003736 0.010752 0.014069 1557.512700 1555.861300 0.018097 38.972120 160.000000 155486 211434.450605 -2261.233048 2795.170225
 -0.013099 -0.009687 -0.007903 1566.912500 1513.290200 0.018108 115.878145 53.000000 44665 211439.526323 -3236.750899 -1570.156341
 0.017023 -0.001041 -0.006196 1502.350000 1516.761000 0.018145 -109.965181 93.999200 78237 211439.811196 3388.245269 -1230.889139
 ...

The definitions follow:

:gx: x-component of scattering vector (along the beam) with all angles at zero (units 1/Angstrom)
:gy: y-component of scattering vector (toward the door) with all angles at zero (units 1/Angstrom)
:gz: z-component of scattering vector (roughly up) with all angles at zero (units 1/Angstrom)
:xc: spatially corrected peak x-position on detector (units: pixels)
:yc: spatially corrected peak y-position on detector (units: pixels)
:ds: 1/d-spacing - modulus of scattering vector: :math:`ds = \frac{\lambda}{ 2 \sin \theta}`.
:eta: azimuthal angle.
:omega: rotation angle of scan
:spotid: an spot identifier to follow the individual spot at any process step
:xl: x-component of scattering vector in the laboratory coordinate system (along the beam) with all angles at zero (units microns)
:yl: y-component of scattering vector in the laboratory coordinate system (toward the door) with all angles at zero (units microns)
:zl: z-component of scattering vector in the laboratory coordinate system (roughly up) with all angles at zero (units microns)

| Written by: ImageD11/guitransformer.py 
| Read by: ImageD11/indexing.py and GrainSpotter and possibly other programs.. 


Input/Output : UBI matrices
===========================

After successful(?) completion of some indexing a set of UBI matrices can be saved from 
the indexing.py script or guiindexing interface. These are just 3x3 matrices separated 
by blank lines::

  0.743414 5.539980 -1.573009
  5.484442 -1.159425 -1.505417
  -1.758983 -1.292534 -5.361137

  -3.422333 4.030841 -2.331791
  -3.311710 -4.157822 -2.286774
  -3.261669 -0.019438 4.770287
  
| It should be the case that (CHECK/TESTCASE) 
| <math> \begin{matrix} h & k & l \end{matrix} \begin{matrix} UBI_{11} & UBI_{12} & 
| UBI_{13} \\ UBI_{21} & UBI_{22} & UBI_{23} \\ UBI_{31} & UBI_{32} & UBI_{33} 
| \end{matrix} \begin{matrix} xr & yr & zr \end{matrix} </math> 
| Written by: ImageD11/indexing.py ImageD11/guiindexer.py 
| Read by: ImageD11/refinegrains.py

Parameter files
===============

This files should all have the format name value.  
For the transformation module the parameters are::

 Unit cell
 cell_a
 cell_b
 cell_c
 cell_alpha
 cell_beta
 cell_gamma
 cell_lattice_[P,A,B,C,I,F,R] F

Diffractometer angles and geometry:
 
:chi: rotation of detector around beam (not tested, probably only in CVS)
:wedge: tilt of axis around y
:distance: sample detector distance (units: millimetres)
:tilt-y: first detector tilt
:tilt-z: second detector tilt
:wavelength: of the incoming x-rays
:y-center: beam centre on the detector (units: pixels)
:y-size: pixel size in the y direction (units: microns)
:z-center: beam centre on the detector (units: pixels)
:z-size: pixel size in the z direction (units: microns)

CHECK that y/z are the right way around and also the hidden detector flip matrix. fit_tolerance is used to decide which spots to assign to hkl rings in fitting the geometrical parameters (units: degrees). 

The indexing parameters are:
FIXME
