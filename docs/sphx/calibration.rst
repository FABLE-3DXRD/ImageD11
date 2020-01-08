===========
Calibration
===========
Calibration of the sample-detector distance, detector tilts and beam centre are one of the fundamental reasons for creating the ImageD11 program and having the sometimes annoying graphical interface. The menu offers the option to read in previously determined parameters from a file::

  cell__a 5.8
  cell__b 5.8
  cell__c 5.8
  cell_alpha 90.0
  cell_beta 90.0
  cell_gamma 90
  cell_lattice_[P,A,B,C,I,F,R] F
  distance 294662.658
  fit_tolerance 0.2
  o11 1
  o12 0
  o21 0
  o22 -1
  omegasign 1.0
  t_x 129.322427382
  t_y 11.8135159025
  t_z 0
  tilt_x 0.00000
  tilt_y -0.013173
  tilt_z 0.002378
  wavelength 0.154
  y-center 1016.328171
  y-size 48.0815
  z-center 984.924425
  z-size 46.77648
  
The meaning of the parameters is mostly self explanatory. Begins with the unit cell parameters (in angstroms/degrees), the lattice centering using the first letter of the space group symbol. The sample detector distance is normally in millimeters. The fit tolerance is the difference in degrees between observed and computed peak positions to use in parameter refinement. The tilts are in radians about the y/z axes in the ID11 co-ordinate system (z is vertical, y is horizontal). The position of the beam centre and pixel sizes are used to transform the spatially corrected peak positions into real space peaks positions. You MUST use the values for pixel size recorded in the spline file if you have made a spatial distortion for things to work out properly. Having read in some data which was output from the peaksearching menu you can make a plot of two theta / eta and then add unit cell peaks.


Perhaps you will need to go to the plotting menu and clear plot before 
making this plot to get a sensible x/y range.

The horizontal axis represents the peak positions in two theta with the eta 
value being the angle around a powder ring measured from the vertical.

If the initial parameters for the transformation are good then you will hopefully find nearly vertical lines made up of the diffraction spots from all of your grains with crosses overlaying them indicating the computed peak positions. If the parameters are not so good then the lines will be deformed and perhaps not match the computed positions. To improve the fit you can change the parameters manually via the edit parameters menu option and then replot. You can also choose to fit the parameters where peaks which are within fit_tolerance of the computed positions are used for a simplex minimisation procedure. Usually this works better if you zoom in on the low angle data first and then repeat the fit after replotting and unzooming to use the full range of data. It works better that way because you are less likely to mis-assign low angle peaks as there is less peak overlap.

Once the parameters are refined and the transformation gives vertical lines you can hope to make a good computation of g-vectors. Note that the wavelength is not normally refined and should be given as correctly as possible to obtain the correct curvature of the Ewald sphere (calibrate it via some other method please). Changing wavelength is eventually equivalent to changing the unit cell parameters so you might get away with a poor initial value.

Before computing the g-vectors (scattering vectors or normals to the diffracting hkl planes) you should check the value of the wedge angle: this is zero when the incident beam is perpendicular to the rotation axis. In cases where a heavy stress rig is used at beamline ID11 together with the Laue optics there is an inclination angle between the rotation axis and the beam - named wedge after the motor in spec for adjusting that angle. You should use the negative value of the two theta angle of the monochromator (eg about -2.8 degrees at 80keV). Once you have clicked on compute g-vectors you can save these vectors into a file. There is also an option for saving a file suitable for introduction into the Graindex software developed by the Riso group (FIXME encourage an expert to explain how to do that).

Making the transformations independent of the gui for scripting
===============================================================

If you have a parameter file from a calibration you can directly make the gve file like this::

 $ python
 from ImageD11 import transformer
 obj = transformer.transformer()
 obj.loadfiltered(  "peaks_t20.flt" )
 obj.loadfileparameters(  "my_fine_parameters.prm" )
 obj.compute_tth_eta( )
 obj.addcellpeaks( )
 obj.computegv( )
 obj.savegv(  'peaks_t20.gve' )
 
... and that's why they call it "object oriented" programming.
