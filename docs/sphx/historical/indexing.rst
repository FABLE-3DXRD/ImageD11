Indexing
========
Once your diffraction spots have been found and converted to g-vectors (scattering vectors) you can attempt to find the orientations of the grains in the sample. ImageD11 only includes rudimentary facilities for indexing suitable for relatively small numbers of grains in the sample (probably <100, although the limit depends on the data and unit cell).

Indexing can be carried out via the indexing menu on the graphical user interface. Begin by reading in the g-vectors which were saved at the and of the transformation procedure. If you have collected data for a single crystal then a useful diagnostic is produced by the plot x/y/z menu option which gives a 3 dimensional representation of the observed data. For a single crystal these should form a nice regular lattice if everything has gone well. For samples containing many grains you get something like a 3d representation of the pole figure measured, but with no intensity information currently.
(FIXME - colorize the plot). 

Indexing then proceeds by assigning hkl values to two peaks and using these indexed peaks together with the unit cell parameters to determine the orientation of a grain. Trial orientations are scored according to how many peaks they index. The algorithm needs some parameters to carry out this procedure. The ds_tol parameter is used to determine which powder style diffraction rings overlap and which peaks to assign hkl values according to the length of the scattering vector. The units are normally 1/angstroms if you input unit cell parameters and wavelength in angstroms. Peaks are considered if they coincide with computed positions within a tolerance of 1/d.

Once the peaks are assigned to powder rings then the program outputs information in the console window about the assignments. Orientations can be generated using two rings selected via the indexing parameters. For each peak (i) found in the first ring the (cosine of the) angle between peak (i) and all of the peaks in the second ring in computed. The two peaks can be in the same ring if it has a non-parallel lattice components (eg (100) and (010) are OK for cubic unit cells, but not for triclinic). For each peak the first ring the program identifies the best agreeing peak in the second ring in terms of the (cosine of the) angle between them, compared to the computed values from the unit cell parameters. A trial orientation is generated if the cosine is within cosine_tol of the theoretical value. Trial orientations can then be scored to see how many peaks they index, within a tolerance of hkl_tol which is applied in terms of the computed hkl indices from the positions the observed peaks. Also a threshold level of minpks is used to ensure that a new orientation matrix indexes at least minpks peaks which were previously unindexed. Usually problems with indexing just need an increase of the hkl_tol value, perhaps combined with a review of the parameters and data quality.

The generated orientation gives a matrix relating hkl indices to g-vectors via the equation::

  h = UBI g
  g = UB h
  
where UBI is the matrix the program writes out and h and g are the hkl indices and reciprocal vector respectively. UBI is the inverse of UB. ImageD11 does NOT introduce a factor of 2pi into the g-vectors as many programs do. Unit cell parameters can then be derived from (UBI).transpose(UBI) giving the metric tensor. The program does this when you write out indexed peaks.

By default the UBI matrices are refined according to the following algorithm::

     From Paciorek et al Acta A55 543 (1999)
     UB = R H-1
        where:
        R = sum_n r_n h_n^t
        H = sum_n h_n h_n^t
        r = g-vectors
        h = hkl indices
		
The values determined depend on the hkl_tol value used and correpondingly the peaks used for refinement (FIXME future, introduce a weighting depending on the goodness of fit of a particular peak). (Dated) sample output from writing out indexing is shown below. Mostly it should be self explanatory::

  Grain: 0   Npeaks=18   <drlv>=0.044535
  UBI:
  [[-1.56672605  2.96622445  4.73147779]
   [-2.91865211  3.75463023 -3.32027434]
   [-4.76097042 -3.27784447  0.47843108]]
  Peak   (  h       k       l      )   drlv             x       y    Omega_obs Omega_calc   Eta_obs Eta_calc   tth_obs tth_calc
  11     ( -1.9730 -0.0505  0.0020 )   0.05733662    757.3  1181.7    -34.3625  -30.4737      143.4444  144.6923      3.0028    3.0430
  12     ( -0.0134  1.9746 -0.0092 )   0.03010606    828.5   717.8    -33.7500  -36.0039     -125.4869 -124.9361      3.0045    3.0430
  [...lots of output snipped...]
  And now listing via peaks which were assigned to rings

  Peak= 2     Ring= 0     gv=[ -0.0814  0.2322 -0.1629 ]   omega=  -15.8973   eta= -123.6402   tth=    2.6044
  Peak not assigned, closest=[  0.0455  1.6506 -0.4516 ] for grain 0

  Total number of peaks was 293
  Peaks assigned to grains 34
  Peaks assigned to rings but remaining unindexed 195
  Peaks not assigned to rings at all 64
  
Nowadays the unit cell parameters should be shown as well. The orientation matrix can be determined from the (UBI) matrix by experts. Non experts should just substitute (1,0,0), then (0,1,0) and finally (0,0,1) into the equations for h and g above to find out where their diffraction planes are pointing in the sample. A recent feature added to the program is a plot of the number of peaks found with gradually increasing hkl_tol value to use for validation and setting of this parameter. A real grain should have a peak close to zero tolerance, and then randomly appearing peaks should increase as something roughly proportional to the cube of the tolerance for a texture free sample.
