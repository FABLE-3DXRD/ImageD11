

Making a new peaksearch.

Input:
    Scan information:
        # Counters (including monitor, etc)
        Detector_1:
            List of images (URL's or counts)
        Detector_2:
            List of images (URL's)
        ...
	# Motors (including Epoch, etc)
	Motor_1:
            List of positions (values or URL's)
	Motor_2:
            List of positions (values or URL's)
        ...
	Trajectories:
	    Lists of connections between scan points
	    Might be actual or virtual (e.g. what is the inner loop)


One big table of _all_ of the images/counters/positions recorded. Assuming
single points in each position.

Links between these would give subsets of points. e.g.
      1D scan on omega
      Move dty
      1D scan on omega
      Move dty
      ... == 2D scan on omega/dty

We can want to see a rocking curve in any scan direction (dty or omega).

Realistic case:
    Continuous rotation on omega
    Sinusoidal movement on dty (prime factor on rotation speed)
    Covers sinogram in diagonal directions
    
"Experiment" is a database of "scan_points".
     Point_id, results (images, motor pos, monitor, time, etc)
	   
Trajectories:
     Lists of links between the scan points (directed for now)
         Link_id(k) : Point_id(i) -> Point_id(j)

Examples:

	sweep_musst type scan:
             images = ["data%04d.edf"%(i) for i in range(npts)]
	     omega  = [ start + step*i    for i in range(npts)]
	     ... deal with missing images on the file system
	interlaced scans:
	     ... iflip or not
	2D scans y/omega:
	     ... sweep_musst style or interlaced


Processing
==========

"Store": an intermediate result that the user should be able
to inspect. It might be written to disk, or the recipe to
make it again on the fly is stored (or linked).


Point by point :
   Corrections (dark, flat, ...)
   Store:
       Point_id/detector_id/corrected_image
   Normalisation (monitor or srcur, ...)
   Store:
       Point_id/detector_id/normalised_image
   Radial integration
   Store:
       Point_id/detector_id/1D_integrations
   Background estimation / removal (can depend on neighbor images or 1D)
   Store:
       Point_id/detector_id/background
   Masking (diamond peaks, ice rings, detector gaps)
   Thresholding (peak versus smooth background)
   Store:
       Point_id/detector_id/mask [bg/peak/static_mask/dynamic_mask]
   Object labelling (group pixels into lists of spots)
   Store:
       Point_id/detector_id/peak_id/lists_of_[pixels/values/weights]
   Object properties (intensity, center of mass, etc)
   Store:
       Point_id/detector_id/peak_id/[list_of_properties]
		                  

Trajectories :
   1D case was scan on omega
   Check overlap of peaks objects from one frame to next
      Score for whether peaks are overlapped could be fraction of pixels
      overlapped or mismatch in position, etc
      image_id(i) -> image_id(j):
        list of peak_id(i) -> peak_id(j) : overlap_score(s)
   Collect "unique" peaks (cliques) or clusters
   Any clique has 1 or more members:
      list of (2X point_id/image_id/peak_id)
   Store:
      

   












