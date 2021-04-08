
"""ImageD11 project files

We want to be able to:
 - Load / save your work
 - Log what has and has not been done

Use cases:

 - Single 3DXRD scan, far field detector for strain / center of mass
    - List of source images, omega angles
    - Files / information to get corrected data (dark, flat, spline, monitor)
    - Background estimation
    - Thresholding
    - Peak locations (per image)
    - Merging peaks in 3D across omega steps
    - Conversion to scattering vectors, geometry calibration
    - Assignment of peaks to grains
    - Refinement of grain sizes, positions, unit cells, intensities per (h,k,l)

 - Near field scan
    - Interface to DCT code !
    - List of source images, omega angles
    - Files / information to get corrected data (dark, flat, spline, monitor)
    - Background estimation
    - Thresholding
    - Peak locations (per image)
    - Merging peaks in 3D across omega steps

 - Scanning 3DXRD scans (difftomo)
    - Same as single 3DXRD but also with a y-step, so:
    - List of source images, omega angles, y-steps
    
 - 3DXRD scans (difftomo) as a function of temperature, time, load, etc
    - As above but with extra "load" column

We need to save data with "parent"-"child" relationships so we can see what 
items "depend-on" other items. This should map onto typical file formats
and work across programming languages. We should try to take inspiration from 
NeXus and ImageCIF (etc) so far as resources allow. Also work with existing 
experimental data. Also get something working in the short term (e.g. for the 
next 5 years).

For now we will use a nested python dictionary to hold the project. This
will be loaded and saved to a yaml (or json, pickle, etc) file.

Initial jobs from this experiment description file:
    - Locate the interesting things in it
    - Validate whether the dataURLs exist (md5 matches if recorded, etc)
    - Assert that all of the dataURLs are unique (no repeats)
    - Load / save a description file


Processing Description:
    - Reference to the raw data description
    - Get corrected images (monitor, dark, flat, mask)
    - Concept of "previous" and "next" images in scan
    - Estimate background (single image method, file series method)

"""