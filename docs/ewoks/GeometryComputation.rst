GeometryComputation
============================

**Description**
This task takes the spatially corrected peak positions on the detector as a columnfile, and a description of the detector geometry in the laboratory frame, and computes the following for each peak:
2theta - radial angle of the peak away from the beam centre (degrees)
eta - azimuthal angle of the peak (degrees)
xl, yl, zl - positions of the spots that hit the detector in the lab frame
gx, gy, gz - reciprocal scattering vectors, in units of 1/(wavelength unit)
ds - reciprocal length of the scattering vectors, in units of 1/(wavelength unit)
These new arrays will be added as columns to the existing columnfile.

**Inputs*

- **TDXRD Geometry File Path** (str): A valid Geometry TDXRD .par file path - describes the detector position, wavelength, etc.

- **spatial_corrected_3d_columnfile** (str): A columnfile file path stored as .h5 file inside the *AnalysisPath*/{sample}/{sample}_{dataset}/spatial_corrected_3d_file.h5

**Outputs**

- **spatial_corrected_3d_columnfile** (str): A columnfile file path stored as .h5 file inside the *AnalysisPath*/{sample}/{sample}_{dataset}/spatial_corrected_3d_file.h5


**Usage**

This task will create a folder called 'par_folder' inside {AnalysisPath}/{sample}/{sample}_{dataset} folder 
The provided geometry .par file will be placed into 'par_folder'.

The spatially corrected peaks will be loaded as an `ImageD11.columnfile.columnfile` class instance using `colfile_from_hdf`.
The geometry .par file will be loaded as key, value pairs, and used to update the *spatial_corrected_3d_columnfile* using the function 
columnfile.updateGeometry(pars=loaded_geometry_key_value).
The new columnfile will be saved to disk at *spatial_corrected_3d_columnfile* with the new columns.

**Note**

We think, it is necessary to keep the geometry .par file inside the analysis folder, so that it maintains the record of 
what recent .par file parameter used for geometry update.

There will be another Ewoks Task to help you create/modify a geometry .par file.