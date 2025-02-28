SpatialCorrection
=================

**Description**
This task takes the raw peak positions on the detector (s_raw, f_raw) as a columnfile, and a spatial correction file describing the spatial distortion of the detector, and performs a spatial correction on *segmented_3d_columnfile* to generate new (sc, fc) columns.

**Inputs**

- **segmented_3d_columnfile** (str): 
  A columnfile file path, located inside *AnalysisPath*/{sample}/{sample}_{dataset}/segmented_3d_file.h5

- **spatial_correction_files** (str):
  Either:
  A string - a single path to a .spline file (Fit2D spatial distortion standard)
  A list/tuple of strings - two paths to .e2dx/.e2dy files

**Note**: 
The *spatial_correction_files* .spline file path, .e2dx/.e2dy files are mutually exclusive.

**Outputs**

- **spatial_corrected_3d_columnfile** (str): A file path stored as .h5 file inside the *AnalysisPath*/{sample}/{sample}_{dataset}/spatial_corrected_3d_file.h5

**Usage**

The 3D peaks will be loaded as an `ImageD11.columnfile.columnfile` class instance using *colfile_from_hdf*.
Either ImageD11.blobcorrector.correct_cf_with_dxdyfiles or ImageD11.blobcorrector.correct_cf_with_spline functions will be used to compute the corrected peak positions.
The new columnfile will be saved to disk at *spatial_corrected_3d_columnfile* with the new columns.