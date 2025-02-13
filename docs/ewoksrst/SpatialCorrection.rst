SpatialCorrection
=================

**Description**
This task does spatial correction on *segmented_3d_columnfile*

**Inputs**

- **segmented_3d_columnfile** (str): 
  A merged 2d peaks columnfile path, it is located inside *AnalysisPath*/{sample}/{sample}_{dataset}/segmented_3d_file.h5

- **spatial_correction_files** (str): 
  A string or list of string, it is farmer if the correction is based on spline file, 
  it is later when the correction is based on E2dxy files.

**Note**: 
The *spatial_correction_files* Spline File path, E2Dx,y files are mutually exclusive.

**Outputs**

- **spatial_corrected_3d_columnfile** (str): A file path stored as .h5 file inside the *AnalysisPath*/{sample}/{sample}_{dataset}/spatial_corrected_3d_file.h5

**Usage**

This task uses *colfile_from_hdf*, and colfile_to_hdf* from columnfile.py python file to extract merged 2D peaks and to save the spatial correction result.
It uses either correct_cf_with_dxdyfiles or correct_cf_with_spline functions to correct the deformed detector plane from the file blobcorrector.py of ImageD11 library.


**Question**

Is there a clever way of dumbing all the configuration of *segmenter parameters*, *folder file* settings inside the 
spatial_corrected_3d_file.h5 or segmented_3d_file.h5 without affecting the columnfile usage inside ImageD11?
