SegmentScan
===========

**Description**
This task segments an entire scan folder by extracting information from *folder_file_config*, *segment_config* dicts.

**Inputs**

- **folder_file_config** (dict): A dictionary output of the *InitFolderFileConfig*
  - *Detector* (str): The detector string to look for in the RAW_DATA masterfile under `{scan}/measurement/{detector}` - e.g. "frelon1", "frelon3", or "eiger".
  - *OmegaMotor* (str): The omega motor string to look for in the RAW_DATA masterfile under `{scan}/instrument/positioners` - e.g. "diffrz".
  - *DtyMotor* (str): The dty motor string (unused for regular 3DXRD, but required) to look for in the RAW_DATA masterfile under `{scan}/instrument/positioners` - e.g. "diffty".
  - *RawScanFolderPath*  (str): Path to the raw scan folder, its parent should contain Masterfile, will be checked by this task as well.
  - *AnalysisPath*  (str): Analysis Path, the result of subsequent value will be stored in AnalysisPath/{sample}/{sample}_{dataset}

- **segment_config** (dict): A dictionary containing the following keys:
  - *BgFilePath* (str): Background Correction File Path a valid path string or None.
  - *MaskFilePath* (str): Mask File Path a valid path string or None.
  - *Threshold* (int): Zeroing Below Intensity
  - *SmoothSigma* (float): Gaussian Filter Sigma
  - *Bgc* (float): Remove Bg per peak
  - *MinPx* (int): Required Min Pixels for a peak
  - *m_offset_thresh* (int): Constantization if less than equal to
  - *m_ratio_thresh* (int): Normalization within the range from greater than equal to this number


**Outputs**

- **segmented_3d_columnfile**: A merged 3d peaks columnfile path, it is located inside the *AnalysisPath*/{sample}/{sample}_{dataset}/segmented_3d_file.h5

**Usage**

This task uses *collect_all_frames_peaks*, *merge_all_frames_peaks*, *do3dmerge* functions from `ImageD11.frelon_peaksearch` and 
*colfile_to_hdf* and "colfile_from_dict" functions from `ImageD11.columnfile`.
