SegmentFrame
============

**Description**
This task segment a single frame from the provided **folder_file_config**  dict information

**Inputs**

- **folder_file_config** (dict): A dictionary output of the *InitFolderFileConfig*
  - *Detector* (str): The detector frelon1, frelon3, or eiger.
  - *OmegaMotor* (str): The omega motor diffrz, diffry, or diffrx.
  - *DtyMotor* (str): The dty motor diffty, difftx, or difftz.
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

- **raw_image** (2D Numpy array) : A 2D array to plot.

- **bg_corrected_image** (2D Numpy array): A 2D array to plot.

- **found_peaks** (List contains two arrays): A highlighted peaks from the selected frame (plot using x-array and y-array positions).

- **folder_file_config** (dict):

- **segment_config** (dict):

**Usage**

This task uses *choose_frame* function from the file *segmenter_gui_utils.py* ImageD11 library.
It uses our saved omega values (from a file inside analysispath), and the values provided in *folder_file_config*, *segment_config* config.
