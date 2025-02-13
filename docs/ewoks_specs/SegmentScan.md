# SegmentScan

**Description** This task segment a scan folder by extracting
information from *folder_file_config*, *segment_config* dicts.

**Inputs**

-   **folder_file_config** (dict): A dictionary output of the
    *InitFolderFileConfig*
    -   *Detector*
    -   *OmegaMotor*
    -   *DtyMotor*
    -   *RawScanFolderPath*
    -   *AnalysisPath*
-   **segment_config** (dict): A dictionary containing the following
    keys:
    -   *BgFilePath* (str): Background Correction File Path a valid path
        string or None.
    -   *MaskFilePath* (str): Mask File Path a valid path string or
        None.
    -   *Threshold* (int): Zeroing Below Intensity
    -   *SmoothSigma* (float): Gaussian Filter Sigma
    -   *Bgc* (float): Remove Bg per peak
    -   *MinPx* (int): Required Min Pixels for a peak
    -   *m_offset_thresh* (int): Constantization if less than equal to
    -   *m_ratio_thresh* (int): Normalization within the range from
        greater than equal to this number

**Outputs**

-   **segmented_3d_columnfile**: A merged 2d peaks columnfile path, it
    is located inside the
    *AnalysisPath*/{sample}/{sample}\_{dataset}/segmented_3d_file.h5

**Usage**

This task uses *collect_all_frames_peaks*, *merge_all_frames_peaks*,
*do3dmerge* functions from frelon_peaksearch.py and *colfile_to_hdf* and
\"colfile_from_dict\" functions from columnfile.py of ImageD11 library.
