# SegmentScan

**Description** This task segment a scan folder by extracting
information from *folder_file_config*, *segment_config* dicts.

**Inputs** - **folder_file_config** (dict): A dictionary output of the
*InitFolderFileConfig* - [Detector]{.title-ref} -
[OmegaMotor]{.title-ref} - [DtyMotor]{.title-ref} -
[RawScanFolderPath]{.title-ref} - [AnalysisPath]{.title-ref}

-   **segment_config** (dict): A dictionary containing the following
    keys:
    -   \`BgFilePath\`: Background Correction File Path a valid path
        string or None.
    -   \`MaskFilePath\`: Mask File Path a valid path string or None.
    -   \`Threshold\`: Zeroing Below Intensity, type: int
    -   \`SmoothSigma\`: Gaussian Filter Sigma, type: float
    -   \`Bgc\`: Remove Bg per peak, type: float
    -   \`MinPx\`: Required Min Pixels for a peak, type: int
    -   \`m_offset_thresh\`: Constantization if less than equal to,
        type: int
    -   \`m_ratio_thresh\`: Normalization within the range from greater
        than equal to this number: int

**Outputs** - **segmented_3d_columnfile**: A merged 2d peaks columnfile
path, it is located inside the
*AnalysisPath*/{sample}/{sample}\_{dataset}/segmented_3d_file.h5

**Usage** This task uses *collect_all_frames_peaks*,
*merge_all_frames_peaks*, *do3dmerge* functions from
frelon_peaksearch.py and *colfile_to_hdf* and \"colfile_from_dict\"
functions from columnfile.py of ImageD11 library.
