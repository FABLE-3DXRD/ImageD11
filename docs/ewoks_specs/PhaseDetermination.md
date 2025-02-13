# PhaseDetermination

**Description** This task provide necessary data for choosing the phase
.par file, and provides the safe mechanism to play around with Poly
(Multi) Phase lattice parameter (.par) files with multiple version.

**Inputs** - **spatial_corrected_3d_columnfile** (str): A file path
stored as .h5 file inside the
*AnalysisPath*/{sample}/{sample}\_{dataset}/spatial_corrected_3d_file.h5 -
**PhaseName** (str): A Phase Name for the provided phase file path. -
**Phase_file_path** (str) (Optional) : A Phase .par file path - **Tag**
(str) (Optional): An optional short string *tag*, so that user can play
with same phase name but with different values, it is optional.

> If *Tag* was not provided, then provided .par file (or .par for the
> given PhaseName) will be shiped to as
> *AnalysisPath*/{sample}/{sample}\_{dataset}/par_folder/{PhaseName}.par
>
> If *Tag* provided then provided .par file (or .par for the given
> PhaseName) will be shiped to as
> *AnalysisPath*/{sample}/{sample}\_{dataset}/par_folder/{Tag}/{PhaseName}.par

**Outputs** - **dstar** (Numby 1D Array): Reciprocal Distance - **eta**
(Numpy 1D Array): The angle, this array and dstar array should have same
length. - **ucell_dstar** (Numpy 1D Array): Provides the dstar value for
the provided phase .par file, Its length depends on the multiplicity
(crystallographic stuff).

-   **spatial_corrected_3d_columnfile** (str): A file path of .h5 file
    inside
    *AnalysisPath*/{sample}/{sample}\_{dataset}/spatial_corrected_3d_file.h5

\- **PhaseName**: Same as input *PhaseName* \` **Tag**: Same as input
*Tag*

**Usage**

This task first will check the presence of \'par_folder\' in
\'{analysis_path}/{sample}/{sample}\_{dataset}\' folder also the
geometry tdxrd.par file inside \'par_folder\', if it is not there, this
task will throws error. May be the correct way is, reading
*spatial_corrected_3d_columnfile*, and check whether it has *two_th,
eta, xl, yl, zl, gx, gy, gz, eta, etc*, to enusre that g-vector
calculation was actually done. Because without g-vector we can\'t do
phase Determination/filtering.
