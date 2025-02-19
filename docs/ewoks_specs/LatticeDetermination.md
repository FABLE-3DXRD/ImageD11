# LatticeDetermination

**Description** Assemble the given lattice parameter and space group
into a parameter (.par ) file. Provides naming (tag on the phase)
mechanism to play safely around with Multiple lattice version of same
Lattice (.par).

**Inputs**

-   **geometry_transformed_3d_columnfile** (str): A file path stored as
    .h5 file inside the
    *AnalysisPath*/{sample}/{sample}\_{dataset}/geometry_transformed_3d_peaks_file.h5

-   **lattice_lengths_angle** (tuple): A tuple of length 6 contains
    (cell\_\_a, cell\_\_b, cell\_\_c, cell_alpha, cell_beta, cell_gamma)

-   **lattice_space_group** (str or integer): if it is integer it should
    be between 1 to 230 inclusive. if it is string, it should be in
    \"P\",\"A\",\"B\",\"C\",\"I\",\"F\",\"R\"

-   **LatticeName** (str): A Lattice Name for the provided lattice
    parameters.

-   **Tag** (str) (Optional): A short string *tag*, so that user can
    play with same lattice name but with multiple version.

    The provided lattice parameter along with space group will be
    aranged as .par file in the path
    *AnalysisPath*/{sample}/{sample}\_{dataset}/par_folder/{LatticeName}\_{Tag}.par

**Outputs**

-   **dstar** (Numby 1D Array): Reciprocal Distance
-   **eta** (Numpy 1D Array): The angle, this array and dstar array
    should have same length.
-   **ucell_dstar** (Numpy 1D Array): Provides the dstar value for the
    provided phase .par file, Its length depends on the multiplicity
    (crystallographic stuff).

**Usage**

First will create a lattice parameter .par file in \'par_folder\' of
\'{analysis_path}/{sample}/{sample}\_{dataset}\' folder It checks the
provided geometry_transformed 3d column file is correct by way is
checking eta, two_th etc. Finally provides you the data, to check the
lattice matches the sample or not. We use unitcell and columnfile
classes here.
