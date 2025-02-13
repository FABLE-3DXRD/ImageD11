# GeometryVectorTransformation

**Description** This task does the transformation from Peaks to the
Geometry Vector, mostly adding few more columns on the columnfile like
2theta, eta, xl, yl, zl (lab coordinates), gx, gy, gz (I think it is a
ray) and ds (reciprocal distance).

**Inputs\* -**TDXRD Geometry File Path\*\* (str): A valid Geometry TDXRD
.par file path - **spatial_corrected_3d_columnfile** (str): A file path
stored as .h5 file inside the
*AnalysisPath*/{sample}/{sample}\_{dataset}/spatial_corrected_3d_file.h5

**Outputs** - **spatial_corrected_3d_columnfile** (str): A file path
stored as .h5 file inside the
*AnalysisPath*/{sample}/{sample}\_{dataset}/spatial_corrected_3d_file.h5

**Usage** This task will create a folder called \'par_folder\' inside
{AnalysisPath}/{sample}/{sample}\_{dataset} folder Ship the provided
geometry tdxrd.par file into \'par_folder\',

It loads all the lines from geometry tdxrd .par file as key, value
pairs, and update the *spatial_corrected_3d_columnfile* using the
function updateGeometry(pars=loaded_geometry_key_value), by using the
columnfile class instance initialized using colfile_from_hdf of
columnfile.py.

**Note** We think, it is necessary to keep the geometry .par file inside
the analysis folder, so that it maintains the record of what recent .par
file parameter used for geometry update.

There will be another Ewoks Task will help you to create/modify a
geometry tdxrd .par file.
