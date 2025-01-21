"""
Python script to automatically end-to-end test our Jupyter notebooks
Currently implemented: nothing (indev)
To run this notebook, you need papermill in your Python environment.
As of 2025/01/21, this is not available in the default Jupyter environment.

I suggest to do the following:
cd to your ImageD11 git checkout folder
$ pip install papermill ansicolors -t . --no-deps

This file will add its parent folder (../) to the system path so it can be imported
So you get local ImageD11 and local papermill
"""
import sys
sys.path.insert(0, '../')

import os


import nbformat
import pytest
import papermill
from nbconvert.preprocessors import ExecutePreprocessor

nb_base_prefix = os.path.join('..', 'ImageD11', 'nbGui')
scan_nb_prefix = os.path.join(nb_base_prefix, 'S3DXRD')

# there are two levels of testing
# does the notebook work without errors?
# does the notebook give you the output you expect?


def notebook_exec(notebook):
    print('Testing notebook', notebook)
    with open(notebook) as f:
        nb = nbformat.read(f, as_version=4)
        ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
        try:
            assert ep.preprocess(nb) is not None, f"Got empty notebook for {notebook}"
        except Exception:
            assert False, f"Failed executing {notebook}"

            
def noteboook_exec_pmill(nb_input_path, nb_output_path, params_dict):
    papermill.execute_notebook(
       nb_input_path,
       nb_output_path,
       parameters=params_dict
    )
    
    
# test tomographic route in order
def tomographic_route():
    scan_nb_names = [
        'import_test_data.ipynb',
        'tomo_1_index.ipynb',
        'tomo_2_map.ipynb',
        'tomo_3_refinement.ipynb',
        '4_visualise.ipynb'
    ]
    scan_nb_paths = [os.path.join(scan_nb_prefix, name) for name in scan_nb_names]
    for notebook in scan_nb_paths:
        notebook_exec(notebook)

# test pbp route in order, after tomo route finished so we have a dataset file...
def pbp_route():
    scan_nb_names = [
        'pbp_1_indexing.ipynb',
        'pbp_2_visualise.ipynb',
        'pbp_3_refinement.ipynb',
        '4_visualise.ipynb'
    ]
    scan_nb_paths = [os.path.join(scan_nb_prefix, name) for name in scan_nb_names]
    for notebook in scan_nb_paths:
        notebook_exec(notebook)

        


def test_simplest_import():
    test_dir = 'import_test'
    scan_nb_names = [
        'import_test_data.ipynb',
    ]
    scan_nb_params = [
        {'download_dir': test_dir,
         'PYTHONPATH': sys.path[0]}
    ]
    scan_nb_paths = [os.path.join(scan_nb_prefix, name) for name in scan_nb_names]
    for notebook_in, scan_params in zip(scan_nb_paths, scan_nb_params):
        notebook_out = notebook_in.replace('.ipynb', '_out.ipynb')
        noteboook_exec_pmill(notebook_in, notebook_out, scan_params)

# def test_routes():
#     tomographic_route()
#     pbp_route()

    
if __name__=='__main__':
    print(ImageD11.__path__)
    print(papermill.__path__)