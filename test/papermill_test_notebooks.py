# flake8: noqa
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

            
def noteboook_exec_pmill(nb_input_path, nb_output_path, params_dict):
    print('Executing notebook', nb_input_path)
    # change output path if it already exists, in case we run the same notebook twice
    if os.path.exists(nb_output_path):
        nb_output_path = nb_output_path.replace('.ipynb', '_2.ipynb')
    papermill.execute_notebook(
       nb_input_path,
       nb_output_path,
       parameters=params_dict
    )

def notebook_route(base_dir, notebook_paths, notebook_param_dicts, notebook_out_dir=None):
    """
    Execute multiple notebooks in the order they are given.
    base_dir: The path to the output folder for the test. Must not already exist.
    notebook_paths: Ordered list of paths to the notebooks, to be deployed one after the other
    notebook_param_dicts: Ordered list of dictionaries of parameters, one dict per notebook to be executed
    """
    if os.path.exists(base_dir):
        raise ValueError('output test directory already exists:', base_dir)
    os.mkdir(base_dir)
    if notebook_out_dir is None:
        notebook_out_dir = os.path.join(base_dir, 'nb_out')
    os.mkdir(notebook_out_dir)
    for notebook_in_path, notebook_param_dict in zip(notebook_paths, notebook_param_dicts):
        notebook_out_path = os.path.join(notebook_out_dir, os.path.split(notebook_in_path)[1].replace('.ipynb', '_out.ipynb'))
        noteboook_exec_pmill(notebook_in_path, notebook_out_path, notebook_param_dict)


# test the full tomographic route from start to finish
def test_tomographic_route():
    tomo_dir = 'tomo_route'                    
    scan_nb_names = [
        'import_test_data.ipynb',
        'tomo_1_index.ipynb',
        'tomo_2_map.ipynb',
        'tomo_3_refinement.ipynb',
        '4_visualise.ipynb'
    ]
    dset_file = os.path.join(tomo_dir, 'processed', 'Si_cube', 'Si_cube_S3DXRD_nt_moves_dty', 'Si_cube_S3DXRD_nt_moves_dty_dataset.h5')
    scan_nb_params = [
        {'download_dir': tomo_dir,  # import_test_data.ipynb
         'PYTHONPATH': sys.path[0]},
        {'PYTHONPATH': sys.path[0],  # tomo_1_index.ipynb
         'dset_file': dset_file,
         'phase_str': 'Si',
         'cf_strong_frac': 0.9939,
         'cf_strong_dsmax': 1.594,
         'cf_strong_dstol': 0.005,
         'indexer_ds_tol': 0.01,
         'rings_for_gen': [0, 1, 3],
         'rings_for_scoring': [0, 1, 2, 3, 4],
         'hkl_tols_seq': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.075],
         'fracs': [0.9, 0.7],
         'max_grains': 1000,
         'peak_assign_tol': 0.05,
        },
        {'PYTHONPATH': sys.path[0],  # tomo_2_map.ipynb
         'dset_file': dset_file,
         'phase_str': 'Si',
         'cf_strong_frac': 0.9939,
         'cf_strong_dstol': 0.005,
         'is_half_scan': False,
         'halfmask_radius': 25,
         'peak_assign_tol': 0.25,
         'manual_threshold': None,
         'hkltol': 0.25,
         'correct_sinos_with_ring_current': False,
         'first_tmap_cutoff_level': 0.4,
         'niter': 500,
         'second_tmap_cutoff_level': 0.05
        },
        {'PYTHONPATH': sys.path[0],  # tomo_3_refinement.ipynb
         'dset_file': dset_file,
         'phase_str': 'Si',
         'default_npks': 20,
         'default_nuniq': 20,
         'hkl_tol_origins': 0.05,
         'hkl_tol_refine': 0.1,
         'hkl_tol_refine_merged': 0.05,
         'ds_tol': 0.004,
         'ifrac': 7e-3,
         'rings_to_refine': None,
         'use_cluster': False
        },
        {'PYTHONPATH': sys.path[0],  # 4_visualise.ipynb
         'dset_file': dset_file,
         'phase_str': 'Si',
         'min_unique': 400,
        }
    ]
    scan_nb_paths = [os.path.join(scan_nb_prefix, name) for name in scan_nb_names]
    notebook_route(tomo_dir, scan_nb_paths, scan_nb_params)
    


# test the full point-by-point route from start to finish
def test_pbp_route():
    tomo_dir = 'pbp_route'                    
    scan_nb_names = [
        'import_test_data.ipynb',
        'pbp_1_indexing.ipynb',
        'pbp_2_visualise.ipynb',
        'pbp_3_refinement.ipynb',
        '4_visualise.ipynb'
    ]
    dset_file = os.path.join(tomo_dir, 'processed', 'Si_cube', 'Si_cube_S3DXRD_nt_moves_dty', 'Si_cube_S3DXRD_nt_moves_dty_dataset.h5')
    scan_nb_params = [
        {'download_dir': tomo_dir,  # import_test_data.ipynb
         'PYTHONPATH': sys.path[0]},
        {'PYTHONPATH': sys.path[0],  # pbp_1_indexing.ipynb
         'dset_file': dset_file,
         'phase_str': 'Si',
         'minpkint': 0,
         'hkl_tol':  0.025,
         'fpks': 0.9,
         'ds_tol': 0.004,
         'etacut': 0.1,
         'ifrac': 5e-3,
         'y0': 24.24,
         'symmetry': "cubic",
         'foridx': [0, 1, 3, 5, 7],
         'forgen': [1, 5, 7],
         'uniqcut': 0.85,
         'use_cluster': False
        },
        {'PYTHONPATH': sys.path[0],  # pbp_2_visualise.ipynb
         'dset_file': dset_file,
         'phase_str': 'Si',
         'min_unique': 20
        },
        {'PYTHONPATH': sys.path[0],  # pbp_3_refinement.ipynb
         'dset_file': dset_file,
         'phase_str': 'Si',
         'min_unique': 20,
         'manual_threshold': None,
         'y0': 24.24,
         'hkl_tol_origins': 0.05,
         'hkl_tol_refine': 0.1,
         'hkl_tol_refine_merged': 0.05,
         'ds_tol': 0.004,
         'ifrac': 7e-3,
         'rings_to_refine': None,
         'set_mask_from_input': False,
         'use_cluster': False
        },
        {'PYTHONPATH': sys.path[0],  # 4_visualise.ipynb
         'dset_file': dset_file,
         'phase_str': 'Si',
         'min_unique': 400,
        }
    ]
    scan_nb_paths = [os.path.join(scan_nb_prefix, name) for name in scan_nb_names]
    notebook_route(tomo_dir, scan_nb_paths, scan_nb_params)


def test_FeAu_JADB_tomo():
    tomo_dir = '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/PROCESSED_DATA/20250123_JADB/tomo_route'                    
    scan_nb_names = [
        '0_segment_and_label.ipynb',
        'tomo_1_index.ipynb',
        'tomo_1_index_minor_phase.ipynb',
        'tomo_2_map.ipynb',
        'tomo_2_map_minor_phase.ipynb',
        'tomo_3_refinement.ipynb',  # for major phase
        'tomo_3_refinement.ipynb',  # for minor phase
        '4_visualise.ipynb',  # for major phase
        '4_visualise.ipynb',  # for minor phase
        '5_combine_phases.ipynb',
        '6_stack_layers.ipynb'
        
    ]
    sample = 'FeAu_0p5_tR_nscope'
    dataset = 'top_200um'  # first of two layers
    dset_file = os.path.join(tomo_dir, sample, f'{sample}_{dataset}', f'{sample}_{dataset}_dataset.h5')
    scan_nb_params = [
        {'PYTHONPATH': sys.path[0],  # 0_segment_and_label.ipynb
         'maskfile': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/mask_with_gaps_E-08-0173.edf',
         'e2dxfile': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/e2dx_E-08-0173_20231127.edf',
         'e2dyfile': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/e2dy_E-08-0173_20231127.edf',
         'detector': 'eiger',
         'omegamotor': 'rot_center',
         'dtymotor': 'dty',
         'options': { 'cut' : 1, 'pixels_in_spot' : 3, 'howmany' : 100000 },
         'dataroot': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/RAW_DATA/',
         'analysisroot': tomo_dir,
         'sample': 'FeAu_0p5_tR_nscope',
         'dataset': 'top_200um',
         'dset_prefix': "top_"
        },
        {'PYTHONPATH': sys.path[0],  # tomo_1_index.ipynb
         'dset_file': dset_file,
         'par_file': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/pars.json',
         'phase_str': 'Fe',
         'cf_strong_frac': 0.9939,
         'cf_strong_dsmax': 1.594,
         'cf_strong_dstol': 0.005,
         'indexer_ds_tol': 0.01,
         'rings_for_gen': [0, 1, 3],
         'rings_for_scoring': [0, 1, 2, 3, 4],
         'hkl_tols_seq': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.075],
         'fracs': [0.9, 0.7],
         'max_grains': 1000,
         'peak_assign_tol': 0.05,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # tomo_1_index_minor_phase.ipynb
         'dset_file': dset_file,
         'major_phase_str': 'Fe',
         'minor_phase_str': 'Au',
         'major_phase_cf_frac': 0.99418,
         'major_phase_cf_dsmax': 1.594,
         'major_phase_cf_dstol': 0.0035,
         'minor_phase_cf_frac': 0.9975,
         'minor_phase_cf_dsmax': 1.594,
         'minor_phase_cf_dstol': 0.0045,
         'indexer_ds_tol': 0.0045,
         'rings_for_gen': [0, 4, 5],
         'rings_for_scoring': [0, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13],
         'hkl_tols_seq': [0.01, 0.02, 0.03, 0.04, 0.05],
         'fracs': [0.9],
         'max_grains': 1000,
         'peak_assign_tol': 0.05,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # tomo_2_map.ipynb
         'dset_file': dset_file,
         'phase_str': 'Fe',
         'cf_strong_frac': 0.9939,
         'cf_strong_dstol': 0.005,
         'is_half_scan': False,
         'halfmask_radius': 25,
         'peak_assign_tol': 0.25,
         'manual_threshold': None,
         'hkltol': 0.25,
         'correct_sinos_with_ring_current': True,
         'first_tmap_cutoff_level': 0.4,
         'niter': 500,
         'second_tmap_cutoff_level': 0.05,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # tomo_2_map_minor_phase.ipynb
         'dset_file': dset_file,
         'major_phase_str': 'Fe',
         'minor_phase_str': 'Au',
         'major_phase_cf_frac': 0.994,
         'major_phase_cf_dstol': 0.005,
         'minor_phase_cf_frac': 0.9975,
         'minor_phase_cf_dstol': 0.005,
         'is_half_scan': False,
         'halfmask_radius': 25,
         'peak_assign_tol': 0.25,
         'hkltol': 0.25,
         'correct_sinos_with_ring_current': True,
         'first_tmap_cutoff_level': 0.4,
         'niter': 500,
         'second_tmap_cutoff_level': 0.5,
         'grain_too_many_px': 10,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # tomo_3_refinement.ipynb - major phase
         'dset_file': dset_file,
         'phase_str': 'Fe',
         'default_npks': 20,
         'default_nuniq': 20,
         'hkl_tol_origins': 0.05,
         'hkl_tol_refine': 0.1,
         'hkl_tol_refine_merged': 0.05,
         'ds_tol': 0.004,
         'ifrac': 7e-3,
         'rings_to_refine': None,
         'use_cluster': False,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # tomo_3_refinement.ipynb - minor phase
         'dset_file': dset_file,
         'phase_str': 'Au',
         'default_npks': 20,
         'default_nuniq': 20,
         'hkl_tol_origins': 0.05,
         'hkl_tol_refine': 0.1,
         'hkl_tol_refine_merged': 0.05,
         'ds_tol': 0.006,
         'ifrac': 1e-3,
         'rings_to_refine': [0, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13],
         'use_cluster': False,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # 4_visualise.ipynb - major phase
         'dset_file': dset_file,
         'phase_str': 'Fe',
         'min_unique': 250,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # 4_visualise.ipynb - minor phase
         'dset_file': dset_file,
         'phase_str': 'Au',
         'min_unique': 0,
         'dset_prefix': "top_",
        },
         {'PYTHONPATH': sys.path[0],  # 5_combine_phases.ipynb
         'dset_file': dset_file,
         'phase_strs': ['Fe', 'Au'],
         'combine_refined': True,
         'dset_prefix': "top_",
        },
         {'PYTHONPATH': sys.path[0],  # 6_stack_layers.ipynb
         'dset_file': dset_file,
         'stack_combined': True,
         'stack_refined': True,
         'zstep': 50.0,
         'dset_prefix': "top_",
        },
    ]
    if len(scan_nb_names) != len(scan_nb_params):
        raise ValueError('Mismatch between number of notebooks and param dicts!')
    scan_nb_paths = [os.path.join(scan_nb_prefix, name) for name in scan_nb_names]
    notebook_route(tomo_dir, scan_nb_paths, scan_nb_params)

    
def test_FeAu_JADB_pbp():
    tomo_dir = '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/PROCESSED_DATA/20250123_JADB/pbp_route'                    
    scan_nb_names = [
        '0_segment_and_label.ipynb',
        'pbp_1_indexing.ipynb',  # for major phase
        'pbp_1_indexing.ipynb',  # for minor phase
        'pbp_2_visualise.ipynb',  # for major phase
        'pbp_2_visualise.ipynb',  # for minor phase
        'pbp_3_refinement.ipynb',  # for major phase
        'pbp_3_refinement.ipynb',  # for minor phase
        '4_visualise.ipynb',  # for major phase
        '4_visualise.ipynb',  # for minor phase
        '5_combine_phases.ipynb',
        '6_stack_layers.ipynb'
        
    ]
    sample = 'FeAu_0p5_tR_nscope'
    dataset = 'top_200um'  # first of two layers
    dset_file = os.path.join(tomo_dir, sample, f'{sample}_{dataset}', f'{sample}_{dataset}_dataset.h5')
    scan_nb_params = [
        {'PYTHONPATH': sys.path[0],  # 0_segment_and_label.ipynb
         'maskfile': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/mask_with_gaps_E-08-0173.edf',
         'e2dxfile': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/e2dx_E-08-0173_20231127.edf',
         'e2dyfile': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/e2dy_E-08-0173_20231127.edf',
         'detector': 'eiger',
         'omegamotor': 'rot_center',
         'dtymotor': 'dty',
         'options': { 'cut' : 1, 'pixels_in_spot' : 3, 'howmany' : 100000 },
         'dataroot': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/RAW_DATA/',
         'analysisroot': tomo_dir,
         'sample': 'FeAu_0p5_tR_nscope',
         'dataset': 'top_200um',
         'dset_prefix': "top_"
        },
        {'PYTHONPATH': sys.path[0],  # pbp_1_indexing.ipynb - major phase
         'dset_file': dset_file,
         'par_file': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/pars.json',
         'phase_str': 'Fe',
         'minpkint': 5,
         'hkl_tol': 0.03,
         'fpks': 30,
         'ds_tol': 0.008,
         'etacut': 0.1,
         'ifrac': 2e-3,
         'y0': -16.0,
         'symmetry': 'cubic',
         'foridx': [0, 1, 3,  5, 7],
         'forgen': [1, 5, 7],
         'uniqcut': 0.85,
         'use_cluster': False,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # pbp_1_indexing.ipynb - minor phase
         'dset_file': dset_file,
         'par_file': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/pars.json',
         'phase_str': 'Au',
         'minpkint': 5,
         'hkl_tol': 0.03,
         'fpks': 30,
         'ds_tol': 0.008,
         'etacut': 0.1,
         'ifrac': 2e-3,
         'y0': -16.0,
         'symmetry': 'cubic',
         'foridx': [0, 1, 3,  5, 7],
         'forgen': [1, 5, 7],
         'uniqcut': 0.85,
         'use_cluster': False,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # pbp_2_visualise.ipynb - major phase
         'dset_file': dset_file,
         'phase_str': 'Fe',
         'min_unique': 20,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # pbp_2_visualise.ipynb - minor phase
         'dset_file': dset_file,
         'phase_str': 'Au',
         'min_unique': 10,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # pbp_3_refinement.ipynb - major phase
         'dset_file': dset_file,
         'phase_str': 'Fe',
         'min_unique': 20,
         'manual_threshold': None,
         'y0': -16.0,
         'hkl_tol_origins': 0.05,
         'hkl_tol_refine': 0.1,
         'hkl_tol_refine_merged': 0.05,
         'ds_tol': 0.004,
         'ifrac': 7e-3,
         'rings_to_refine': None,
         'set_mask_from_input': True,
         'use_cluster': False,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # pbp_3_refinement.ipynb - minor phase
         'dset_file': dset_file,
         'phase_str': 'Au',
         'min_unique': 10,
         'manual_threshold': None,
         'y0': -16.0,
         'hkl_tol_origins': 0.05,
         'hkl_tol_refine': 0.1,
         'hkl_tol_refine_merged': 0.05,
         'ds_tol': 0.004,
         'ifrac': 7e-3,
         'rings_to_refine': None,
         'set_mask_from_input': True,
         'use_cluster': False,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # 4_visualise.ipynb - major phase
         'dset_file': dset_file,
         'phase_str': 'Fe',
         'min_unique': 250,
         'dset_prefix': "top_",
        },
        {'PYTHONPATH': sys.path[0],  # 4_visualise.ipynb - minor phase
         'dset_file': dset_file,
         'phase_str': 'Au',
         'min_unique': 100,
         'dset_prefix': "top_",
        },
         {'PYTHONPATH': sys.path[0],  # 5_combine_phases.ipynb
         'dset_file': dset_file,
         'phase_strs': ['Fe', 'Au'],
         'combine_refined': True,
         'dset_prefix': "top_",
        },
         {'PYTHONPATH': sys.path[0],  # 6_stack_layers.ipynb
         'dset_file': dset_file,
         'stack_combined': True,
         'stack_refined': True,
         'zstep': 50.0,
         'dset_prefix': "top_",
        },
    ]
    if len(scan_nb_names) != len(scan_nb_params):
        raise ValueError('Mismatch between number of notebooks and param dicts!')
    scan_nb_paths = [os.path.join(scan_nb_prefix, name) for name in scan_nb_names]
    notebook_route(tomo_dir, scan_nb_paths, scan_nb_params)
    

if __name__=='__main__':
    print(papermill.__path__)
    test_FeAu_JADB_tomo()
    test_FeAu_JADB_pbp()