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

os.environ['PYDEVD_DISABLE_FILE_VALIDATION'] = '1'  # ignore papermill debugger warnings

import papermill

from ImageD11.nbGui.nb_utils import find_datasets_to_process, prepare_notebooks_for_datasets, notebook_exec_pmill

nb_base_prefix = os.path.join('..', 'ImageD11', 'nbGui')
scan_nb_prefix = os.path.join(nb_base_prefix, 'S3DXRD')
bb_nb_prefix = os.path.join(nb_base_prefix, 'TDXRD')


def notebook_route(base_dir, notebook_paths, notebook_param_dicts, notebook_out_dir=None, skip_dir_check=False):
    """
    Execute multiple notebooks in the order they are given
    base_dir: The path to the output folder for the test. Must not already exist.
    notebook_paths: Ordered list of paths to the notebooks, to be deployed one after the other
    notebook_param_dicts: Ordered list of dictionaries of parameters, one dict per notebook to be executed
    """
    if len(notebook_paths) != len(notebook_param_dicts):
        raise ValueError('Mismatch between number of notebooks and param dicts!')
    if os.path.exists(base_dir) and not skip_dir_check:
        raise ValueError('output test directory already exists:', base_dir)
    if not os.path.exists(base_dir):
        os.mkdir(base_dir)
    if notebook_out_dir is None:
        notebook_out_dir = os.path.join(base_dir, 'nb_out')
    if not os.path.exists(notebook_out_dir):
        os.mkdir(notebook_out_dir)
    for notebook_in_path, notebook_param_dict in zip(notebook_paths, notebook_param_dicts):
        notebook_out_path = os.path.join(notebook_out_dir,
                                         os.path.split(notebook_in_path)[1].replace('.ipynb', '_out.ipynb'))
        notebook_exec_pmill(notebook_in_path, notebook_out_path, notebook_param_dict, rename_colliding=True)


# there are two levels of testing
# does the notebook work without errors?
# does the notebook give you the output you expect?


# test the full tomographic route from start to finish
def test_tomographic_route():
    tomo_dir = 'tomo_route'
    dataroot = os.path.join(tomo_dir, 'raw')
    analysisroot = os.path.join(tomo_dir, 'processed')

    PYTHONPATH = sys.path[0]
    sample = 'Si_cube'
    dataset = 'S3DXRD_nt_moves_dty'
    samples_dict = {sample: [dataset]}
    # first, run the import_test_data.ipynb notebook to set up the file structure
    notebook_route(tomo_dir, [os.path.join(scan_nb_prefix, 'import_test_data.ipynb')],
                   [{'download_dir': tomo_dir, 'PYTHONPATH': PYTHONPATH}])

    nb_params = [
        ('tomo_1_index.ipynb',
         {'phase_str': 'Si',
          'min_frames_per_peak': 0,
          'cf_strong_frac': 0.9939,
          'cf_strong_dsmax': 1.594,
          'cf_strong_dstol': 0.005,
          'rings_for_gen': [0, 1, 3],
          'rings_for_scoring': [0, 1, 2, 3, 4],
          'hkl_tols_seq': [0.01, 0.02, 0.03, 0.04],
          'fracs': [0.9, 0.7],
          'max_grains': 1000,
          'peak_assign_tol': 0.025,
          }
         ),
        ('tomo_2_map.ipynb',
         {'phase_str': 'Si',
          'cf_strong_frac': 0.9939,
          'cf_strong_dstol': 0.005,
          'is_half_scan': False,
          'halfmask_radius': 25,
          'peak_assign_tol': 0.25,
          'draw_mask_interactive': False,
          'manual_threshold': None,
          'hkltol': 0.25,
          'correct_sinos_with_ring_current': False,
          'first_tmap_cutoff_level': 0.4,
          'niter': 500,
          'second_tmap_cutoff_level': 0.05
          }
         ),
        ('tomo_3_refinement.ipynb',
         {'phase_str': 'Si',
          'default_npks': 20,
          'default_nuniq': 20,
          'hkl_tol_origins': 0.05,
          'hkl_tol_refine': 0.1,
          'hkl_tol_refine_merged': 0.05,
          'ds_tol': 0.004,
          'ifrac': 7e-3,
          'rings_to_refine': None,
          'use_cluster': False
          }
         ),
        ('4_visualise.ipynb',
         {'phase_str': 'Si',
          'min_unique': 400,
          }
         )
    ]

    notebooks_to_execute = prepare_notebooks_for_datasets(samples_dict,
                                                          nb_params,
                                                          dataroot,
                                                          analysisroot,
                                                          PYTHONPATH=PYTHONPATH,
                                                          notebook_parent_dir=scan_nb_prefix)

    for nb_path in notebooks_to_execute:
        notebook_exec_pmill(nb_path, nb_path, None)


# test the full point-by-point route from start to finish
def test_pbp_route():
    tomo_dir = 'pbp_route'
    dataroot = os.path.join(tomo_dir, 'raw')
    analysisroot = os.path.join(tomo_dir, 'processed')

    PYTHONPATH = sys.path[0]
    sample = 'Si_cube'
    dataset = 'S3DXRD_nt_moves_dty'
    samples_dict = {sample: [dataset]}
    # first, run the import_test_data.ipynb notebook to set up the file structure
    notebook_route(tomo_dir, [os.path.join(scan_nb_prefix, 'import_test_data.ipynb')],
                   [{'download_dir': tomo_dir, 'PYTHONPATH': PYTHONPATH}])

    nb_params = [
        ('pbp_1_indexing.ipynb',
         {'phase_str': 'Si',
          'minpkint': 0,
          'hkl_tol': 0.025,
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
          }
         ),
        ('pbp_2_visualise.ipynb',
         {'phase_str': 'Si',
          'min_unique': 20
          }
         ),
        ('pbp_3_refinement.ipynb',
         {'phase_str': 'Si',
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
          }
         ),
        ('4_visualise.ipynb',
         {'phase_str': 'Si',
          'min_unique': 400,
          }
         )
    ]

    notebooks_to_execute = prepare_notebooks_for_datasets(samples_dict,
                                                          nb_params,
                                                          dataroot,
                                                          analysisroot,
                                                          PYTHONPATH=PYTHONPATH,
                                                          notebook_parent_dir=scan_nb_prefix)

    for nb_path in notebooks_to_execute:
        notebook_exec_pmill(nb_path, nb_path, None)


def test_FeAu_JADB_tomo():
    # where is the data?
    dataroot = '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/RAW_DATA'
    analysisroot = '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/PROCESSED_DATA/20250228_JADB/tomo_route'
    PYTHONPATH = sys.path[0]
    # find layers to process
    sample = 'FeAu_0p5_tR_nscope'
    first_dataset = 'top_200um'
    dset_prefix = "top"
    skips_dict = {sample: []}
    sample_list = [sample]
    samples_dict = find_datasets_to_process(dataroot, skips_dict, dset_prefix, sample_list)

    nb_params = [
        ('0_segment_and_label.ipynb',
         {'maskfile': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/mask_with_gaps_E-08-0173.edf',
          'e2dxfile': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/e2dx_E-08-0173_20231127.edf',
          'e2dyfile': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/e2dy_E-08-0173_20231127.edf',
          'detector': 'eiger',
          'omegamotor': 'rot_center',
          'dtymotor': 'dty',
          'options': {'cut': 1, 'pixels_in_spot': 3, 'howmany': 100000},
          'normalise_intensities_to_monitor': True,
          'monitor_name': 'fpico6'
          },
         ),
        ('tomo_1_index.ipynb',
         {'par_file': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/pars.json',
          'phase_str': 'Fe',
          'min_frames_per_peak': 0,
          'cf_strong_frac': 0.9939,
          'cf_strong_dsmax': 1.594,
          'cf_strong_dstol': 0.005,
          'rings_for_gen': [0, 1, 3],
          'rings_for_scoring': [0, 1, 2, 3, 4],
          'hkl_tols_seq': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.075],
          'fracs': [0.9, 0.7],
          'max_grains': 1000,
          'peak_assign_tol': 0.05,
          }
         ),
        ('tomo_1_index_minor_phase.ipynb',
         {'major_phase_strs': ['Fe'],
          'minor_phase_str': 'Au',
          'remove_major_phase_peaks': True,
          'min_frames_per_peak': 0,
          'major_phase_cf_dstol': 0.0035,
          'minor_phase_cf_frac': 0.9,
          'minor_phase_cf_dsmax': 1.594,
          'minor_phase_cf_dstol': 0.0045,
          'rings_for_gen': [0, 4, 5],
          'rings_for_scoring': [0, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13],
          'hkl_tols_seq': [0.01, 0.02, 0.03, 0.04, 0.05],
          'fracs': [0.9, 0.7],
          'max_grains': 1000,
          'peak_assign_tol': 0.05,
          }
         ),
        ('tomo_2_map.ipynb',
         {'phase_str': 'Fe',
          'cf_strong_frac': 0.9975,
          'cf_strong_dstol': 0.005,
          'is_half_scan': False,
          'halfmask_radius': 25,
          'peak_assign_tol': 0.05,
          'draw_mask_interactive': False,
          'manual_threshold': None,
          'hkltol': 0.25,
          'correct_sinos_with_ring_current': True,
          'first_tmap_cutoff_level': 0.4,
          'niter': 500,
          'second_tmap_cutoff_level': 0.05,
          }
         ),
        ('tomo_2_map_minor_phase.ipynb',
         {'major_phase_strs': ['Fe'],
          'minor_phase_str': 'Au',
          'remove_major_phase_peaks': True,
          'major_phase_cf_dstol': 0.005,
          'minor_phase_cf_frac': 0.9975,
          'minor_phase_cf_dstol': 0.005,
          'is_half_scan': False,
          'halfmask_radius': 25,
          'peak_assign_tol': 0.05,
          'hkltol': 0.25,
          'correct_sinos_with_ring_current': True,
          'first_tmap_cutoff_level': 0.4,
          'niter': 500,
          'second_tmap_cutoff_level': 0.5,
          'grain_too_many_px': 10,
          }
         ),
        ('tomo_3_refinement.ipynb',
         {'phase_str': 'Fe',
          'default_npks': 20,
          'default_nuniq': 20,
          'hkl_tol_origins': 0.05,
          'hkl_tol_refine': 0.1,
          'hkl_tol_refine_merged': 0.05,
          'ds_tol': 0.004,
          'ifrac': 7e-3,
          'rings_to_refine': None,
          'use_cluster': False,
          }
         ),
        ('tomo_3_refinement.ipynb',
         {'phase_str': 'Au',
          'default_npks': 20,
          'default_nuniq': 20,
          'hkl_tol_origins': 0.05,
          'hkl_tol_refine': 0.1,
          'hkl_tol_refine_merged': 0.05,
          'ds_tol': 0.006,
          'ifrac': 1e-3,
          'rings_to_refine': [0, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13],
          'use_cluster': False,
          }
         ),
        ('4_visualise.ipynb',
         {'phase_str': 'Fe',
          'min_unique': 250,
          }
         ),
        ('4_visualise.ipynb',
         {'phase_str': 'Au',
          'min_unique': 0,
          }
         ),
        ('5_combine_phases.ipynb',
         {'phase_strs': ['Fe', 'Au'],
          'combine_refined': True,
          }
         ),
    ]

    notebooks_to_execute = prepare_notebooks_for_datasets(samples_dict,
                                                          nb_params,
                                                          dataroot,
                                                          analysisroot,
                                                          PYTHONPATH=PYTHONPATH,
                                                          notebook_parent_dir=scan_nb_prefix)

    for nb_path in notebooks_to_execute:
        notebook_exec_pmill(nb_path, nb_path, None)

    # now run the final notebook to merge slices together
    # work out the path to the first dataset
    dset_path = os.path.join(analysisroot, sample, f'{sample}_{first_dataset}', f'{sample}_{first_dataset}_dataset.h5')
    nb_param = {'PYTHONPATH': PYTHONPATH,  # 7_stack_layers.ipynb
                'dset_path': dset_path,
                'dset_prefix': dset_prefix,
                'stack_combined': True,
                'stack_refined': True,
                'zstep': 50.0,
                }

    nb_path = os.path.join(scan_nb_prefix, '7_stack_layers.ipynb')
    notebook_route(analysisroot, [nb_path], [nb_param], skip_dir_check=True)


def test_FeAu_JADB_pbp():
    # where is the data?
    dataroot = '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/RAW_DATA'
    analysisroot = '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/PROCESSED_DATA/20250228_JADB/pbp_route'
    PYTHONPATH = sys.path[0]
    # find layers to process
    sample = 'FeAu_0p5_tR_nscope'
    first_dataset = 'top_200um'
    dset_prefix = "top"
    skips_dict = {sample: []}
    sample_list = [sample]
    samples_dict = find_datasets_to_process(dataroot, skips_dict, dset_prefix, sample_list)

    nb_params = [
        ('0_segment_and_label.ipynb',
         {'maskfile': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/mask_with_gaps_E-08-0173.edf',
          'e2dxfile': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/e2dx_E-08-0173_20231127.edf',
          'e2dyfile': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/e2dy_E-08-0173_20231127.edf',
          'detector': 'eiger',
          'omegamotor': 'rot_center',
          'dtymotor': 'dty',
          'options': {'cut': 1, 'pixels_in_spot': 3, 'howmany': 100000},
          'normalise_intensities_to_monitor': True,
          'monitor_name': 'fpico6'
          },
         ),
        ('pbp_1_indexing.ipynb',
         {'par_file': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/pars.json',
          'phase_str': 'Fe',
          'minpkint': 5,
          'hkl_tol': 0.03,
          'fpks': 30,
          'ds_tol': 0.008,
          'etacut': 0.1,
          'ifrac': 2e-3,
          'y0': -16.0,
          'symmetry': 'cubic',
          'foridx': [0, 1, 3, 5, 7],
          'forgen': [1, 5, 7],
          'uniqcut': 0.85,
          'use_cluster': False,
          }
         ),
        ('pbp_1_indexing.ipynb',
         {'par_file': '/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu/pars/pars.json',
          'phase_str': 'Au',
          'minpkint': 5,
          'hkl_tol': 0.03,
          'fpks': 30,
          'ds_tol': 0.008,
          'etacut': 0.1,
          'ifrac': 2e-3,
          'y0': -16.0,
          'symmetry': 'cubic',
          'foridx': [0, 3, 5, 7],
          'forgen': [0, 5, 7],
          'uniqcut': 0.85,
          'use_cluster': False,
          }
         ),
        ('pbp_2_visualise.ipynb',
         {'phase_str': 'Fe',
          'min_unique': 20,
          }
         ),
        ('pbp_2_visualise.ipynb',
         {'phase_str': 'Au',
          'min_unique': 10,
          }
         ),
        ('pbp_3_refinement.ipynb',
         {'phase_str': 'Fe',
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
          }
         ),
        ('pbp_3_refinement.ipynb',
         {'phase_str': 'Au',
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
          }
         ),
        ('4_visualise.ipynb',
         {'phase_str': 'Fe',
          'min_unique': 250,
          }
         ),
        ('4_visualise.ipynb',
         {'phase_str': 'Au',
          'min_unique': 100,
          }
         ),
        ('5_combine_phases.ipynb',
         {'phase_strs': ['Fe', 'Au'],
          'combine_refined': True,
          }
         ),
    ]

    notebooks_to_execute = prepare_notebooks_for_datasets(samples_dict,
                                                          nb_params,
                                                          dataroot,
                                                          analysisroot,
                                                          PYTHONPATH=PYTHONPATH,
                                                          notebook_parent_dir=scan_nb_prefix)

    for nb_path in notebooks_to_execute:
        notebook_exec_pmill(nb_path, nb_path, None)

    # now run the final notebook to merge slices together
    # work out the path to the first dataset
    dset_path = os.path.join(analysisroot, sample, f'{sample}_{first_dataset}', f'{sample}_{first_dataset}_dataset.h5')
    nb_param = {'PYTHONPATH': PYTHONPATH,  # 7_stack_layers.ipynb
                'dset_path': dset_path,
                'dset_prefix': dset_prefix,
                'stack_combined': True,
                'stack_refined': True,
                'zstep': 50.0,
                }

    nb_path = os.path.join(scan_nb_prefix, '7_stack_layers.ipynb')
    notebook_route(analysisroot, [nb_path], [nb_param], skip_dir_check=True)

    
def test_FeAu_f2scan_JADB_pbp():
    # where is the data?
    dataroot = "/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu_f2scan/RAW_DATA"
    analysisroot = "/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu_f2scan/PROCESSED_DATA/20250306_JADB"
    PYTHONPATH = sys.path[0]
    # find layers to process
    sample = "FeAu_No1_190um"
    first_dataset = "2um_redo_z_0"
    dset_prefix = "2um_redo_z"
    skips_dict = {sample: []}
    sample_list = [sample]
    samples_dict = find_datasets_to_process(dataroot, skips_dict, dset_prefix, sample_list)

    nb_params = [
        ('0_segment_and_label.ipynb',
         {'maskfile': "/data/id11/nanoscope/Eiger/eiger_mask_E-08-0144_20240205.edf",
          'e2dxfile': "/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu_f2scan/pars/e2dx_E-08-0144_20240205.edf",
          'e2dyfile': "/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu_f2scan/pars/e2dy_E-08-0144_20240205.edf",
          'detector': 'eiger',
          'omegamotor': 'owisRz_cen360',
          'dtymotor': 'diffty',
          'options': {'cut': 1, 'pixels_in_spot': 3, 'howmany': 100000},
          'normalise_intensities_to_monitor': True,
          'monitor_name': 'fpico3'
          },
         ),
        ('pbp_1_indexing.ipynb',
         {'par_file': "/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu_f2scan/pars/pars.json",
          'phase_str': 'Fe_bcc',
          'minpkint': 5,
          'hkl_tol': 0.05,
          'fpks': 0.9,
          'ds_tol': 0.004,
          'etacut': 0.1,
          'ifrac': 5e-3,
          'y0': 14.30168621868912,
          'symmetry': 'cubic',
          'foridx': [0, 1, 3, 5, 7],
          'forgen': [1, 5, 7],
          'uniqcut': 0.85,
          'use_cluster': False,
          }
         ),
        ('pbp_1_indexing.ipynb',
         {'par_file': "/data/id11/inhouse2/test_data_3DXRD/S3DXRD/FeAu_f2scan/pars/pars.json",
          'phase_str': 'Au_fcc',
          'minpkint': 5,
          'hkl_tol': 0.05,
          'fpks': 0.9,
          'ds_tol': 0.004,
          'etacut': 0.1,
          'ifrac': 5e-3,
          'y0': 14.30168621868912,
          'symmetry': 'cubic',
          'foridx': [0, 1, 2, 3, 4],
          'forgen': [0, 1, 4],
          'uniqcut': 0.85,
          'use_cluster': False,
          }
         ),
        ('pbp_2_visualise.ipynb',
         {'phase_str': 'Fe_bcc',
          'min_unique': 25,
          }
         ),
        ('pbp_2_visualise.ipynb',
         {'phase_str': 'Au_fcc',
          'min_unique': 25,
          }
         ),
        ('pbp_3_refinement.ipynb',
         {'phase_str': 'Fe_bcc',
          'min_unique': 22,
          'manual_threshold': None,
          'y0': 14.30168621868912,
          'hkl_tol_origins': 0.075,
          'hkl_tol_refine': 0.15,
          'hkl_tol_refine_merged': 0.075,
          'ds_tol': 0.006,
          'ifrac': 0,
          'rings_to_refine': None,
          'set_mask_from_input': True,
          'use_cluster': False,
          }
         ),
        ('pbp_3_refinement.ipynb',
         {'phase_str': 'Au_fcc',
          'min_unique': 22,
          'manual_threshold': None,
          'y0': 14.30168621868912,
          'hkl_tol_origins': 0.075,
          'hkl_tol_refine': 0.125,
          'hkl_tol_refine_merged': 0.075,
          'ds_tol': 0.006,
          'ifrac': 0,
          'rings_to_refine': [0, 2, 3, 4, 5, 6, 7, 8],
          'set_mask_from_input': True,
          'use_cluster': False,
          }
         ),
        ('4_visualise.ipynb',
         {'phase_str': 'Fe_bcc',
          'min_unique': 120,
          }
         ),
        ('4_visualise.ipynb',
         {'phase_str': 'Au_fcc',
          'min_unique': 120,
          }
         ),
        ('5_combine_phases.ipynb',
         {'phase_strs': ['Fe_bcc', 'Au_fcc'],
          'combine_refined': True,
          }
         ),
    ]

    notebooks_to_execute = prepare_notebooks_for_datasets(samples_dict,
                                                          nb_params,
                                                          dataroot,
                                                          analysisroot,
                                                          PYTHONPATH=PYTHONPATH,
                                                          notebook_parent_dir=scan_nb_prefix)

    for nb_path in notebooks_to_execute:
        notebook_exec_pmill(nb_path, nb_path, None)

    # now run the final notebook to merge slices together
    # work out the path to the first dataset
    dset_path = os.path.join(analysisroot, sample, f'{sample}_{first_dataset}', f'{sample}_{first_dataset}_dataset.h5')
    nb_param = {'PYTHONPATH': PYTHONPATH,  # 7_stack_layers.ipynb
                'dset_path': dset_path,
                'dset_prefix': dset_prefix,
                'stack_combined': True,
                'stack_refined': True,
                'zstep': 50.0,
                }

    nb_path = os.path.join(scan_nb_prefix, '7_stack_layers.ipynb')
    notebook_route(analysisroot, [nb_path], [nb_param], skip_dir_check=True)
    
    

def test_FeAu_JADB_bb():
    # where is the data?
    dataroot = '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/RAW_DATA/'
    analysisroot = '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/PROCESSED_DATA/20250304_JADB/default'
    PYTHONPATH = sys.path[0]
    # find layers to process
    sample = 'FeAu_0p5_tR'
    first_dataset = 'ff1'
    dset_prefix = "ff"
    skips_dict = {sample: []}
    sample_list = [sample]
    samples_dict = find_datasets_to_process(dataroot, skips_dict, dset_prefix, sample_list)

    # some extra stuff for notebook 0
    bgfile = None
    darkfile = None
    flatfile = None
    maskfile = '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/pars/mask.edf'

    nb_params = [
        ('0_segment_frelon.ipynb',
         {
             'splinefile': ['/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/pars/frelon36_spline_20240604_dx.edf',
                            '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/pars/frelon36_spline_20240604_dy.edf'],
             'bgfile': bgfile,
             'maskfile': maskfile,
             'darkfile': darkfile,
             'flatfile': flatfile,
             'detector': 'frelon3',
             'omegamotor': 'diffrz',
             'dtymotor': 'diffty',
             'options': {
                 "bgfile": bgfile,
                 "maskfile": maskfile,
                 "darkfile": darkfile,
                 "flatfile": flatfile,
                 "threshold": 70,
                 "smoothsigma": 1.0,
                 "bgc": 0.9,
                 "minpx": 3,
                 "m_offset_thresh": 100,
                 "m_ratio_thresh": 150,
             },
             'normalise_intensities_to_monitor': True,
             'monitor_name': 'fpico4'
         }  # end this dict
         ),  # end this tuple for this notebook
        ('1_index_default.ipynb',
         {
             'phase_str': 'Fe',
             'parfile': '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/pars/pars_tdxrd.json',
             'cf_strong_frac': 0.9837,
             'cf_strong_dsmax': 1.01,
             'cf_strong_dstol': 0.01,
             'rings_for_gen': [0, 1],
             'rings_for_scoring': [0, 1, 2, 3],
             'hkl_tols_seq': [0.01, 0.02, 0.03, 0.04],
             'fracs': [0.9, 9.75],
             'max_grains': 1000,
             'makemap_hkl_tol_seq': [0.05, 0.025, 0.01],
             'symmetry': 'cubic',
             'absolute_minpks': 120,
             'dset_prefix': "ff"
         }
         )
    ]

    notebooks_to_execute = prepare_notebooks_for_datasets(samples_dict,
                                                          nb_params,
                                                          dataroot,
                                                          analysisroot,
                                                          PYTHONPATH=PYTHONPATH,
                                                          notebook_parent_dir=bb_nb_prefix)

    for nb_path in notebooks_to_execute:
        notebook_exec_pmill(nb_path, nb_path, None)

    # now run the final notebook to merge slices together
    # work out the path to the first dataset
    dset_path = os.path.join(analysisroot, sample, f'{sample}_{first_dataset}', f'{sample}_{first_dataset}_dataset.h5')
    nb_param = {  # 3_merge_slices.ipynb
        'PYTHONPATH': PYTHONPATH,
        'dset_path': dset_path,
        'phase_str': 'Fe',
        'z_translation_motor': 'samtz',
        'dset_prefix': "ff"
    }

    nb_path = os.path.join(bb_nb_prefix, '3_merge_slices.ipynb')
    notebook_route(analysisroot, [nb_path], [nb_param], skip_dir_check=True)


def test_FeAu_JADB_bb_grid():
    # where is the data?
    dataroot = '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/RAW_DATA/'
    analysisroot = '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/PROCESSED_DATA/20250304_JADB/grid'
    PYTHONPATH = sys.path[0]
    # find layers to process
    sample = 'FeAu_0p5_tR'
    first_dataset = 'ff1'
    dset_prefix = "ff"
    skips_dict = {sample: []}
    sample_list = [sample]
    samples_dict = find_datasets_to_process(dataroot, skips_dict, dset_prefix, sample_list)

    # some extra stuff for notebook 0
    bgfile = None
    darkfile = None
    flatfile = None
    maskfile = '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/pars/mask.edf'

    nb_params = [
        ('0_segment_frelon.ipynb',
         {
             'splinefile': ['/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/pars/frelon36_spline_20240604_dx.edf',
                            '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/pars/frelon36_spline_20240604_dy.edf'],
             'bgfile': bgfile,
             'maskfile': maskfile,
             'darkfile': darkfile,
             'flatfile': flatfile,
             'detector': 'frelon3',
             'omegamotor': 'diffrz',
             'dtymotor': 'diffty',
             'options': {
                 "bgfile": bgfile,
                 "maskfile": maskfile,
                 "darkfile": darkfile,
                 "flatfile": flatfile,
                 "threshold": 70,
                 "smoothsigma": 1.0,
                 "bgc": 0.9,
                 "minpx": 3,
                 "m_offset_thresh": 100,
                 "m_ratio_thresh": 150,
             },
             'normalise_intensities_to_monitor': True,
             'monitor_name': 'fpico4'
         }  # end this dict
         ),  # end this tuple for this notebook
        ('1_index_grid.ipynb',
         {
             'phase_str': 'Fe',
             'parfile': '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/pars/pars_tdxrd.json',
             'cf_strong_frac': 0.9837,
             'cf_strong_dsmax': 1.01,
             'cf_strong_dstol': 0.01,
             'rings_to_use': [0, 1, 3],
             'symmetry': 'cubic',
             'makemap_tol_seq': [0.02, 0.015, 0.01],
             'gridpars': {
                 'DSTOL': 0.004,
                 'RING1': [1, 0, ],
                 'RING2': [0, ],
                 'NUL': True,
                 'FITPOS': True,
                 'tolangle': 0.50,
                 'toldist': 100.,
                 'NTHREAD': 1,
             },
             'grid_xlim': 600,
             'grid_ylim': 600,
             'grid_zlim': 200,
             'grid_step': 100,
             'frac': 0.85,
             'absolute_minpks': 56,
         }
         )
    ]

    notebooks_to_execute = prepare_notebooks_for_datasets(samples_dict,
                                                          nb_params,
                                                          dataroot,
                                                          analysisroot,
                                                          PYTHONPATH=PYTHONPATH,
                                                          notebook_parent_dir=bb_nb_prefix)

    for nb_path in notebooks_to_execute:
        notebook_exec_pmill(nb_path, nb_path, None)

    # now run the final notebook to merge slices together
    # work out the path to the first dataset
    dset_path = os.path.join(analysisroot, sample, f'{sample}_{first_dataset}', f'{sample}_{first_dataset}_dataset.h5')
    nb_param = {  # 3_merge_slices.ipynb
        'PYTHONPATH': PYTHONPATH,
        'dset_path': dset_path,
        'phase_str': 'Fe',
        'z_translation_motor': 'samtz',
        'dset_prefix': "ff"
    }

    nb_path = os.path.join(bb_nb_prefix, '3_merge_slices.ipynb')
    notebook_route(analysisroot, [nb_path], [nb_param], skip_dir_check=True)


def test_FeAu_JADB_bb_friedel():
    # where is the data?
    dataroot = '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/RAW_DATA/'
    analysisroot = '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/PROCESSED_DATA/20250304_JADB/friedel'
    PYTHONPATH = sys.path[0]
    # find layers to process
    sample = 'FeAu_0p5_tR'
    first_dataset = 'ff1'
    dset_prefix = "ff"
    skips_dict = {sample: []}
    sample_list = [sample]
    samples_dict = find_datasets_to_process(dataroot, skips_dict, dset_prefix, sample_list)

    # some extra stuff for notebook 0
    bgfile = None
    darkfile = None
    flatfile = None
    maskfile = '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/pars/mask.edf'

    nb_params = [
         ('0_segment_frelon.ipynb',
         {
             'splinefile': ['/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/pars/frelon36_spline_20240604_dx.edf',
                            '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/pars/frelon36_spline_20240604_dy.edf'],
             'bgfile': bgfile,
             'maskfile': maskfile,
             'darkfile': darkfile,
             'flatfile': flatfile,
             'detector': 'frelon3',
             'omegamotor': 'diffrz',
             'dtymotor': 'diffty',
             'options': {
                 "bgfile": bgfile,
                 "maskfile": maskfile,
                 "darkfile": darkfile,
                 "flatfile": flatfile,
                 "threshold": 70,
                 "smoothsigma": 1.0,
                 "bgc": 0.9,
                 "minpx": 3,
                 "m_offset_thresh": 100,
                 "m_ratio_thresh": 150,
             },
             'normalise_intensities_to_monitor': True,
             'monitor_name': 'fpico4'
         }  # end this dict
         ),  # end this tuple for this notebook
        ('1_index_friedel.ipynb',
         {'phase_str': 'Fe',
          'parfile': '/data/id11/inhouse2/test_data_3DXRD/TDXRD/FeAu/pars/pars_tdxrd.json',
          'cf_strong_frac': 0.991,
          'cf_strong_dsmax': 1.01,
          'cf_strong_dstol': 0.01,
          'womega': 1.0,
          'weta': 1.0,
          'wtth': 1.5,
          'wI': 0.5,
          'indexer_ds_tol': 0.003,
          'rings_for_gen': [1, 3],
          'rings_for_scoring': [0, 1, 2, 3],
          'hkl_tols_seq': [0.01, 0.02],
          'fracs': [0.9, 0.6],
          'max_grains': 1000,
          'symmetry': 'cubic',
          'gridpars': {
              'DSTOL': 0.004,
              'NUL': True,
              'FITPOS': True,
              'tolangle': 0.25,
              'toldist': 100.,
              'NTHREAD': 1,
              'NPKS': 25
          },
          'absolute_minpks': 25,
          },
         )
    ]

    notebooks_to_execute = prepare_notebooks_for_datasets(samples_dict,
                                                          nb_params,
                                                          dataroot,
                                                          analysisroot,
                                                          PYTHONPATH=PYTHONPATH,
                                                          notebook_parent_dir=bb_nb_prefix)

    for nb_path in notebooks_to_execute:
        notebook_exec_pmill(nb_path, nb_path, None)

    # now run the final notebook to merge slices together
    # work out the path to the first dataset
    dset_path = os.path.join(analysisroot, sample, f'{sample}_{first_dataset}', f'{sample}_{first_dataset}_dataset.h5')
    nb_param = {  # 3_merge_slices.ipynb
        'PYTHONPATH': PYTHONPATH,
        'dset_path': dset_path,
        'phase_str': 'Fe',
        'z_translation_motor': 'samtz',
        'dset_prefix': "ff"
    }

    nb_path = os.path.join(bb_nb_prefix, '3_merge_slices.ipynb')
    notebook_route(analysisroot, [nb_path], [nb_param], skip_dir_check=True)


if __name__ == '__main__':
    print(papermill.__path__)
    # test_tomographic_route()
    # test_pbp_route()
    # test_FeAu_JADB_tomo()
    # test_FeAu_JADB_pbp()
    test_FeAu_f2scan_JADB_pbp()
    test_FeAu_JADB_bb()
    # test_FeAu_JADB_bb_grid()
    # test_FeAu_JADB_bb_friedel()
