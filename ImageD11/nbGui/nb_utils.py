import os
import subprocess
import time
from collections import OrderedDict

import h5py
import numpy as np
import fabio
from matplotlib import pyplot as plt
from tqdm.auto import tqdm
from scipy.optimize import curve_fit
from skimage.feature import blob_log

import ImageD11.cImageD11
import ImageD11.columnfile
import ImageD11.grain
import ImageD11.indexing
import ImageD11.refinegrains
import ImageD11.unitcell
import ImageD11.sinograms.roi_iradon
import ImageD11.sinograms.properties
from ImageD11.peakselect import select_ring_peaks_by_intensity


### General utilities (for all notebooks)


def is_notebook_executed(nb_path):
    import nbformat
    with open(nb_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)
    
    for cell in nb['cells']:
        if cell.cell_type == 'code' and 'execution_count' in cell and cell.execution_count is not None:
            return True  # At least one cell has been executed
    return False  # No executed cells found


def clear_notebook(nb_path):
    """Clears outputs of a Jupyter notebook."""
    import nbformat
    from nbconvert.preprocessors import ClearOutputPreprocessor
    with open(nb_path, "r", encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)

    # Use ClearOutputPreprocessor to remove outputs
    preprocessor = ClearOutputPreprocessor()
    preprocessor.preprocess(nb, {})

    # Save cleared notebook
    with open(nb_path, "w", encoding="utf-8") as f:
        nbformat.write(nb, f)


def notebook_prepare_pmill(nb_input_path, nb_output_path, params_dict, rename_colliding=False):
    """
    Prepare, but not execute, a notebook using papermill
    """
    import papermill
    if os.path.exists(nb_output_path):
        if rename_colliding:
            nb_output_path = nb_output_path.replace('.ipynb', '_2.ipynb')
        else:
            raise ValueError('Notebook already present', nb_output_path)
    papermill.execute_notebook(
       nb_input_path,
       nb_output_path,
       parameters=params_dict,
       prepare_only=True  # don't execute, just prepare
    )
    # clear outputs of the notebook
    clear_notebook(nb_output_path)
    return nb_output_path


def notebook_exec_pmill(nb_input_path, nb_output_path, params_dict, rename_colliding=False):
    import papermill
    # change output path if it already exists, in case we run the same notebook twice
    if os.path.exists(nb_output_path) and rename_colliding:
        nb_output_path = nb_output_path.replace('.ipynb', '_2.ipynb')
    print('Executing notebook', nb_output_path)
    papermill.execute_notebook(
       nb_input_path,
       nb_output_path,
       parameters=params_dict
    )


def prepare_notebooks_for_datasets(samples_dict, notebooks, dataroot, analysisroot, PYTHONPATH=None, notebook_parent_dir=None):
    """
    Prepare, but not execute, a series of notebooks for each dataset in samples_dict.
    Places the prepared notebooks for each dataset like PROCESSED_DATA/sample/sample_dataset/foo.ipynb
    Returns a list of absolute paths of notebooks to execute.
    
    samples_dict: dict of {sample1: [ds1, ds2, ds3], sample2: [ds1, ds2, ds3]} etc.
    notebooks: list of tuples of [(notebook_filename.ipynb, {params_for_notebook_1.ipynb})] etc. Param dicts should not contain dataroot, analysisroot, sample, dataset or dsfile information - those are intsead prepared by this function.
    dataroot: path to raw data folder
    analysisroot: path to root of analysis folder (usually PROCESSED_DATA)
    PYTHONPATH: Python path
    notebook_parent_dir: path to parent directory of input notebooks. Default: current working directory
    """
    if notebook_parent_dir is None:
        notebook_parent_dir = os.path.abspath('./')
    
    notebooks_to_execute = []
    for sample, datasets in samples_dict.items():
        for dataset in datasets:
            print("Preparing notebooks for " + sample + ":" + dataset)
            # Make a dataset so we know file paths
            ds = ImageD11.sinograms.dataset.DataSet(dataroot=dataroot,
                                                    analysisroot=analysisroot,
                                                    sample=sample,
                                                    dset=dataset)
            # if the analyispath doesn't exist, make it
            if not os.path.exists(ds.analysispath):
                os.makedirs(ds.analysispath)
            
            for (nb_name, nb_params) in notebooks:
                nb_in = os.path.join(notebook_parent_dir, nb_name)  # use the notebook from the current folder
                nb_out = os.path.join(ds.analysispath, nb_name)
                # prepare parameters for this notebook
                if PYTHONPATH is not None:
                    nb_params['PYTHONPATH'] = PYTHONPATH
                if nb_name.startswith('0'):
                    # the first notebook, segmentation, so we don't have a dataset name yet
                    nb_params['dataroot'] = ds.dataroot
                    nb_params['analysisroot'] = ds.analysisroot
                    nb_params['sample'] = sample
                    nb_params['dataset'] = dataset
                else:
                    # a later notebook, so all we need is the dataset path
                    nb_params['dset_path'] = ds.dsfile
                try:
                    nb_out = notebook_prepare_pmill(nb_in, nb_out, nb_params, rename_colliding=True)
                    notebooks_to_execute.append(nb_out)
                    print('Made notebook ' + nb_name + ' in ' + sample + ':' + dataset)
                except ValueError:  # we already found a notebook with this name in the folder
                    # has it been executed already? If yes, skip it
                    if is_notebook_executed(nb_out):
                        print('Already found executed notebook ' + nb_name + ' in ' + sample + ':' + dataset + ', skipping')
                        continue
                    else:
                        print('Found existing unexecuted notebook ' + nb_name + ' in ' + sample + ':' + dataset + ', will execute')
                        notebooks_to_execute.append(nb_out)
    
    return notebooks_to_execute


## Cluster related stuff (GOTO ImageD11.futures)

def slurm_submit_and_wait(bash_script_path, wait_time_sec=60):
    if not os.path.exists(bash_script_path):
        raise IOError("Bash script not found!")
    submit_command = "sbatch {}".format(bash_script_path)
    sbatch_submit_result = subprocess.run(submit_command, capture_output=True, shell=True).stdout.decode("utf-8")

    print(sbatch_submit_result.replace("\n", ""))

    slurm_job_number = None

    if sbatch_submit_result.startswith("Submitted"):
        slurm_job_number = sbatch_submit_result.replace("\n", "").split("job ")[1]

    print(slurm_job_number)

    assert slurm_job_number is not None

    slurm_job_finished = False

    while not slurm_job_finished:
        squeue_results = subprocess.run("squeue -u $USER", capture_output=True, shell=True).stdout.decode("utf-8")

        if slurm_job_number not in squeue_results:
            print("Slurm job finished!")
            slurm_job_finished = True
        else:
            print("Slurm job not finished! Waiting {} seconds...".format(wait_time_sec))
            time.sleep(wait_time_sec)


def slurm_submit_many_and_wait(bash_script_paths, wait_time_sec=60):
    for bash_script_path in bash_script_paths:
        if not os.path.exists(bash_script_path):
            raise IOError("Bash script not found!")

    slurm_job_numbers = []
    for bash_script_path in bash_script_paths:
        submit_command = "sbatch {}".format(bash_script_path)
        sbatch_submit_result = subprocess.run(submit_command, capture_output=True, shell=True).stdout.decode("utf-8")

        print(sbatch_submit_result.replace("\n", ""))

        slurm_job_number = None

        if sbatch_submit_result.startswith("Submitted"):
            slurm_job_number = sbatch_submit_result.replace("\n", "").split("job ")[1]

        # print(slurm_job_number)

        assert slurm_job_number is not None

        slurm_job_numbers.append(slurm_job_number)

    slurm_job_finished = False

    while not slurm_job_finished:
        squeue_results = subprocess.run("squeue -u $USER", capture_output=True, shell=True).stdout.decode("utf-8")

        jobs_still_running = False
        for slurm_job_number in slurm_job_numbers:
            if slurm_job_number in squeue_results:
                jobs_still_running = True

        if jobs_still_running:
            print("Slurm jobs not finished! Waiting {} seconds...".format(wait_time_sec))
            time.sleep(wait_time_sec)
        else:
            print("Slurm jobs all finished!")
            slurm_job_finished = True


def prepare_mlem_bash(ds, grains, id11_code_path, n_simultaneous_jobs=50, cores_per_task=8):
    slurm_mlem_path = os.path.join(ds.analysispath, "slurm_mlem")

    if os.path.exists(slurm_mlem_path):
        if len(os.listdir(slurm_mlem_path)) > 0:
            raise OSError("Slurm MLEM logs folder exists and is not empty!")
    else:
        os.mkdir(slurm_mlem_path)

    recons_path = os.path.join(ds.analysispath, "mlem_recons")

    if os.path.exists(recons_path):
        if len(os.listdir(recons_path)) > 0:
            raise OSError("MLEM recons folder exists and is not empty!")
    else:
        os.mkdir(recons_path)

    bash_script_path = os.path.join(slurm_mlem_path, ds.dsname + '_mlem_recon_slurm.sh')
    python_script_path = os.path.join(id11_code_path, "ImageD11/nbGui/S3DXRD/run_mlem_recon.py")
    outfile_path = os.path.join(slurm_mlem_path, ds.dsname + '_mlem_recon_slurm_%A_%a.out')
    errfile_path = os.path.join(slurm_mlem_path, ds.dsname + '_mlem_recon_slurm_%A_%a.err')
    log_path = os.path.join(slurm_mlem_path,
                            ds.dsname + '_mlem_recon_slurm_$SLURM_ARRAY_JOB_ID_$SLURM_ARRAY_TASK_ID.log')

    reconfile = os.path.join(recons_path, ds.dsname + "_mlem_recon_$SLURM_ARRAY_TASK_ID.txt")

    # python 3 version (de-indent whole below comment):

    #     bash_script_string = f"""#!/bin/bash
    # #SBATCH --job-name=mlem-recon
    # #SBATCH --output={outfile_path}
    # #SBATCH --error={errfile_path}
    # #SBATCH --array=0-{len(grains) - 1}%{n_simultaneous_jobs}
    # #SBATCH --time=02:00:00
    # # define memory needs and number of tasks for each array job
    # #SBATCH --ntasks=1
    # #SBATCH --cpus-per-task={cores_per_task}
    # #
    # date
    # echo python3 {python_script_path} {ds.grainsfile} $SLURM_ARRAY_TASK_ID {reconfile} {pad} {niter} {dohm} {mask_cen} > {log_path} 2>&1
    # python3 {python_script_path} {ds.grainsfile} $SLURM_ARRAY_TASK_ID {reconfile} {pad} {niter} {dohm} {mask_cen} > {log_path} 2>&1
    # date
    #     """

    # python 2 version
    bash_script_string = """#!/bin/bash
#SBATCH --job-name=mlem-recon
#SBATCH --output={outfile_path}
#SBATCH --error={errfile_path}
#SBATCH --array=0-{njobs}%{n_simultaneous_jobs}
#SBATCH --time=02:00:00
# define memory needs and number of tasks for each array job
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cores_per_task}
#SBATCH --mem-per-cpu=20G
#
date
echo PYTHONPATH={id11_code_path} python3 {python_script_path} {grainsfile} $SLURM_ARRAY_TASK_ID {dsfile} {reconfile} {cores_per_task} > {log_path} 2>&1
PYTHONPATH={id11_code_path} python3 {python_script_path} {grainsfile} $SLURM_ARRAY_TASK_ID {dsfile} {reconfile} {cores_per_task} > {log_path} 2>&1
date
    """.format(outfile_path=outfile_path,
               errfile_path=errfile_path,
               njobs=len(grains) - 1,
               n_simultaneous_jobs=n_simultaneous_jobs,
               cores_per_task=cores_per_task,
               python_script_path=python_script_path,
               id11_code_path=id11_code_path,
               grainsfile=ds.grainsfile,
               reconfile=reconfile,
               dsfile=ds.dsfile,
               log_path=log_path)

    with open(bash_script_path, "w") as bashscriptfile:
        bashscriptfile.writelines(bash_script_string)

    return bash_script_path, recons_path


def prepare_astra_bash(ds, grainsfile, id11_code_path, group_name='grains', memory=150):
    slurm_astra_path = os.path.join(ds.analysispath, "slurm_astra")

    if not os.path.exists(slurm_astra_path):
        os.mkdir(slurm_astra_path)

    bash_script_path = os.path.join(slurm_astra_path, ds.dsname + '_astra_recon_slurm.sh')
    python_script_path = os.path.join(id11_code_path, "ImageD11/nbGui/S3DXRD/run_astra_recon.py")
    outfile_path = os.path.join(slurm_astra_path, ds.dsname + '_astra_recon_slurm_%A.out')
    errfile_path = os.path.join(slurm_astra_path, ds.dsname + '_astra_recon_slurm_%A.err')
    log_path = os.path.join(slurm_astra_path,
                            ds.dsname + '_astra_recon_slurm_$SLURM_ARRAY_JOB_ID.log')

    # python 2 version
    bash_script_string = """#!/bin/bash
#SBATCH --job-name=astra-recon
#SBATCH --output={outfile_path}
#SBATCH --error={errfile_path}
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
# define memory needs and number of tasks for each array job
#SBATCH --ntasks=1
#SBATCH --mem={memory}G
date
module load cuda
echo PYTHONPATH={id11_code_path} python3 {python_script_path} {grainsfile} {dsfile} {group_name} > {log_path} 2>&1
PYTHONPATH={id11_code_path} python3 {python_script_path} {grainsfile} {dsfile} {group_name} > {log_path} 2>&1
date
    """.format(outfile_path=outfile_path,
               errfile_path=errfile_path,
               python_script_path=python_script_path,
               id11_code_path=id11_code_path,
               grainsfile=grainsfile,
               dsfile=ds.dsfile,
               group_name=group_name,
               memory=memory,
               log_path=log_path)

    with open(bash_script_path, "w") as bashscriptfile:
        bashscriptfile.writelines(bash_script_string)

    return bash_script_path


def prepare_pbp_bash(pbp_object, id11_code_path, minpkint):
    ds = pbp_object.dset

    slurm_pbp_path = os.path.join(ds.analysispath, "slurm_pbp")

    if not os.path.exists(slurm_pbp_path):
        os.mkdir(slurm_pbp_path)

    bash_script_path = os.path.join(slurm_pbp_path, ds.dsname + '_pbp_recon_slurm.sh')
    python_script_path = os.path.join(id11_code_path, "ImageD11/nbGui/S3DXRD/run_pbp_recon.py")
    outfile_path = os.path.join(slurm_pbp_path, ds.dsname + '_pbp_recon_slurm_%A.out')
    errfile_path = os.path.join(slurm_pbp_path, ds.dsname + '_pbp_recon_slurm_%A.err')
    log_path = os.path.join(slurm_pbp_path,
                            ds.dsname + '_pbp_recon_slurm_$SLURM_JOB_ID.log')

    # python 2 version
    bash_script_string = """#!/bin/bash
#SBATCH --job-name=pbp_scanning
#SBATCH --output={outfile_path}
#SBATCH --error={errfile_path}
#SBATCH --time=48:00:00
#SBATCH --partition=nice-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64

#
date
source /cvmfs/hpc.esrf.fr/software/packages/linux/x86_64/jupyter-slurm/latest/envs/jupyter-slurm/bin/activate
echo OMP_NUM_THREADS=1 PYTHONPATH={id11_code_path} python3 {python_script_path} {dsfile} {hkltol} {fpks} {dstol} {etacut} {ifrac} {costol} {y0} {symmetry} {foridx} {forgen} {uniqcut} {phase_name} {minpkint} > {log_path} 2>&1
OMP_NUM_THREADS=1 PYTHONPATH={id11_code_path} python3 {python_script_path} {dsfile} {hkltol} {fpks} {dstol} {etacut} {ifrac} {costol} {y0} {symmetry} {foridx} {forgen} {uniqcut} {phase_name} {minpkint} > {log_path} 2>&1
date
        """.format(outfile_path=outfile_path,
                   errfile_path=errfile_path,
                   python_script_path=python_script_path,
                   id11_code_path=id11_code_path,
                   dsfile=ds.dsfile,
                   hkltol=pbp_object.hkl_tol,
                   fpks=pbp_object.fpks,
                   dstol=pbp_object.ds_tol,
                   etacut=pbp_object.etacut,
                   ifrac=pbp_object.ifrac,
                   costol=pbp_object.cosine_tol,
                   y0=pbp_object.y0,
                   symmetry=pbp_object.symmetry,
                   foridx=str(pbp_object.foridx).replace(" ", ""),
                   forgen=str(pbp_object.forgen).replace(" ", ""),
                   uniqcut=pbp_object.uniqcut,
                   phase_name=str(pbp_object.phase_name),
                   minpkint=minpkint,
                   log_path=log_path)

    with open(bash_script_path, "w") as bashscriptfile:
        bashscriptfile.writelines(bash_script_string)

    return bash_script_path





## IO related stuff

# Helper funtions

# GOTO Silx.IO? is silx a dep already? might as well be if not
# or save grain map output
def save_array(grp, name, ary):
    cmp = {'compression': 'gzip',
           'compression_opts': 2,
           'shuffle': True}

    hds = grp.require_dataset(name,
                              shape=ary.shape,
                              dtype=ary.dtype,
                              **cmp)
    hds[:] = ary
    return hds


def find_datasets_to_process(rawdata_path, skips_dict, dset_prefix, sample_list):
    # check rawdata path exists, return empty dict otherwise
    if not os.path.exists(rawdata_path):
        return {}

    samples_dict = {}

    for sample in sample_list:
        sample_path = os.path.join(rawdata_path, sample)
        if os.path.exists(sample_path):
            all_dset_folders_for_sample = os.listdir(sample_path)
            dsets_list = []
            for folder in all_dset_folders_for_sample:
                if dset_prefix in folder:
                    dset_name = folder.split(sample + "_")[1]
                    if sample in skips_dict.keys():
                        if dset_name not in skips_dict[sample]:
                            dsets_list.append(dset_name)
                    else:
                        dsets_list.append(dset_name)

            samples_dict[sample] = sorted(dsets_list)

    return samples_dict

def save_ubi_map(ds, ubi_map, eps_map, misorientation_map, ipf_x_col_map, ipf_y_col_map, ipf_z_col_map):
    raise ValueError('This function is deprecated" Use ImageD11.sinograms.tensor_map.TensorMap instead')
    # with h5py.File(ds.pbpubifile, 'w') as hout:
    #     grp = hout.create_group('arrays')
    #     save_array(grp, 'ubi_map', ubi_map).attrs['description'] = 'Refined UBI values at each pixel'
    #     save_array(grp, 'eps_map', eps_map).attrs['description'] = 'Strain matrices (sample ref) at each pixel'
    #     save_array(grp, 'misorientation_map', misorientation_map).attrs[
    #         'description'] = 'Misorientation to grain avg at each pixel'
    #     ipfxdset = save_array(grp, 'ipf_x_col_map', ipf_x_col_map)
    #     ipfxdset.attrs['description'] = 'IPF X color at each pixel'
    #     ipfxdset.attrs['CLASS'] = 'IMAGE'
    #     ipfydset = save_array(grp, 'ipf_y_col_map', ipf_y_col_map)
    #     ipfydset.attrs['description'] = 'IPF Y color at each pixel'
    #     ipfydset.attrs['CLASS'] = 'IMAGE'
    #     ipfzdset = save_array(grp, 'ipf_z_col_map', ipf_z_col_map)
    #     ipfzdset.attrs['description'] = 'IPF Z color at each pixel'
    #     ipfzdset.attrs['CLASS'] = 'IMAGE'


### Sinogram stuff


# GOTO should be fixed by monitor in assemble_label
# should just be a sinogram numpy array and a monitor spectrum
# def correct_sinogram_rows_with_ring_current(grain, ds):
#     grain.ssino = grain.ssino / ds.ring_currents_per_scan_scaled[:, None]


### Peak manipulation

# GOTO class method for Peaks2Grain class
def assign_peaks_to_grains(grains, cf, tol):
    """Assigns peaks to the best fitting grain"""
    # assign peaks to grains

    # column to store the grain labels
    labels = np.zeros(cf.nrows, 'i')
    # get all g-vectors from columnfile (updateGeometry)
    # should we instead calculate considering grain translations? (probably!)
    gv = np.transpose((cf.gx, cf.gy, cf.gz)).astype(float)
    # column to store drlv2 (error in hkl)
    drlv2 = np.ones(cf.nrows, 'd')
    # iterate over all grains
    print("Scoring and assigning {} grains".format(len(grains)))
    for inc, g in enumerate(tqdm(grains)):
        n = ImageD11.cImageD11.score_and_assign(g.ubi, gv, tol, drlv2, labels, inc)

    # add the labels column to the columnfile
    cf.addcolumn(labels, 'grain_id')
    cf.addcolumn(drlv2, 'drlv2')


### Plotting

def plot_index_results(ind, colfile, title):
    # Generate a histogram of |drlv| for a ubi matrix
    ind.histogram_drlv_fit()
    # indexer.fight_over_peaks()

    fig, axs = plt.subplots(3, 2, layout="constrained", figsize=(9, 12))
    axs_flat = axs.ravel()

    # For each grain, plot the error in hkl vs the number of peaks with that error

    for grh in ind.histogram:
        axs_flat[0].plot(ind.bins[1:-1], grh[:-1], "-")

    axs_flat[0].set(ylabel="number of peaks",
                    xlabel="error in hkl (e.g. hkl versus integer)",
                    title=title)

    # set a mask of all non-assigned g-vectors

    # m = ind.ga == -1
    m = colfile.grain_id == -1

    # plot the assigned g-vectors omega vs dty (sinograms)

    axs_flat[1].scatter(colfile.omega[~m],
                        colfile.dty[~m],
                        c=colfile.grain_id[~m],
                        s=2,
                        cmap='tab20')

    axs_flat[1].set(title='Sinograms of {} grains'.format(colfile.grain_id.max() + 1),
                    xlabel=r'$\omega~(\degree)$',
                    ylabel='dty')

    # Define weak peaks as all non-assigned peaks with intensity 1e-4 of max
    if m.sum() > 0:
        cut = colfile.sum_intensity[m].max() * 1e-4
        weak = colfile.sum_intensity[m] < cut

        # Plot unassigned peaks in omega vs dty
        axs_flat[2].scatter(colfile.omega[m][weak], colfile.dty[m][weak], s=2, label='weak')
        axs_flat[2].scatter(colfile.omega[m][~weak], colfile.dty[m][~weak], s=2, label='not weak')

    axs_flat[2].set(title='Sinograms of unassigned peaks',
                    xlabel=r'$\omega~(\degree)$',
                    ylabel='dty')
    axs_flat[2].legend()

    # Plot d-star vs intensity for all assigned peaks

    axs_flat[3].scatter(colfile.ds[~m], colfile.sum_intensity[~m], s=2)
    axs_flat[3].set(title='Intensity of all assigned peaks',
                    xlabel=r'$d^{*}~(\AA^{-1})$',
                    ylabel='Intensity',
                    yscale='log')

    # Plot d-star vs intensity for all unassigned peaks
    if m.sum() > 0:
        axs_flat[4].scatter(colfile.ds[m][weak], colfile.sum_intensity[m][weak], s=2, label='weak')
        axs_flat[4].scatter(colfile.ds[m][~weak], colfile.sum_intensity[m][~weak], s=2, label='not weak')

    axs_flat[4].set(title='Intensity of all unassigned peaks',
                    xlabel=r'$d^{*}~(\AA^{-1})$',
                    ylabel='Intensity',
                    yscale='log')
    axs_flat[4].legend()

    # Get the number of peaks per grain

    npks = [(colfile.grain_id == i).sum() for i in range(len(ind.ubis))]

    # Plot histogram of number of peaks per grain

    axs_flat[5].hist(npks, bins=32)
    axs_flat[5].set(title='Hist of peaks per grain',
                    xlabel='Number of peaks',
                    ylabel='Number of grains')

    for ax in axs_flat:
        ax.set_box_aspect(0.7)

    plt.show()


def plot_grain_sinograms(grains, cf, n_grains_to_plot=None):
    if n_grains_to_plot is None:
        n_grains_to_plot = len(grains)

    grains_step = len(grains) // n_grains_to_plot

    grid_size = np.ceil(np.sqrt(len(grains[::grains_step]))).astype(int)
    nrows = (len(grains[::grains_step]) + grid_size - 1) // grid_size

    fig, axs = plt.subplots(grid_size, nrows, figsize=(10, 10), layout="constrained", sharex=True, sharey=True)
    if grid_size == 1 & nrows == 1:
        # only 1 grain
        g = grains[0]
        m = cf.grain_id == g.gid
        axs.scatter(cf.omega[m], cf.dty[m], c=cf.sum_intensity[m], s=2)
        axs.set_title(g.gid)
    else:
        for i, ax in enumerate(axs.ravel()):
            if i < len(grains[::grains_step]):
                # get corresponding grain for this axis
                g = grains[::grains_step][i]
                m = cf.grain_id == g.gid
                ax.scatter(cf.omega[m], cf.dty[m], c=cf.sum_intensity[m], s=2)
                ax.set_title('Grain ' + str(g.gid))

    fig.supxlabel(r'$\omega~(\degree)$')
    fig.supylabel("dty")

    plt.show()


def get_rgbs_for_grains(grains):
    # get the UB matrices for each grain
    UBs = np.array([g.UB for g in grains])

    # get the reference unit cell of one of the grains (should be the same for all)
    ref_ucell = grains[0].ref_unitcell

    # get a meta orientation for all the grains
    meta_ori = ref_ucell.get_orix_orien(UBs)

    rgb_x_all = ref_ucell.get_ipf_colour_from_orix_orien(meta_ori, axis=np.array([1., 0, 0]))
    rgb_y_all = ref_ucell.get_ipf_colour_from_orix_orien(meta_ori, axis=np.array([0., 1, 0]))
    rgb_z_all = ref_ucell.get_ipf_colour_from_orix_orien(meta_ori, axis=np.array([0., 0, 1]))

    for grain, rgb_x, rgb_y, rgb_z in zip(grains, rgb_x_all, rgb_y_all, rgb_z_all):
        grain.rgb_x = rgb_x
        grain.rgb_y = rgb_y
        grain.rgb_z = rgb_z


def plot_inverse_pole_figure_from_meta_orien(meta_orien, ref_ucell, axis=np.array([0., 0, 1]), **plot_kwargs):
    try:
        from orix.vector.vector3d import Vector3d
    except ImportError:
        raise ImportError("Missing diffpy and/or orix, can't compute orix phase!")

    ipf_direction = Vector3d(axis)

    # get the RGB colours
    rgb = ref_ucell.get_ipf_colour_from_orix_orien(meta_orien, axis=ipf_direction)

    # scatter the meta orientation using the colours
    meta_orien.scatter("ipf", c=rgb, direction=ipf_direction, **plot_kwargs)


def plot_all_ipfs_from_meta_orien(meta_orien, ref_ucell, **plot_kwargs):
    plot_inverse_pole_figure_from_meta_orien(meta_orien, ref_ucell, axis=np.array([1., 0., 0.]), **plot_kwargs)
    plot_inverse_pole_figure_from_meta_orien(meta_orien, ref_ucell, axis=np.array([0., 1., 0.]), **plot_kwargs)
    plot_inverse_pole_figure_from_meta_orien(meta_orien, ref_ucell, axis=np.array([0., 0., 1.]), **plot_kwargs)


def plot_inverse_pole_figure(grains, axis=np.array([0., 0, 1]), **plot_kwargs):
    # get the UB matrices for each grain
    UBs = np.array([g.UB for g in grains])

    # get the reference unit cell of one of the grains (should be the same for all)
    ref_ucell = grains[0].ref_unitcell

    # get a meta orientation for all the grains
    meta_orien = ref_ucell.get_orix_orien(UBs)

    plot_inverse_pole_figure_from_meta_orien(meta_orien, ref_ucell, axis, **plot_kwargs)


def plot_direct_pole_figure(grains, uvw=np.array([1., 0., 0.]), **plot_kwargs):
    # get the UB matrices for each grain
    UBs = np.array([g.UB for g in grains])

    # get the reference unit cell of one of the grains (should be the same for all)
    ref_ucell = grains[0].ref_unitcell

    # make a combined orientation from them (makes plot much faster)
    meta_orien = ref_ucell.get_orix_orien(UBs)

    try:
        from orix.vector import Miller
    except ImportError:
        raise ImportError("Missing orix, can't compute pole figure!")

    # make Miller object from uvw
    m1 = Miller(uvw=uvw, phase=ref_ucell.orix_phase).symmetrise(unique=True)

    # get outer product of all orientations with the crystal direction we're interested in
    uvw_all = (~meta_orien).outer(m1)

    uvw_all.scatter(hemisphere="both", axes_labels=["X", "Y"], **plot_kwargs)


def plot_all_ipfs(grains, **plot_kwargs):
    plot_inverse_pole_figure(grains, axis=np.array([1., 0, 0]), **plot_kwargs)
    plot_inverse_pole_figure(grains, axis=np.array([0., 1, 0]), **plot_kwargs)
    plot_inverse_pole_figure(grains, axis=np.array([0., 0, 1]), **plot_kwargs)



def plot_grain_positions(grains, colour='npks', centre_plot=False, size_scaling=0.5):
    """
    colour: choose from 'npks' or one of 'x', 'y', 'z' for IPF scaling
    centre_plot: choose whether to centre the plot horizontally (x and y)
    size_scaling: we only know relative grain sizes, adjust this to scale the diameter of the points on the plot
    """
    if colour.lower() not in ['npks', 'x', 'y', 'z']:
        raise ValueError("colour should be one of ['npks', 'x', 'y', 'z']")
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(projection='3d', proj_type="ortho")
    xx = [grain.translation[0] for grain in grains]
    yy = [grain.translation[1] for grain in grains]
    zz = [grain.translation[2] for grain in grains]
    if colour == 'npks':
        ax.set_title("Grain centre-of-mass positions coloured by number of peaks indexed")
        col = [float(grain.npks) for grain in grains]
    elif colour.lower() in ['x', 'y', 'z']:
        rgbattr = 'rgb_' + colour.lower()
        try:
            col = [getattr(grain, rgbattr) for grain in grains]  # IPF colour
        except AttributeError:
            # couldn't get the IPF attributes
            # try to compute it first
            # will still fail if we don't have reference unitcells
            get_rgbs_for_grains(grains)
            col = [getattr(grain, rgbattr) for grain in grains]  # IPF colour
        ax.set_title("Grain centre-of-mass positions coloured by IPF " +  colour.lower())
    # sizes in MPL 3D scale the area of the plot
    # intensity info is proportional to volume
    # decrease to radius then scale to area with power(x, 2/3)
    sizes = [size_scaling*np.power((float(grain.intensity_info.split("mean = ")[1].split(" , ")[0].replace("'", ""))), 2/3) for grain in grains]
    if centre_plot:
        scatterplot = ax.scatter(xx-np.mean(xx), yy-np.mean(yy), zz, c=col, s=sizes)
    else:
        scatterplot = ax.scatter(xx, yy, zz, c=col, s=sizes)
    if colour == 'npks':
        plt.colorbar(scatterplot)

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.show()
    

    
def plot_grain_histograms(fltfile, ubifile, parfile, OmSlop, OmFloat=True, nbins=30, tol=0.05):
    o=ImageD11.refinegrains.refinegrains(OmFloat=OmFloat, OmSlop=OmSlop, )

    o.loadparameters(parfile)
    o.readubis(ubifile)
    o.loadfiltered(fltfile)
    o.tolerance=tol
    o.generate_grains()
    o.assignlabels(quiet=True)

    # indexed peaks only
    d = o.scandata[fltfile]
    d.filter(d.labels >= 0)

    drlv_bins = np.linspace( 0, tol, nbins )
    ng = int(d.labels.max())+1
    drlv = np.sqrt( d.drlv2 )
    dp5 = [ drlv[d.labels==i] for i in range(ng)]
    hl = [ np.histogram(dpi, drlv_bins)[0] for dpi in dp5 ]

    if drlv_bins.shape[0] != hl[0].shape[0]:
        plotbins = (drlv_bins[1:] + drlv_bins[:-1])/2
    
    fig, axs = plt.subplots(2, 1, layout='constrained', sharex=True, figsize=(10, 7))
    for i in range(ng): 
        axs[0].plot(plotbins,hl[i],label=str(i))
    
    hist = axs[1].hist2d( drlv, d.labels, (drlv_bins, np.arange(-0.5,ng,1.)),vmin=0.5)
    fig.colorbar(hist[-1], ax=axs[1])
    #fig.supylabel("Grain")
    #fig.supxlabel("drlv")
    axs[0].set_ylabel('N peaks')
    axs[1].set_ylabel('Grain ID')
    fig.supxlabel('HKL Error')
    plt.show()


# backwards compatible
do_index = ImageD11.indexing.do_index


# GOTO refinegrains somewhere?
def refine_grain_positions(cf_3d, ds, grains, parfile, symmetry="cubic", cf_frac=0.85, cf_dstol=0.01,
                           hkl_tols=(0.05, 0.025, 0.01)):
    sample = ds.sample
    dataset = ds.dset
    cf_strong_allrings = select_ring_peaks_by_intensity(cf_3d, frac=cf_frac, dsmax=cf_3d.ds.max(), doplot=None,
                                                        dstol=cf_dstol)
    print("Got {} strong peaks for makemap".format(cf_strong_allrings.nrows))
    cf_strong_allrings_path = '{}_{}_3d_peaks_strong_all_rings.flt'.format(sample, dataset)
    cf_strong_allrings.writefile(cf_strong_allrings_path)

    tmp_ubi_path = '{}_{}_grains.ubi'.format(sample, dataset)
    tmp_map_path = '{}_{}_grains.map'.format(sample, dataset)

    new_flt_path = '{}_{}_3d_peaks_strong_all_rings.flt.new'.format(sample,
                                                                    dataset)  # flt file containing assignments from makemap
    unindexed_flt_path = '{}_{}_3d_peaks_strong_all_rings.flt.unindexed'.format(sample,
                                                                                dataset)  # remaining unassigned peaks from makemap

    ImageD11.grain.write_grain_file(tmp_ubi_path, grains)

    omegas_sorted = np.sort(ds.omega)[0]
    omega_slop = np.round(np.diff(omegas_sorted).mean(), 3)

    for inc, makemap_tol in enumerate(hkl_tols):
        print("Running makemap {}/{}".format(inc + 1, len(hkl_tols)))
        if inc == 0:  # ubi into map
            makemap_command = "makemap.py -p {} -u {} -U {} -f {} -F {} -s {} -t {} --omega_slop={} --no_sort".format(
                parfile, tmp_ubi_path, tmp_map_path, cf_strong_allrings_path, unindexed_flt_path, symmetry,
                hkl_tols[inc], omega_slop)
            makemap_output = subprocess.run(makemap_command, capture_output=True, shell=True).stdout.decode("utf-8")

            # makemap_output = !makemap.py -p {parfile} -u {tmp_ubi_path} -U {tmp_map_path} -f {cf_strong_allrings_path} -F {unindexed_flt_path} -s {symmetry} -t {hkl_tols[inc]} --omega_slop={omega_slop} --no_sort
        else:  # map into map
            makemap_command = "makemap.py -p {} -u {} -U {} -f {} -F {} -s {} -t {} --omega_slop={} --no_sort".format(
                parfile, tmp_map_path, tmp_map_path, cf_strong_allrings_path, unindexed_flt_path, symmetry,
                hkl_tols[inc], omega_slop)
            makemap_output = subprocess.run(makemap_command, capture_output=True, shell=True).stdout.decode("utf-8")
            # makemap_output = !makemap.py -p {parfile} -u {tmp_map_path} -U {tmp_map_path} -f {cf_strong_allrings_path} -F {unindexed_flt_path} -s {symmetry} -t {hkl_tols[inc]} --omega_slop={omega_slop} --no_sort

    grains2 = ImageD11.grain.read_grain_file(tmp_map_path)

    return grains2

def stereo( v, ix, iy, iz ):
    """
    Stereographic projection

    v = vector [3,n]
    ix = XX axis choice
    iy = YY axis choice
    iz = ZZ axis choice (the normal to the plot)

    Reflects to the positive hemisphere.
        e.g.: if v[iz] < 0 we use -v
    Normalises and plots the projection onto the plane to the point at 0,0,-1 (e.g. 1+z further away)
    """
    v = np.asarray( v )
    n = v / np.linalg.norm(v, axis=0 )
    X = n[ix] * np.sign( n[iz] ) / (1 + abs(n[iz]))
    Y = n[iy] * np.sign( n[iz] ) / (1 + abs(n[iz]))
    return X, Y


def plot_grains_polefig( hkls, cf, grains,
                        dstol=0.006, Imin=0, hbins=128, powder_max_pts=1e6,
                        cmapname='jet'
                       ):
    """
    Plot a plot figure of hkl and show where the grains are

    hkls = list of hkls to use. Might be ucell.ringhkls[ ucell.ringds[ j ] ]
    cf   = columnfile for plotting. Probably a cf_2d
    grains = list of grains for plotting
    dstol = which peaks from cf.ds to be plotted on pole figure
    Imin = cutoff to filter cf and remove noise
    pbins = binning for drawing the pole figure. Depends a bit on your screen resolution.
    cmapname = color scheme for your grains
    """
    # Pick one of the hkls for the plot titles:
    hklT = np.transpose( hkls )
    ilabel = np.argmax( np.sum( hklT, axis=0) )
    title = tuple(hklT[:,ilabel])
    # Generate the computed spot positions for all of our grains:
    gcalc = []
    for g in grains:
        gc = g.UB.dot( hklT )
        gcalc.append(gc)
    gcalc = np.concatenate(gcalc, axis=0)
    dsvals = np.linalg.norm( gcalc, axis=0 )
    # average and std of computed peaks
    m , s = dsvals.mean(), dsvals.std()
    if (dsvals.max()  > m + dstol) or (dsvals.min() < m - dstol ):
        import warnings
        warnings.warn('Your dstol is lower than the spread of your grains hkl peak positions')
    # The peaks for the pole figure
    pks = ( abs(cf.ds - m) < dstol ) & (cf.sum_intensity > Imin )
    # plot the observed data
    # ...and the grains
    # x,y   y,z   x,z
    gve = cf.gx[pks], cf.gy[pks], cf.gz[pks]
    I = cf.sum_intensity[ pks ]
    # powderskip :
    end = cf.nrows-1
    if powder_max_pts > 0:
        end = min( powder_max_pts, end )
    # Now the plotting:
    f, a = plt.subplots( 2,2, figsize=(8,6), constrained_layout=True)
    a  = a.ravel()
    a[0].plot( cf.ds[:end], cf.sum_intensity[:end], ',', label='powder')
    a[0].plot( cf.ds[pks], cf.sum_intensity[pks], '.', label='on pole figure')
    a[0].set( yscale='log', xlabel='dstar',)
    a[0].legend()
    ka = 1 # which axis to plot
    rng = -1.1, 1.1 # which range for gve
    cmap = plt.matplotlib.colormaps[ cmapname ]
    colors = cmap(np.linspace(0, 1, len(grains)))

    for ix,iy,iz in ( 0, 1, 2), (1,2,0), (0, 2, 1) : # the x/y/z axis choices
        xx, yy = stereo( gve, ix, iy, iz)
        a[ka].hist2d( xx, yy, weights=np.sqrt(I), bins=hbins, norm='log', cmap='gray_r')
        a[ka].set( xlabel='xyz'[ix], ylabel='xyz'[iy], aspect='equal', xlim=rng, ylim=rng,
                 title=str(title)+' : '+'xyz'[iz])
        for ig, g in enumerate(grains):
            xx, yy = stereo( g.UB.dot( hklT ), ix, iy, iz)
            a[ka].scatter( xx , yy, s=50, facecolors='none', edgecolors=colors[ig])
        ka += 1
    # Doesn't return the figure.... should it?

def plot_eta_vs_omega_error( cf, title=None ):
    """
    Diagnostic plot for a wedge error on a .flt.new columnfile coming out of makemap.py

    cf = columnfile .flt.new from makemap.py
    title = title for the plot
    """
    f, a = plt.subplots(1,1,constrained_layout=True)
    m = cf.tth_per_grain > 0
    a.plot( cf.eta_per_grain[m], (cf.omegacalc_per_grain -cf.omega)[m], ".")
    a.set(title=title, xlabel='eta', ylabel='omega error')


def plot_strain_errors( cf, grains, maskfile=None, wavelength=None ):
    """
    Diagnostic plot for a strain errors on a .flt.new columnfile coming out of makemap.py

    cf = columnfile
    grains = grains to plot the peaks (matching cf.labels)
    maskfile = for the f_raw, s_raw plot
    wavelength = needed for tth to ds if not in colf.parameters
    """
    f, ax = plt.subplots(2,3,constrained_layout=True, figsize=(10,6))
    ax = ax.ravel()
    if wavelength is None:
        wavelength = cf.parameters.get('wavelength')
    ds_per_grain = 2 * np.sin( np.radians( cf.tth_per_grain/2 ) ) / wavelength
    for i, g in enumerate( grains ):
        gcalc = g.UB.dot( (cf.h, cf.k, cf.l) )
        dscalc = np.linalg.norm( gcalc, axis=0 )
        m = cf.labels == i
        strain = ( ds_per_grain[m] - dscalc[m] ) / dscalc[m]
        for j, name in enumerate(( 'eta_per_grain', 'omega', 'tth_per_grain', 'Number_of_pixels', 'sum_intensity' )):
            ax[j].plot( cf[name][m], strain, '.')
            ax[j].set( xlabel=name, ylabel='strain')
        if maskfile is not None:
            ax[5].imshow( fabio.open( maskfile ).data, vmin=0, vmax=5, cmap='gray_r', origin='lower' )
        scat = ax[5].scatter( cf.f_raw[m], cf.s_raw[m], c=strain )
    ax[4].set( xscale='log' )  # intensity
    ax[5].set( aspect='equal')
    f.colorbar(scat, ax=ax[5])


### (hopefully) no longer used

# GOTO follow wherever we put scipy curve fit

# def calcy(cos_omega, sin_omega, sol):
#     return sol[0] + cos_omega * sol[1] + sin_omega * sol[2]


# def fity(y, cos_omega, sin_omega, wt=1):
#     """
#     Fit a sinogram to get a grain centroid
#     # calc = d0 + x*co + y*so
#     # dc/dpar : d0 = 1
#     #         :  x = co
#     #         :  y = so
#     # gradients
#     # General linear least squares
#     # Solution by the normal equation
#     # wt is weights (1/sig? or 1/sig^2?) 
#     # 
#     """
#     g = [wt * np.ones(y.shape, float), wt * cos_omega, wt * sin_omega]  # gradient
#     nv = len(g)
#     m = np.zeros((nv, nv), float)
#     r = np.zeros(nv, float)
#     for i in range(nv):
#         r[i] = np.dot(g[i], wt * y)  # A^T . b
#         for j in range(i, nv):
#             m[i, j] = np.dot(g[i], g[j])  # (A^T . A) . a = A^T . b
#             m[j, i] = m[i, j]
#     sol = np.dot(np.linalg.inv(m), r)
#     return sol


# def fity_robust(dty, co, so, nsigma=5, doplot=False):
#     cen, dx, dy = fity(dty, co, so)
#     calc2 = calc1 = calcy(co, so, (cen, dx, dy))
#     # mask for columnfile, we're selecting specific 4D peaks
#     # that come from the right place in y
#     selected = np.ones(co.shape, bool)
#     for i in range(3):
#         err = dty - calc2
#         estd = max(err[selected].std(), 1.0)  # 1 micron
#         # print(i,estd)
#         es = estd * nsigma
#         selected = abs(err) < es
#         cen, dx, dy = fity(dty, co, so, selected.astype(float))
#         calc2 = calcy(co, so, (cen, dx, dy))
#     # bad peaks are > 5 sigma
#     if doplot:
#         f, a = plt.subplots(1, 2)
#         theta = np.arctan2(so, co)
#         a[0].plot(theta, calc1, ',')
#         a[0].plot(theta, calc2, ',')
#         a[0].plot(theta[selected], dty[selected], "o")
#         a[0].plot(theta[~selected], dty[~selected], 'x')
#         a[1].plot(theta[selected], (calc2 - dty)[selected], 'o')
#         a[1].plot(theta[~selected], (calc2 - dty)[~selected], 'x')
#         a[1].set(ylim=(-es, es))
#         plt.show()
#     return selected, cen, dx, dy


# def graincen(gid, colf, doplot=True, nsigma=5):
#     # Get peaks beloging to this grain ID
#     m = colf.grain_id == gid
#     # Get omega values of peaks in radians
#     romega = np.radians(colf.omega[m])
#     # Calculate cos and sin of omega
#     co = np.cos(romega)
#     so = np.sin(romega)
#     # Get dty values of peaks
#     dty = colf.dty[m]
#     selected, cen, dx, dy = fity_robust(dty, co, so, nsigma=nsigma, doplot=doplot)
#     return selected, cen, dx, dy
