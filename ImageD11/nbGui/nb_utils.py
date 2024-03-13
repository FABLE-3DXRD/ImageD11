import os
import subprocess
import time

import h5py
import numba
import numpy as np
from matplotlib import pyplot as plt
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map
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

from ImageD11.blobcorrector import eiger_spatial


### General utilities (for all notebooks)

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


def prepare_mlem_bash(ds, grains, pad, is_half_scan, id11_code_path, n_simultaneous_jobs=50, cores_per_task=8,
                      niter=50):
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

    if is_half_scan:
        dohm = "Yes"
        mask_cen = "Yes"
    else:
        dohm = "No"
        mask_cen = "No"

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
#
date
echo python3 {python_script_path} {id11_code_path} {grainsfile} $SLURM_ARRAY_TASK_ID {reconfile} {pad} {niter} {dohm} {mask_cen} > {log_path} 2>&1
python3 {python_script_path} {id11_code_path} {grainsfile} $SLURM_ARRAY_TASK_ID {reconfile} {pad} {niter} {dohm} {mask_cen} > {log_path} 2>&1
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
               pad=pad,
               niter=niter,
               dohm=dohm,
               mask_cen=mask_cen,
               log_path=log_path)

    with open(bash_script_path, "w") as bashscriptfile:
        bashscriptfile.writelines(bash_script_string)

    return bash_script_path, recons_path


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
    samples_dict = {}

    for sample in sample_list:
        all_dset_folders_for_sample = os.listdir(os.path.join(rawdata_path, sample))
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


## Grain IO (needs its own class really)

# Box beam
# GOTO
# Make higher level HDF5/Nexus Generic IO class
# Check if Silx already has this
# Make dataset and GrainsIO both inherit from this
def save_3dxrd_grains(grains, ds):
    with h5py.File(ds.grainsfile, 'w') as hout:
        grn = hout.create_group('grains')
        for g in tqdm(grains):
            gg = grn.create_group(str(g.gid))
            save_array(gg, 'peaks_3d_indexing', g.peaks_3d).attrs[
                'description'] = "Strong 3D peaks that were assigned to this grain during indexing"
            gg.attrs.update({'ubi': g.ubi,
                             'translation': g.translation})


def read_3dxrd_grains(ds):
    with h5py.File(ds.grainsfile, 'r') as hin:
        grains_group = 'grains'

        grains = []
        for gid_string in tqdm(sorted(hin[grains_group].keys(), key=lambda x: int(x))):
            gg = hin[grains_group][gid_string]
            ubi = gg.attrs['ubi'][:]
            translation = gg.attrs['translation'][:]
            g = ImageD11.grain.grain(ubi, translation=translation)
            g.gid = int(gid_string)
            g.peaks_3d = gg['peaks_3d_indexing'][:]
            grains.append(g)

    return grains


# S3DXRD

def save_s3dxrd_grains_after_indexing(grains, ds):
    with h5py.File(ds.grainsfile, 'w') as hout:
        grn = hout.create_group('grains')
        for g in tqdm(grains):
            gg = grn.create_group(str(g.gid))
            save_array(gg, 'peaks_4d_indexing', g.peaks_4d).attrs[
                'description'] = "Strong 4D peaks that were assigned to this grain during indexing"
            gg.attrs.update({'ubi': g.ubi})


def save_s3dxrd_grains_minor_phase_after_indexing(grains, ds, phase_name=None):
    if not hasattr(ds, "grainsfile_minor_phase"):
        if phase_name is None:
            raise ValueError(
                "DataSet has no grainsfile_minor_phase attribute and you didn't provide a phase_name to generate one!")
        else:
            ds.grainsfile_minor_phase = os.path.join(ds.analysispath, ds.dsname + '_grains_' + phase_name + '.h5')
    with h5py.File(ds.grainsfile_minor_phase, 'w') as hout:
        grn = hout.create_group('grains')
        for g in tqdm(grains):
            gg = grn.create_group(str(g.gid))
            save_array(gg, 'peaks_4d_indexing', g.peaks_4d).attrs[
                'description'] = "Strong 4D peaks that were assigned to this grain during indexing"
            gg.attrs.update({'ubi': g.ubi})


def read_s3dxrd_grains_for_recon(ds):
    with h5py.File(ds.grainsfile, 'r') as hin:
        grains_group = 'grains'

        grains = []
        for gid_string in tqdm(sorted(hin[grains_group].keys(), key=lambda x: int(x))):
            gg = hin[grains_group][gid_string]
            ubi = gg.attrs['ubi'][:]
            g = ImageD11.grain.grain(ubi)
            g.gid = int(gid_string)
            grains.append(g)

    return grains


def read_s3dxrd_grains_minor_phase_for_recon(ds, phase_name=None):
    if not hasattr(ds, "grainsfile_minor_phase"):
        if phase_name is None:
            raise ValueError(
                "DataSet has no grainsfile_minor_phase attribute and you didn't provide a phase_name to generate one!")
        else:
            ds.grainsfile_minor_phase = os.path.join(ds.analysispath, ds.dsname + '_grains_' + phase_name + '.h5')
    with h5py.File(ds.grainsfile_minor_phase, 'r') as hin:
        grains_group = 'grains'

        grains = []
        for gid_string in tqdm(sorted(hin[grains_group].keys(), key=lambda x: int(x))):
            gg = hin[grains_group][gid_string]
            ubi = gg.attrs['ubi'][:]
            g = ImageD11.grain.grain(ubi)
            g.gid = int(gid_string)
            grains.append(g)

    return grains


def save_s3dxrd_grains_for_mlem(grains, ds, gord, inds, whole_sample_mask, y0):
    with h5py.File(ds.grainsfile, 'r+') as hout:
        try:
            grp = hout.create_group('peak_assignments')
        except ValueError:
            grp = hout['peak_assignments']

        ds_gord = save_array(grp, 'gord', gord)
        ds_gord.attrs['description'] = 'Grain ordering: g[i].pks = gord[ inds[i] : inds[i+1] ]'
        ds_inds = save_array(grp, 'inds', inds)
        ds_inds.attrs['description'] = 'Grain indices: g[i].pks = gord[ inds[i] : inds[i+1] ]'

        grains_group = 'grains'
        for g in tqdm(grains):
            gg = hout[grains_group][str(g.gid)]
            # save stuff for sinograms

            save_array(gg, 'ssino', g.ssino).attrs['description'] = 'Sinogram of peak intensities sorted by omega'
            save_array(gg, 'sinoangles', g.sinoangles).attrs['description'] = 'Projection angles for sinogram'
            save_array(gg, 'og_recon', g.og_recon).attrs['description'] = 'Original ID11 iRadon reconstruction'
            save_array(gg, 'circle_mask', whole_sample_mask).attrs[
                'description'] = 'Reconstruction mask to use for MLEM'

            # might as well save peaks stuff while we're here
            save_array(gg, 'translation', g.translation).attrs['description'] = 'Grain translation in lab frame'
            save_array(gg, 'peaks_2d_sinograms', g.peaks_2d).attrs[
                'description'] = "2D peaks from strong 4D peaks that were assigned to this grain for sinograms"
            save_array(gg, 'peaks_4d_sinograms', g.peaks_4d).attrs[
                'description'] = "Strong 4D peaks that were assigned to this grain for sinograms"

            gg.attrs['cen'] = g.cen
            gg.attrs['y0'] = y0


def save_s3dxrd_grains_after_recon(grains, ds, raw_intensity_array, grain_labels_array, rgb_x_array, rgb_y_array,
                                   rgb_z_array):
    with h5py.File(ds.grainsfile, 'r+') as hout:
        try:
            grp = hout.create_group('slice_recon')
        except ValueError:
            grp = hout['slice_recon']
        save_array(grp, 'intensity', raw_intensity_array).attrs['description'] = 'Raw intensity array for all grains'
        save_array(grp, 'labels', grain_labels_array).attrs['description'] = 'Grain labels array for all grains'

        ipfxdset = save_array(grp, 'ipf_x_col_map', rgb_x_array)
        ipfxdset.attrs['description'] = 'IPF X color at each pixel'
        ipfxdset.attrs['CLASS'] = 'IMAGE'
        ipfydset = save_array(grp, 'ipf_y_col_map', rgb_y_array)
        ipfydset.attrs['description'] = 'IPF Y color at each pixel'
        ipfydset.attrs['CLASS'] = 'IMAGE'
        ipfzdset = save_array(grp, 'ipf_z_col_map', rgb_z_array)
        ipfzdset.attrs['description'] = 'IPF Z color at each pixel'
        ipfzdset.attrs['CLASS'] = 'IMAGE'

        grains_group = 'grains'

        for g in tqdm(grains):
            gg = hout[grains_group][str(g.gid)]

            save_array(gg, 'recon', g.recon).attrs['description'] = 'Final reconstruction'


def save_s3dxrd_grains_minor_phase_after_recon(grains, ds, raw_intensity_array, grain_labels_array, rgb_x_array,
                                               rgb_y_array,
                                               rgb_z_array, phase_name=None):
    if not hasattr(ds, "grainsfile_minor_phase"):
        if phase_name is None:
            raise ValueError(
                "DataSet has no grainsfile_minor_phase attribute and you didn't provide a phase_name to generate one!")
        else:
            ds.grainsfile_minor_phase = os.path.join(ds.analysispath, ds.dsname + '_grains_' + phase_name + '.h5')

    # delete existing file, because our grain numbers have changed
    if os.path.exists(ds.grainsfile_minor_phase):
        os.remove(ds.grainsfile_minor_phase)

    with h5py.File(ds.grainsfile_minor_phase, 'w-') as hout:  # fail if exists
        #         try:
        #             grp = hout.create_group('peak_assignments')
        #         except ValueError:
        #             grp = hout['peak_assignments']

        #         ds_gord = save_array( grp, 'gord', gord )
        #         ds_gord.attrs['description'] = 'Grain ordering: g[i].pks = gord[ inds[i] : inds[i+1] ]'
        #         ds_inds = save_array( grp, 'inds', inds )
        #         ds_inds.attrs['description'] = 'Grain indices: g[i].pks = gord[ inds[i] : inds[i+1] ]'

        try:
            grp = hout.create_group('slice_recon')
        except ValueError:
            grp = hout['slice_recon']
        save_array(grp, 'intensity', raw_intensity_array).attrs['description'] = 'Raw intensity array for all grains'
        save_array(grp, 'labels', grain_labels_array).attrs['description'] = 'Grain labels array for all grains'

        ipfxdset = save_array(grp, 'ipf_x_col_map', rgb_x_array)
        ipfxdset.attrs['description'] = 'IPF X color at each pixel'
        ipfxdset.attrs['CLASS'] = 'IMAGE'
        ipfydset = save_array(grp, 'ipf_y_col_map', rgb_y_array)
        ipfydset.attrs['description'] = 'IPF Y color at each pixel'
        ipfydset.attrs['CLASS'] = 'IMAGE'
        ipfzdset = save_array(grp, 'ipf_z_col_map', rgb_z_array)
        ipfzdset.attrs['description'] = 'IPF Z color at each pixel'
        ipfzdset.attrs['CLASS'] = 'IMAGE'

        grains_group = hout.create_group('grains')
        for g in tqdm(grains):
            gg = grains_group.create_group(str(g.gid))
            # save stuff for sinograms

            gg.attrs.update({'ubi': g.ubi})

            save_array(gg, 'peaks_4d_indexing', g.peaks_4d).attrs[
                'description'] = "Strong 4D peaks that were assigned to this grain during indexing"

            save_array(gg, 'ssino', g.ssino).attrs['description'] = 'Sinogram of peak intensities sorted by omega'
            save_array(gg, 'sinoangles', g.sinoangles).attrs['description'] = 'Projection angles for sinogram'
            save_array(gg, 'og_recon', g.recon).attrs['description'] = 'Original ID11 iRadon reconstruction'
            save_array(gg, 'recon', g.recon).attrs['description'] = 'Final reconstruction'

            # might as well save peaks stuff while we're here
            save_array(gg, 'translation', g.translation).attrs['description'] = 'Grain translation in lab frame'
            save_array(gg, 'peaks_2d_sinograms', g.peaks_2d).attrs[
                'description'] = "2D peaks from strong 4D peaks that were assigned to this grain for sinograms"
            save_array(gg, 'peaks_4d_sinograms', g.peaks_4d).attrs[
                'description'] = "Strong 4D peaks that were assigned to this grain for sinograms"

            gg.attrs['cen'] = g.cen


def read_s3dxrd_grains_after_recon(ds):
    with h5py.File(ds.grainsfile, 'r') as hin:
        grp = hin['slice_recon']

        raw_intensity_array = grp['intensity'][:]
        grain_labels_array = grp['labels'][:]
        rgb_x_array = grp['ipf_x_col_map'][:]
        rgb_y_array = grp['ipf_y_col_map'][:]
        rgb_z_array = grp['ipf_z_col_map'][:]

        grains_group = 'grains'

        grains = []
        for gid_string in tqdm(sorted(hin[grains_group].keys(), key=lambda x: int(x))):
            gg = hin[grains_group][gid_string]
            ubi = gg.attrs['ubi'][:]

            g = ImageD11.grain.grain(ubi)

            # general grain properties

            g.gid = int(gid_string)
            g.translation = gg['translation'][:]
            g.cen = gg.attrs['cen']
            g.y0 = gg.attrs['y0'][()]
            g.sample_mask = gg['circle_mask'][:]

            # sinogram stuff
            g.ssino = gg['ssino'][:]
            g.sinoangles = gg['sinoangles'][:]

            # reconstructions
            g.og_recon = gg['og_recon'][:]
            g.recon = gg['recon'][:]

            grains.append(g)

    return grains, raw_intensity_array, grain_labels_array, rgb_x_array, rgb_y_array, rgb_z_array


def read_s3dxrd_grains_minor_phase_after_recon(ds, phase_name=None):
    if not hasattr(ds, "grainsfile_minor_phase"):
        if phase_name is None:
            raise ValueError(
                "DataSet has no grainsfile_minor_phase attribute and you didn't provide a phase_name to generate one!")
        else:
            ds.grainsfile_minor_phase = os.path.join(ds.analysispath, ds.dsname + '_grains_' + phase_name + '.h5')

    with h5py.File(ds.grainsfile_minor_phase, 'r') as hin:
        grp = hin['slice_recon']

        raw_intensity_array = grp['intensity'][:]
        grain_labels_array = grp['labels'][:]
        rgb_x_array = grp['ipf_x_col_map'][:]
        rgb_y_array = grp['ipf_y_col_map'][:]
        rgb_z_array = grp['ipf_z_col_map'][:]

        grains_group = 'grains'

        grains = []
        for gid_string in tqdm(sorted(hin[grains_group].keys(), key=lambda x: int(x))):
            gg = hin[grains_group][gid_string]
            ubi = gg.attrs['ubi'][:]

            g = ImageD11.grain.grain(ubi)

            # general grain properties

            g.gid = int(gid_string)
            g.translation = gg['translation'][:]
            g.cen = gg.attrs['cen']

            # sinogram stuff
            g.ssino = gg['ssino'][:]
            g.sinoangles = gg['sinoangles'][:]

            # reconstructions
            g.og_recon = gg['og_recon'][:]
            g.recon = gg['recon'][:]

            grains.append(g)

    return grains, raw_intensity_array, grain_labels_array, rgb_x_array, rgb_y_array, rgb_z_array


def save_ubi_map(ds, ubi_map, eps_map, misorientation_map, ipf_x_col_map, ipf_y_col_map, ipf_z_col_map):
    with h5py.File(ds.pbpubifile, 'w') as hout:
        grp = hout.create_group('arrays')
        save_array(grp, 'ubi_map', ubi_map).attrs['description'] = 'Refined UBI values at each pixel'
        save_array(grp, 'eps_map', eps_map).attrs['description'] = 'Strain matrices (sample ref) at each pixel'
        save_array(grp, 'misorientation_map', misorientation_map).attrs[
            'description'] = 'Misorientation to grain avg at each pixel'
        ipfxdset = save_array(grp, 'ipf_x_col_map', ipf_x_col_map)
        ipfxdset.attrs['description'] = 'IPF X color at each pixel'
        ipfxdset.attrs['CLASS'] = 'IMAGE'
        ipfydset = save_array(grp, 'ipf_y_col_map', ipf_y_col_map)
        ipfydset.attrs['description'] = 'IPF Y color at each pixel'
        ipfydset.attrs['CLASS'] = 'IMAGE'
        ipfzdset = save_array(grp, 'ipf_z_col_map', ipf_z_col_map)
        ipfzdset.attrs['description'] = 'IPF Z color at each pixel'
        ipfzdset.attrs['CLASS'] = 'IMAGE'


# Other

# GOTO Move to DataSet
def correct_half_scan(ds):
    """Pads the dataset to become bigger and symmetric"""
    # NOTE: We need to keep track of which bins are "real" and measured and which aren't
    # Sinomask
    c0 = 0
    # check / fix the centre of rotation
    # get value of bin closest to c0
    central_bin = np.argmin(abs(ds.ybincens - c0))
    # get centre dty value of this vin
    central_value = ds.ybincens[central_bin]

    lo_side = ds.ybincens[:central_bin + 1]
    hi_side = ds.ybincens[central_bin:]

    # get the hi/lo side which is widest
    # i.e if you go from -130 to +20, it selects -130
    yrange = max(hi_side[-1] - hi_side[0], lo_side[-1] - lo_side[0])

    # round to nearest multiple of ds.ystep
    yrange = np.ceil(yrange / ds.ystep) * ds.ystep

    # make new ymin and ymax that are symmetric around central_value
    ds.ymin = central_value - yrange
    ds.ymax = central_value + yrange

    new_yrange = ds.ymax - ds.ymin

    # determine new number of y bins
    ny = int(new_yrange // ds.ystep) + 1

    ds.ybincens = np.linspace(ds.ymin, ds.ymax, ny)
    ds.ybinedges = np.linspace(ds.ymin - ds.ystep / 2, ds.ymax + ds.ystep / 2, ny + 1)

    print(len(ds.ybincens))
    print(ds.ybincens)
    print(ds.ystep)
    print(yrange)
    print(ny)


### IPF Colour stuff
# GOTO New file, Orix interface, taking a grain instance (or a U) as an argument
# Check sym_u inside ImageD11

def grain_to_rgb(g, ax=(0, 0, 1)):
    return hkl_to_color_cubic(crystal_direction_cubic(g.ubi, ax))


def crystal_direction_cubic(ubi, axis):
    hkl = np.dot(ubi, axis)
    # cubic symmetry implies:
    #      24 permutations of h,k,l
    #      one has abs(h) <= abs(k) <= abs(l)
    hkl = abs(hkl)
    hkl.sort()
    return hkl


def hkl_to_color_cubic(hkl):
    """
    https://mathematica.stackexchange.com/questions/47492/how-to-create-an-inverse-pole-figure-color-map
        [x,y,z]=u⋅[0,0,1]+v⋅[0,1,1]+w⋅[1,1,1].
            These are:
                u=z−y, v=y−x, w=x
                This triple is used to assign each direction inside the standard triangle
                
    makeColor[{x_, y_, z_}] := 
         RGBColor @@ ({z - y, y - x, x}/Max@{z - y, y - x, x})                
    """
    x, y, z = hkl
    assert x <= y <= z
    assert z >= 0
    u, v, w = z - y, y - x, x
    m = max(u, v, w)
    r, g, b = u / m, v / m, w / m
    return (r, g, b)


def hkl_to_pf_cubic(hkl):
    x, y, z = hkl
    assert x <= y <= z
    assert z >= 0
    m = np.sqrt((hkl ** 2).sum())
    return x / (z + m), y / (z + m)


def get_rgbs_for_grains(grains):
    for grain in grains:
        grain.rgb_z = grain_to_rgb(grain, ax=(0, 0, 1), )  # symmetry = Symmetry.cubic)
        grain.rgb_y = grain_to_rgb(grain, ax=(0, 1, 0), )  # symmetry = Symmetry.cubic)
        grain.rgb_x = grain_to_rgb(grain, ax=(1, 0, 0), )  # symmetry = Symmetry.cubic)


### Segmentation

# GOTO Inside blobcorrector

def correct_pixel(pixel, spline_file):
    sr, fr = pixel
    sc, fc = ImageD11.blobcorrector.correctorclass(spline_file).correct(sr, fr)
    return sc, fc


def apply_spatial(cf, spline_file, workers):
    # sc = np.zeros(cf.nrows)
    # fc = np.zeros(cf.nrows)

    print("Spatial correction...")

    raw_pixels = np.vstack((cf['s_raw'], cf['f_raw'])).T

    corrected_pixels = process_map(correct_pixel, raw_pixels, [spline_file] * len(raw_pixels), max_workers=workers,
                                   chunksize=len(raw_pixels) // workers)

    sc, fc = [list(t) for t in zip(*corrected_pixels)]

    cf.addcolumn(sc, "sc")
    cf.addcolumn(fc, "fc")

    return cf


def apply_spatial_lut(cf, spline_file):
    # sc = np.zeros(cf.nrows)
    # fc = np.zeros(cf.nrows)

    print("Spatial correction...")

    corrector = ImageD11.blobcorrector.correctorclass(spline_file)
    corrector.correct_px_lut(cf)

    return cf


### Sinogram stuff

# GOTO sinograms/geometry

def sine_function(x, offset, a, b):
    return b * np.sin(np.radians(x)) + a * np.cos(np.radians(x)) + offset


def fit_sine_wave(x_data, y_data, initial_guess):
    # Fit the sine function to the data
    popt, _ = curve_fit(sine_function, x_data, y_data, p0=initial_guess, method='trf', loss='soft_l1', max_nfev=10000)

    offset, a, b = popt

    return offset, a, b


def fit_grain_position_from_sino(grain, cf_strong):
    initial_guess = (0, 0.5, 0.5)

    offset, a, b = fit_sine_wave(cf_strong.omega[grain.mask_4d], cf_strong.dty[grain.mask_4d], initial_guess)

    grain.cen = offset

    grain.dx = -b
    grain.dy = -a


def fit_grain_position_from_recon(grain, ds, y0):
    grain.bad_recon = False
    blobs = blob_log(grain.recon, min_sigma=1, max_sigma=10, num_sigma=10, threshold=.01)
    blobs_sorted = sorted(blobs, key=lambda x: x[2], reverse=True)
    try:
        largest_blob = blobs_sorted[0]

        # we now have the blob position in recon space
        # we need to go back to microns

        # first axis (vertical) is x
        # second axis (horizontal) is y

        x_recon_space = largest_blob[0]
        y_recon_space = largest_blob[1]

        # centre of the recon image is centre of space

        # the below should be independent, tested, inside sinograms/geometry
        x_microns = (x_recon_space - grain.recon.shape[0] // 2) * ds.ystep + y0
        y_microns = -(y_recon_space - grain.recon.shape[1] // 2) * ds.ystep + y0

        grain.x_blob = x_microns
        grain.y_blob = y_microns
    except IndexError:
        # didn't find any blobs
        # for small grains like these, if we didn't find a blob, normally indicates recon is bad
        # we will exclude it from maps and export
        grain.bad_recon = True


# should be fixed by monitor in assemble_label
# should just be a sinogram numpy array and a monitor spectrum
def correct_sinogram_rows_with_ring_current(grain, ds):
    grain.ssino = grain.ssino / ds.ring_currents_per_scan_scaled[:, None]


def get_ring_current_per_scan(ds):
    """Gets the ring current for each scan (i.e rotation/y-step)
       Stores it inside ds.ring_currents_per_scan and a scaled version inside ds.ring_currents_per_scan_scaled"""
    if not hasattr(ds, "ring_currents_per_scan"):
        ring_currents = []
        with h5py.File(ds.masterfile, "r") as h5in:
            for scan in ds.scans:
                ring_current = float(h5in[scan]["instrument/machine/current"][()])
                ring_currents.append(ring_current)

        ds.ring_currents_per_scan = np.array(ring_currents)
        ds.ring_currents_per_scan_scaled = np.array(ring_currents / np.max(ring_currents))


### Peak manipulation
# GOTO

@numba.njit(parallel=True)
def pmax(ary):
    """ Find the min/max of an array in parallel """
    mx = ary.flat[0]
    mn = ary.flat[0]
    for i in numba.prange(1, ary.size):
        mx = max(ary.flat[i], mx)
        mn = min(ary.flat[i], mn)
    return mn, mx


@numba.njit(parallel=True)
def palloc(shape, dtype):
    """ Allocate and fill an array with zeros in parallel """
    ary = np.empty(shape, dtype=dtype)
    for i in numba.prange(ary.size):
        ary.flat[i] = 0
    return ary


# counting sort by grain_id
@numba.njit
def counting_sort(ary, maxval=None, minval=None):
    """ Radix sort for integer array. Single threaded. O(n)
    Numpy should be doing this...
    """
    if maxval is None:
        assert minval is None
        minval, maxval = pmax(ary)  # find with a first pass
    maxval = int(maxval)
    minval = int(minval)
    histogram = palloc((maxval - minval + 1,), np.int64)
    indices = palloc((maxval - minval + 2,), np.int64)
    result = palloc(ary.shape, np.int64)
    for gid in ary:
        histogram[gid - minval] += 1
    indices[0] = 0
    for i in range(len(histogram)):
        indices[i + 1] = indices[i] + histogram[i]
    i = 0
    for gid in ary:
        j = gid - minval
        result[indices[j]] = i
        indices[j] += 1
        i += 1
    return result, histogram


@numba.njit(parallel=True)
def find_grain_id(spot3d_id, grain_id, spot2d_label, grain_label, order, nthreads=20):
    """
    Assignment grain labels into the peaks 2d array
    spot3d_id = the 3d spot labels that are merged and indexed
    grain_id = the grains assigned to the 3D merged peaks
    spot2d_label = the 3d label for each 2d peak
    grain_label => output, which grain is this peak
    order = the order to traverse spot2d_label sorted
    """
    assert spot3d_id.shape == grain_id.shape
    assert spot2d_label.shape == grain_label.shape
    assert spot2d_label.shape == order.shape
    T = nthreads
    print("Using", T, "threads")
    for tid in numba.prange(T):
        pcf = 0  # thread local I hope?
        for i in order[tid::T]:
            grain_label[i] = -1
            pkid = spot2d_label[i]
            while spot3d_id[pcf] < pkid:
                pcf += 1
            if spot3d_id[pcf] == pkid:
                grain_label[i] = grain_id[pcf]


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
    for g in tqdm(grains):
        n = ImageD11.cImageD11.score_and_assign(g.ubi, gv, tol, drlv2, labels, g.gid)

    # add the labels column to the columnfile
    cf.addcolumn(labels, 'grain_id')


def get_2d_peaks_from_4d_peaks(ds, cf):
    # Big scary block
    # Must understand what this does!

    # Ensure cf is sorted by spot3d_id
    # NOTE: spot3d_id should be spot4d_id, because we have merged into 4D?
    assert (np.argsort(cf.spot3d_id) == np.arange(cf.nrows)).all()

    # load the 2d peak labelling output
    pks = ImageD11.sinograms.properties.pks_table.load(ds.pksfile)

    # Grab the 2d peak centroids
    p2d = pks.pk2d(ds.omega, ds.dty)

    # NOTE: These are not spatially corrected?!

    numba_order, numba_histo = counting_sort(p2d['spot3d_id'])

    grain_2d_id = palloc(p2d['spot3d_id'].shape, np.dtype(int))

    cleanid = cf.grain_id.copy()

    find_grain_id(cf.spot3d_id, cleanid, p2d['spot3d_id'], grain_2d_id, numba_order)

    gord, counts = counting_sort(grain_2d_id)

    inds = np.concatenate(((0,), np.cumsum(counts)))

    # I think what we end up with is:
    # inds
    # this is an array which tells you which 2D spots each grain owns
    # the 2D spots are sorted by spot ID
    # inds tells you for each grain were you can find its associated 2D spots

    return gord, inds, p2d


# GOTO dataset method
# add optional peaks mask
# read only the bits of pkd that we need?
def tocolf(pkd, ds):
    """ Converts a dictionary of peaks into an ImageD11 columnfile
    adds on the geometric computations (tth, eta, gvector, etc) """
    spat = eiger_spatial(dxfile=ds.e2dxfile, dyfile=ds.e2dyfile)
    cf = ImageD11.columnfile.colfile_from_dict(spat(pkd))
    cf.parameters.loadparameters(ds.parfile)
    cf.updateGeometry()
    return cf


# hkl assignment needs grains2peaks class
# GOTO: part of a little collection of functions that put stuff on the Grain object
def do_sinos(g, p2d, ds, hkltol=0.25):
    # g.peaks_2d are 2d peaks that were merged into the 4d peaks that were assigned to this grain (obviously!)
    flt = tocolf({p: p2d[p][g.peaks_2d] for p in p2d}, ds)  # convert it to a columnfile and spatially correct

    hkl_real = np.dot(g.ubi, (flt.gx, flt.gy, flt.gz))  # calculate hkl of all assigned peaks
    hkl_int = np.round(hkl_real).astype(int)  # round to nearest integer
    dh = ((hkl_real - hkl_int) ** 2).sum(axis=0)  # calculate square of difference

    # g.dherrall = dh.mean()  # mean hkl error across all assigned peaks
    # g.npksall = flt.nrows  # total number of assigned peaks
    flt.filter(dh < hkltol * hkltol)  # filter all assigned peaks to be less than hkltol squared
    hkl_real = np.dot(g.ubi, (flt.gx, flt.gy, flt.gz))  # recalculate error after filtration
    hkl_int = np.round(hkl_real).astype(int)
    # dh = ((hkl_real - hkl_int) ** 2).sum(axis=0)
    # g.dherr = dh.mean()  # dherr is mean hkl error across assigned peaks after hkltol filtering
    # g.npks = flt.nrows  # total number of assigned peaks after hkltol filtering
    g.etasigns_2d_strong = np.sign(flt.eta)
    g.hkl_2d_strong = hkl_int  # integer hkl of assigned peaks after hkltol filtering
    g.sinoangles, g.ssino, g.hits = map_grain_from_peaks(g, flt, ds)
    return g


# GOTO GrainSinogram class - has grain object and dataset object as attributes
# has IO methods
def map_grain_from_peaks(g, flt, ds):
    """
    Computes sinogram
    flt is already the peaks for this grain
    Returns angles, sino
    """
    NY = len(ds.ybincens)  # number of y translations
    iy = np.round((flt.dty - ds.ybincens[0]) / (ds.ybincens[1] - ds.ybincens[0])).astype(
        int)  # flt column for y translation index

    # The problem is to assign each spot to a place in the sinogram
    hklmin = g.hkl_2d_strong.min(axis=1)  # Get minimum integer hkl (e.g -10, -9, -10)
    dh = g.hkl_2d_strong - hklmin[:, np.newaxis]  # subtract minimum hkl from all integer hkls
    de = (g.etasigns_2d_strong.astype(int) + 1) // 2  # something signs related
    #   4D array of h,k,l,+/-
    # pkmsk is whether a peak has been observed with this HKL or not
    pkmsk = np.zeros(list(dh.max(axis=1) + 1) + [2, ],
                     int)  # make zeros-array the size of (max dh +1) and add another axis of length 2
    pkmsk[dh[0], dh[1], dh[2], de] = 1  # we found these HKLs for this grain
    #   sinogram row to hit
    pkrow = np.cumsum(pkmsk.ravel()).reshape(pkmsk.shape) - 1  #
    # counting where we hit an HKL position with a found peak
    # e.g (-10, -9, -10) didn't get hit, but the next one did, so increment

    npks = pkmsk.sum()
    destRow = pkrow[dh[0], dh[1], dh[2], de]
    sino = np.zeros((npks, NY), 'f')
    hits = np.zeros((npks, NY), 'f')
    angs = np.zeros((npks, NY), 'f')
    adr = destRow * NY + iy
    # Just accumulate 
    sig = flt.sum_intensity
    ImageD11.cImageD11.put_incr64(sino, adr, sig)
    ImageD11.cImageD11.put_incr64(hits, adr, np.ones(len(de), dtype='f'))
    ImageD11.cImageD11.put_incr64(angs, adr, flt.omega)

    sinoangles = angs.sum(axis=1) / hits.sum(axis=1)
    # Normalise:
    sino = (sino.T / sino.max(axis=1)).T
    # Sort (cosmetic):
    order = np.lexsort((np.arange(npks), sinoangles))
    sinoangles = sinoangles[order]
    ssino = sino[order].T

    return sinoangles, ssino, hits[order].T


# GOTO apply_halfmask function to existing sinogram
# roi_iradon
def run_iradon_id11(sino, angles, pad=20, y0=0, workers=1, sample_mask=None, apply_halfmask=False,
                    mask_central_zingers=False):
    outsize = sino.shape[0] + pad

    if apply_halfmask:
        halfmask = np.zeros_like(sino)

        halfmask[:len(halfmask) // 2 - 1, :] = 1
        halfmask[len(halfmask) // 2 - 1, :] = 0.5

        sino_to_recon = sino * halfmask
    else:
        sino_to_recon = sino

    # # pad the sample mask
    # sample_mask_padded = np.pad(sample_mask, pad//2)

    # Perform iradon transform of grain sinogram, store result (reconstructed grain shape) in g.recon
    recon = ImageD11.sinograms.roi_iradon.iradon(sino_to_recon,
                                                 theta=angles,
                                                 mask=sample_mask,
                                                 output_size=outsize,
                                                 projection_shifts=np.full(sino.shape, -y0),
                                                 filter_name='hamming',
                                                 interpolation='linear',
                                                 workers=workers)

    if mask_central_zingers:
        grs = recon.shape[0]
        xpr, ypr = -grs // 2 + np.mgrid[:grs, :grs]
        inner_mask_radius = 25
        outer_mask_radius = inner_mask_radius + 2

        inner_circle_mask = (xpr ** 2 + ypr ** 2) < inner_mask_radius ** 2
        outer_circle_mask = (xpr ** 2 + ypr ** 2) < outer_mask_radius ** 2

        mask_ring = inner_circle_mask & outer_circle_mask
        # we now have a mask to apply
        fill_value = np.median(recon[mask_ring])
        recon[inner_circle_mask] = fill_value

    return recon


def iradon_grain(grain, pad=20, y0=0, workers=1, sample_mask=None, apply_halfmask=False, mask_central_zingers=False):
    sino = grain.ssino
    angles = grain.sinoangles
    recon = run_iradon_id11(sino, angles, pad, y0, workers, sample_mask, apply_halfmask, mask_central_zingers)
    grain.recon = recon

    return grain


# GOTO peakselect module
# holds many useful peak filtering functions
#
# merge sandbox/ringselect
def unitcell_peaks_mask(cf, dstol, dsmax):
    cell = ImageD11.unitcell.unitcell_from_parameters(cf.parameters)
    cell.makerings(dsmax)
    m = np.zeros(cf.nrows, bool)
    for v in cell.ringds:
        if v < dsmax:
            m |= (abs(cf.ds - v) < dstol)

    return m


# GOTO columnfile method
def strongest_peaks(colf, uself=True, frac=0.995, B=0.2, doplot=None):
    # correct intensities for structure factor (decreases with 2theta)
    cor_intensity = colf.sum_intensity * (np.exp(colf.ds * colf.ds * B))
    if uself:
        lf = ImageD11.refinegrains.lf(colf.tth, colf.eta)
        cor_intensity *= lf
    order = np.argsort(cor_intensity)[::-1]  # sort the peaks by intensity
    sortedpks = cor_intensity[order]
    cums = np.cumsum(sortedpks)
    cums /= cums[-1]
    # anything above this should be calculated one and added to the CF
    # check if the column exists before calculating
    # slider for frac?
    # all above only needs calculating once
    enough = np.searchsorted(cums, frac)
    # Aim is to select the strongest peaks for indexing.
    cutoff = sortedpks[enough]
    mask = cor_intensity > cutoff
    if doplot is not None:
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))
        axs[0].plot(cums / cums[-1], ',')
        axs[0].set(xlabel='npks', ylabel='fractional intensity')
        axs[0].plot([mask.sum(), ], [frac, ], "o")
        axs[1].plot(cums / cums[-1], ',')
        axs[1].set(xlabel='npks logscale', ylabel='fractional intensity', xscale='log', ylim=(doplot, 1.),
                   xlim=(np.searchsorted(cums, doplot), len(cums)))
        axs[1].plot([mask.sum(), ], [frac, ], "o")
        plt.show()
    return mask


def selectpeaks(cf, dstol=0.005, dsmax=10, frac=0.99, doplot=None):
    m = unitcell_peaks_mask(cf, dstol=dstol, dsmax=dsmax)
    cfc = cf.copy()
    cfc.filter(m)
    ms = strongest_peaks(cfc, frac=frac, doplot=doplot)
    cfc.filter(ms)
    return cfc


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

    m = ind.ga == -1
    # m = colfile.grain_id == -1

    # plot the assigned g-vectors omega vs dty (sinograms)

    axs_flat[1].scatter(colfile.omega[~m],
                        colfile.dty[~m],
                        c=ind.ga[~m],
                        s=2,
                        cmap='tab20')

    axs_flat[1].set(title='Sinograms of {} grains'.format(ind.ga.max() + 1),
                    xlabel='Omega/deg',
                    ylabel='dty/um')

    # Define weak peaks as all non-assigned peaks with intensity 1e-4 of max
    cut = colfile.sum_intensity[m].max() * 1e-4
    weak = colfile.sum_intensity[m] < cut

    # Plot unassigned peaks in omega vs dty

    axs_flat[2].scatter(colfile.omega[m][weak], colfile.dty[m][weak], s=2, label='weak')
    axs_flat[2].scatter(colfile.omega[m][~weak], colfile.dty[m][~weak], s=2, label='not weak')

    axs_flat[2].set(title='Sinograms of unassigned peaks',
                    xlabel='Omega/deg',
                    ylabel='dty/um')
    axs_flat[2].legend()

    # Plot d-star vs intensity for all assigned peaks

    axs_flat[3].scatter(colfile.ds[~m], colfile.sum_intensity[~m], s=2)
    axs_flat[3].set(title='Intensity of all assigned peaks',
                    xlabel='d-star',
                    ylabel='Intensity',
                    yscale='log')

    # Plot d-star vs intensity for all unassigned peaks

    axs_flat[4].scatter(colfile.ds[m][weak], colfile.sum_intensity[m][weak], s=2, label='weak')
    axs_flat[4].scatter(colfile.ds[m][~weak], colfile.sum_intensity[m][~weak], s=2, label='not weak')

    axs_flat[4].set(title='Intensity of all unassigned peaks',
                    xlabel='d-star',
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
    for i, ax in enumerate(axs.ravel()):
        if i < len(grains[::grains_step]):
            # get corresponding grain for this axis
            g = grains[::grains_step][i]
            m = cf.grain_id == g.gid
            ax.scatter(cf.omega[m], cf.dty[m], c=cf.sum_intensity[m], s=2)
            ax.set_title(g.gid)

    fig.supxlabel("Omega")
    fig.supylabel("Y translation (um)")

    plt.show()


# GOTO follow orix colouring stuff
def triangle():
    """ compute a series of point on the edge of the triangle """
    xy = [np.array(v) for v in ((0, 1, 1), (0, 0, 1), (1, 1, 1))]
    xy += [xy[2] * (1 - t) + xy[0] * t for t in np.linspace(0.1, 1, 5)]
    return np.array([hkl_to_pf_cubic(np.array(p)) for p in xy])


# GOTO follow orix colouring stuff
def plot_ipfs(grains):
    f, a = plt.subplots(1, 3, figsize=(15, 5))
    ty, tx = triangle().T
    for i, title in enumerate('xyz'):
        ax = np.zeros(3)
        ax[i] = 1.
        hkl = [crystal_direction_cubic(g.ubi, ax) for g in grains]
        xy = np.array([hkl_to_pf_cubic(h) for h in hkl])
        rgb = np.array([hkl_to_color_cubic(h) for h in hkl])
        for j in range(len(grains)):
            grains[j].rgb = rgb[j]
        a[i].scatter(xy[:, 1], xy[:, 0],
                     c=rgb)  # Note the "x" axis of the plot is the 'k' direction and 'y' is h (smaller)
        a[i].set(title=title, aspect='equal', facecolor='silver', xticks=[], yticks=[])
        a[i].plot(tx, ty, 'k-', lw=1)


# GOTO GrainMap class
def build_slice_arrays(grains, cutoff_level=0.0):
    grain_labels_array = np.zeros_like(grains[0].recon) - 1

    redx = np.zeros_like(grains[0].recon)
    grnx = np.zeros_like(grains[0].recon)
    blux = np.zeros_like(grains[0].recon)

    redy = np.zeros_like(grains[0].recon)
    grny = np.zeros_like(grains[0].recon)
    bluy = np.zeros_like(grains[0].recon)

    redz = np.zeros_like(grains[0].recon)
    grnz = np.zeros_like(grains[0].recon)
    bluz = np.zeros_like(grains[0].recon)

    raw_intensity_array = np.zeros_like(grains[0].recon)

    raw_intensity_array.fill(cutoff_level)

    def norm(r):
        m = r > r.max() * 0.2
        return (r / r[m].mean()).clip(0, 1)

    for g in tqdm(grains):
        i = g.gid

        g_raw_intensity = norm(g.recon)

        g_raw_intensity_mask = g_raw_intensity > raw_intensity_array

        g_raw_intensity_map = g_raw_intensity[g_raw_intensity_mask]

        raw_intensity_array[g_raw_intensity_mask] = g_raw_intensity_map

        redx[g_raw_intensity_mask] = g_raw_intensity_map * g.rgb_x[0]
        grnx[g_raw_intensity_mask] = g_raw_intensity_map * g.rgb_x[1]
        blux[g_raw_intensity_mask] = g_raw_intensity_map * g.rgb_x[2]

        redy[g_raw_intensity_mask] = g_raw_intensity_map * g.rgb_y[0]
        grny[g_raw_intensity_mask] = g_raw_intensity_map * g.rgb_y[1]
        bluy[g_raw_intensity_mask] = g_raw_intensity_map * g.rgb_y[2]

        redz[g_raw_intensity_mask] = g_raw_intensity_map * g.rgb_z[0]
        grnz[g_raw_intensity_mask] = g_raw_intensity_map * g.rgb_z[1]
        bluz[g_raw_intensity_mask] = g_raw_intensity_map * g.rgb_z[2]

        grain_labels_array[g_raw_intensity_mask] = i

    raw_intensity_array[raw_intensity_array == cutoff_level] = 0

    rgb_x_array = np.transpose((redx, grnx, blux), axes=(1, 2, 0))
    rgb_y_array = np.transpose((redy, grny, bluy), axes=(1, 2, 0))
    rgb_z_array = np.transpose((redz, grnz, bluz), axes=(1, 2, 0))

    return rgb_x_array, rgb_y_array, rgb_z_array, grain_labels_array, raw_intensity_array


### Indexing


# GOTO
def do_index(cf,
             dstol=0.05,
             hkl_tols=(0.01, 0.02, 0.03, 0.04, 0.05, 0.1),
             fracs=(0.9, 0.8, 0.7, 0.6, 0.5),
             cosine_tol=np.cos(np.radians(90 - 0.25)),
             max_grains=1000,
             forgen=(),
             foridx=()):
    print("Indexing {} peaks".format(cf.nrows))
    # replace Fe with something else
    Fe = ImageD11.unitcell.unitcell_from_parameters(cf.parameters)
    Fe.makerings(cf.ds.max())
    indexer = ImageD11.indexing.indexer_from_colfile(cf)

    ImageD11.indexing.loglevel = 3

    indexer.ds_tol = dstol
    indexer.assigntorings()
    indexer.max_grains = max_grains

    ImageD11.cImageD11.cimaged11_omp_set_num_threads(2)

    for ringid in forgen:
        if ringid not in foridx:
            raise ValueError("All rings in forgen must be in foridx!")

    n_peaks_expected = 0
    rings = []
    for i, dstar in enumerate(indexer.unitcell.ringds):
        multiplicity = len(indexer.unitcell.ringhkls[indexer.unitcell.ringds[i]])
        counts_on_this_ring = (indexer.ra == i).sum()
        # is this ring going to be used for indexing?
        if i in foridx:
            n_peaks_expected += multiplicity  # we expect peaks from this ring
            if i in forgen:  # we are generating orientations from this ring
                rings.append((counts_on_this_ring, multiplicity, i))

    rings.sort()

    print("{} peaks expected".format(n_peaks_expected))
    print("Trying these rings (counts, multiplicity, ring number): {}".format(rings))
    indexer.cosine_tol = np.abs(cosine_tol)

    for frac in fracs:
        for tol in hkl_tols:
            indexer.minpks = n_peaks_expected * frac
            indexer.hkl_tol = tol
            for i in range(len(rings)):
                for j in range(i, len(rings)):
                    indexer.ring_1 = rings[i][2]
                    indexer.ring_2 = rings[j][2]

                    indexer.find()
                    indexer.scorethem()

            print(frac, tol, len(indexer.ubis))

    grains = [ImageD11.grain.grain(ubi, translation=np.array([0., 0., 0.])) for ubi in indexer.ubis]
    print("Found {} grains".format(len(grains)))

    return grains, indexer


# GOTO refinegrains somewhere?
def refine_grain_positions(cf_3d, ds, grains, parfile, symmetry="cubic", cf_frac=0.85, cf_dstol=0.01,
                           hkl_tols=(0.05, 0.025, 0.01)):
    sample = ds.sample
    dataset = ds.dset
    cf_strong_allrings = selectpeaks(cf_3d, frac=cf_frac, dsmax=cf_3d.ds.max(), doplot=None, dstol=cf_dstol)
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
