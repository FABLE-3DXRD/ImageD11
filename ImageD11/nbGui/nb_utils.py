import os
import subprocess
import time

import h5py
import numpy as np
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
#
date
echo python3 {python_script_path} {id11_code_path} {grainsfile} $SLURM_ARRAY_TASK_ID {dsfile} {reconfile} {cores_per_task} > {log_path} 2>&1
python3 {python_script_path} {id11_code_path} {grainsfile} $SLURM_ARRAY_TASK_ID {dsfile} {reconfile} {cores_per_task} > {log_path} 2>&1
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

    grains = [ImageD11.grain.grain(ubi) for ubi in indexer.ubis]
    print("Found {} grains".format(len(grains)))

    return grains, indexer


# GOTO refinegrains somewhere?
def refine_grain_positions(cf_3d, ds, grains, parfile, symmetry="cubic", cf_frac=0.85, cf_dstol=0.01,
                           hkl_tols=(0.05, 0.025, 0.01)):
    sample = ds.sample
    dataset = ds.dset
    cf_strong_allrings = select_ring_peaks_by_intensity(cf_3d, frac=cf_frac, dsmax=cf_3d.ds.max(), doplot=None, dstol=cf_dstol)
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
