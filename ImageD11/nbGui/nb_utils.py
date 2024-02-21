import os
import subprocess
import time

import numba
import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

import ImageD11.cImageD11
import ImageD11.columnfile
import ImageD11.grain
import ImageD11.indexing
import ImageD11.refinegrains
import ImageD11.unitcell

from ImageD11.blobcorrector import eiger_spatial


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


def triangle():
    """ compute a series of point on the edge of the triangle """
    xy = [np.array(v) for v in ((0, 1, 1), (0, 0, 1), (1, 1, 1))]
    xy += [xy[2] * (1 - t) + xy[0] * t for t in np.linspace(0.1, 1, 5)]
    return np.array([hkl_to_pf_cubic(np.array(p)) for p in xy])


def calcy(cos_omega, sin_omega, sol):
    return sol[0] + cos_omega * sol[1] + sin_omega * sol[2]


def fity(y, cos_omega, sin_omega, wt=1):
    """
    Fit a sinogram to get a grain centroid
    # calc = d0 + x*co + y*so
    # dc/dpar : d0 = 1
    #         :  x = co
    #         :  y = so
    # gradients
    # What method is being used here???????????
    """
    g = [wt * np.ones(y.shape, float), wt * cos_omega, wt * sin_omega]
    nv = len(g)
    m = np.zeros((nv, nv), float)
    r = np.zeros(nv, float)
    for i in range(nv):
        r[i] = np.dot(g[i], wt * y)
        for j in range(i, nv):
            m[i, j] = np.dot(g[i], g[j])
            m[j, i] = m[i, j]
    sol = np.dot(np.linalg.inv(m), r)
    return sol


def fity_robust(dty, co, so, nsigma=5, doplot=False):
    # NEEDS COMMENTING
    cen, dx, dy = fity(dty, co, so)
    calc2 = calc1 = calcy(co, so, (cen, dx, dy))
    # mask for columnfile, we're selecting specific 4D peaks
    # that come from the right place in y, I think?
    selected = np.ones(co.shape, bool)
    for i in range(3):
        err = dty - calc2
        estd = max(err[selected].std(), 1.0)  # 1 micron
        # print(i,estd)
        es = estd * nsigma
        selected = abs(err) < es
        cen, dx, dy = fity(dty, co, so, selected.astype(float))
        calc2 = calcy(co, so, (cen, dx, dy))
    # bad peaks are > 5 sigma
    if doplot:
        f, a = plt.subplots(1, 2)
        theta = np.arctan2(so, co)
        a[0].plot(theta, calc1, ',')
        a[0].plot(theta, calc2, ',')
        a[0].plot(theta[selected], dty[selected], "o")
        a[0].plot(theta[~selected], dty[~selected], 'x')
        a[1].plot(theta[selected], (calc2 - dty)[selected], 'o')
        a[1].plot(theta[~selected], (calc2 - dty)[~selected], 'x')
        a[1].set(ylim=(-es, es))
        plt.show()
    return selected, cen, dx, dy


def graincen(gid, colf, doplot=True):
    # Get peaks beloging to this grain ID
    m = colf.grain_id == gid
    # Get omega values of peaks in radians
    romega = np.radians(colf.omega[m])
    # Calculate cos and sin of omega
    co = np.cos(romega)
    so = np.sin(romega)
    # Get dty values of peaks
    dty = colf.dty[m]
    selected, cen, dx, dy = fity_robust(dty, co, so, doplot=doplot)
    return selected, cen, dx, dy


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


def tocolf(pkd, parfile, dxfile, dyfile):
    """ Converts a dictionary of peaks into an ImageD11 columnfile
    adds on the geometric computations (tth, eta, gvector, etc) """
    spat = eiger_spatial(dxfile=dxfile, dyfile=dyfile)
    cf = ImageD11.columnfile.colfile_from_dict(spat(pkd))
    cf.parameters.loadparameters(parfile)
    cf.updateGeometry()
    return cf


def unitcell_peaks_mask(cf, dstol, dsmax):
    cell = ImageD11.unitcell.unitcell_from_parameters(cf.parameters)
    cell.makerings(dsmax)
    m = np.zeros(cf.nrows, bool)
    for v in cell.ringds:
        if v < dsmax:
            m |= (abs(cf.ds - v) < dstol)

    return m


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


def assign_peaks_to_grains(grains, cf, tol):
    # assign peaks to grains

    # column to store the grain labels
    labels = np.zeros(cf.nrows, 'i')
    # get all g-vectors from columnfile
    gv = np.transpose((cf.gx, cf.gy, cf.gz)).astype(float)
    # column to store drlv2 (error in hkl)
    drlv2 = np.ones(cf.nrows, 'd')
    # iterate over all grains
    print("Scoring and assigning {} grains".format(len(grains)))
    for g in tqdm(grains):
        n = ImageD11.cImageD11.score_and_assign(g.ubi, gv, tol, drlv2, labels, g.gid)

    # add the labels column to the columnfile
    cf.addcolumn(labels, 'grain_id')


def do_index(cf,
             dstol=0.05,
             max_mult=13,
             min_ring_count=0,
             hkl_tols=(0.01, 0.02, 0.03, 0.04, 0.05, 0.1),
             fracs=(0.9, 0.8, 0.7, 0.6, 0.5),
             cosine_tol=np.cos(np.radians(90.25)),
             max_grains=1000):
    print("Indexing {} peaks".format(cf.nrows))
    Fe = ImageD11.unitcell.unitcell_from_parameters(cf.parameters)
    Fe.makerings(cf.ds.max())
    indexer = ImageD11.indexing.indexer_from_colfile(cf)

    ImageD11.indexing.loglevel = 3

    indexer.ds_tol = dstol
    indexer.assigntorings()
    indexer.max_grains = max_grains

    n_peaks_expected = 0
    rings = []
    for i, dstar in enumerate(indexer.unitcell.ringds):
        multiplicity = len(indexer.unitcell.ringhkls[indexer.unitcell.ringds[i]])
        counts_on_this_ring = (indexer.ra == i).sum()
        if counts_on_this_ring > min_ring_count:
            n_peaks_expected += multiplicity
            if multiplicity < max_mult:
                rings.append((counts_on_this_ring, multiplicity, i))

    rings.sort()

    print("{} peaks expected".format(n_peaks_expected))
    print("Trying these rings (counts, multiplicity, ring number): {}".format(rings))
    indexer.cosine_tol = cosine_tol

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


def build_slice_arrays(grains, cutoff_level=0.0):
    grain_labels_array = np.zeros_like(grains[0].recon) - 1
    red = np.zeros_like(grains[0].recon)
    grn = np.zeros_like(grains[0].recon)
    blu = np.zeros_like(grains[0].recon)

    raw_intensity_array = np.zeros_like(grains[0].recon)

    raw_intensity_array.fill(cutoff_level)

    def norm(r):
        m = r > r.max() * 0.2
        return (r / r[m].mean()).clip(0, 1)

    for g in tqdm(grains):
        i = g.gid

        g_raw_intensity = norm(g.recon)

        g_raw_intensity_mask = g_raw_intensity > raw_intensity_array

        g_raw_intenstiy_map = g_raw_intensity[g_raw_intensity_mask]

        raw_intensity_array[g_raw_intensity_mask] = g_raw_intenstiy_map

        red[g_raw_intensity_mask] = g_raw_intenstiy_map * g.rgb_z[0]
        grn[g_raw_intensity_mask] = g_raw_intenstiy_map * g.rgb_z[1]
        blu[g_raw_intensity_mask] = g_raw_intenstiy_map * g.rgb_z[2]

        grain_labels_array[g_raw_intensity_mask] = i

    raw_intensity_array[raw_intensity_array == cutoff_level] = 0

    rgb_array = np.transpose((red, grn, blu), axes=(1, 2, 0))

    return rgb_array, grain_labels_array, raw_intensity_array


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


def correct_half_scan(ds):
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
