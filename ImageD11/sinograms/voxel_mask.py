"""Written by Axel Henningsson, 2026.

Updated and slightly mangled by James Ball, 2026.

A voxel at (sx,sy) is illuminated when:
| y0 - sx*sin(omega) - sy*cos(omega) - dty | <= ystep
for one voxel that traces a sinusoid in (omega, dty) space.
The thing we want to avoid is finding the peaks for a given (sx,sy) by
testing this predicate against every peak in the sinogram, which is very slow.

We do this in the following way (see VoxelSinoMasker.partition() for details):
1. We bin the peaks in omega and dty, either by supplied omega bin centres or by a heuristic that guesses a good omega bin size. dty bins are fixed by ystep.
2. We sort the peaks by omega bin indices (primary) and raw dty values (secondary) to determine a peak ordering. We use raw dty because it changes within an omega bin.
3. We make partition arrays, which track the start of each omega and dty bin in the sorted arrays
omega_partitions[i] is the index of the first (sorted) peak in omega bin i
dty_partitions[i,j] is the offset, within omega bin i, of the first (sorted) peak in dty bin j
4. We precompute sin(omega) and cos(omega) for each peak for speed
5. We precompute sin(omega_bin) and cos(omega_bin) for each omega bin for speed

Then, to find the peaks that illuminate a voxel at (sx,sy) (see fill_voxel_idx() for details):
1. We walk the omega bins, and for each bin compute the dty value of the voxel at that omega
2. We determine the corresponding dty value for the -next- omega bin
3. We convert the difference in dty values between omega bins into a dty bin padding for that omega bin
4. This tells us, as a function of omega, how "wide" our search space is in dty bins
5. We walk the omega bins again, computing the y coordinate of the voxel at that omega, and the corresponding dty bin
6. Including the dty bin padding from the previous step, we now know the dty partitions to visit for this omega bin
7. We get the start and end indices from dty_partitions and omega_partitions
8. Within an omega partition, peaks are monotonically increasing in dty, so we can just do peaks[start:end]
9. For any peak in that range, we compute the exact y distance to the voxel and check if it is within ystep. If so, we keep it.
10. We fill the idx_buffer and ydist_buffer with the (sorted!) indices and distances of the peaks that illuminate the voxel.
We also keep track of how many peaks we found (m), and return that number. If the buffers are too small, we return -1.
This whole process means that we never miss peaks. The only purpose is to reduce the number of peaks we have to test the predicate for.

We now have an idx_buffer and ydist_buffer that contain the indices and distances of the peaks that illuminate the voxel at (sx,sy).
By just returning the first m elements of the buffers, we can avoid any memory allocation and reuse the buffers for all voxels in a grid.
These buffers can be optimally sized for a given (sx,sy) grid - see max_candidates for details.
"""

import os
import sys

import h5py
import numba
import numpy as np

# -------- CACHING -------- #
# It is handy to cache the results of a partition to disk, so that we can reuse it later without recomputing it.
# This is useful for memory-mapping between workers (pbp indexing)
# and also going between a local machine and a cluster machine.

_CACHE_VERSION = 7


## Save and load caches
_PEAK_H5GROUP = "PBPPeakIndex"


def default_index_filename(manager):
    """Generate a default filename for the peak + gve index.

    Parameters
    ----------
    manager
        PBP or PBPRefine object

    Returns
    -------
        Filename to save/load peak index/gve caches

    Raises
    ------
    ValueError
        If no icolf_filename on the manager object
    """
    base = getattr(manager, "icolf_filename", None)
    if base is None:
        raise ValueError("no icolf_filename to derive an index cache path from")
    return os.path.splitext(base)[0] + "_pbpindex.h5"


def save_peak_index_cache(idx, filename):
    """Save peak index to an on-disk H5 cache, with versioning.

    Parameters
    ----------
    idx
        Peak index to cache
    filename
        Filename to save to

    """
    with h5py.File(filename, "a") as hout:
        if _PEAK_H5GROUP in hout:
            del hout[_PEAK_H5GROUP]
        g = hout.create_group(_PEAK_H5GROUP)
        g.attrs["version"] = _CACHE_VERSION
        for k in ("nom", "ndty"):
            g.attrs[k] = int(idx[k])
        # optional: only the indexer pre-computes it
        if idx.get("maxlocal") is not None:
            g.attrs["maxlocal"] = int(idx["maxlocal"])
        for k in ("ymin", "ray_margin", "dev"):
            g.attrs[k] = float(idx[k])
        # no chunks, no compression -> contiguous -> id.get_offset() works,
        # which is what makes the mmap path below possible. Do not add
        # compression here without also disabling that path.
        for k in ("order", "omega_partitions", "dty_partitions", "usin", "ucos"):
            g.create_dataset(k, data=idx[k])


def load_peak_index_cache(filename, mmap=False, verbose=True):
    """Load peak index from an on-disk H5 cache, with optional mem-map.

    Use it from multiprocessing workers: the arrays are identical in every
    worker and never written, so one mapping shared through the page cache
    replaces one heap copy per process. Falls back to a read for any dataset
    that is not contiguous.

    Parameters
    ----------
    filename
        Filename to load from
    mmap, optional
        returns read-only np.memmap views instead of copies. by default False
    verbose, optional
        Print info, by default True

    Returns
    -------
        Dictionary representing a peak index
    """

    if filename is None:
        return None
    if not os.path.exists(filename):
        if verbose:
            print("index: %s does not exist" % filename, file=sys.stderr)
        return None
    try:
        hin = h5py.File(filename, "r")
    except OSError as e:
        # HDF5 file locking on GPFS/NFS lands here. Do not report it as a
        # stale cache -- set HDF5_USE_FILE_LOCKING=FALSE.
        print("index: cannot open %s: %s" % (filename, e), file=sys.stderr)
        return None
    with hin:
        if _PEAK_H5GROUP not in hin:
            print("index: no %s group" % _PEAK_H5GROUP, file=sys.stderr)
            return None
        g = hin[_PEAK_H5GROUP]
        v = g.attrs.get("version", -1)
        if v != _CACHE_VERSION:
            print("index: version %s, expected %d -- rerun setpeaks()"
                  % (v, _CACHE_VERSION), file=sys.stderr)
            return None
 
        def get(k):
            d = g[k]
            if not mmap:
                return d[:]
            off = d.id.get_offset()
            if off is None:          # chunked/compressed: cannot map
                return d[:]
            return np.memmap(filename, dtype=d.dtype, mode="r",
                             offset=off, shape=d.shape)
 
        return dict(
            nom=int(g.attrs["nom"]), ndty=int(g.attrs["ndty"]),
            ymin=float(g.attrs["ymin"]),
            ray_margin=float(g.attrs["ray_margin"]),
            dev=float(g.attrs["dev"]),
            maxlocal=(int(g.attrs["maxlocal"])
                      if "maxlocal" in g.attrs else None),
            order=get("order"),
            omega_partitions=get("omega_partitions"),
            dty_partitions=get("dty_partitions"),
            usin=get("usin"), ucos=get("ucos"),
        )

# -------- OMEGA BINNING -------- #
# A voxel at (sx,sy) is illuminated when:
# | y0 - sx*sin(omega) - sy*cos(omega) - dty | <= ystep
# for one voxel that traces a sinusoid in (omega, dty) space.
# The thing we want to avoid is finding the peaks for a given (sx,sy) by
# testing this predicate against every peak in the sinogram, which is very slow.
# 
# Instead:
# 1. Put peaks into (omega, dty) bins
# 2. Given a (sx,sy) voxel, walk over omega bins
# 3. For each omega bin, compute the dty range that satisfies the predicate
# 4. Use the dty partition to find the peaks in that range, and test the predicate on them.
# 
# This is basically a heirarchical spatial hash: the omega bins are the first level,
# and the dty bins are the second level.
# 
# We can choose different ways to do this based on what information we have:
# 1. Frame number
# If you segmented recently, you'll have frame number available in cf_2d
# It's a direct lookup from frame number to omega bin index
# 2. dset.obinedges
# If you segmented before frame number was available, you can use the dset.obinedges to bin omega
# 3. Heuristic
# If you don't have either of the above, we can use a heuristic to bin omega based on the dty bins

def _omega_bin_index_from_frm(manager, omega, verbose=True):
    """Omega bin index via the per-peak frame number.

    frm says which frame each peak came from. We do NOT assume the
    frame-in-scan index tracks omega: a zigzag scan, or a rotation that runs
    continuously across dty rows, breaks that and there is nothing in the file
    that promises otherwise. Instead we bin each *frame's* omega once, straight
    out of dset.omega, and index that table with frm. Exact for any scan
    ordering, and it digitizes nscans*nframes values rather than one per peak.

    Parameters
    ----------
    manager
        PBP or PBPRefine object
    omega
        Columnfile column
    verbose, optional
        Print diagnostics, by default True

    Returns
    -------
    iomega
        Omega bin indices per peak
    cens
        Omega bin centers
    dev
        Maximum deviation from bin centers
    """

    icolf = manager.icolf
    if "frm" not in icolf.titles:
        return None
    dset = getattr(manager, "dset", None)
    om_img = getattr(dset, "omega", None)
    edges = getattr(dset, "obinedges", None)
    cens = getattr(dset, "obincens", None)
    if om_img is None or edges is None or cens is None:
        if verbose:
            print("frm present but dset.omega/obinedges missing; "
                  "falling back to omega binning")
        return None
    om_img = np.asarray(om_img, dtype=np.float64)
    edges = np.asarray(edges, dtype=np.float64)
    cens = np.asarray(cens, dtype=np.float64)
    frm = np.round(np.asarray(icolf.frm, dtype=np.float64)).astype(np.int64)
    if frm.size == 0 or frm.min() < 0 or frm.max() >= om_img.size:
        if verbose:
            print("frm values fall outside dset.omega (%d frames); "
                  "falling back to omega binning" % om_img.size)
        return None

    flat = om_img.ravel()
    if edges[0] >= -1e-9 and flat.min() < -1e-9:
        flat = flat % 360.0
    iom_of_frame = np.digitize(flat, edges) - 1
    np.clip(iom_of_frame, 0, len(cens) - 1, out=iom_of_frame)
    iomega = iom_of_frame[frm]
    dev = float(np.abs(omega - cens[iomega]).max())

    if verbose:
        live = int((np.bincount(iomega, minlength=len(cens)) > 0).sum())
        note = ""
        if om_img.ndim == 2 and "dtyi" in icolf.titles:
            # does the scan row implied by frm agree with the dty bin?
            iscan = frm // om_img.shape[1]
            db = np.asarray(icolf.dtyi, dtype=np.int64)
            note = (", scan row agrees with dty bin for %.1f%% of peaks"
                    % (100.0 * float((iscan == db).mean())))
        print("omega bins via frm: %d bins (%d with peaks), %.1f peaks/bin, "
              "|omega - bin centre| max %.5f deg%s"
              % (len(cens), live, frm.size / max(live, 1), dev, note))
    return iomega, cens, dev

def _omega_bin_index(manager, omega, omega_step=None, verbose=True):
    """Map each peak onto an omega frame without needing a frm column.
 
    For data segmented before frm existed. Bins on dset.obinedges, which is
    one entry per omega position, so the result is the same shape as the frm
    decode.
 
    Returns (iomega, omega_centres), or None if there is nothing better than
    guessing, in which case the caller should use the heuristic.

    Parameters
    ----------
    manager
        PBP or PBPRefine object
    omega
        Columnfile column
    omega_step, optional
        _description_, by default None
    verbose, optional
        Print diagnostics, by default True

    Returns
    -------
    iomega
        Omega bin indices per peak
    cens
        Omega bin centers
    """

    dset = getattr(manager, "dset", None)
 
    if omega_step is not None:
        o0 = float(omega.min())
        iomega = np.round((omega - o0) / float(omega_step)).astype(np.int64)
        iomega -= iomega.min()
        cens = o0 + np.arange(iomega.max() + 1) * float(omega_step)
        if verbose:
            print("omega bins: %d, from omega_step=%g" % (len(cens), omega_step))
        return iomega, cens
 
    edges = getattr(dset, "obinedges", None) if dset is not None else None
    cens = getattr(dset, "obincens", None) if dset is not None else None
    if edges is None or cens is None:
        if verbose:
            print("omega bins: no frm and no dset.obinedges -- using the "
                  "heuristic, which only guesses at the sample radius")
        return None
 
    edges = np.asarray(edges, dtype=np.float64)
    cens = np.asarray(cens, dtype=np.float64)
    om = omega % 360.0 if edges[0] >= -1e-9 and omega.min() < -1e-9 else omega
    iomega = np.digitize(om, edges) - 1
    np.clip(iomega, 0, len(cens) - 1, out=iomega)
    iomega = iomega.astype(np.int64)
    dev = float(np.abs(omega - cens[iomega]).max())
    if verbose:
        live = int((np.bincount(iomega, minlength=len(cens)) > 0).sum())
        print("omega bins via obinedges (no frm): %d bins (%d with peaks), "
              "|omega - bin centre| max %.5f deg" % (len(cens), live, dev))
    return iomega, cens


def choose_omega_bins(manager, omega, omega_binsize="auto", omega_step=None,
                      verbose=True):
    """Decide how to bin omega, best method first.

    1. frm column decoded through dset.omega / obinedges  -- exact
    2. dset.obinedges alone                               -- pre-frm data
    3. VoxelSinoMasker's heuristic                        -- last resort

    An explicit omega_step or omega_binsize overrides all three.

    manager needs .icolf and .dset; both PBPRefine and PBP have them.

    Returns (kw, iomega) where kw goes straight to partition(**kw) and
    iomega is the per-peak bin index, or None when the heuristic is used and
    the caller must reconstruct it. Callers that need `dev` should use
    iomega rather than the masker's permuted omega, which keep_columns=False
    does not build.

    Parameters
    ----------
    manager
        PBP or PBPRefine object
    omega
        Columnfile column
    omega_binsize, optional
        _description_, by default "auto"
    omega_step, optional
        _description_, by default None
    verbose, optional
        _description_, by default True

    Returns
    -------
    (kw, iomega)
        iomega: Like dtyi but for omega. Integer omega_steps away from omega min.
        kw: dict of the stuff you need to build
    """

    if omega_step is not None:
        if verbose:
            print("omega bins: omega_step=%g (explicit)" % omega_step)
        return dict(omega_binsize=float(omega_step)), None
 
    if omega_binsize not in (None, "auto"):
        if verbose:
            print("omega bins: omega_binsize=%g (explicit)" % omega_binsize)
        return dict(omega_binsize=float(omega_binsize)), None
 
    got = _omega_bin_index_from_frm(manager, omega, verbose)          # 1
    if got is None:  # frm not available, fall back to obinedges
        got = _omega_bin_index(manager, omega, verbose=verbose)  # 2
    if got is not None:  # frm or obinedges passed
        # _frame_index_from_frm returns 3 values, _omega_frame_index 2;
        # the first two are the same in both
        return dict(omega_bin_indices=got[0], omega_bins=got[1]), got[0]
 
    return dict(omega_binsize=None), None                      # 3


@numba.njit(cache=True)
def fill_voxel_idx(
    sx,
    sy,
    y0,
    ystep,
    ymin,
    omega_partitions,
    dty_partitions,
    dty_sorted,
    sin_omega_sorted,
    cos_omega_sorted,
    sin_omega_bins,
    cos_omega_bins,
    idx_buffer,    # | you will want to reuse these for all voxels!
    ydist_buffer,  # |
):
    """Get the voxel indices and y-distances for a given voxel centroid position.

    This tracks the voxel in sinogram space and returns all peaks the voxel possibly intersects.

    The idea is to use a partitioning scheme that allows for retrival of these peaks without touching
    the entire sinogram.

    Args:
        sx (:obj:`float`): The x-coordinate of the voxel position in sample coordinates.
        sy (:obj:`float`): The y-coordinate of the voxel position in sample coordinates.
        y0 (:obj:`float`): Rotaiton axis offset.
        ystep (:obj:`float`): The y-step size. Voxels with centroids within ystep from the beam
            center are masked.
        ymin (:obj:`float`): The minimum y-coordinate expected dty_sorted.
        omega_partitions (:obj:`np.ndarray`): The omega partitions. See VoxelSinoMasker.partition() for details.
        dty_partitions (:obj:`np.ndarray`): The dty partitions. See VoxelSinoMasker.partition() for details.
        dty_sorted (:obj:`np.ndarray`): The sorted dty values. See VoxelSinoMasker.partition() for details.
        sin_omega_sorted (:obj:`np.ndarray`): The sorted sin_omega values. See VoxelSinoMasker.partition() for details.
        cos_omega_sorted (:obj:`np.ndarray`): The sorted cos_omega values. See VoxelSinoMasker.partition() for details.
        sin_omega_bins (:obj:`np.ndarray`): The sin_omega bins. See VoxelSinoMasker.partition() for details.
        cos_omega_bins (:obj:`np.ndarray`): The cos_omega bins. See VoxelSinoMasker.partition() for details.
        idx_buffer (:obj:`np.ndarray`): The buffer for the voxel indices. This is intended to be reused when masking
            voxels over a grid. The buffer size is user-defined allowing for optimizations. Setting it to dty_sorted.size
            is 100% safe, but likely very very overkill.
        idx_buffer (:obj:`np.ndarray`): The buffer for the voxel indices. This is intended to be reused when masking
            voxels over a grid. The buffer size is user-defined allowing for optimizations. Setting it to dty_sorted.size
            is 100% safe, but likely very very overkill.

    Returns:
        :obj:`int`: the number of peaks written into the buffers, or -1 if the
        buffers were too small.

    Note:
        This returns a sentinel rather than raising because it is meant to be
        callable from inside a numba prange, where an exception is not
        reliably propagated: the thread stops writing and every item it still
        owned is silently left unfilled. get_voxel_idx() wraps this and raises,
        for Python-level callers.
    """
    n_bins = sin_omega_bins.size

    if n_bins == 0:  # this should never happen, but just to be clean...
        return 0

    # figure out a safe dty padding for each omega bin such that no peak is forgotten.
    y_bins_diff = np.empty(n_bins, dtype=np.int64)
    inv_ystep = 1.0 / ystep

    # First dty value for the first omega bin
    y_prev = y0 - sx * sin_omega_bins[0] - sy * cos_omega_bins[0]

    # Loop over omega bins
    # This lets us work out how many dty bins we need to consider for each omega bin, so that we don't miss any peaks.
    for i in range(n_bins - 1):  # every omega partition bin...
        y_next = y0 - sx * sin_omega_bins[i + 1] - sy * cos_omega_bins[i + 1]
        # Difference in y between the next bin and the previous bin
        diff = (
            y_next - y_prev
        )  # TODO: this is a safe upper bound, but it could likely be tightened with a second order Taylor or similar.
        if diff < 0.0:
            diff = -diff
        y_bins_diff[i] = int(diff * inv_ystep + 0.5) + 1
        y_prev = y_next
    if n_bins > 1:
        y_bins_diff[n_bins - 1] = y_bins_diff[n_bins - 2]
    else:
        y_bins_diff[0] = 1
    # end of dty padding for each omega bin

    m = 0

    row_len = dty_partitions.shape[1]
    row_len_minus_1 = row_len - 1

    buffer_size = idx_buffer.size
    for i in range(n_bins):  # for every omega bin...
        # compute where we are in y
        y = y0 - sx * sin_omega_bins[i] - sy * cos_omega_bins[i]

        # convert to integer index
        # floor(x + 0.5), not int(x + 0.5): int() truncates toward zero and so
        # rounds the wrong way for negative indices, which occur whenever the
        # voxel sinusoid swings below ymin.
        y_index = int(np.floor((y - ymin) * inv_ystep + 0.5))

        # get the dty padding for this omega bin, i.e we will not
        # just grab the y_index in the bin but some amount of neighboring
        # dty values on either side of the y_index must be considered.
        dty_partition_padding = y_bins_diff[i]

        padded_low = y_index - dty_partition_padding
        # clip to the valid range of dty partitions, which is [0, row_len-1]
        if padded_low < 0:
            padded_low = 0
        padded_high = y_index + dty_partition_padding + 1
        if padded_high > row_len_minus_1:
            padded_high = row_len_minus_1

        # if there is simply no scan close to the voxel for this omega bin, then skip it completely.
        if padded_high <= 0:
            continue

        if padded_low >= row_len_minus_1:
            continue

        # if we got here, for a given (sx,sy) voxel, at this omega,
        # there is a dty sub-partition that we need to consider.

        # get the safe start and end of the dty sub-partition for this omega bin
        dty_partition_start = dty_partitions[i, padded_low]
        dty_partition_end = dty_partitions[i, padded_high]

        # get the start of the omega partition
        omega_partition_start = omega_partitions[i]

        # iterate over the feasible dty chunk in global peak indices
        for j in range(
            omega_partition_start + dty_partition_start,
            omega_partition_start + dty_partition_end,
        ):  # for each peak in the dty-sub-partition...
            # Now we compute the exact y-distance to the voxel centroid
            dty_val = dty_sorted[j]
            s = sin_omega_sorted[j]
            c = cos_omega_sorted[j]
            ydist = abs(y0 - sx * s - sy * c - dty_val)

            # Check if the peak is within ystep of the voxel centroid
            # if so, then keep it.
            if ydist <= ystep:
                if m >= buffer_size:  # buffer is too small, return sentinel
                    return -1
                idx_buffer[m] = j  # these now refer to the sorted peaks!
                ydist_buffer[m] = ydist
                m += 1

    return m


@numba.njit(cache=True)
def get_voxel_idx(
    sx,
    sy,
    y0,
    ystep,
    ymin,
    omega_partitions,
    dty_partitions,
    dty_sorted,
    sin_omega_sorted,
    cos_omega_sorted,
    sin_omega_bins,
    cos_omega_bins,
    idx_buffer,
    ydist_buffer,
):
    """Slice-returning wrapper around fill_voxel_idx(). See that for details.

    Raises:
        :obj:`ValueError`: If the buffer arrays are too small. Then simply try to set your buffer size larger.
    """
    m = fill_voxel_idx(
        sx, sy, y0, ystep, ymin,
        omega_partitions, dty_partitions, dty_sorted,
        sin_omega_sorted, cos_omega_sorted,
        sin_omega_bins, cos_omega_bins,
        idx_buffer, ydist_buffer,
    )
    if m < 0:
        raise ValueError("Buffer is too small")
    return idx_buffer[:m], ydist_buffer[:m]


@numba.njit(cache=True, parallel=True)
def max_candidates(
    sx_arr,
    sy_arr,
    y0,
    ystep,
    ymin,
    dty_partitions,
    sin_omega_bins,
    cos_omega_bins,
):
    """Upper bound on how many peaks can be found per voxel.

    Works very similarly to fill_voxel_idx().
    Here, we only read the partition offsets to find out how many peaks are in the dty sub-partition for each omega bin.
    We never compute the predicate on a peak level, so this remains fast.
    This is used to determine the appropriate buffer size for the idx_buffer and ydist_buffer that are used in fill_voxel_idx().

    Parameters
    ----------
    sx_arr
        Array of x-coordinates of the voxel positions in sample coordinates.
    sy_arr
        Array of y-coordinates of the voxel positions in sample coordinates.
    y0
        The rotation axis offset.
    ystep
        The y-step size.
    ymin
        The minimum y-coordinate.
    dty_partitions
        The dty partitions.
    sin_omega_bins
        The sin_omega bins.
    cos_omega_bins
        The cos_omega bins.

    Returns
    -------
        The maximum number of candidates across all supplied voxels
    """

    n_bins = sin_omega_bins.size
    row_len_minus_1 = dty_partitions.shape[1] - 1
    inv_ystep = 1.0 / ystep
    best = np.zeros(sx_arr.size, dtype=np.int64)
    # loop over voxels in parallel
    for v in numba.prange(sx_arr.size):  # ty: ignore[not-iterable]
        sx = sx_arr[v]
        sy = sy_arr[v]
        total = 0
        # compute the y-value for the first omega bin
        y_prev = y0 - sx * sin_omega_bins[0] - sy * cos_omega_bins[0]
        pad = 1
        # loop over other omega bins
        for i in range(n_bins):
            if i + 1 < n_bins:
                y_next = y0 - sx * sin_omega_bins[i + 1] - sy * cos_omega_bins[i + 1]
                diff = y_next - y_prev
                if diff < 0.0:
                    diff = -diff
                pad = int(diff * inv_ystep + 0.5) + 1
            else:
                y_next = y_prev
            y_index = int(np.floor((y_prev - ymin) * inv_ystep + 0.5))
            padded_low = y_index - pad
            if padded_low < 0:
                padded_low = 0
            padded_high = y_index + pad + 1
            if padded_high > row_len_minus_1:
                padded_high = row_len_minus_1
            if padded_high > 0 and padded_low < row_len_minus_1:
                # Count the number of peaks in the dty sub-partition for this omega bin
                total += dty_partitions[i, padded_high] - dty_partitions[i, padded_low]
            y_prev = y_next
        best[v] = total
    return best.max()


class VoxelSinoMasker:
    """Masker for voxels in sinogram space.

    This class is used to mask voxels in sinogram space. It is designed to be used with a regular grid of voxels.

    A partition scheme is used such that peaks are sorted by omega (primary) and dty (secondary).
    This allows for fast retrieval of peaks that are within a voxel by only looking at the relevant bins.

    The partition scheme is built on the assumption that the scan defines a regular grid in dty. Omega values can be
    arbitrary, (i.e merged 3d peaks should work as well as single 2d peaks).

    omega_partitions[i] is the index of the first (sorted) peak in omega bin i
    dty_partitions[i,j] is the offset, within omega bin i, of the first (sorted) peak in dty bin j

    Adaptable omega bin sizes are used to chunk the data over angles. For a call to mask() the algorithm then
    oterates over the omega bins and collects partitions of dty values that are candidates for the voxel.

    Args:
        omega (:obj:`np.ndarray`): The omega values, one per peak, shape=(n_peaks,). (degrees)
        dty (:obj:`np.ndarray`): The dty values, one per peak, shape=(n_peaks,).
        dty_stepsize (:obj:`float`): The dty step size, assumed to be constant such that the
            scan defines a regular grid in dty.
        ymin: ds.ybincens.min()

    Attributes:
        omega_binsize (:obj:`float`): The omega bin size.
        peak_ordering (:obj:`np.ndarray`): The peak ordering.
        omega_partitions (:obj:`np.ndarray`): The omega partitions.
        omega_bins (:obj:`np.ndarray`): The omega bins.
        dty_partitions (:obj:`np.ndarray`): The dty partitions.
        ymin (:obj:`float`): The minimum y-coordinate.
        n_peaks (:obj:`int`): The number of peaks.
        omega (:obj:`np.ndarray`): The omega values.
        sinomega (:obj:`np.ndarray`): The sin_omega values.
        cosomega (:obj:`np.ndarray`): The cos_omega values.
        dty (:obj:`np.ndarray`): The dty values.
        dty_stepsize (:obj:`float`): The dty step size.
    """

    def __init__(self, omega, dty, dty_stepsize, ymin=None):
        self.uniq_dtys = np.array([np.abs(dty).max()])
        self.ymin = dty.min() if ymin is None else float(ymin)

        # dty need NOT lie on a grid. The partition is a spatial hash and the
        # exact predicate is applied per peak with that peak's own dty, so a
        # continuously-moving dty motor is fine -- the buckets just fill
        # evenly rather than clustering. What must hold is that dty_stepsize
        # is a sensible bucket width for the range being covered.
        span = float(dty.max() - self.ymin)
        nbin = span / dty_stepsize
        if not np.isfinite(nbin) or nbin < 1 or nbin > 1e7:
            raise ValueError(
                "dty spans %.4g with dty_stepsize %.4g -- that is %.3g bins. "
                "Check dty_stepsize." % (span, dty_stepsize, nbin))
        if dty.min() < self.ymin - 0.5 * dty_stepsize:
            raise ValueError(
                "dty.min() = %.6g is below ymin = %.6g by more than half a "
                "step, so peaks would bin to negative rows."
                % (dty.min(), self.ymin))

        self.n_peaks = len(dty)
        self.omega = omega
        self.sinomega = None
        self.cosomega = None
        self.dty = dty
        self.dty_stepsize = dty_stepsize

        self.omega_binsize = None
        self.peak_ordering = None
        self.omega_partitions = None
        self.omega_bins = None
        self.dty_partitions = None
    
    def __str__(self):
        return (f"VoxelSinoMasker(n_peaks={self.n_peaks}, "
                f"omega_range=[{self.omega.min():.2f}, {self.omega.max():.2f}], "
                f"dty_range=[{self.dty.min():.6g}, {self.dty.max():.6g}], "
                f"dty_stepsize={self.dty_stepsize:.6g}, "
                f"omega_binsize={self.omega_binsize})")

    def _heuristic_omega_binsize(self):
        # find the omega binsize such that for any voxel on a regular grid
        # the movement of the voxel over one bin step is about 2 * dty_stepsize.
        # this is a heuristic compromise between the required dty padding per
        # omega bin collection and the number of omega bins that are required.
        # more omega bins is slower, and more dty padding is slower, so we
        # are looking for the optimal compromise that will be the fastest.
        # it is also possible to do actual test runs to get a final truth..
        omega_binsize = np.inf
        for y in self.uniq_dtys:
            dydw_rad = np.sqrt(2.0) * y  # change in y per radian of omega
            dydw_deg = (np.pi / 180.0) * dydw_rad  # change in y per degree of omega
            if dydw_deg == 0.0:
                continue
            dom = np.abs(
                self.dty_stepsize / dydw_deg
            )  # degrees such that change in y becomes dty_stepsize
            if dom < omega_binsize:
                omega_binsize = dom  # we track the bound
        return (
            2 * omega_binsize
        )  # a factor of 2 here seems to work nicely from empricial testing

    def partition(self, omega_binsize=None, inplace=False,
                  omega_bin_indices=None, omega_bins=None,
                  keep_columns=True):

        """Partition the peaks for fast peak masking.

        A partition scheme is used such that peaks are sorted by omega (primary) and dty (secondary).
        This allows for fast retrieval of peaks that are within a voxel by only looking at the relevant bins.

        The partition scheme is built on the assumption that the scan defines a regular grid in dty. Omega values can be
        arbitrary, (i.e merged 3d peaks should work as well as single 2d peaks).

        Adaptable omega bin sizes are used to chunk the data over angles. For a call to mask() the algorithm then
        oterates over the omega bins and collects partitions of dty values that are candidates for the voxel.

        For a finer omega_binsize the algorithm can decrease the size of the dty partitions it needs to visit, however,
        more omega bins are then required. For any dataset there exists an optimal omega_binsize that minimizes the number
        of flops. This can also be machine dependent, so for absolutely best performance you should test different values.
        However, the default heuristic should be very performant for most datasets.

        Args:
            omega_binsize (:obj:`float`): The omega bin size. If None, a heuristic will be used to find a good value. (degrees)
                (note that this is not the same as the experimental stepsize of the omega motor.)
            inplace (:obj:`bool`): If True, reorder self.omega and self.dty in place. These are the
                same arrays the caller passed to __init__, so in place means the caller's columns are
                reordered too and every other column they hold is then misaligned against them.
                Defaults to False, which costs one copy of each.
            omega_bin_indices (:obj:`np.ndarray`): Optional precomputed bin index per peak, e.g. decoded
                from a frame number. Overrides omega_binsize.
            omega_bins (:obj:`np.ndarray`): Bin centres in degrees matching omega_bin_indices.
            keep_columns (:obj:`bool`): If False, the permuted dty/omega/
                sinomega/cosomega copies and the mask() buffers are not
                built. Saves four full columns of memory. mask() then raises;
                use this when you only want peak_ordering and the partition
                arrays and will permute your own columns.

        Raises:
            :obj:`ValueError`: If the omega binsize is not set and the heuristic fails to find a good value.
        """
        # Prepare the omega bins
        # Depends on whether the caller has supplied a frm column/obincens or not. If they have, we can use that to get the exact binning.
        # Otherwise, we need to use the heuristic to find a good binning.
        if omega_bin_indices is not None:
            # caller supplied the binning, e.g. exact frames decoded from a
            # frame number column
            self.omega_binsize = None
            omega_bin_indices = np.ascontiguousarray(omega_bin_indices, dtype=np.int64)
            omega_bins = np.asarray(omega_bins, dtype=np.float64)[
                : omega_bin_indices.max() + 1
            ]
        else:
            if omega_binsize is None:
                self.omega_binsize = self._heuristic_omega_binsize()
            else:
                self.omega_binsize = omega_binsize

            omin = self.omega.min()

            omega_bin_indices = np.round(
                (self.omega - omin) / self.omega_binsize
            ).astype(np.int64)
            max_idx = omega_bin_indices.max()
            omega_bins = omin + self.omega_binsize * np.arange(max_idx + 1)

        # Lexsort by dty first, then omega, so that the omega partitions are contiguous and the dty partitions are contiguous within each omega partition.
        # omega_bin_indices is npeaks long and could be unsorted
        peak_ordering = np.lexsort((self.dty, omega_bin_indices))
        omega_bin_indices_sorted = omega_bin_indices[peak_ordering]
        # work out the length of each omega partition and the start of each partition in the sorted array
        omega_partition_lengths = np.bincount(omega_bin_indices_sorted, minlength=omega_bins.size)
        omega_partitions = np.zeros(omega_bins.size + 1, dtype=np.int64)
        omega_partitions[1:] = np.cumsum(omega_partition_lengths)

        # get also the dty partitions
        # sort them in the same way as omega
        dty_sorted = self.dty[peak_ordering]
        # this is dtyi
        dty_bin_indices = np.round((dty_sorted - self.ymin) / self.dty_stepsize).astype(
            np.int64
        )
        ymax_index = dty_bin_indices.max()

        nbytes = omega_bins.size * (ymax_index + 2) * 8
        if nbytes > 2e9:
            raise ValueError(
                "the partition would need %.1f GB (%d omega bins x %d dty "
                "bins). Pass a coarser omega_binsize -- the dty window per "
                "bin widens to compensate, and the optimum is fairly flat."
                % (nbytes / 1e9, omega_bins.size, ymax_index + 2))
        if nbytes > 2e8:
            print("note: partition is %.0f MB (%d x %d)"
                  % (nbytes / 1e6, omega_bins.size, ymax_index + 2))

        dty_partitions = np.empty((omega_bins.size, ymax_index + 2), dtype=np.int64)
        for i in range(omega_bins.size):
            dty_indices_for_bin = dty_bin_indices[
                omega_partitions[i] : omega_partitions[i + 1]
            ]
            dty_sub_partition_lengths = np.bincount(
                dty_indices_for_bin, minlength=ymax_index + 1
            )

            dty_sub_partitions = np.zeros(
                dty_sub_partition_lengths.size + 1, dtype=np.int64
            )
            dty_sub_partitions[1:] = np.cumsum(dty_sub_partition_lengths)
            dty_partitions[i] = dty_sub_partitions

        self.peak_ordering = peak_ordering
        self.omega_partitions = omega_partitions
        self.omega_bins = omega_bins
        self.dty_partitions = dty_partitions

        if not keep_columns:
            # partition arrays and peak_ordering only -- callers that permute
            # their own columns do not need these, and they are four full
            # columns of memory.
            self.dty = None
            self.omega = None
            self.sinomega = None
            self.cosomega = None
            self.sinomega_bins = np.sin(np.radians(omega_bins))
            self.cosomega_bins = np.cos(np.radians(omega_bins))
            self.buffer_size = 0
            self.idx_buffer = None
            self.ydist_buffer = None
            return
 
        if inplace:
            # saves memory, but overwrites the arrays the caller handed us
            self.dty[:] = dty_sorted
            self.omega[:] = self.omega[peak_ordering]
        else:
            self.dty = dty_sorted
            self.omega = self.omega[peak_ordering]
        self.sinomega = np.sin(np.radians(self.omega))
        self.cosomega = np.cos(np.radians(self.omega))
 
        self.sinomega_bins = np.sin(np.radians(omega_bins))
        self.cosomega_bins = np.cos(np.radians(omega_bins))
 
        # we default to a 5% buffer, this should be very safe, altough
        # we have fallbacks in case of a sharp corner.
        # this can be optimised via max_candidates() if the caller knows the voxel grid in advance
        self._build_buffers(self.n_peaks // 20)

    def max_candidates(self, xi0s, yi0s, ystep, y0):
        """Largest number of peaks any of these voxels can pull out of the partition.

        If you size the buffers with this, fill_voxel_idx() can never overflow, which
        matters when using Numba to avoid getting invalid results.
        """
        return int(
            max_candidates(
                np.ascontiguousarray(xi0s, dtype=np.float64),
                np.ascontiguousarray(yi0s, dtype=np.float64),
                y0,
                ystep,
                self.ymin,
                self.dty_partitions,
                self.sinomega_bins,
                self.cosomega_bins,
            )
        )

    def _build_buffers(self, buffer_size):
        self.buffer_size = buffer_size
        self.idx_buffer = np.empty(buffer_size, dtype=np.int64)
        self.ydist_buffer = np.empty(buffer_size, dtype=np.float64)

    def _mask(self, xi0, yi0, ystep, y0):
        if self.dty is None:
            raise ValueError(
                "partition(keep_columns=False) was used, so the permuted "
                "columns and mask buffers were not built. Re-partition with "
                "keep_columns=True to use mask().")

        return get_voxel_idx(
            xi0,
            yi0,
            y0,
            ystep,
            self.ymin,
            self.omega_partitions,
            self.dty_partitions,
            self.dty,
            self.sinomega,
            self.cosomega,
            self.sinomega_bins,
            self.cosomega_bins,
            self.idx_buffer,
            self.ydist_buffer,
        )

    def mask(self, xi0, yi0, ystep, y0):
        while self.buffer_size < self.n_peaks:
            try:
                idx, ydist = self._mask(xi0, yi0, ystep, y0)
                return idx, ydist
            except ValueError:
                # Increasing buffer size dynamically, this is overhead on first call,
                # if the buffer is too small.
                self._build_buffers(2 * self.buffer_size)
                continue
        idx, ydist = self._mask(xi0, yi0, ystep, y0)
        return idx, ydist

if __name__ == "__main__":
    dty_stepsize = 0.003  # 3 microns scan step
    omega_stepsize = 0.05
    frames_omega = np.arange(0, 360, omega_stepsize)
    unique_dty = np.arange(-0.4, 0.4, dty_stepsize)

    print("Number of unique dty values: {}".format(len(unique_dty)))
    print("Number of frames per scan step: {}".format(len(frames_omega)))

    omega = []
    dty = []

    for i in range(len(unique_dty)):
        # make sure we can handle 3d merged peaks
        # and not just the 2D merged case.
        n = np.random.randint(
            2 * frames_omega.size, 16 * frames_omega.size
        )  # 2- 16 peaks per frame
        omega.extend(list(np.random.uniform(frames_omega.min(), frames_omega.max(), n)))
        dty.extend([unique_dty[i]] * n)

    omega = np.array(omega)
    dty = np.array(dty)

    n_peaks = len(dty)

    permute = np.random.permutation(n_peaks)
    omega = omega[permute]
    dty = dty[permute]

    # permute floating point motor positions to make sure we can handle noise in the motor positions
    omega += np.random.uniform(-omega_stepsize * 0.001, omega_stepsize * 0.001, n_peaks)
    dty += np.random.uniform(-dty_stepsize * 0.001, dty_stepsize * 0.001, n_peaks)

    xi0, yi0 = 0.22, 0.27
    y0 = -0.023
    ystep = dty_stepsize

    def brute_force_mask(xi0, yi0, y0, ystep, sinomega, cosomega, dty):
        """Reference: test the predicate against every peak in the sinogram.

        This is what the partitioning scheme has to reproduce exactly, and
        what it is avoiding -- it touches all n_peaks for every voxel.
        """
        ydist = np.abs(y0 - xi0 * sinomega - yi0 * cosomega - dty)
        idx = np.where(ydist <= ystep)[0]
        return idx, ydist[idx]

    peak_selector = VoxelSinoMasker(omega, dty, dty_stepsize)
    peak_selector.partition()

    print("Number of peaks: {}".format(n_peaks))
    print("Omega binsize: {}".format(peak_selector.omega_binsize))
    print("Number of omega bins: {}".format(peak_selector.omega_bins.size))

    idx, ydist = peak_selector.mask(xi0, yi0, ystep, y0)

    # What we expect to get is this:
    expected_idx, expected_ydist = brute_force_mask(
        xi0,
        yi0,
        y0,
        ystep,
        peak_selector.sinomega,
        peak_selector.cosomega,
        peak_selector.dty,
    )
    assert np.allclose(np.sort(idx), np.sort(expected_idx))
    assert np.allclose(np.sort(ydist), np.sort(expected_ydist))

    # this should look like a sinogram of a point-source
    import matplotlib.pyplot as plt

    fontsize = 22
    ticksize = 22
    plt.rcParams["font.size"] = fontsize
    plt.rcParams["xtick.labelsize"] = ticksize
    plt.rcParams["ytick.labelsize"] = ticksize
    plt.rcParams["font.family"] = "Times New Roman"
    plt.style.use("dark_background")
    fig, ax = plt.subplots(1, 1, figsize=(13, 7))
    ax.scatter(
        peak_selector.omega[idx],
        peak_selector.dty[idx],
        s=4,
        alpha=0.5,
        label="voxel peaks",
    )
    ax.set_ylim(dty.min(), dty.max())
    ax.axhline(y0, color="red", linestyle="--", linewidth=2, label="y0")
    ax.legend(fontsize=fontsize)
    ax.set_xlim(0, 360)
    ax.set_xlabel("omega")
    ax.set_ylabel("dty")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.grid(True, alpha=0.25)

    # and now for some benchmark testing...

    ys = unique_dty[0 :: len(unique_dty) // 20]
    import time

    # early warmup
    idx, ydist = peak_selector.mask(xi0, yi0, ystep, y0)
    t1 = time.perf_counter()
    for i in range(len(ys)):
        for j in range(len(ys)):
            idx, ydist = peak_selector.mask(ys[i], ys[j], ystep, y0)
    t2 = time.perf_counter()
    time_per_call_voxel_mask = (t2 - t1) / (len(ys) * len(ys))
    print(
        "Time per voxel_mask.get_voxel_idx() call: {}".format(time_per_call_voxel_mask)
    )

    # early warmup
    idx, ydist = brute_force_mask(
        xi0,
        yi0,
        y0,
        ystep,
        peak_selector.sinomega,
        peak_selector.cosomega,
        peak_selector.dty,
    )
    t1 = time.perf_counter()
    for i in range(len(ys)):
        for j in range(len(ys)):
            idx, ydist = brute_force_mask(
                ys[i],
                ys[j],
                y0,
                ystep,
                peak_selector.sinomega,
                peak_selector.cosomega,
                peak_selector.dty,
            )
    t2 = time.perf_counter()

    time_per_call_brute_force = (t2 - t1) / (len(ys) * len(ys))
    print(
        "Time per full-scan reference mask call: {}".format(time_per_call_brute_force)
    )

    print(
        "Speedup: {:.1f} x".format(time_per_call_brute_force / time_per_call_voxel_mask)
    )

    plt.show()