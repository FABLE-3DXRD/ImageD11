"""Written by Axel Henningsson, 2026."""

import numba
import numpy as np


@numba.njit(cache=True)
def fill_voxel_idx(
    xi0,
    yi0,
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
    idx_buffer,  # you will want to reuse these for all voxels!
    ydist_buffer,  # you will want to reuse these for all voxels
):
    """Get the voxel indices and y-distances for a given voxel centroid position.

    This tracks the voxel in sinogram space and returns all peaks the voxel possibly intersects.

    The idea is to use a partitioning scheme that allows for retrival of these peaks without touching
    the entire sinogram.

    Args:
        xi0 (:obj:`float`): The x-coordinate of the voxel position in sample coordinates.
        yi0 (:obj:`float`): The y-coordinate of the voxel position in sample coordinates.
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
        :obj:`tuple`: A tuple containing the voxel indices and y distances.
        The voxel indices are the indices of the peaks that are within the voxela and
        refers to peaks that follow the same ordering as dty_sorted.
        The y distances are the absolute distances of the peaks from the voxel centroid position. These are guaranteed to be
        less than or equal to ystep.

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

    y_prev = y0 - xi0 * sin_omega_bins[0] - yi0 * cos_omega_bins[0]

    for i in range(n_bins - 1):  # every omega partition bin...
        y_next = y0 - xi0 * sin_omega_bins[i + 1] - yi0 * cos_omega_bins[i + 1]
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
        y = y0 - xi0 * sin_omega_bins[i] - yi0 * cos_omega_bins[i]

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
        if padded_low < 0:
            padded_low = 0
        padded_high = y_index + dty_partition_padding + 1
        if padded_high > row_len_minus_1:
            padded_high = row_len_minus_1

        # if there is simply no scan close to the voxel for this omega bin, then skip it completely.
        # this corre
        if padded_high <= 0:
            continue

        if padded_low >= row_len_minus_1:
            continue

        # get the safe start and end of the dty sub-partition for this omega bin
        dty_partition_start = dty_partitions[i, padded_low]
        dty_partition_end = dty_partitions[i, padded_high]

        # get the start of the omega partition
        omega_partition_start = omega_partitions[i]

        # iterate over the feasible dty chunk in global indices
        for j in range(
            omega_partition_start + dty_partition_start,
            omega_partition_start + dty_partition_end,
        ):  # for each peak in the dty-sub-partition...
            # Now we compute the exact y-distance to the voxel centroid
            dty_val = dty_sorted[j]
            s = sin_omega_sorted[j]
            c = cos_omega_sorted[j]
            ydist = abs(y0 - xi0 * s - yi0 * c - dty_val)

            # Check if the peak is within ystep of the voxel centroid
            # if so, then keep it.
            if ydist <= ystep:
                if m >= buffer_size:
                    return -1
                idx_buffer[m] = j  # these now refer to the sorted peaks!
                ydist_buffer[m] = ydist
                m += 1

    return m


@numba.njit(cache=True)
def get_voxel_idx(
    xi0,
    yi0,
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
        xi0, yi0, y0, ystep, ymin,
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
    xi0s,
    yi0s,
    y0,
    ystep,
    ymin,
    omega_partitions,
    dty_partitions,
    sin_omega_bins,
    cos_omega_bins,
):
    """Upper bound on how many peaks any of these voxels can pull out.

    Walks the same omega bins and dty sub-partitions as fill_voxel_idx() but
    reads only the partition offsets, never the peak data. Use it to size the
    buffers up front so fill_voxel_idx() can never signal an overflow, which
    is what you want when calling it from a prange.
    """
    n_bins = sin_omega_bins.size
    row_len_minus_1 = dty_partitions.shape[1] - 1
    inv_ystep = 1.0 / ystep
    best = np.zeros(xi0s.size, dtype=np.int64)
    for v in numba.prange(xi0s.size):
        xi0 = xi0s[v]
        yi0 = yi0s[v]
        total = 0
        y_prev = y0 - xi0 * sin_omega_bins[0] - yi0 * cos_omega_bins[0]
        pad = 1
        for i in range(n_bins):
            if i + 1 < n_bins:
                y_next = y0 - xi0 * sin_omega_bins[i + 1] - yi0 * cos_omega_bins[i + 1]
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

    Adaptable omega bin sizes are used to chunk the data over angles. For a call to mask() the algorithm then
    oterates over the omega bins and collects partitions of dty values that are candidates for the voxel.

    Args:
        omega (:obj:`np.ndarray`): The omega values, one per peak, shape=(n_peaks,). (degrees)
        dty (:obj:`np.ndarray`): The dty values, one per peak, shape=(n_peaks,).
        dty_stepsize (:obj:`float`): The dty step size, assumed to be constant such that the
            scan defines a regular grid in dty.

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

    def __init__(self, omega, dty, dty_stepsize):
        self.uniq_dtys = np.unique(dty)
        self.ymin = dty.min()

        float_yindex = (dty - self.ymin) / dty_stepsize
        delta_yindex_stepsize = np.abs(float_yindex - np.round(float_yindex))
        if not np.less_equal(delta_yindex_stepsize, 0.01).all():
            raise ValueError(
                "dty_stepsize must be a divisor of the unique dty values, however it \
                appears you have scanned an irregular grid in dty since there is \
                more than 1 percent drift in the dty values - this is not supported!"
            )

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
                  omega_bin_indices=None, omega_bins=None):
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
        However, the defualt heuristic should be very performant for most datasets.

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

        Raises:
            :obj:`ValueError`: If the omega binsize is not set and the heuristic fails to find a good value.
        """
        if omega_bin_indices is not None:
            # caller supplied the binning, e.g. exact frames decoded from a
            # frame number column
            self.omega_binsize = None
            omega_bin_indices = np.ascontiguousarray(omega_bin_indices, dtype=np.int64)
            omega_bin_indices = omega_bin_indices - omega_bin_indices.min()
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

        peak_ordering = np.lexsort((self.dty, omega_bin_indices))
        omega_bin_indices_sorted = omega_bin_indices[peak_ordering]
        omega_partition_lengths = np.bincount(omega_bin_indices_sorted)
        omega_partitions = np.zeros(omega_bins.size + 1, dtype=np.int64)
        omega_partitions[1:] = np.cumsum(omega_partition_lengths)

        # get also the dty partitions
        dty_sorted = self.dty[peak_ordering]
        dty_bin_indices = np.round((dty_sorted - self.ymin) / self.dty_stepsize).astype(
            np.int64
        )
        ymax_index = dty_bin_indices.max()

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

        if inplace:
            # saves memory, but overwrites the arrays the caller handed us
            self.dty[:] = self.dty[peak_ordering]
            self.omega[:] = self.omega[peak_ordering]
        else:
            self.dty = self.dty[peak_ordering]
            self.omega = self.omega[peak_ordering]
        self.sinomega = np.sin(np.radians(self.omega))
        self.cosomega = np.cos(np.radians(self.omega))

        self.sinomega_bins = np.sin(np.radians(omega_bins))
        self.cosomega_bins = np.cos(np.radians(omega_bins))

        # we default to a 5% buffer, this should be very safe, altough
        # we have fallbacks in case of a sharp corner.
        self._build_buffers(self.n_peaks // 20)

    def sort_by_partitions(self, peaks):
        """In place sorting of peak columns, such as sc, fc, sum_intensity, etc. by partitions

        This method is to be called once after a desired partition has been set. The peak
        order is then compatible with the partition indices, and can be used to mask peaks
        fast with the class mask method for single thread test and with get_voxel_idx()
        integrated code.

        Args:
            peaks (:obj:`np.ndarray` | :obj:`list` | :obj:`dict`): Peak columns to sort, such as sc, fc, sum_intensity, etc. shape=(n_peaks,).
            these can be either a single column (numpy array) or a multi-column (list or dict) of columns.
        """
        if self.peak_ordering is None:
            raise ValueError(
                "Peak ordering is not set, please call partition() first to set the partition indices."
            )
        if isinstance(peaks, np.ndarray):
            peaks[:] = peaks[self.peak_ordering]
        elif isinstance(peaks, list):
            for value in peaks:
                self.sort_by_partitions(value)
        elif isinstance(peaks, dict):
            for key, value in peaks.items():
                self.sort_by_partitions(value)
        else:
            raise ValueError("Unsupported type: {}".format(type(peaks)))

    def max_candidates(self, xi0s, yi0s, ystep, y0):
        """Largest number of peaks any of these voxels can pull out of the partition.

        Size the buffers with this and fill_voxel_idx() can never overflow, which
        matters when calling it from a prange where the -1 return is the only
        signal you get.
        """
        return int(
            max_candidates(
                np.ascontiguousarray(xi0s, dtype=np.float64),
                np.ascontiguousarray(yi0s, dtype=np.float64),
                y0,
                ystep,
                self.ymin,
                self.omega_partitions,
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

    peak_selector = VoxelSinoMasker(omega, dty, dty_stepsize)
    peak_selector.partition()

    print("Number of peaks: {}".format(n_peaks))
    print("Omega binsize: {}".format(peak_selector.omega_binsize))
    print("Number of omega bins: {}".format(peak_selector.omega_bins))

    idx, ydist = peak_selector.mask(xi0, yi0, ystep, y0)

    # What we expect to get is this:
    import ImageD11.sinograms.point_by_point as pbp

    expected_idx, expected_ydist = pbp.get_voxel_idx(
        y0,
        xi0,
        yi0,
        peak_selector.sinomega,
        peak_selector.cosomega,
        peak_selector.dty,
        ystep,
    )
    expected_ydist = expected_ydist[expected_idx]
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
    idx, ydist = pbp.get_voxel_idx(
        y0,
        xi0,
        yi0,
        peak_selector.sinomega,
        peak_selector.cosomega,
        peak_selector.dty,
        ystep,
    )
    t1 = time.perf_counter()
    for i in range(len(ys)):
        for j in range(len(ys)):
            idx, ydist = pbp.get_voxel_idx(
                y0,
                xi0,
                yi0,
                peak_selector.sinomega,
                peak_selector.cosomega,
                peak_selector.dty,
                ystep,
            )
    t2 = time.perf_counter()

    time_per_call_pbp = (t2 - t1) / (len(ys) * len(ys))
    print(
        "Time per ImageD11.sinograms.point_by_point.get_voxel_idx call: {}".format(
            time_per_call_pbp
        )
    )

    print("Speedup: {:.1f} x".format(time_per_call_pbp / time_per_call_voxel_mask))

    plt.show()