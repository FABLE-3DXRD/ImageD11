import numba
import numpy as np
from ImageD11 import sparseframe

@numba.njit
def select(img, mask, row, col, val, cut):
    # TODO: This is in now cImageD11.tosparse_{u16|f32}
    # Choose the pixels that are > cut and put into sparse arrays
    k = 0
    for s in range(img.shape[0]):
        for f in range(img.shape[1]):
            if img[s, f] * mask[s, f] > cut:
                row[k] = s
                col[k] = f
                val[k] = img[s, f]
                k += 1
    return k


@numba.njit
def top_pixels(nnz, row, col, val, howmany, thresholds):
    """
    selects the strongest pixels from a sparse collection
    - THRESHOLDS should be a sorted array of potential cutoff values to try
    that are higher than the original cutoff used to select data
    - howmany is the maximum number of pixels to return
    """
    # quick return if there are already few enough pixels
    if nnz <= howmany:
        return nnz
    # histogram of how many pixels are above each threshold
    h = np.zeros(len(thresholds), dtype=np.uint32)
    for k in range(nnz):
        for i, t in enumerate(thresholds):
            if val[k] > t:
                h[i] += 1
            else:
                break
    # choose the one to use. This is the first that is lower than howmany
    tcut = thresholds[-1]
    for n, t in zip(h, thresholds):
        if n < howmany:
            tcut = t
            break
    # now we filter the pixels
    n = 0
    for k in range(nnz):
        if val[k] > tcut:
            row[n] = row[k]
            col[n] = col[k]
            val[n] = val[k]
            n += 1
            if n >= howmany:
                break
    return n


class frmtosparse:
    def __init__(self, mask, dtype):
        # cache the mallocs on this function. Should be one per process
        self.row = np.empty(mask.size, np.uint16)
        self.col = np.empty(mask.size, np.uint16)
        self.val = np.empty(mask.size, dtype)
        self.mask = mask

    def __call__(self, frm, cut):
        nnz = select(frm, self.mask, self.row, self.col, self.val, cut)
        return nnz, self.row[:nnz], self.col[:nnz], self.val[:nnz]


def clean_utils(nnz, row, col, val, howmany, thresholds, mask_shape, pixels_in_spot):
    #flake8: global OPTIONS
    if nnz == 0:
        return None
    if nnz > howmany:
        nnz = top_pixels(nnz, row, col, val, howmany, thresholds)
        # Now get rid of the single pixel 'peaks'
        #   (for the mallocs, data is copied here)
        s = sparseframe.sparse_frame(
            row[:nnz].copy(), col[:nnz].copy(), mask_shape
        )
        s.set_pixels("intensity", val[:nnz].copy())
    else:
        s = sparseframe.sparse_frame(row, col, mask_shape)
        s.set_pixels("intensity", val)
    if pixels_in_spot <= 1:
        return s
    # label them according to the connected objects
    s.set_pixels("f32", s.pixels["intensity"].astype(np.float32))
    npk = sparseframe.sparse_connected_pixels(
        s, threshold=0, data_name="f32", label_name="cp"
    )
    # only keep spots with more than 3 pixels ...
    #    mom = sparseframe.sparse_moments( s,
    #                                     intensity_name="f32",
    #                                     labels_name="cp" )
    #    npx = mom[:, cImageD11.s2D_1]
    npx = np.bincount(s.pixels["cp"], minlength=npk)
    pxcounts = npx[s.pixels["cp"]]
    pxmsk = pxcounts >= pixels_in_spot
    if pxmsk.sum() == 0:
        return None
    sf = s.mask(pxmsk)
    return sf
