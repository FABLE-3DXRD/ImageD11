"""Sinogram-related geometric functions pulled from various places"""
# TODO: Needs testing/validating!


import numpy as np
from scipy.optimize import curve_fit
from skimage.feature import blob_log


def dtymask(i, j, cosomega, sinomega, dtyi, step, y0):
    """
    Selects the peaks using integer bins in dty.
    Needs a UI / plot to see what is happening (map the i,j to the sinogram and recon)
    """
    # TODO: Not sure if this is correct (what if sinogram is padded?)
    dty_icalc = np.round((step * (i * cosomega + j * sinomega) + y0) / step).astype(int)
    return dtyi == dty_icalc


def recon_space_to_real_space(i, j, recon_shape, ystep, y0):
    """Convert pixel position in reconstruction space (i, j) to real space
    i, j are indices to the reconstruction image
    recon_shape is a tuple of the shape of the reconstruction array
    ystep is the dty step (microns or mm per pixel, dty motor unit)"""
    x = (i - recon_shape[0] // 2) * ystep + y0
    y = -(j - recon_shape[1] // 2) * ystep + y0

    return x, y


def real_space_to_recon_space(x, y, recon_shape, ystep, y0):
    """Convert real space (x, y) position to reconstruction space.
    Should be the exact opposite of recon_space_to_real_space()
       Accounts for the shape of the reconstruction because the origin in reconstruction space is in the corner of the image."""
    i = recon_shape[0] // 2 + (x - y0) / ystep
    j = -(recon_shape[1] // 2 + (y - y0) / ystep)

    return i, j


def sine_function(omega, offset, a, b):
    """Phased sine function with vertical offset using linear combination of sin and cos"""
    return b * np.sin(np.radians(omega)) + a * np.cos(np.radians(omega)) + offset


def fit_sine_wave(omega, dty, initial_guess):
    """Fits a sine wave to omega and dty data
    Returns sine wave coefficients and offset
    To convert coefficients into lab positions, use sine_coeffs_to_lab_position()"""
    popt, _ = curve_fit(
        sine_function,
        omega,
        dty,
        p0=initial_guess,
        method="trf",
        loss="soft_l1",
        max_nfev=10000,
    )

    offset, a, b = popt

    return offset, a, b


def sine_coeffs_to_lab_position(offset, a, b):
    """Converts coefficients used in sine_function() to lab coordinates"""
    x = -b
    y = -a
    cen = offset

    return cen, x, y


def lab_position_to_sine_coeffs(cen, x, y):
    """The opposite of sine_coeffs_to_lab_position"""
    b = -x
    a = -y
    offset = cen

    return offset, a, b


def fit_lab_position_from_peaks(omega, dty):
    """Fits omega and dty data to a sine wave, converts result into lab coordinates"""
    initial_guess = (0, 0.5, 0.5)

    offset, a, b = fit_sine_wave(omega, dty, initial_guess)

    cen, x, y = sine_coeffs_to_lab_position(offset, a, b)

    return cen, x, y


def fit_lab_position_from_recon(recon, ystep, y0):
    """Fits grain position by doing a LoG search with skimage on the reconstruction image
    Useful if grain is very small (close to pointlike)
    Returns position in lab frame (unit matches ystep)
    Returns None if we couldn't find any blobs"""
    blobs = blob_log(recon, min_sigma=1, max_sigma=10, num_sigma=10, threshold=0.01)
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

        x, y = recon_space_to_real_space(
            x_recon_space, y_recon_space, recon.shape, ystep, y0
        )

        return x, y
    except IndexError:
        # didn't find any blobs
        # if we didn't find a blob, normally indicates recon is bad
        return None
