"""Sinogram-related geometric functions pulled from various places.
In this file there are three reference frames:

1. The sample frame

2. Step space:
   step from the rotation axis center,  coords are (si, sj)


3. Reconstruction space:
   units are (ystep), origin in corner (matches iradon output) when plotted with origin="lower", coords are (ri, rj)


Diagrams are below, S indicates the centre of the sample

sample frame:

         ^
         |
         | x
         |
<------- S (0, 0)
    y

Step space:

                +ve
                /\
                |
              i |
                |
-ve <--------   S (0, 0)--------> +ve
                |           j
                |
                |
               \/
              -ve


Reconstruction space (iradon output when plotted with origin="lower"):

   ^
   |
 i |      S
   |
(0, 0) ------->
         j


"""

import numpy as np
from scipy.optimize import curve_fit
from skimage.feature import blob_log


def sample_to_step(x, y, ystep):
    """Converts sample (x, y) position to step space (si, sj)"""
    si = x / ystep
    sj = -y / ystep

    return si, sj


def step_to_sample(si, sj, ystep):
    """Converts step space (si, sj) to sample position (x, y)"""
    x = si * ystep
    y = -sj * ystep

    return x, y


def step_to_recon(si, sj, recon_shape):
    """Converts step space (si, sj) to reconstruction space (ri, rj)"""
    ri = si + (recon_shape[0] // 2)
    rj = sj + (recon_shape[1] // 2)

    return ri, rj


def recon_to_step(ri, rj, recon_shape):
    """Converts reconstruction space (ri, rj) to step space (si, sj)"""
    si = ri - (recon_shape[0] // 2)
    sj = rj - (recon_shape[1] // 2)

    return si, sj


def sample_to_recon(x, y, recon_shape, ystep):
    """Converts sample space (x, y) to reconstruction space (ri, rj)"""
    si, sj = sample_to_step(x, y, ystep)
    ri, rj = step_to_recon(si, sj, recon_shape)

    return ri, rj


def recon_to_sample(ri, rj, recon_shape, ystep):
    """Converts reconstruction space (ri, rj) to sample space (x, y)"""
    si, sj = recon_to_step(ri, rj, recon_shape)
    x, y = step_to_sample(si, sj, ystep)

    return x, y


def x_y_y0_omega_to_dty(omega, x, y, y0):
    """This is the dty position where the grain should be in the middle of the beam
       For a given omega value
       A grain on the centre of rotation has a = b = 0
       To get the grain into the beam, dty should equal y0
       So the offset equals y0.
       Phased sine function with vertical offset using linear combination of sin and cos
       A grain at (0, +y), you need to move dty negative to compensate
       At omega = 0, cos(0) = 1
       So cosine term is multiplied by -y
       A grain at (+x, 0) when omega = +90, sin(90) = 1
       So b goes to -x"""
    return y0 - x * np.sin(np.radians(omega)) - y * np.cos(np.radians(omega))


def fit_sine_wave(omega, dty, initial_guess):
    """Fits a sine wave to omega and dty data
    Returns directly x, y, y0 in sample frame"""
    popt, _ = curve_fit(x_y_y0_omega_to_dty,
                        omega,
                        dty,
                        p0=initial_guess,
                        method="trf",
                        loss="soft_l1",
                        max_nfev=10000,
                        )

    x, y, y0 = popt

    return x, y, y0


def dty_omega_to_x_y_y0(dty, omega):
    """Fits sine wave to dty vs omega plot, extracts x, y, y0"""
    initial_guess = (0.5, 0.5, 0)

    x, y, y0 = fit_sine_wave(omega, dty, initial_guess)

    return x, y, y0


def dty_to_dtyi(dty, ystep):
    """Converts dty value (lab frame) to step space, then rounds to integer
       Note that this will invert the sign of the values"""
    _, dty_step = sample_to_step(0, dty, ystep)
    dtyi = np.round(dty_step).astype(int)
    return dtyi


def dty_to_dtyi_for_sinogram(dty, ystep, ymin):
    """Converts dty value (lab frame) to something like scan number"""
    dtyi = np.round((dty - ymin) / ystep).astype(int)
    return dtyi


def step_omega_to_dty(si, sj, omega, y0, ystep):
    """Convert step value (si, sj) to (x, y)
       Determine corresponding dty values from (x, y) given omega"""
    x, y = step_to_sample(si, sj, ystep)
    dty = x_y_y0_omega_to_dty(omega, x, y, y0)
    return dty


def step_omega_to_dtyi(si, sj, omega, y0, ystep):
    """Convert step value (si, sj) to (x, y)
       Determine ccorresponding dty values from (x, y) given omega
       Compute dtyi using dty_to_dtyi then return"""
    x, y = step_to_sample(si, sj, ystep)
    dty = x_y_y0_omega_to_dty(omega, x, y, y0)
    dtyi = dty_to_dtyi(dty, ystep)
    return dtyi


def recon_omega_to_dty(ri, rj, omega, y0, recon_shape, ystep):
    """Convert recon value (ri, rj) to (x, y)
       Determine ccorresponding dty values from (x, y) given omega"""
    x, y = recon_to_sample(ri, rj, recon_shape, ystep)
    dty = x_y_y0_omega_to_dty(omega, x, y, y0)
    return dty


def recon_omega_to_dtyi(ri, rj, omega, y0, recon_shape, ystep):
    """Convert recon value (ri, rj) to (x, y)
       Determine ccorresponding dty values from (x, y) given omega
       Compute dtyi using dty_to_dtyi then return"""
    x, y = recon_to_sample(ri, rj, recon_shape, ystep)
    dty = x_y_y0_omega_to_dty(omega, x, y, y0)
    dtyi = dty_to_dtyi(dty, ystep)
    return dtyi


def dtyimask_from_step(si, sj, omega, dtyi, y0, ystep):
    """Convert step value (si, sj) to (x, y)
       Determine ccorresponding dty values from (x, y) given omega
       Compute dtyi_calc using dty_to_dtyi
       Compare dtyi_calc to supplied dtyi and produce a mask"""
    dtyi_calc = step_omega_to_dtyi(si, sj, omega, y0, ystep)
    return dtyi == dtyi_calc


def dtyimask_from_recon(ri, rj, omega, dtyi, y0, ystep, recon_shape):
    """Convert recon value (ri, rj) to (x, y)
       Determine ccorresponding dty values from (x, y) given omega
       Compute dtyi_calc using dty_to_dtyi
       Compare dtyi_calc to supplied dtyi and produce a mask"""
    dtyi_calc = recon_omega_to_dtyi(ri, rj, omega, y0, recon_shape, ystep)
    return dtyi == dtyi_calc


def fit_sample_position_from_recon(recon, ystep):
    """Fits grain position by doing a LoG search with skimage on the reconstruction image
    Useful if grain is very small (close to pointlike)
    Returns position in sample frame (unit matches ystep)
    Returns None if we couldn't find any blobs"""
    blobs = blob_log(recon, min_sigma=1, max_sigma=10, num_sigma=10, threshold=0.01)
    blobs_sorted = sorted(blobs, key=lambda x: x[2], reverse=True)
    try:
        largest_blob = blobs_sorted[0]

        ri, rj, sigma = largest_blob
        # we now have the blob position in recon space
        # we need to go back to sample space
        x, y = recon_to_sample(ri, rj, recon.shape, ystep)

        return x, y
    except IndexError:
        # didn't find any blobs
        # if we didn't find a blob, normally indicates recon is bad
        return None
