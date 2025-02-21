"""Sinogram-related geometric functions pulled from various places.
The choice of reference frames is tricky.
In regular 3DXRD, the origin of the lab reference frame is defined as the intersection between the middle of the beam
and the rotation axis.
However, because the rotation axis is moving, this gets more difficult for scanning 3DXRD.
I have therefore defined the reference frames in the following way:

1. Static lab frame - vectors defined by the beam vector and the dty translation axis vector.
   Its origin is defined directly as (0,0,0):
   x: The intersection of the beam vector and the dty translation axis vector
   y: The centre of the beam horizontally
   z: The centre of the beam vertically
   dtyi is a discretisation of the dty motor values.

2. Sample frame
   Origin is the rotation axis
   The rotation axis is translated by dty.
   The sample frame rotates with omega
   This translates with dty and rotates with omega, CW about the rotation axis (looking top down).
   This has the rotation axis at its origin.
   y0 is the true value of dty when the rotation axis intersects the beam.

3. Step space
   This is simply an integer discretisation of the sample frame.


4. Reconstruction space:
   units are (ystep), origin in corner (matches iradon output) when plotted with origin="lower", coords are (ri, rj)


Diagrams are below, S indicates the rotation axis

lab frame:

         ^
         |
         | x
         |
<------- O (0, 0) centre of beam
    y


sample frame (could be rotated about omega then translated along lab y by dty):

         ^
         |
         | x
         |
<------- S (0, 0)
    y

Step space:

                +ve
                ^
                |
              i |
                |
-ve <--------   S (0, 0)--------> +ve
                |           j
                |
                |
                v
              -ve


Reconstruction space (iradon output when plotted with origin="lower"):

   ^
   |
 i |      S
   |
(0, 0) ------->
         j


"""
from functools import partial

import numpy as np
from scipy.optimize import curve_fit
from skimage.feature import blob_log


def sample_to_lab_sincos(sx, sy, y0, dty, sinomega, cosomega):
    """
    Converts position in sample frame (sx, sy) to position in lab frame (lx, ly).
    The units of sx, sy, y0, and dty must agree.

    :param sx: X-coordinate in sample reference frame
    :type sx: (float, np.ndarray)
    :param sy: Y-coordinate in sample reference frame
    :type sy: (float, np.ndarray)
    :param y0: the true value of dty when the rotation axis intersects the beam
    :type y0: float
    :param dty: the dty motor underneath the rotation axis
    :type dty: (float, np.ndarray)
    :param sinomega: the sine of the omega value of the rotation axis
    :type sinomega: (float, np.ndarray)
    :param cosomega: the cosine of the omega value of the rotation axis
    :type cosomega: (float, np.ndarray)
    :return: lx, ly: the x and y coordinates in the lab reference frame
    :rtype: ((float, np.ndarray), (float, np.ndarray))
    """
    # first, rotate sx and sy by omega about its origin (rotation axis)
    sxr = sx * cosomega - sy * sinomega
    syr = sx * sinomega + sy * cosomega

    # then translate about the rotated frame
    lx = sxr
    ly = syr + dty - y0

    return lx, ly


def sample_to_lab(sx, sy, y0, dty, omega):
    """Calls sample_to_lab_sincos after converting omega (degrees) into sinomega, cosomega"""

    omega_rad = np.radians(omega)
    sinomega = np.sin(omega_rad)
    cosomega = np.cos(omega_rad)
    lx, ly = sample_to_lab_sincos(sx, sy, y0, dty, sinomega, cosomega)

    return lx, ly


def lab_to_sample_sincos(lx, ly, y0, dty, sinomega, cosomega):
    """
    Converts position in lab frame (lx, ly) to position in lab frame (sx, sy).
    The units of sx, sy, y0, and dty must agree.

    :param lx: X-coordinate in lab reference frame
    :type lx: (float, np.ndarray)
    :param ly: Y-coordinate in lab reference frame
    :type ly: (float, np.ndarray)
    :param y0: the true value of dty when the rotation axis intersects the beam
    :type y0: float
    :param dty: the dty motor underneath the rotation axis
    :type dty: (float, np.ndarray)
    :param sinomega: the sine of the omega value of the rotation axis
    :type sinomega: (float, np.ndarray)
    :param cosomega: the cosine of the omega value of the rotation axis
    :type cosomega: (float, np.ndarray)
    :return: sx, su: the x and y coordinates in the sample reference frame
    :rtype: ((float, np.ndarray), (float, np.ndarray))
    """
    # first, translate
    sxr = lx
    syr = ly - dty + y0

    # then unrotate
    sx = sxr * cosomega + syr * sinomega
    sy = -sxr * sinomega + syr * cosomega

    return sx, sy


def lab_to_sample(lx, ly, y0, dty, omega):
    """Calls lab_to_sample_sincos after converting omega (degrees) into sinomega, cosomega"""

    omega_rad = np.radians(omega)
    sinomega = np.sin(omega_rad)
    cosomega = np.cos(omega_rad)
    sx, sy = lab_to_sample_sincos(lx, ly, y0, dty, sinomega, cosomega)

    return sx, sy


def sample_to_step(sx, sy, ystep):
    """Converts sample (sx, sy) position to step space (si, sj)"""
    si = sx / ystep
    sj = -sy / ystep
    return si, sj


def step_to_sample(si, sj, ystep):
    """Converts step space (si, sj) to sample position (sx, sy)"""
    sx = si * ystep
    sy = -sj * ystep
    return sx, sy


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


def sample_to_recon(sx, sy, recon_shape, ystep):
    """Converts sample space (sx, sy) to reconstruction space (ri, rj)"""
    si, sj = sample_to_step(sx, sy, ystep)
    ri, rj = step_to_recon(si, sj, recon_shape)
    return ri, rj


def recon_to_sample(ri, rj, recon_shape, ystep):
    """Converts reconstruction space (ri, rj) to sample space (sx, sy)"""
    si, sj = recon_to_step(ri, rj, recon_shape)
    sx, sy = step_to_sample(si, sj, ystep)
    return sx, sy


def lab_to_step(lx, ly, y0, dty, omega, ystep):
    """Converts lab space (lx, ly) to step space (si, sj)"""
    sx, sy = lab_to_sample(lx, ly, y0, dty, omega)
    si, sj = sample_to_step(sx, sy, ystep)
    return si, sj


def step_to_lab(si, sj, y0, dty, omega, ystep):
    """Converts step space (si, sj) to lab space (lx, ly)"""
    sx, sy = step_to_sample(si, sj, ystep)
    lx, ly = sample_to_lab(sx, sy, y0, dty, omega)
    return lx, ly


def lab_to_recon(lx, ly, y0, dty, omega, recon_shape, ystep):
    """Converts lab space (lx, ly) to recon space (ri, rj)"""
    si, sj = lab_to_step(lx, ly, y0, dty, omega, ystep)
    ri, rj = step_to_recon(si, sj, recon_shape)
    return ri, rj


def recon_to_lab(ri, rj, y0, dty, omega, recon_shape, ystep):
    """Converts recon space (ri, rj) to lab space (lx, ly)"""
    si, sj = recon_to_step(ri, rj, recon_shape)
    lx, ly = step_to_lab(si, sj, y0, dty, omega, ystep)
    return lx, ly


def dty_values_grain_in_beam_sincos(sx, sy, y0, sinomega, cosomega):
    """
    Take a grain positioned at (sx, sy) in the sample reference frame
    Determine the dty values needed to bring the grain into the beam as you rotate by omega.
    Equivalently, solve the sample<->lab conversion where ly = 0 (grain is in-beam).
    From sample_to_lab_sincos:

    sxr = sx * cosomega - sy * sinomega
    syr = sx * sinomega + sy * cosomega
    
    lx = sxr
    ly = syr + dty - y0
    
    Therefore:
    ly = sx * sinomega + sy * cosomega + dty - y0
    
    Solving for ly = 0:
    0 = sx * sinomega + sy * cosomega + dty - y0
    dty = y0 - sx * sinomega - sy * cosomega
    """
    dty = y0 - sx * sinomega - sy * cosomega
    return dty


def dty_values_grain_in_beam(sx, sy, y0, omega):
    """Calls dty_values_grain_in_beam_sincos after converting omega (degrees) into sinomega, cosomega"""
    omega_rad = np.radians(omega)
    sinomega = np.sin(omega_rad)
    cosomega = np.cos(omega_rad)
    return dty_values_grain_in_beam_sincos(sx, sy, y0, sinomega, cosomega)


def x_y_y0_omega_to_dty(omega, x, y, y0):
    return dty_values_grain_in_beam(sx=x, sy=y, y0=y0, omega=omega)


dtycalc = x_y_y0_omega_to_dty


def fit_sine_wave(omega, dty, initial_guess, weights=None):
    """
    Take a series of (omega, dty) values that originate from the same point in the sample.
    Determine the point in sample coordinates (sx, sy) by fitting a sine wave to the data.
    Includes a possible shift of the rotation axis in y (y0).
    In practice, we perform a curve_fit on dty_values_grain_in_beam to determine sx, sy, y0
    """
    partial_function = partial(lambda _om, _sx, _sy, _y0: dty_values_grain_in_beam(_sx, _sy, _y0, _om))

    popt, _ = curve_fit(
        partial_function,
        omega,
        dty,
        p0=initial_guess,
        method="trf",
        loss="soft_l1",
        max_nfev=10000,
        sigma=weights,
    )

    sx_out, sy_out, y0_out = popt

    return sx_out, sy_out, y0_out


def sx_sy_y0_from_dty_omega(dty, omega, weights=None):
    """Fits sine wave to dty vs omega plot, extracts sx, sy, y0"""
    initial_guess = (0.5, 0.5, 0)

    sx, sy, y0 = fit_sine_wave(omega, dty, initial_guess, weights=weights)

    return sx, sy, y0


def dty_to_dtyi(dty, ystep, ymin):
    """
    We take continuous dty space and discretise it
    We do this by counting the number of ysteps away from ymin
    ymin should probably be ds.ymin
    This should be equivalent to determining the closest bin in ds.ybincens
    But this is neater as it accounts for values outside the ds.ybincens range

    dtyi = (dty - ymin) / ystep
    """
    dtyi = np.round((dty - ymin) / ystep).astype(int)
    return dtyi


def dtyi_to_dty(dtyi, ystep, ymin):
    """
    Simply:
    dty = (dtyi * ystep) + ymin
    """
    dty = dtyi * ystep + ymin
    return dty


def step_omega_to_dty(si, sj, omega, y0, ystep):
    """
    Convert step space (si, sj) to sample space (sx, sy)
    Then get corresponding dty values which puts (sx, sy) into the beam given omega
    """
    sx, sy = step_to_sample(si, sj, ystep)
    dty = dty_values_grain_in_beam(sx, sy, y0, omega)
    return dty


def step_omega_to_dtyi(si, sj, omega, y0, ystep, ymin):
    """
    Convert step space (si, sj) to sample space (sx, sy)
    Then get corresponding dty values which puts (sx, sy) into the beam given omega
    Then converts dty to dtyi
    """
    sx, sy = step_to_sample(si, sj, ystep)
    dty = dty_values_grain_in_beam(sx, sy, y0, omega)
    dtyi = dty_to_dtyi(dty, ystep, ymin)
    return dtyi


def recon_omega_to_dty(ri, rj, omega, y0, recon_shape, ystep):
    """
    Convert recon space (ri, rj) to step space (si, sj)
    Then get corresponding dty values which puts (si, sj) into the beam given omega
    """
    sx, sy = recon_to_sample(ri, rj, recon_shape, ystep)
    dty = dty_values_grain_in_beam(sx, sy, y0, omega)
    return dty


def recon_omega_to_dtyi(ri, rj, omega, y0, recon_shape, ystep, ymin):
    """
    Convert recon space (ri, rj) to step space (si, sj)
    Then get corresponding dtyi values which puts (si, sj) into the beam given omega
    Then converts dty to dtyi
    """
    sx, sy = recon_to_sample(ri, rj, recon_shape, ystep)
    dty = dty_values_grain_in_beam(sx, sy, y0, omega)
    dtyi = dty_to_dtyi(dty, ystep, ymin)
    return dtyi


def dtyimask_from_sample(sx, sy, omega, dtyi, y0, ystep, ymin):
    """
    Given a position in the sample (sx, sy) and arrays of (omega, dtyi) values
    Mask each (omega, dtyi) value if they put the position (sx, sy) into the beam
    """
    dty_calc = dty_values_grain_in_beam(sx, sy, y0, omega)
    dtyi_calc = dty_to_dtyi(dty_calc, ystep, ymin)
    mask = dtyi_calc == dtyi
    return mask


def dtyimask_from_sample_sincos(sx, sy, sinomega, cosomega, dtyi, y0, ystep, ymin):
    """
    Given a position in the sample (sx, sy) and arrays of (sinomega, cosomega, dtyi) values
    Mask each (sinomega, cosomega, dtyi) value if they put the position (sx, sy) into the beam
    """
    dty_calc = dty_values_grain_in_beam_sincos(sx, sy, y0, sinomega, cosomega)
    dtyi_calc = dty_to_dtyi(dty_calc, ystep, ymin)
    mask = dtyi_calc == dtyi
    return mask


def dtyimask_from_step(si, sj, omega, dtyi, y0, ystep, ymin):
    """
    Given a position in step space (si, sj) and arrays of (omega, dtyi) values
    Convert step space (si, sj) to sample space (sx, sy)
    Mask each (omega, dtyi) value if they put the position (sx, sy) into the beam
    """
    sx, sy = step_to_sample(si, sj, ystep)
    return dtyimask_from_sample(sx, sy, omega, dtyi, y0, ystep, ymin)


def dtyimask_from_step_sincos(si, sj, sinomega, cosomega, dtyi, y0, ystep, ymin):
    """
    Given a position in step space (si, sj) and arrays of (sinomega, cosomega, dtyi) values
    Convert step space (si, sj) to sample space (sx, sy)
    Mask each (sinomega, cosomega, dtyi) value if they put the position (sx, sy) into the beam
    """
    sx, sy = step_to_sample(si, sj, ystep)
    return dtyimask_from_sample_sincos(sx, sy, sinomega, cosomega, dtyi, y0, ystep, ymin)


def dtyimask_from_recon(ri, rj, omega, dtyi, y0, ystep, ymin, recon_shape):
    """
    Given a position in recon space (ri, rj) and arrays of (omega, dtyi) values
    Convert recon space (ri, rj) to sample space (sx, sy)
    Mask each (omega, dtyi) value if they put the position (sx, sy) into the beam
    """
    sx, sy = recon_to_sample(ri, rj, recon_shape, ystep)
    return dtyimask_from_sample(sx, sy, omega, dtyi, y0, ystep, ymin)


def dtyimask_from_recon_sincos(ri, rj, sinomega, cosomega, dtyi, y0, ystep, ymin, recon_shape):
    """
    Given a position in recon space (ri, rj) and arrays of (sinomega, cosomega, dtyi) values
    Convert recon space (ri, rj) to sample space (sx, sy)
    Mask each (sinomega, cosomega, dtyi) value if they put the position (sx, sy) into the beam
    """
    sx, sy = recon_to_sample(ri, rj, recon_shape, ystep)
    return dtyimask_from_sample_sincos(sx, sy, sinomega, cosomega, dtyi, y0, ystep, ymin)


def fit_sample_position_from_recon(recon, ystep):
    """
    Fits grain position in sample space by doing a LoG search with skimage on the reconstruction image
    Useful if grain is very small (close to pointlike)
    Returns position in sample frame (unit matches ystep)
    Returns None if we couldn't find any blobs
    """
    blobs = blob_log(recon, min_sigma=1, max_sigma=10, num_sigma=10, threshold=0.01)
    blobs_sorted = sorted(blobs, key=lambda x: x[2], reverse=False)
    try:
        largest_blob = blobs_sorted[0]

        ri, rj, sigma = largest_blob
        # we now have the blob position in recon space
        # we need to go back to sample space
        sx, sy = recon_to_sample(ri, rj, recon.shape, ystep)

        return sx, sy
    except IndexError:
        # didn't find any blobs
        # if we didn't find a blob, normally indicates recon is bad
        return None


def sino_shift_and_pad(y0, ny, ymin, ystep):
    """
    Determine the difference in sinogram pixels between
    the central dtyi value of the sinogram and the dtyi value of the rotation axis (from y0)
    To do this, we convert y0 into dtyi without rounding
    Also determine the minimum pad to get the whole sample in the frame
    Which should be double the shift
    """
    # get the middle row of the sinogram in pixel space
    ymid_px = ny/2
    # get the value of y0 in pixel space
    # this is the same as dty_to_dtyi without rounding
    y0_px = (y0 - ymin)/ystep
    shift = ymid_px - y0_px
    pad = np.ceil(np.abs(shift) * 2).astype(int) + 1
    return shift, pad


def step_grid_from_ybincens(ybincens, step_size, gridstep, y0):
    """
    Produce a grid in (si, sj) step space to reconstruct on.
    Constructs a symmetric grid with the origin on y0
    Therefore the step space is truly aligned to the rotation axis.
    The size of the grid is determined by the extremes of ybincens.
    The grid will always be symmetric such that (0, 0) is the centre of an image formed from the grid
    :param ybincens: An evenly sampled dty space
    :type ybincens: np.ndarray
    :param step_size: the step size of the grid (same units as ybincens and y0)
    :type step_size: float
    :param gridstep: the sampling resolution of the grid (e.g 2 = [0, 2, 4, 6, ...])
    :type gridstep: int
    :return: (si, sj) for si in ints for sj in ints
    :rtype: tuple(int, int)
    """
    # centre y range on y0
    y_relative = ybincens - y0
    # find minimum and maximum values of y
    y_min = y_relative.min()
    y_max = y_relative.max()
    # work out which is larger, so our grid is always symmetric
    # i.e never go from (-65, 65) to (100, 100)
    # always go from (-100, -100) to (100, 100)
    y_largest = np.max(np.abs([y_min, y_max]))
    y_min = -y_largest
    y_max = y_largest
    y_min_int = np.floor(y_min/step_size).astype(int)
    y_max_int = np.ceil(y_max/step_size).astype(int)
    ints = range(y_min_int, y_max_int + 1, gridstep)
    step_points = [(si, sj) for si in ints for sj in ints]
    return step_points
