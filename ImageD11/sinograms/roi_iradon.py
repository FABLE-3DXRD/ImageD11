"""This code came from skimage.transform and then was modified

- Added mask ROI for back projection. Optimisation for small grains
  in big maps. Should allow threading (one thread per tile)
  
- Added projection shifts (sinogram shaped array of dty shifts). 
  Should work the same way for both forward and back projection. 
  Probably only for small offsets
  
"""

"""
License: BSD-3-Clause

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of the University nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE HOLDERS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

from scipy.fft import fft, ifft, fftfreq, fftshift
from scipy.interpolate import interp1d
import numpy as np
import concurrent.futures, os
import skimage.transform.radon_transform
import numba
from ImageD11 import cImageD11


def _sinogram_pad(n, o=None):
    if o is None:
        diagonal = int(np.ceil(np.sqrt(2) * n))
    else:
        diagonal = int(np.ceil(np.sqrt(2) * o))
    pad = diagonal - n
    old_center = n // 2
    new_center = diagonal // 2
    pad_before = new_center - old_center
    pad_width = ((pad_before, pad - pad_before), (0, 0))
    return pad_width


def _get_fourier_filter(size, filter_name):
    """Construct the Fourier filter.
    """
    n = np.concatenate((np.arange(1, size / 2 + 1, 2, dtype=int),
                        np.arange(size / 2 - 1, 0, -2, dtype=int)))
    f = np.zeros(size)
    f[0] = 0.25
    f[1::2] = -1 / (np.pi * n) ** 2

    # Computing the ramp filter from the fourier transform of its
    # frequency domain representation lessens artifacts and removes a
    # small bias as explained in [1], Chap 3. Equation 61
    fourier_filter = 2 * np.real(fft(f))  # ramp filter
    if filter_name == "ramp":
        pass
    elif filter_name == "shepp-logan":
        # Start from first element to avoid divide by zero
        omega = np.pi * fftfreq(size)[1:]
        fourier_filter[1:] *= np.sin(omega) / omega
    elif filter_name == "cosine":
        freq = np.linspace(0, np.pi, size, endpoint=False)
        cosine_filter = fftshift(np.sin(freq))
        fourier_filter *= cosine_filter
    elif filter_name == "hamming":
        fourier_filter *= fftshift(np.hamming(size))
    elif filter_name == "hann":
        fourier_filter *= fftshift(np.hanning(size))
    elif filter_name is None:
        fourier_filter[:] = 1
    return fourier_filter[:, np.newaxis]


# skimage.transform.iradon
def iradon(radon_image,
           theta,
           output_size=None,
           filter_name="ramp",
           interpolation="linear",
           projection_shifts=None,
           mask=None,
           workers=1):
    """Inverse radon transform. From skimage.transform. Simplified then ruined.

    - allow projection offset/shifts to be used
        1D = constant offset for this projection (why?)
        2D = offset versus dty for this projection
    - allow a mask to be used for roi reconstruction
    - run over several threads (each one doing a different roi)
    """
    angles_count = len(theta)
    if angles_count != radon_image.shape[1]:
        raise ValueError("The given ``theta`` does not match the number of "
                         "projections in ``radon_image``.")
    if output_size is None:
        output_size = radon_image.shape[0]
    to_pad = _sinogram_pad(radon_image.shape[0], output_size)
    if projection_shifts is not None:
        assert projection_shifts.shape == radon_image.shape
        projection_shifts = np.pad(projection_shifts, to_pad,
                                   mode='constant', constant_values=0)
    radon_image = np.pad(radon_image, to_pad,
                         mode='constant', constant_values=0)
    img_shape = radon_image.shape[0]

    # Resize image to next power of two (but no less than 64) for
    # Fourier analysis; speeds up Fourier and lessens artifacts
    projection_size_padded = max(64, int(2 ** np.ceil(np.log2(2 * img_shape))))
    pad_width = ((0, projection_size_padded - img_shape), (0, 0))
    img = np.pad(radon_image, pad_width, mode='constant', constant_values=0)
    # return img
    # Apply filter in Fourier domain
    if filter_name is not None:
        fourier_filter = _get_fourier_filter(projection_size_padded, filter_name)
        projection = fft(img, axis=0, workers=workers) * fourier_filter
        radon_filtered = np.real(ifft(projection, axis=0, workers=workers)[:img_shape, :])
    else:
        radon_filtered = radon_image

    # Reconstruct image by interpolation
    reconstructed = np.zeros((output_size, output_size),
                             dtype=radon_image.dtype)
    radius = output_size // 2
    xpr, ypr = np.mgrid[:output_size, :output_size] - radius

    if mask is not None:
        xpr = xpr[mask]
        ypr = ypr[mask]
        recm = reconstructed[mask]
    else:
        recm = reconstructed

    x = np.arange(img_shape) - img_shape // 2
    rtheta = np.deg2rad(theta)

    # interpolate data in parallel
    # if we are linear interpolating, np.interp is much faster

    if interpolation == "linear":
        def run_interp(idx):
            result = 0
            for i in idx:
                t = ypr * np.cos(rtheta[i]) - xpr * np.sin(rtheta[i])
                if projection_shifts is not None:
                    xi = x + projection_shifts.T[i]
                else:
                    xi = x
                result += np.interp(t, xi, radon_filtered[:, i], left=0, right=0)
            return result
    else:
        def run_interp(idx):
            result = 0
            for i in idx:
                t = ypr * np.cos(rtheta[i]) - xpr * np.sin(rtheta[i])
                if projection_shifts is not None:
                    xi = x + projection_shifts.T[i]
                else:
                    xi = x
                interpolant = interp1d(xi, radon_filtered[:, i],
                                       kind=interpolation,
                                       copy=False,
                                       assume_sorted=True,
                                       bounds_error=False,
                                       fill_value=0)
                result += interpolant(t)
            return result

    if workers == 1:
        recm += run_interp(range(angles_count))
    else:
        if workers is None or workers < 1:
            workers = cImageD11.cores_available()

        # apply interpolants in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            todo = list(range(angles_count))
            jobs = [ todo[j::workers] for j in range(workers) ]
            for interp_result in pool.map(run_interp, jobs):
                recm += interp_result

    recm *= np.pi / (2 * angles_count)

    if mask is not None:
        reconstructed[mask] = recm
    return reconstructed


# TODO : fixme
#   It would be 'nice' for the radon transform to be closer to iradon
#  in using the same xpr/ypr/shifts formula. That will give pixel co-ords
#  for the 1D projections. The problem is to then 'cut up' the pixels on 
#  a 1D projection (so trapezium integrations). Doing bilinear interpolation
#  gave nasty artifacts at 45 degrees angle. This needs a bit of thought.
#
#  4 corners. 
#    Sort into order(?) or just abs(sin), abs(cos)?
#    First triangle. Middle Part. Last triangle.
#

def radon(image, theta,
          output_size=None,  # sinogram width
          projection_shifts=None,
          mask=None,
          workers=1,
          ):
    """
    From skimage.transform. Modified to have projection shifts and roi
    to match the masked iradon here. And spread theta over threads.
    
    Calculates the radon transform of an image given specified
    projection angles.

    Parameters
    ----------
    image : array_like
        Input image. The rotation axis will be located in the pixel with
        indices ``(image.shape[0] // 2, image.shape[1] // 2)``.
    theta : array_like
        Projection angles (in degrees).
    output_size : shape of the output sinogram. Matches [:, len(theta)]

    Returns
    -------
    radon_image : ndarray
        Radon transform (sinogram).  The tomography rotation axis will lie
        at the pixel index ``radon_image.shape[0] // 2`` along the 0th
        dimension of ``radon_image``.

    References
    ----------
    .. [1] AC Kak, M Slaney, "Principles of Computerized Tomographic
           Imaging", IEEE Press 1988.
    .. [2] B.R. Ramesh, N. Srinivasa, K. Rajgopal, "An Algorithm for Computing
           the Discrete Radon Transform With Some Applications", Proceedings of
           the Fourth IEEE Region 10 International Conference, TENCON '89, 1989

    Notes
    -----
    Based on code of Justin K. Romberg
    (https://www.clear.rice.edu/elec431/projects96/DSP/bpanalysis.html)

    """
    if image.ndim != 2:
        raise ValueError('The input image must be 2-D')

    assert image.dtype == np.float32
    assert len(image.shape) == 2
    assert image.shape[0] == image.shape[1]

    if output_size is None:
        output_size = image.shape[0]

    if projection_shifts is not None:
        assert projection_shifts.shape[1] == len(theta)
        assert projection_shifts.shape[0] == output_size

    radius = output_size // 2
    # padding the image. Shall we bother? Apparently yes.
    pad = [int(np.ceil(output_size - s)) for s in image.shape]
    new_center = [(s + p) // 2 for s, p in zip(image.shape, pad)]
    old_center = [s // 2 for s in image.shape]
    pad_before = [nc - oc for oc, nc in zip(old_center, new_center)]
    pad_width = [(pb, p - pb) for pb, p in zip(pad_before, pad)]
    padded_image = np.pad(image, pad_width, mode='constant',
                          constant_values=0)
    # padded_image is always square
    if padded_image.shape[0] != padded_image.shape[1]:
        raise ValueError('padded_image must be a square')
    center = padded_image.shape[0] // 2
    radon_image = np.zeros((len(theta), padded_image.shape[0]),
                           dtype=image.dtype).T

    angles_count = len(theta)
    rtheta = np.deg2rad(theta)

    def roti(i):
        if projection_shifts is not None:
            dx = projection_shifts.T[i]  # measured positions are shifted
        else:
            dx = None
        rotated = skimage.transform.radon_transform.warp(padded_image, fxyrot,
                                                         map_args={'angle': rtheta[i],
                                                                   'center': center,
                                                                   'projection_shifts': dx},
                                                         clip=False)
        radon_image[:, i] = rotated.sum(0)

    slices = list(range(angles_count))
    if workers == 1:
        for i in slices:
            roti(i)
        return radon_image
    if workers is None or workers < 1:
        workers = cImageD11.cores_available()
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for _ in pool.map(roti, slices):
            pass
    return radon_image


def fxyrot(colrow, angle=0, center=0, projection_shifts=None):
    # used in radon above
    # apply the projection shifts in reverse
    # t = ypr * np.cos(rtheta[i]) - xpr * np.sin(rtheta[i])
    # if projection_shifts is not None:
    #    xi = x + projection_shifts.T[i] # measured positions are shifted
    # else:
    #    xi = x
    col, row = colrow.T - center
    n = int(np.sqrt(col.shape[0]))
    assert n * n == col.shape[0]
    col.shape = n, n
    row.shape = n, n
    cos_a, sin_a = np.cos(angle), np.sin(angle)
    if projection_shifts is not None:
        ct = col.T
        ct += projection_shifts
    x = cos_a * col + sin_a * row
    y = -sin_a * col + cos_a * row
    colrow[:, 0] = x.ravel()
    colrow[:, 1] = y.ravel()
    return colrow + center


@numba.njit(cache=True)
def recon_cens(omega, dty, ybins, imsize, wt, y0=0.0):
    """ Back project the peak centers into a map
    
    omega, dty = peak co-ordinates in the sinogram (degrees)
    ybins = spatial binning. Usually dataset.ybincens
    imsize = probably len(ybins)+1
    wt = intensity, probably ones()

    see test/silicon/check_centers.ipynb
    """
    r = np.zeros((imsize, imsize), dtype=np.dtype(np.float32))
    rc = imsize // 2
    for i in range(len(omega)):
        o = np.radians(omega[i])
        s, c = np.sin(o), np.cos(o)
        # Offset of this peak in dty in units of ystep
        yv = (dty[i] - y0) / (ybins[1] - ybins[0])
        if abs(c) > abs(s):
            # cos(omega) ~ 1 means x axis vertical in the image
            # going down image, beam is along x
            for p in range(imsize): # walk up/down image
                k = p - rc  # -rc -> rc, distance to center
                v = ((-yv - k * s) / c) + rc # needs a picture
                j = int(np.floor(v))
                f = (j + 1) - v
                if (j > 0) & (j <= imsize):
                    r[p, imsize - j] += wt[i] * f
                f = v - j
                if ((j + 1) > 0) & ((j + 1) <= imsize):
                    r[p, imsize-(j + 1)] += wt[i] * f
        else:
            # going along image
            for p in range(1,imsize+1):
                j = p - rc  # -rc -> rc
                v = ((-yv - j * c) / s) + rc
                k = int(np.floor(v))
                f = (k + 1) - v
                if (k >= 0) & (k < imsize):
                    r[k, imsize - p] += wt[i] * f
                f = v - k
                if ((k + 1) >= 0) & ((k + 1) < imsize):
                    r[k + 1, imsize - p] += wt[i] * f
    return r


def mlem(sino,
         theta,
         startvalue=1,
         projection_shifts=None,
         output_size=None,
         mask=None,
         workers=1,
         niter=50,
         divtol=1e-5):
    """
    # Also called "MART" for Multiplicative ART
    # This keeps a positivity constraint for both the data and reconstruction
    #
    # This implementation was inspired from from:
    # https://www.youtube.com/watch?v=IhETD4nSJec
    # by Andrew Reader
    #
    # ToDo : implement a mask
    # ToDo : check about the corners / circle=False aspects
    #

    An "MLEM" algorithm from XRDUA was used in this paper:

    "Impurity precipitation in atomized particles evidenced by nano x-ray diffraction computed tomography"
    Anne Bonnin; Jonathan P. Wright; Rémi Tucoulou; Hervé Palancher
    Appl. Phys. Lett. 105, 084103 (2014) https://doi.org/10.1063/1.4894009

    This python code implements something similar based on a youtube video (https://www.youtube.com/watch?v=IhETD4nSJec)

    There are lots of papers from mathematicians in the literature about MART (multiplicative ART).
    The conversion of latex algebra back and forth into computer code seems to be a bit of a
    problem for me (Jon Wright - Nov 2023).
    """

    #
    # Also called "MART" for Multiplicative ART
    # This keeps a positivity constraint for both the data and reconstruction
    #
    # This implementation was inspired from from:
    # https://www.youtube.com/watch?v=IhETD4nSJec
    # by Andrew Reader
    #

    # output_size here is the shape of the output in real space

    if output_size is None:
        output_size = sino.shape[0]

    sensitivity_image = iradon(np.ones_like(sino),
                               theta,
                               output_size=output_size,
                               mask=None,
                               filter_name=None,
                               workers=workers,
                               )

    recip_sensitivity_image = 1. / sensitivity_image
    # The image reconstruction:
    mlem_rec = np.empty(sensitivity_image.shape, np.float32)
    mlem_rec[:] = startvalue

    # modify projection_shifts for radon
    # only need to pad [0]

    if projection_shifts is not None:
        pad = [int(np.ceil(output_size - s)) for s in projection_shifts.shape]
        new_center = [(s + p) // 2 for s, p in zip(projection_shifts.shape, pad)]
        old_center = [s // 2 for s in projection_shifts.shape]
        pad_before = [nc - oc for oc, nc in zip(old_center, new_center)]
        pad_width = [(pb, p - pb) for pb, p in zip(pad_before, pad)]
        pad_width[1] = (0, 0)

        proj_shifts_padded = np.pad(projection_shifts, pad_width, mode='constant', constant_values=0)

    else:
        proj_shifts_padded = None

    for i in range(niter):
        calc_sino = radon(mlem_rec,
                          theta,
                          mask=None,
                          projection_shifts=proj_shifts_padded,
                          workers=workers
                          )

        pad = [int(np.ceil(output_size - s)) for s in sino.shape]
        new_center = [(s + p) // 2 for s, p in zip(sino.shape, pad)]
        old_center = [s // 2 for s in sino.shape]
        pad_before = [nc - oc for oc, nc in zip(old_center, new_center)]
        pad_width = [(pb, p - pb) for pb, p in zip(pad_before, pad)]
        pad_width[1] = (0, 0)

        sino_padded = np.pad(sino, pad_width, mode='constant', constant_values=0)

        ratio = sino_padded / (calc_sino + divtol)

        correction = recip_sensitivity_image * iradon(ratio,
                                                      theta,
                                                      output_size=output_size,
                                                      mask=None,
                                                      projection_shifts=proj_shifts_padded,
                                                      filter_name=None,
                                                      workers=workers,
                                                      )

        mlem_rec *= correction

    if mask is not None:
        mlem_rec[~mask] = 0.

    return mlem_rec


def apply_halfmask_to_sino(sino):
    """Applies halfmask correction to sinogram"""
    halfmask = np.zeros_like(sino)

    halfmask[:len(halfmask) // 2 - 1, :] = 1
    halfmask[len(halfmask) // 2 - 1, :] = 0.5

    sino_halfmasked = sino.copy() * halfmask

    return sino_halfmasked


def correct_recon_central_zingers(recon, radius=25):
    recon_corrected = recon.copy()
    grs = recon.shape[0]
    xpr, ypr = -grs // 2 + np.mgrid[:grs, :grs]
    inner_mask_radius = radius
    outer_mask_radius = inner_mask_radius + 2

    inner_circle_mask = (xpr ** 2 + ypr ** 2) < inner_mask_radius ** 2
    outer_circle_mask = (xpr ** 2 + ypr ** 2) < outer_mask_radius ** 2

    mask_ring = inner_circle_mask & outer_circle_mask
    # we now have a mask to apply
    fill_value = np.median(recon_corrected[mask_ring])
    recon_corrected[inner_circle_mask] = fill_value

    return recon_corrected


def run_iradon(sino, angles, pad=20, shift=0,
               workers=1,
               mask=None,
               apply_halfmask=False,
               mask_central_zingers=False,
               central_mask_radius=25,
               filter_name='hamming',
               ):
    """Applies an iradon to a sinogram, with an optional pad
       Calculates scaled-up output size from pad value
       Applies projection shifts of shift
       Optionally mask the returned reconstruction with sample_mask
       Optionally apply halfmask to correct artifacts with half-scans (from -ny to 0 rather than -ny/2 to +ny/2
       Optionally apply a circular median filter mask to the center of the reconstruction given a radius in central_mask_radius"""

    # TODO Move apply_halfmask method to sinogram class

    outsize = sino.shape[0] + pad

    if apply_halfmask:
        sino_to_recon = apply_halfmask_to_sino(sino)
    else:
        sino_to_recon = sino

    # Perform iradon transform of grain sinogram, store result (reconstructed grain shape) in g.recon
    recon = iradon(sino_to_recon,
                   theta=angles,
                   mask=mask,
                   output_size=outsize,
                   projection_shifts=np.full(sino.shape, shift),
                   filter_name=filter_name,
                   interpolation='linear',
                   workers=workers)

    if mask_central_zingers:
        recon = correct_recon_central_zingers(recon, radius=central_mask_radius)

    return recon


def run_mlem(sino, angles, mask=None, pad=20, shift=0, workers=1, niter=20, apply_halfmask=False,
             mask_central_zingers=False, central_mask_radius=25):
    """Applies an MLEM reconstruction to a sinogram, with an optional pad
       Calculates scaled-up output size from pad value
       Applies projection shifts of shift
       Optionally mask the returned reconstruction with sample_mask
       Optionally apply halfmask to correct artifacts with half-scans (from -ny to 0 rather than -ny/2 to +ny/2
       Optionally apply a circular median filter mask to the center of the reconstruction given a radius in central_mask_radius"""

    outsize = sino.shape[0] + pad

    if apply_halfmask:
        sino_to_recon = apply_halfmask_to_sino(sino)
    else:
        sino_to_recon = sino

    # Perform iradon transform of grain sinogram, store result (reconstructed grain shape) in g.recon
    recon = mlem(sino_to_recon,
                 theta=angles,
                 mask=mask,
                 workers=workers,
                 output_size=outsize,
                 projection_shifts=np.full(sino_to_recon.shape, shift),
                 niter=niter)

    if mask_central_zingers:
        recon = correct_recon_central_zingers(recon, radius=central_mask_radius)

    return recon
