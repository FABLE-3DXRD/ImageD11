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
        def run_interp(i):
            t = ypr * np.cos(rtheta[i]) - xpr * np.sin(rtheta[i])
            if projection_shifts is not None:
                xi = x + projection_shifts.T[i]
            else:
                xi = x

            result = np.interp(t, xi, radon_filtered[:, i], left=0, right=0)
            return result
    else:
        def run_interp(i):
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
            result = interpolant(t)
            return result

    if workers == 1:
        for angle in range(angles_count):
            recm += run_interp(angle)
    else:
        if workers is None or workers < 1:
            workers = cImageD11.cores_available()

        # apply interpolants in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            for interp_result in pool.map(run_interp, range(angles_count)):
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


@numba.njit(boundscheck=False)
def recon_cens(omega, dty, ybins, imsize, wt, y0=0.0):
    """ Back project the peak centers into a map 
    
    omega, dty = peak co-ordinates in the sinogram
    ybins = spatial binning
    imsize = probably len(ybins)+1
    wt = intensity, probably ones()
    """
    r = np.zeros((imsize, imsize), dtype=np.float32)
    rc = imsize // 2
    for i in range(len(omega)):
        s, c = np.sin(omega[i]), np.cos(omega[i])
        yv = (dty[i] - y0) / (ybins[1] - ybins[0])
        if abs(c) > abs(s):
            # going across image, beam is along x
            for p in range(imsize):
                k = p - rc  # -rc -> rc
                v = ((-yv - k * s) / c) + rc
                j = int(np.floor(v))
                f = (j + 1) - v
                if (j >= 0) & (j < imsize):
                    r[j, p] += wt[i] * f
                f = v - j
                if ((j + 1) >= 0) & ((j + 1) < imsize):
                    r[j + 1, p] += wt[i] * f
        else:
            # going along image
            for p in range(imsize):
                j = p - rc  # -rc -> rc
                v = ((-yv - j * c) / s) + rc
                k = int(np.floor(v))
                f = (k + 1) - v
                if (k >= 0) & (k < imsize):
                    r[p, k] += wt[i] * f
                f = v - k
                if ((k + 1) >= 0) & ((k + 1) < imsize):
                    r[p, k + 1] += wt[i] * f
    return r


def mlem(sino,
         theta,
         startvalue=1,
         projection_shifts=None,
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

    def backproject(sino, theta, mask=mask, projection_shifts=None):
        """ project the sinogram into the sample """
        return iradon(sino,
                      theta,
                      filter_name=None,
                      mask=mask,
                      workers=workers,
                      projection_shifts=projection_shifts)

    def forwardproject(sample, theta, projection_shifts=None):
        """ project the sample into the experiment (sinogram) """
        return radon(sample,
                     theta,
                     workers=workers,
                     projection_shifts=projection_shifts)

    #
    # Also called "MART" for Multiplicative ART
    # This keeps a positivity constraint for both the data and reconstruction
    #
    # This implementation was inspired from from:
    # https://www.youtube.com/watch?v=IhETD4nSJec
    # by Andrew Reader
    #

    # ToDo : check about the corners / circle=False aspects
    #
    # Number of pixels hitting each output in the sample:
    sensitivity_image = backproject(np.ones_like(sino),
                                    theta,
                                    mask=None,
                                    projection_shifts=projection_shifts)
    recip_sensitivity_image = 1. / sensitivity_image
    # The image reconstruction:
    mlem_rec = np.empty(sensitivity_image.shape, np.float32)
    mlem_rec[:] = startvalue
    for i in range(niter):
        calc_sino = forwardproject(mlem_rec,
                                   theta,
                                   projection_shifts=projection_shifts)
        ratio = sino / (calc_sino + divtol)
        correction = recip_sensitivity_image * backproject(ratio,
                                                           theta,
                                                           projection_shifts=projection_shifts)
        mlem_rec *= correction

    if mask is not None:
        mlem_rec[~mask] = 0.

    return mlem_rec
