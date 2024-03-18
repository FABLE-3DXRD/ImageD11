import numpy as np

import ImageD11.grain
import ImageD11.sinograms.dataset
import ImageD11.cImageD11
import ImageD11.sinograms.roi_iradon


class GrainSinogram:
    """Class to build, hold and reconstruct sinograms of grains for Scanning 3DXRD data"""

    def __init__(self, grain_obj, dataset):
        if not isinstance(grain_obj, ImageD11.grain.grain):
            raise TypeError("grain_obj must be an ImageD11.grain.grain instance")
        if not isinstance(dataset, ImageD11.sinograms.dataset.DataSet):
            raise TypeError("dataset must be an ImageD11.sinograms.dataset.DataSet instance")

        self.grain = grain_obj
        self.ds = dataset

        # Peak information to build the sinogram
        self.etasigns_2d_strong = None
        self.hkl_2d_strong = None

        # Sinogram information
        self.sinoangles = None
        self.sino = None
        self.ssino = None  # sorted by angle

        # Reconstruction information
        # self.recons is a dict that stores reconstructions by method
        # i.e. self.recons["iradon"] is the iradon reconstruction
        self.recons = {}
        self.recon_mask = None
        self.recon_pad = None
        self.recon_y0 = None

    def prepare_peaks_for_sinogram(self, peak_indices_2d, hkltol=0.25):
        """Filters the 2D peaks assigned to the grain based on the HKL tolerance"""

        # Filter the peaks table to the peaks for this grain only
        p2d = {p: self.ds.pk2d[p][peak_indices_2d] for p in self.ds.pk2d}

        # Make a spatially corrected columnfile from the filtered peaks table
        flt = self.ds.get_colfile_from_peaks_dict(peaks_dict=p2d)

        hkl_real = np.dot(self.grain.ubi, (flt.gx, flt.gy, flt.gz))  # calculate hkl of all assigned peaks
        hkl_int = np.round(hkl_real).astype(int)  # round to nearest integer
        dh = ((hkl_real - hkl_int) ** 2).sum(axis=0)  # calculate square of difference

        flt.filter(dh < hkltol * hkltol)  # filter all assigned peaks to be less than hkltol squared
        hkl_real = np.dot(self.grain.ubi, (flt.gx, flt.gy, flt.gz))  # recalculate error after filtration
        hkl_int = np.round(hkl_real).astype(int)

        self.etasigns_2d_strong = np.sign(flt.eta)
        self.hkl_2d_strong = hkl_int  # integer hkl of assigned peaks after hkltol filtering
        self.flt = flt

    def build_sinogram(self):
        """
        Computes sinogram from peaks data
        Flt is the filtered assigned 2D peaks for this grain

        """
        NY = len(self.ds.ybincens)  # number of y translations
        iy = np.round((self.flt.dty - self.ds.ybincens[0]) / (self.ds.ybincens[1] - self.ds.ybincens[0])).astype(
            int)  # flt column for y translation index

        # The problem is to assign each spot to a place in the sinogram
        hklmin = self.hkl_2d_strong.min(axis=1)  # Get minimum integer hkl (e.g -10, -9, -10)
        dh = self.hkl_2d_strong - hklmin[:, np.newaxis]  # subtract minimum hkl from all integer hkls
        de = (self.etasigns_2d_strong.astype(int) + 1) // 2  # something signs related
        #   4D array of h,k,l,+/-
        # pkmsk is whether a peak has been observed with this HKL or not
        pkmsk = np.zeros(list(dh.max(axis=1) + 1) + [2, ],
                         int)  # make zeros-array the size of (max dh +1) and add another axis of length 2
        pkmsk[dh[0], dh[1], dh[2], de] = 1  # we found these HKLs for this grain
        #   sinogram row to hit
        pkrow = np.cumsum(pkmsk.ravel()).reshape(pkmsk.shape) - 1  #
        # counting where we hit an HKL position with a found peak
        # e.g (-10, -9, -10) didn't get hit, but the next one did, so increment

        npks = pkmsk.sum()
        destRow = pkrow[dh[0], dh[1], dh[2], de]
        sino = np.zeros((npks, NY), 'f')
        hits = np.zeros((npks, NY), 'f')
        angs = np.zeros((npks, NY), 'f')
        adr = destRow * NY + iy
        # Just accumulate
        sig = self.flt.sum_intensity
        ImageD11.cImageD11.put_incr64(sino, adr, sig)
        ImageD11.cImageD11.put_incr64(hits, adr, np.ones(len(de), dtype='f'))
        ImageD11.cImageD11.put_incr64(angs, adr, self.flt.omega)

        sinoangles = angs.sum(axis=1) / hits.sum(axis=1)
        # Normalise:
        self.sino = (sino.T / sino.max(axis=1)).T
        # Sort (cosmetic):
        order = np.lexsort((np.arange(npks), sinoangles))
        self.sinoangles = sinoangles[order]
        self.ssino = self.sino[order].T

    def correct_halfmask(self):
        """Applies halfmask correction to sinogram"""
        self.ssino = ImageD11.sinograms.roi_iradon.apply_halfmask_to_sino(self.ssino)

    def correct_ring_current(self):
        """Corrects each row of the sinogram to the ring current of the corresponding scan"""
        self.ssino = self.ssino / self.ds.ring_currents_per_scan_scaled[:, None]

    def update_recon_parameters(self, pad, y0, mask):
        """Update all reconstruction parameters in one go"""
        self.recon_pad = pad
        self.recon_y0 = y0
        self.recon_mask = mask

    def recon(self, method="iradon", workers=1):
        """Performs reconstruction given reconstruction method"""
        if method not in ["iradon", "mlem"]:
            raise ValueError("Unsupported method!")

        sino = self.ssino
        angles = self.sinoangles

        if method == "iradon":
            recon_function = ImageD11.sinograms.roi_iradon.run_iradon
        elif method == "mlem":
            recon_function = ImageD11.sinograms.roi_iradon.run_mlem

        recon = recon_function(sino=sino,
                               angles=angles,
                               pad=self.recon_pad,
                               y0=self.recon_y0,
                               workers=workers,
                               mask=self.recon_mask)

        self.recons[method] = recon

    def mask_central_zingers(self, method="iradon", radius=25):
        self.recons[method] = ImageD11.sinograms.roi_iradon.correct_recon_central_zingers(self.recons[method],
                                                                                          radius=radius)

    def to_h5py_group(self, parent_group, group_name):
        """Creates a H5Py group for this GrainSinogram.
           parent_group is the parent H5py Group
           group_name is the name of the H5py Group for this name
           Very useful for saving lists of GrainSinograms to an H5 file"""

        # create a group for this specific GrainSinogram
        grain_group = parent_group.require_group(group_name)

        # save peak information
        peak_info_group = grain_group.require_group("peak_info")
        for peak_info_attr in ["etasigns_2d_strong", "hkl_2d_strong"]:
            peak_info_var = getattr(self, peak_info_attr)
            if peak_info_var is not None:
                save_array(peak_info_group, peak_info_attr, peak_info_var)

        # save sinograms
        sinogram_group = grain_group.require_group("sinograms")

        for sino_attr in ["sino", "ssino", "sinoangles"]:
            sino_var = getattr(self, sino_attr)
            if sino_var is not None:
                save_array(sinogram_group, sino_attr, sino_var)

        # save reconstruction parameters

        recon_par_group = grain_group.require_group("recon_parameters")

        for recon_par_attr in ["recon_mask", "recon_pad", "recon_y0"]:
            recon_par_var = getattr(self, recon_par_attr)
            if recon_par_var is not None:
                recon_par_group[recon_par_attr] = recon_par_var

        # save reconstructions

        recon_group = grain_group.require_group("recons")

        for recon_attr in self.recons.keys():
            save_array(recon_group, recon_attr, self.recons[recon_attr])

        return grain_group

    @classmethod
    def from_h5py_group(cls, group, ds, grain):
        """Creates a GrainSinogram object from an h5py group, dataset and grain object"""
        grainsino_obj = GrainSinogram(grain_obj=grain, dataset=ds)

        if "peak_info" in group.keys():
            for peak_info_attr in ["etasigns_2d_strong", "hkl_2d_strong"]:
                peak_info_var = group["peak_info"].get(peak_info_attr)[:]
                setattr(grainsino_obj, peak_info_attr, peak_info_var)

        if "recon_parameters" in group.keys():
            for recon_par_attr in ["recon_mask", "recon_pad", "recon_y0"]:
                recon_par_var = group["recon_parameters"].get(recon_par_attr)[()]
                setattr(grainsino_obj, recon_par_attr, recon_par_var)

        if "recons" in group.keys():
            for recon_attr in group["recons"].keys():
                recon_var = group["recons"].get(recon_attr)[:]
                grainsino_obj.recons[recon_attr] = recon_var

        return grainsino_obj


# TODO: Use Silx Nexus IO instead?
# This is a temporary sensible middle ground!
def save_array(grp, name, ary):
    # TODO: Move this helper function somewhere else
    cmp = {'compression': 'gzip',
           'compression_opts': 2,
           'shuffle': True}

    hds = grp.require_dataset(name,
                              shape=ary.shape,
                              dtype=ary.dtype,
                              **cmp)
    hds[:] = ary
    return hds
