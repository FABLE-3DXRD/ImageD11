from functools import partial

import h5py
import numba
import numpy as np

import ImageD11.grain
import ImageD11.cImageD11
import ImageD11.columnfile
import ImageD11.sinograms.dataset
import ImageD11.sinograms.properties
import ImageD11.sinograms.roi_iradon
import ImageD11.sinograms.geometry


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
        self.cf_for_sino = None

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
        self.recon_shift = None
        self.recon_niter = None
        self.recon_y0 = None

    def prepare_peaks_from_2d(self, cf_2d, grain_label, hkltol=0.25):
        """Prepare peaks used for sinograms from 2D peaks data
           Performs greedy assignment of all 2d peaks to each grain"""

        # get all g-vectors from columnfile
        gv = np.transpose((cf_2d.gx, cf_2d.gy, cf_2d.gz)).astype(float)

        # column to store the grain labels
        labels = np.zeros(cf_2d.nrows, 'i')

        # column to store drlv2 (error in hkl)
        drlv2 = np.ones(cf_2d.nrows, 'd')

        _ = ImageD11.cImageD11.score_and_assign(self.grain.ubi, gv, hkltol, drlv2, labels, grain_label)

        mask_2d = labels == grain_label

        # get needed columns in the flt file from the mask (without creating a full new FLT for this grain)
        dty = cf_2d.dty[mask_2d]
        omega = cf_2d.omega[mask_2d]
        gx = cf_2d.gx[mask_2d]
        gy = cf_2d.gy[mask_2d]
        gz = cf_2d.gz[mask_2d]
        eta = cf_2d.eta[mask_2d]
        sum_intensity = cf_2d.sum_intensity[mask_2d]

        cf_dict = {"dty": dty,
                   "omega": omega,
                   "gx": gx,
                   "gy": gy,
                   "gz": gz,
                   "eta": eta,
                   "sum_intensity": sum_intensity}

        grain_flt = ImageD11.columnfile.colfile_from_dict(cf_dict)

        self.cf_for_sino = grain_flt

    def prepare_peaks_from_4d(self, cf_4d, gord, inds, grain_label, hkltol=0.25):
        """Prepares peaks used for sinograms from 4D peaks data.
           cf_4d should contain grain assignments
           grain_label should be the label for this grain"""
        # fail if not already assigned
        if 'grain_id' not in cf_4d.titles:
            raise ValueError("cf_4d does not contain grain assignments!")

        grain_peaks_2d = gord[inds[grain_label + 1]:inds[grain_label + 2]]

        # Filter the peaks table to the peaks for this grain only
        p2d = {p: self.ds.pk2d[p][grain_peaks_2d] for p in self.ds.pk2d}

        # Make a spatially corrected columnfile from the filtered peaks table
        flt = self.ds.get_colfile_from_peaks_dict(peaks_dict=p2d)
        # Make sure the columnfile has the correct geometry and phase
        flt.parameters = cf_4d.parameters
        flt.updateGeometry()

        # Filter the columnfile by hkl tol
        hkl_real = np.dot(self.grain.ubi, (flt.gx, flt.gy, flt.gz))  # calculate hkl of all assigned peaks
        hkl_int = np.round(hkl_real).astype(int)  # round to nearest integer
        dh = ((hkl_real - hkl_int) ** 2).sum(axis=0)  # calculate square of difference

        flt.filter(dh < hkltol * hkltol)  # filter all assigned peaks to be less than hkltol squared

        self.cf_for_sino = flt

    def build_sinogram(self, columns=('omega',)):
        """
        Computes sinogram for this grain using all peaks in self.cf_for_sino

        columns = list of columns to produce sinograms of intensity*value
                for example: "ds", "eta", "omega", etc
        """
        for a in columns:
            assert a in self.cf_for_sino.titles
        assert 'omega' in columns

        NY = len(self.ds.ybincens)  # number of y translations
        iy = np.round(
            (self.cf_for_sino.dty - self.ds.ybincens[0]) / (self.ds.ybincens[1] - self.ds.ybincens[0])).astype(
            int)  # flt column for y translation index

        hkl = np.round(np.dot(self.grain.ubi, (self.cf_for_sino.gx, self.cf_for_sino.gy, self.cf_for_sino.gz))).astype(
            int)
        etasigns = np.sign(self.cf_for_sino.eta)

        # The problem is to assign each spot to a place in the sinogram
        hklmin = hkl.min(axis=1)  # Get minimum integer hkl (e.g -10, -9, -10)
        dh = hkl - hklmin[:, np.newaxis]  # subtract minimum hkl from all integer hkls
        de = (etasigns.astype(int) + 1) // 2  # something signs related
        #   4D array of h,k,l,+/-
        # pkmsk is whether a peak has been observed with this HKL or not
        pkmsk = np.zeros(list(dh.max(axis=1) + 1) + [2, ],
                         int)  # make zeros-array the size of (max dh +1) and add another axis of length 2
        pkmsk[dh[0], dh[1], dh[2], de] = 1  # we found these HKLs for this grain
        #   sinogram row to hit
        pkrow = np.cumsum(pkmsk.ravel()).reshape(pkmsk.shape) - 1  #
        # counting where we hit an HKL position with a found peak
        # e.g (-10, -9, -10) didn't get hit, but the next one did, so increment
        #
        # found peak indexes into the pkmsk array
        pkhkle = np.arange( np.prod( pkmsk.shape ), dtype=int )[ pkmsk.flat == 1 ]
        # hkl indices (transpose for 3,N versus N,3)
        pkindices = np.array( np.unravel_index(pkhkle, pkmsk.shape) )
        pkindices[:3] += hklmin[:,np.newaxis]
        npks = pkmsk.sum()
        destRow = pkrow[dh[0], dh[1], dh[2], de]
        sino = np.zeros((npks, NY), 'f')
        hits = np.zeros((npks, NY), 'f')
        adr = destRow * NY + iy
        # Just accumulate
        sig = self.cf_for_sino.sum_intensity
        ImageD11.cImageD11.put_incr64(sino, adr, sig)
        ImageD11.cImageD11.put_incr64(hits, adr, np.ones(len(de), dtype='f'))
        # intensity weighted sums
        angs = {}
        for name in columns:
            if name == 'sum_intensity':
                continue
            angs[name] = np.zeros((npks, NY), 'f')
            ImageD11.cImageD11.put_incr64(angs[name], adr, self.cf_for_sino[name] * sig )
        sinoangles = angs['omega'] .sum(axis=1) / sino.sum(axis=1)
        # Normalise:
        self.proj_scale = sino.max(axis=1)
        self.sino = sino / self.proj_scale[:, np.newaxis]
        # Sort (cosmetic):
        order = np.lexsort((np.arange(npks), sinoangles))
        self.sinoangles = sinoangles[order]
        self.ssino = self.sino[order].T
        self.proj_scale = self.proj_scale[order]
        if len(columns)>1:
            self.angle_wt_sinos = { name : angs[name][order].T
                                    for name in columns
                                    if name != 'sum_intensity'
                                   }
            if 'sum_intensity' in columns:
                self.angle_wt_sinos['sum_intensity'] = sino[order].T
        self.hits = hits
        self.hkle = pkindices[:, order]  # dims are [ (h,k,l,sign(eta)) , nprojections ]

    def update_lab_position_from_peaks(self, cf_4d, grain_label):
        """Updates translation of self.grain using peaks in assigned 4D colfile.
           Also updates self.recon_y0 with centre"""
        mask_4d = cf_4d.grain_id == grain_label
        omega = cf_4d.omega[mask_4d]
        dty = cf_4d.dty[mask_4d]
        sx, sy, y0 = ImageD11.sinograms.geometry.sx_sy_y0_from_dty_omega(dty, omega)
        self.grain.translation = np.array([sx, sy, 0])
        self.recon_y0 = y0

    def update_lab_position_from_recon(self, method="iradon"):
        """Updates translation of self.grain by finding centre-of-mass of reconstruction
           Only really valid for very small grains.
           Does not update translation if no fitting result found."""
        fitting_result = ImageD11.sinograms.geometry.fit_sample_position_from_recon(self.recons[method], self.ds.ystep)
        if fitting_result is None:
            # indicate that the grain has a dodgy reconstruction
            # useful for filtering out bad grains
            self.bad_recon = True
        else:
            x, y = fitting_result
            self.grain.translation = np.array([x, y, 0])

    def correct_halfmask(self):
        """Applies halfmask correction to sinogram"""
        self.ssino = ImageD11.sinograms.roi_iradon.apply_halfmask_to_sino(self.ssino)

    def correct_ring_current(self, is_half_scan=False, min_ring_current_frac=0.5):
        """Corrects each row of the sinogram to the ring current of the corresponding scan"""

        # ignore ring current values below min_ring_current
        mean_ring_current = np.mean(self.ds.ring_currents_per_scan_scaled)
        min_ring_current = min_ring_current_frac * np.max(self.ds.ring_currents_per_scan_scaled)
        ring_current_to_use = np.where(self.ds.ring_currents_per_scan_scaled < min_ring_current, mean_ring_current, self.ds.ring_currents_per_scan_scaled)

        if is_half_scan:
            correction = ring_current_to_use
            addition_length = len(self.ds.ybincens) - len(ring_current_to_use)
            correction_halfmask_addition = np.zeros(addition_length)
            correction_halfmask_addition.fill(correction.max())
            correction = np.concatenate((correction, correction_halfmask_addition))
        else:
            correction = ring_current_to_use
        self.ssino = self.ssino / correction[:, None]

    def update_recon_parameters(self, pad=None, shift=None, mask=None, niter=None, y0=None):
        """Update some or all of the reconstruction parameters in one go"""

        if pad is not None:
            self.recon_pad = pad
        if shift is not None:
            self.recon_shift = shift
        if mask is not None:
            self.recon_mask = mask
        if niter is not None:
            self.recon_niter = niter
        if y0 is not None:
            self.recon_y0 = y0

    def recon(self,
              method="iradon",
              workers=1,
              projections = None,
              **extra_args):
        """Performs reconstruction given reconstruction method

        method = "iradon", "mlem", "astra"
        workers = threads to use - passed to the method

        projections = which projections to use for reconstruction using a subset of data (twins)
                iradon( self.sino[ :, projections ], theta = self.sinoangles[ projections ], ... )
                effectively a sinogram mask - perhaps this should move to iradon code.

        extra_args = passed to the projection function

           iradon -> ImageD11.sinograms.roi_iradon.run_iradon
           { pad=20, shift=0, workers=1, mask=None,
                      apply_halfmask=False,
                      mask_central_zingers=False,
                      central_mask_radius=25,
                      filter_name='hamming'}
           mlem : ImageD11.sinograms.roi_iradon.run_mlem
              { mask=None, pad=20, shift=0, workers=1, niter=20, apply_halfmask=False,
             mask_central_zingers=False, central_mask_radius=25 }
           astra : run_astra in this file
               { shift=0, pad=0, mask=None, niter=100, astra_method='SIRT_CUDA', workers=None }
        """
        if method not in ["iradon", "mlem", "astra"]:
            raise ValueError("Unsupported method!")
        if self.ssino is None or self.sinoangles is None:
            raise ValueError("Sorted sino or sinoangles are missing/have not been computed, unable to reconstruct.")
        if projections is None:
            sino = self.ssino
            angles = self.sinoangles
        else:
            sino = self.ssino[:, projections]
            angles = self.sinoangles[ projections ]
            print("Using subset",sino.shape,"from", self.ssino.shape)
        if method == "iradon":
            recon_function = ImageD11.sinograms.roi_iradon.run_iradon
        elif method == "mlem":
            # MLEM has niter as an extra argument
            recon_function = ImageD11.sinograms.roi_iradon.run_mlem
            # Overwrite the default argument if self.recon_niter is set
            if self.recon_niter is not None and "niter" not in extra_args:
                extra_args[ "niter" ] = self.recon_niter
            # TODO: We should think about default reconstruction arguments in more detail
            # At the moment, we could pass None for pad or shift etc to recon_function if they are not manually set, which seems dangerous
            # Do we check for None at the start of this function?
            # Or do we initialise them to sensible default values inside __init__?
        elif method == "astra":
            recon_function = run_astra
        recon = recon_function(sino=sino,
                               angles=angles,
                               pad=self.recon_pad,
                               shift=self.recon_shift,
                               workers=workers,
                               mask=self.recon_mask,
                               **extra_args)

        self.recons[method] = recon
        return recon

    def get_shape_mask(self, method="iradon", cutoff=None):
        """Gets a boolean mask representing the grain shape from the grain reconstruction.
           Performs a binary threshold of the reconstruction.
           Uses manual threshold value if cutoff value provided.
           Performs global Otsu threshold if no value provided."""
        image = self.recons[method]
        if cutoff is None:
            threshold = threshold_otsu(image)
        else:
            threshold = cutoff
        binary_image = image > threshold
        return binary_image

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
        # peak_info_group = grain_group.require_group("peak_info")
        # for peak_info_attr in ["etasigns_2d_strong", "hkl_2d_strong"]:
        #     peak_info_var = getattr(self, peak_info_attr)
        #     if peak_info_var is not None:
        #         save_array(peak_info_group, peak_info_attr, peak_info_var)

        # save sinograms
        sinogram_group = grain_group.require_group("sinograms")

        for sino_attr in ["sino", "ssino", "sinoangles"]:
            sino_var = getattr(self, sino_attr)
            if sino_var is not None:
                save_array(sinogram_group, sino_attr, sino_var)

        # save reconstruction parameters as attributes

        recon_par_group = grain_group.require_group("recon_parameters")

        for recon_par_attr in ["recon_pad", "recon_shift", "recon_niter", "recon_y0"]:
            recon_par_var = getattr(self, recon_par_attr)
            if recon_par_var is not None:
                recon_par_group.attrs[recon_par_attr] = recon_par_var

        if self.recon_mask is not None:
            save_array(recon_par_group, "recon_mask", self.recon_mask)

        # save reconstructions

        recon_group = grain_group.require_group("recons")

        for recon_attr in self.recons.keys():
            save_array(recon_group, recon_attr, self.recons[recon_attr])

        return grain_group

    @classmethod
    def from_h5py_group(cls, group, ds, grain):
        """Creates a GrainSinogram object from an h5py group, dataset and grain object"""
        grainsino_obj = GrainSinogram(grain_obj=grain, dataset=ds)

        if "sinograms" in group.keys():
            for sino_attr in group["sinograms"].keys():
                sino_var = group["sinograms"].get(sino_attr)[:]
                setattr(grainsino_obj, sino_attr, sino_var)

        if "recon_parameters" in group.keys():
            for recon_par_attr in ["recon_pad", "recon_shift", "recon_niter", "recon_y0"]:
                if group["recon_parameters"].attrs.get(recon_par_attr) is not None:
                    recon_par_var = group["recon_parameters"].attrs.get(recon_par_attr)[()]
                    setattr(grainsino_obj, recon_par_attr, recon_par_var)

            grainsino_obj.recon_mask = group["recon_parameters"].get("recon_mask")[:]

        if "recons" in group.keys():
            for recon_attr in group["recons"].keys():
                recon_var = group["recons"].get(recon_attr)[:]
                grainsino_obj.recons[recon_attr] = recon_var

        return grainsino_obj


def run_astra(sino, angles, shift=0, pad=0, mask=None, niter=100, astra_method='SIRT_CUDA', workers=None):
    import astra
    angles = np.radians(angles)
    allowed_methods = ['BP', 'SIRT', 'BP_CUDA', 'FBP_CUDA', 'SIRT_CUDA', 'SART_CUDA', 'CGLS_CUDA', 'EM_CUDA']
    if astra_method not in allowed_methods:
        raise ValueError("Unsupported method!")
    manual_mask = None
    if astra_method == 'EM_CUDA' and mask is not None:
        # print("Can't use mask with EM_CUDA method!")
        manual_mask = mask.copy()
        mask = None
    
    vol_geom = astra.create_vol_geom((sino.shape[0]+pad, sino.shape[0]+pad))
    proj_geom = astra.create_proj_geom('parallel', 1.0, sino.shape[0], angles)
    if shift != 0:
        proj_geom = astra.functions.geom_postalignment(proj_geom, shift)
    proj_id = astra.create_projector('linear', proj_geom, vol_geom)
    proj_data_id = astra.data2d.create('-sino', proj_geom, data=sino.T)
    if astra_method == 'EM_CUDA':
        # initialise with ones
        rec_id = astra.data2d.create('-vol', vol_geom, 1.0)
    else:
        rec_id = astra.data2d.create('-vol', vol_geom, 0.0)

    cfg = astra.creators.astra_dict(astra_method)
    
    if 'CUDA' not in astra_method:
        cfg['ProjectorId'] = proj_id
    cfg['ProjectionDataId'] = proj_data_id
    cfg['ReconstructionDataId'] = rec_id
    cfg['option'] = {}
    if astra_method != 'EM_CUDA':
        cfg['option']['MinConstraint'] = 0
        cfg['option']['MaxConstraint'] = 1
    
    if mask is not None:
        mask_id = astra.data2d.create('-vol', vol_geom, mask)
        cfg['option']['ReconstructionMaskId'] = mask_id
    
    alg_id = astra.algorithm.create(cfg)
    astra.algorithm.run(alg_id, iterations=niter)
    recon = astra.data2d.get(rec_id)
     
    # Clean up.
    astra.algorithm.delete(alg_id)
    astra.data2d.delete(rec_id)
    astra.projector.delete(proj_id)
    astra.projector.delete(proj_data_id)
    
    if mask is not None:
        astra.data2d.delete(mask_id)

    if astra_method == 'EM_CUDA' and manual_mask is not None:
        # manually mask
        recon = np.where(manual_mask, recon, 0.0)
    
    return recon
    

def write_h5(filename, list_of_sinos, overwrite_grains=False, group_name='grains'):
    """Write list of GrainSinogram objects to H5Py file
       If overwrite_grains is True, will replace grains in group_name"""
    with h5py.File(filename, "a") as hout:
        grains_group = hout.require_group(group_name)

        for gsinc, gs in enumerate(list_of_sinos):
            group_name = str(gsinc)
            gs.to_h5py_group(parent_group=grains_group, group_name=group_name)
            if overwrite_grains:
                gs.grain.to_h5py_group(parent_group=grains_group, group_name=group_name)


def read_h5(filename, ds, group_name='grains'):
    """Read list of GrainSinogram objects from H5Py file
       Will also create self.grain objects from the H5Py file
       Because GrainSinogram objects can't exist without corresponding grain objects"""
    with h5py.File(filename, "r") as hin:
        grains_group = hin[group_name]
        gs_objects = []

        # take all the keys in the grains group, sort them by integer value, iterate
        for gid_string in sorted(grains_group.keys(), key=lambda x: int(x)):
            grain_group = grains_group[gid_string]

            # import the grain object
            g = ImageD11.grain.grain.from_h5py_group(grain_group)

            # create the GrainSinogram object
            gs = GrainSinogram.from_h5py_group(grain_group, ds, g)

            gs_objects.append(gs)

    return gs_objects


def write_slice_recon(filename, slice_arrays):
    rgb_x_array, rgb_y_array, rgb_z_array, grain_labels_array, raw_intensity_array = slice_arrays
    with h5py.File(filename, "a") as hout:
        slice_group = hout.require_group('slice_recon')
        save_array(slice_group, 'intensity', raw_intensity_array)
        save_array(slice_group, 'labels', grain_labels_array)

        save_array(slice_group, 'ipf_x_col_map', rgb_x_array).attrs['CLASS'] = 'IMAGE'
        save_array(slice_group, 'ipf_y_col_map', rgb_y_array).attrs['CLASS'] = 'IMAGE'
        save_array(slice_group, 'ipf_z_col_map', rgb_z_array).attrs['CLASS'] = 'IMAGE'


def write_pbp_strain(filename, pbp_strain_slice_arrays):
    ubi_map, eps_map, misorientation_map, ipf_x_col_map, ipf_y_col_map, ipf_z_col_map = pbp_strain_slice_arrays
    with h5py.File(filename, "a") as hout:
        slice_group = hout.require_group('pbp_strain_recon')
        save_array(slice_group, 'ubi', ubi_map)
        save_array(slice_group, 'eps', eps_map)
        save_array(slice_group, 'misorientation', misorientation_map)

        save_array(slice_group, 'ipf_x_col_map', ipf_x_col_map).attrs['CLASS'] = 'IMAGE'
        save_array(slice_group, 'ipf_y_col_map', ipf_y_col_map).attrs['CLASS'] = 'IMAGE'
        save_array(slice_group, 'ipf_z_col_map', ipf_z_col_map).attrs['CLASS'] = 'IMAGE'


def read_slice_recon(filename):
    with h5py.File(filename, "r") as hin:
        slice_group = hin["slice_recon"]
        intensity = slice_group["intensity"][:]
        labels = slice_group["labels"][:]
        ipf_x = slice_group["ipf_x_col_map"][:]
        ipf_y = slice_group["ipf_y_col_map"][:]
        ipf_z = slice_group["ipf_z_col_map"][:]

    return ipf_x, ipf_y, ipf_z, labels, intensity


def build_slice_arrays(grainsinos, cutoff_level=0.0, method="iradon", grain_labels=None):
    """Build grain maps from individual grain reonstructions
       Optionally provide a different list of grain labels to label the grains"""

    if grain_labels is not None:
        assert len(grainsinos) == len(grain_labels)

    first_recon = grainsinos[0].recons[method]
    grain_labels_array = np.zeros_like(first_recon) - 1

    redx = np.zeros_like(first_recon)
    grnx = np.zeros_like(first_recon)
    blux = np.zeros_like(first_recon)

    redy = np.zeros_like(first_recon)
    grny = np.zeros_like(first_recon)
    bluy = np.zeros_like(first_recon)

    redz = np.zeros_like(first_recon)
    grnz = np.zeros_like(first_recon)
    bluz = np.zeros_like(first_recon)

    raw_intensity_array = np.zeros_like(first_recon)

    raw_intensity_array.fill(cutoff_level)

    def norm(r):
        m = r > r.max() * 0.2
        return (r / r[m].mean()).clip(0, 1)

    for i, gs in enumerate(grainsinos):

        if grain_labels is not None:
            label = grain_labels[i]
        else:
            label = i

        g_raw_intensity = norm(gs.recons[method])

        g_raw_intensity_mask = g_raw_intensity > raw_intensity_array

        g_raw_intensity_map = g_raw_intensity[g_raw_intensity_mask]

        raw_intensity_array[g_raw_intensity_mask] = g_raw_intensity_map

        redx[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_x[0]
        grnx[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_x[1]
        blux[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_x[2]

        redy[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_y[0]
        grny[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_y[1]
        bluy[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_y[2]

        redz[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_z[0]
        grnz[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_z[1]
        bluz[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_z[2]

        grain_labels_array[g_raw_intensity_mask] = label

    raw_intensity_array[raw_intensity_array == cutoff_level] = 0

    # redx has shape (i, j)
    # (redx, grnx, blux) has shape (3, i, j)
    # transpose changes this to (i, j, 3) needed for mpl imshow
    # crucially: (i, j) unaffected

    rgb_x_array = np.transpose((redx, grnx, blux), axes=(1, 2, 0))
    rgb_y_array = np.transpose((redy, grny, bluy), axes=(1, 2, 0))
    rgb_z_array = np.transpose((redz, grnz, bluz), axes=(1, 2, 0))

    return rgb_x_array, rgb_y_array, rgb_z_array, grain_labels_array, raw_intensity_array



from skimage.filters import threshold_otsu
from skimage.filters import threshold_otsu, threshold_local, threshold_yen, threshold_li
from skimage.morphology import disk, white_tophat, remove_small_holes
from skimage.segmentation import watershed
from skimage.measure import label

def norm(r):
    m = r > r.max() * 0.2
    return (r / r[m].mean()).clip(0, 1)

def filter_normalise_recon(recon):
    """Function to filter and normalise reconstructions"""
    # determine threshold of grain vs background
    thresh = threshold_yen(recon)
    
    # determine binary mask of grain
    binary = recon > thresh

    # this still leaves some hot pixels behind
    footprint = disk(8)
    res = white_tophat(binary, footprint)
    binary_clean = binary ^ res
    
    # binary_clean might contain holes
    binary_clean_noholes = remove_small_holes(binary_clean)
    
    # fig, ax = plt.subplots()
    # ax.imshow(binary_clean_noholes, origin="lower")
    # plt.show()
    
    # determine masked recon image
    recon_masked = np.where(binary_clean_noholes, recon, 0)
    
    # to normalise, we first need to segment into contiguous blobs
    # normalise each of them separetely
    # label image regions
    label_image = label(binary_clean_noholes)
    
    recon_masked_normed = np.zeros_like(recon)
    
    
    for region_id in np.unique(label_image.ravel())[1:]:
        recon_this_region_only = np.where(label_image == region_id, recon_masked, 0)
        recon_this_region_only_normed = norm(recon_this_region_only)
        recon_masked_normed = np.where(label_image == region_id, recon_this_region_only_normed, recon_masked_normed)
    
    return recon_masked_normed



def build_slice_arrays_clean(grainsinos, method="iradon", grain_labels=None):
    """Build grain maps from individual grain reonstructions
       Cleans and normalises reconstructions before building slice arrays
       Saves cleaned and normalised reconstruction to gs.recons[method_norm]
       Optionally provide a different list of grain labels to label the grains"""

    if grain_labels is not None:
        assert len(grainsinos) == len(grain_labels)

    first_recon = grainsinos[0].recons[method]
    grain_labels_array = np.zeros_like(first_recon) - 1

    redx = np.zeros_like(first_recon)
    grnx = np.zeros_like(first_recon)
    blux = np.zeros_like(first_recon)

    redy = np.zeros_like(first_recon)
    grny = np.zeros_like(first_recon)
    bluy = np.zeros_like(first_recon)

    redz = np.zeros_like(first_recon)
    grnz = np.zeros_like(first_recon)
    bluz = np.zeros_like(first_recon)

    raw_intensity_array = np.zeros_like(first_recon)

    for i, gs in enumerate(grainsinos):

        if grain_labels is not None:
            label = grain_labels[i]
        else:
            label = i
        
        # the grain should win where its intensity is greater than what's there before
        g_raw_intensity_map = filter_normalise_recon(gs.recons[method])
        gs.recons[method + "_norm"] = g_raw_intensity_map
        
        g_raw_intensity_mask = g_raw_intensity_map > raw_intensity_array
        
        raw_intensity_array = np.where(g_raw_intensity_mask, g_raw_intensity_map, raw_intensity_array)
        
        redx = np.where(g_raw_intensity_mask, g_raw_intensity_map * gs.grain.rgb_x[0], redx)
        redy = np.where(g_raw_intensity_mask, g_raw_intensity_map * gs.grain.rgb_y[0], redy)
        redz = np.where(g_raw_intensity_mask, g_raw_intensity_map * gs.grain.rgb_z[0], redz)
        
        grnx = np.where(g_raw_intensity_mask, g_raw_intensity_map * gs.grain.rgb_x[1], grnx)
        grny = np.where(g_raw_intensity_mask, g_raw_intensity_map * gs.grain.rgb_y[1], grny)
        grnz = np.where(g_raw_intensity_mask, g_raw_intensity_map * gs.grain.rgb_z[1], grnz)
        
        blux = np.where(g_raw_intensity_mask, g_raw_intensity_map * gs.grain.rgb_x[2], blux)
        bluy = np.where(g_raw_intensity_mask, g_raw_intensity_map * gs.grain.rgb_y[2], bluy)
        bluz = np.where(g_raw_intensity_mask, g_raw_intensity_map * gs.grain.rgb_z[2], bluz)
        
        grain_labels_array[g_raw_intensity_mask] = label

    # redx has shape (i, j)
    # (redx, grnx, blux) has shape (3, i, j)
    # transpose changes this to (i, j, 3) needed for mpl imshow
    # crucially: (i, j) unaffected

    rgb_x_array = np.transpose((redx, grnx, blux), axes=(1, 2, 0))
    rgb_y_array = np.transpose((redy, grny, bluy), axes=(1, 2, 0))
    rgb_z_array = np.transpose((redz, grnz, bluz), axes=(1, 2, 0))

    return rgb_x_array, rgb_y_array, rgb_z_array, grain_labels_array, raw_intensity_array


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


# TODO go to some sort of Numba zone or C code?
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
    # print("Using", T, "threads")
    for tid in numba.prange(T):
        pcf = 0  # thread local I hope?
        for i in order[tid::T]:
            grain_label[i] = -1
            pkid = spot2d_label[i]
            while spot3d_id[pcf] < pkid:
                pcf += 1
            if spot3d_id[pcf] == pkid:
                grain_label[i] = grain_id[pcf]


def get_2d_peaks_from_4d_peaks(p2d, cf):
    """
    inds is an array which tells you which 2D spots each grain owns
    the 2D spots are sorted by spot ID
    inds tells you for each grain were you can find its associated 2D spots"""
    # Big scary block

    # Ensure cf is sorted by spot3d_id
    # NOTE: spot3d_id should be spot4d_id, because we have merged into 4D?
    assert (np.argsort(cf.spot3d_id) == np.arange(cf.nrows)).all()

    numba_order, numba_histo = counting_sort(p2d['spot3d_id'])

    grain_2d_id = palloc(p2d['spot3d_id'].shape, np.dtype(int))

    cleanid = cf.grain_id.copy()

    find_grain_id(cf.spot3d_id, cleanid, p2d['spot3d_id'], grain_2d_id, numba_order)

    gord, counts = counting_sort(grain_2d_id)

    inds = np.concatenate(((0,), np.cumsum(counts)))

    return gord, inds
