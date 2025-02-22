# forward_projector.py
# Providing beam, sample and forward projector to compute forward projected peaks, thus creating forward projected diffraction images, therefore cf_2d, cf_3d or cf_4d peaks
# Beam shape and profile, sample structure factor (to be better computed), absorption, Lorentz, Polarization and detector point spread have all been considered
# Suitable for both s3DXRD using pencil beam and box-beam 3DXRD
# Haixing Fang, haixing.fang@esrf.fr
# Jan 2nd, 2025
# updated on January 30, 2025

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from numba import njit
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from sklearn.cluster import DBSCAN
import h5py

import os, sys
import time
import logging
from tqdm import tqdm
from joblib import Parallel, delayed
import inspect
import multiprocessing
import tempfile

from ImageD11.forward_model import forward_model
from ImageD11.forward_model import grainmaps
from ImageD11.forward_model import pars_conversion
from ImageD11.forward_model import io

import ImageD11.parameters
import ImageD11.unitcell
from ImageD11.sinograms import lima_segmenter
from ImageD11.columnfile import colfile_to_hdf
from ImageD11.columnfile import colfile_from_hdf
from ImageD11.nbGui.nb_utils import slurm_submit_and_wait

econst = 12.3984
logging.basicConfig(level=logging.INFO, force=True)
os.environ["JOBLIB_TEMP_FOLDER"] = "/tmp_14_days/slurm/log"

"""
The right-handed laboratory coodinate system is defined as follows:
X - along the X-ray beam
Y - transversely perpendicular to the X-ray beam and rotation axis, pointing outwards the synchrotron
Z - vetical up
omega - rotation about Z axis as positive direction defined as counter clock-wise when looking top-down on the diffractometer
Note that length unit may be used as mm or um (for dty, y0_offset, beam size, voxel size etc.), see the comments for the variables
"""

class beam:
    """
    Beam class to define X-ray beam size, source-to-sample distance, wavelength, bandwidth and profile
    Possible to create parallel box-beam, cone beam or pencil beam with different combinations
    example:
    ss = beam(Lss = 0, FWHM = [0.001, 0.001], energy = 43.56)   # pencil beam 0.001*0.001 mm^2
    ss = beam(Lss = 97e3, FWHM = [0.5, 0.2], energy = 43.56)    # parallel box beam 0.5*0.2 mm^2
    ss = beam(Lss = 10, FWHM = [0.001, 0.001], energy = 43.56)  # cone beam 0.001*0.001 mm^2
    """
    
    def __init__(self, Lss = 0.0, FWHM = [0.001, 0.001], energy = 43.56, propagation_direction = [1, 0, 0], bandwidth = 0.001,
                       polarization = [0, 0.95, 0], flux = 5e14, beam_profile = "gaussian", cutoff = 0.1):

        self.Lss = Lss       # source-to-rotation axis distance [mm]
        self.FWHM = FWHM     # FWHM dimensions in width x height [mm^2]
        self.energy = energy # X-ray energy [keV]
        self.wavelength = econst / energy # X-ray wavelength [angstrom]
        self.propagation_direction = np.array(propagation_direction) # X-ray direction, along X
        self.bandwidth = bandwidth
        self.polarization = polarization      # along Y axis by default
        self.flux = flux / (FWHM[0]*FWHM[1])  # flux at sample position [photons/s/mm^2]
        self.beam_profile = beam_profile      # name of beam profile function, "gaussian", "lorentzian", "pseudo voigt"
        self.cutoff = cutoff                  # cutoff value for considering intensity/weight
        
        self.wave_vector = ((2 * np.pi / self.wavelength) * self.propagation_direction / np.linalg.norm(self.propagation_direction))
        
        if self.Lss > 5e3:
            self.name = 'parallel_beam'     # including box and line-beam
        elif self.Lss >= 0.5 and self.Lss < 100 and self.FWHM[0] < 0.01 and self.FWHM[1] < 0.01:
            self.name = 'cone_beam'
        elif self.Lss < 0.02 and self.FWHM[0] < 0.01 and self.FWHM[1] < 0.01:
            self.name = 'pencil_beam'
            if self.FWHM[0] >= 0.5 or self.FWHM[1] >= 0.5:
                self.FWHM[0] /= 1000.0   # guess the input FWHM was wrongly input as in um, now it is corrected to be in [mm]
                self.FWHM[1] /= 1000.0
            self.set_beam_shape()
        else:
            self.name = None
            print('Beam character cannot be identified; Only supporting parallel, cone and pencil beam for now !')
    
    def set_beam_spectrum(self, plot_flag = False):
        """
        set beam spectrum by considering the bandwidth
        NOTE: this has not been considered in the forward computation yet
        """
        Energy_spectra = np.linspace(self.energy * (1 - self.bandwidth), self.energy * (1 + self.bandwidth), 1000)  # Photon energy [keV]
        Int_spectra = self.flux * 1 / (2 * np.pi * self.bandwidth) * np.exp(-(Energy_spectra - self.energy) ** 2 / (2 * self.bandwidth ** 2))
        print(Energy_spectra.shape)
        print(type(Int_spectra))
        self.spectrum = np.vstack([Energy_spectra, Int_spectra]).T
        if plot_flag:
            plt.figure()
            plt.plot(self.spectrum[:,0], self.spectrum[:,1],'r.', linewidth=1, markersize=6)
            plt.xlabel('X-ray energy (keV)', fontsize = 16)
            plt.ylabel('Intensity (a.u.)', fontsize = 16)
            plt.tick_params(axis='both', which='major', labelsize=14, width=1.5, length=4)
            plt.show()

    def set_beam_shape(self, plot_flag = False):
        """
        beam shape is used as weights for scaling the voxel contribution as a function of its distance to the ray center and source intensity scaling factor, i.e. two effects
        Gaussian, Lorentzian and pseudo_voigt profiles are available
        """
        x = np.arange(0, 5*np.mean(self.FWHM), np.mean(self.FWHM)/3.0) # 5 times of the beam FWHM is considered as potential effective width for interacting with the sample
        if self.beam_profile.lower() in "gaussian":
            intensities = self.gaussian(x)
        elif self.beam_profile.lower() in "lorentzian":
            intensities = self.lorentzian(x)
        elif self.beam_profile.lower() in ["pseudo", "pseudo_voigt", "pseudo voigt"]:
            intensities = self.pseudo_voigt(x, frac = 0.5)
        else:
            logging.info('{} not identified; I will use gaussian beam profile'.format(self.beam_profile))
            self.beam_profile = "gaussian"
            intensities = self.gaussian(x)
        
        indices = np.where(intensities > self.cutoff)[0]
        self.weight_pos = x[indices] / np.mean(self.FWHM) + 0.5 # consider extra half pixel for the fact the pixel-center-to-ray-center distance can't tell smaller half distance
        self.weight = intensities[indices]
        if plot_flag:
            plt.figure()
            plt.plot(self.weight_pos, self.weight,'r.', linewidth=1, markersize=6)
            plt.xlabel('Distance to the ray center / beam FWHM (-)', fontsize = 16)
            plt.ylabel('Weight (a.u.)', fontsize = 16)
            plt.tick_params(axis='both', which='major', labelsize=14, width=1.5, length=4)
            plt.show()
        
    def gaussian(self, x):
        """
        Gaussian function for describing the beam profile, transformed to weights as a function of distance to the ray center.
        Args:
            x (array): distances, e.g. x = np.arange(0, 5, 20)
        Returns:
            intensities (array)
        """
        # back calculate the sigma value when we know the beam FWHM
        sigma_squared = - np.mean(self.FWHM)**2 / (2*np.log(0.5))   # sigma^2
        return np.exp(-((x - 0.0) ** 2) / (2 * sigma_squared))

    def lorentzian(self, x):
        """
        Lorentzian function for describing the beam profile, translated to weights as a function of distance.
        Args:
            x (array): distances, e.g. x = np.arange(0, 5, 50)
        Returns:
            intensities (array)
        """
        gamma = np.mean(self.FWHM) / 2.0
        return (gamma ** 2) / ((x - 0.0) ** 2 + gamma ** 2)

    def pseudo_voigt(self, x, frac = 0.5):
        """
        Pseudo-Voigt function (combination of Gaussian and Lorentzian).
        Args:
            x (array): distances, e.g. x = np.arange(0, 5, 50)
            frac (float): fraction of lorentzian part
        Returns:
            intensities (array)
        """
        assert frac >= 0 and frac <= 1, "frac must be in [0, 1]"
        intensities_gaussian = self.gaussian(x)
        intensities_lorentzian = self.lorentzian(x)
        return frac * intensities_lorentzian + (1 - frac) * intensities_gaussian   


class sample:
    """
    Sample class to define the sample illuminated by the X-ray beam
    It can be obtained from experimentally reconstructed sample (available for now) or virtually created sample (TODO)
    """
    def __init__(self, filename = None, outname = None, min_misori = 3, crystal_system = 'cubic', dis_tol = np.sqrt(3), remove_small_grains = True, min_vol = 3):
        self.filename = filename
        self.outname = outname
        if self.outname is None and self.filename is not None:
            self.outname = os.path.join(os.path.split(self.filename)[0], 'DS.h5')
        self.min_misori = min_misori
        self.crystal_system = crystal_system
        self.dis_tol = dis_tol
        self.remove_small_grains = remove_small_grains
        self.min_vol = min_vol
        print('****************************************** Parameters for operating the grain map: ')
        print('Output file name: {}'.format(self.outname))
        print('min_misori = {}'.format(self.min_misori))
        print('crystal_system: {}'.format(self.crystal_system))
        print('dis_tol: {}'.format(self.dis_tol))
        print('remove_small_grains = {}'.format(self.remove_small_grains))
        print('min_vol = {} voxels'.format(self.min_vol))
        if self.filename is not None:
            self.read_grainmap()
            self.get_Rsample()
        else:
            self.DS = None
           
    def set_rou(self, rou = 2.7):
        'set sample density [g/cm^3]'
        self.rou = rou # default for Al
    
    def set_mass_abs(self, mass_abs = 0.5685):
        'set mass-energy absorption coefficient [cm^2/g]'
        self.mass_abs = mass_abs # default for Al at 40 keV
    
    def get_Rsample(self):
        'get sample radius by assuming a cylinder-like sample shape from the grain map dimensions [mm]'
        try:
            ind_mask = np.argwhere(self.DS['labels'] > -1)
            Rsample = np.linalg.norm( ind_mask[-1,:] - ind_mask[0,:] ) * np.mean(self.DS['voxel_size'][1:3]) / 2000.0  # [mm]           
            self.Rsample = Rsample
        except KeyError as e:
            print("Key error: {}".format(e))
        except Exception as e:
            print("An unexpected error occurred: {}".format(e))
    
    def read_grainmap(self):
        # read grain map from different types of files: pbp_tensormap_refined.h5, xxx_grains.h5 or DS.h5
        if "DS" in self.filename:
            # first try with DS
            try:
                self.read_DS()
                logging.info("DS file '{}' loaded successfully".format(self.filename))
            except Exception as e1:
                logging.info("An unexpected error occurred: {}".format(e1))
                try:
                    # attempt to read the input of a tensor map from pbp_tensormap_refined.h5 or xxx_grains.h5
                    self.read_tensormap()
                    logging.info("Tensor map '{}' loaded successfully".format(self.filename))
                except Exception as e2:
                    logging.info("Failed to load tensor map: {}".format(e2))
        else:
            try:
                # attempt to read the input of a tensor map from pbp_tensormap_refined.h5 or xxx_grains.h5
                self.read_tensormap()
                logging.info("Tensor map '{}' loaded successfully".format(self.filename))
            except Exception as e1:
                logging.info("Failed to load tensor map: {}".format(e1))
                try:
                    # attempt to read the input of a DS like file, e.g. DS.h5
                    self.read_DS()
                    logging.info("DS file '{}' loaded successfully".format(self.filename))
                except Exception as e2:
                    logging.info("Failed to load DS file: {}".format(e2))
                    raise RuntimeError("Unable to load input as either a tensor map or a DS file.")   
    
    def read_DS(self):
        self.DS = io.read_h5_file(self.filename)
        
    def read_tensormap(self):
        # support s3DXRD outputs such as pbp_tensormap_refined.h5, xxx_grains.h5
        gm = grainmaps.grainmap(filename = self.filename, min_misori = self.min_misori,
                                crystal_system = self.crystal_system, remove_small_grains = self.remove_small_grains, min_vol = self.min_vol)
        gm.merge_and_identify_grains(dis_tol = self.dis_tol)
        self.DS = gm.DS
        
    def write_DS(self):
        """
        write DS to an h5 file with filename defined as self.outname'
        """
        if not self.outname.endswith(('.h5', '.hdf5')):
            self.outname = self.outname + '.h5'
        
        with h5py.File(self.outname, 'w') as hout:
            for key, value in self.DS.items():
                hout.create_dataset(key, data = value)
        print('Done with saving DS to {}'.format(self.outname))
        
    """
    TODO
    Virtually creating microstructure as DS-like data
    """


class forward_projector:
    """
    forward_projector class to compute forward diffraction patterns for a given X-ray beam illuminating a sample
    """
    def __init__(self, sample_filename, pars_filename, phase_name, output_folder = None, gids_mask = None,
                 verbose = 1, auto_set = True, detector_mask = None, to_sparse = False, read_fwd_peaks_from_file = True, **kwargs):
        """
        Providing sample_filename, pars_filename, phase_name etc to set the forward_projector.
        if auto_set is True, it will automatically set the inputs of the beam and sample, as well as acquisition parameters
        
        For computing fwd_peaks for selected grains, use gids_mask for the grains of interest;
        For accounting the detector mask, input detector_mask (2D array), where 1 for active and 0 for inactive
        """
        # get essential input arguments
        self.sample_filename = sample_filename       # sample file name, e.g. pbp_tensormap_refined.h5, xxx_grains.h5 or DS.h5
        self.pars_filename = pars_filename           # ImageD11 parameter file, normally pars.json
        self.phase_name = phase_name                 # phase name specified in the pars.json, e.g. "Al"
        if output_folder is None:
            output_folder = os.path.join(os.getcwd(), "output_fwd")
        self.output_folder = output_folder
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        self.gids_mask = gids_mask                   # grain IDs for masking the sample
        self.verbose = verbose                       # logging level, 0, 1 or 2
        self.auto_set = auto_set
        self.to_sparse = to_sparse
        self.detector_mask = detector_mask
        self.read_fwd_peaks_from_file = read_fwd_peaks_from_file  # if the pwd_peaks file already exists, it will read directly from the file without the need to compute again
        self.cf_2d_file = os.path.join(self.output_folder, 'fwd_cf_2d.h5')
        self.cf_3d_file = os.path.join(self.output_folder, 'fwd_cf_3d.h5')
        self.cf_4d_file = os.path.join(self.output_folder, 'fwd_cf_4d.h5')
        
        # get other optional input arguments
        args = {
        "energy": kwargs.get("energy", 43.56),                 # [keV]
        "beam_size": kwargs.get("beam_size", [1e-3, 1e-3]),    # [mm]
        "beam_profile": kwargs.get("beam_profile", "gaussian"),# [-]
        "flux": kwargs.get("flux", 5e14),                      # [photons/s]
        "Lss": kwargs.get("Lss", 0.0),                         # [mm]
        "min_misori": kwargs.get("min_misori", 3.0),           # [deg]
        "crystal_system": kwargs.get("crystal_system", 'cubic'),
        "remove_small_grains": kwargs.get("remove_small_grains", True),
        "min_vol": kwargs.get("min_vol", 3),                   # [voxel]
        "rou": kwargs.get("rou", 2.7),                         # [g/cm^3]
        "mass_abs": kwargs.get("mass_abs", 0.56685),           # [cm^2/g]
        "y0_offset": kwargs.get("y0_offset", 0.0),             # [um]
        "exp_time": kwargs.get("exp_time", 0.002),             # [s]
        "rot_start": kwargs.get("rot_start", -89.975),         # [deg]
        "rot_end": kwargs.get("rot_end", 90.9668),             # [deg]
        "rot_step": kwargs.get("rot_step", 0.05),              # [deg]
        "sparse_omega": kwargs.get("sparse_omega", True),
        "halfy": kwargs.get("halfy", 182.0),                   # [um]
        "dty_step": kwargs.get("dty_step", 1.0),               # [um]
        "ds_max": kwargs.get("ds_max", 1.2),                   # [1/angstrom]
        "plot_peaks": kwargs.get("plot_peaks", False),
        "plot_flag": kwargs.get("plot_flag", False),
        "detector": kwargs.get("detector", "eiger"),
        "int_factors": kwargs.get("int_factors", (0.1065, 0.7807, 0.1065)),
        "dety_merge_range": kwargs.get("dety_merge_range", 10), # within this range 2D peaks will be merged if the signal comes from the same reflection [pixel]
        "detz_merge_range": kwargs.get("detz_merge_range", 10), # within this range 2D peaks will be merged if the signal comes from the same reflection [pixel]
        "use_cluster": kwargs.get("use_cluster", True),
        "slurm_folder": kwargs.get("slurm_folder", "slurm_fwd_proj"),
        "num_cores": kwargs.get("num_cores", 28)
        }
        
        self.args = args
        if self.verbose >= 1:
            logging.info("Got the following optional arguments: {}".format(self.args))
            logging.info("Output folder: {}".format(self.output_folder))
        if self.auto_set:
            self.set_all()
            # update energy and the beam source
            self.args["energy"] = econst / self.pars.parameters['wavelength']
            self.set_beam()
            if self.verbose >= 1:
                logging.info("Updated X-ray energy = {} keV".format(self.args["energy"]))
    
    def set_all(self):
        """
        set all the inputs for the projector: beam, sample, pars, mask and acquisition.
        """
        self.set_beam()
        self.set_sample()
        self.set_pars()
        self.set_sample_mask()
        self.set_acquisition()
        self.set_opts_seg()
        self.set_image_size()
    
    def set_beam(self):
        """
        set the X-ray beam source, including beam energy, FWHM, flux, shape, source-to-sample distance, name etc.
        """
        self.beam_input = beam(energy = self.args['energy'], FWHM = self.args['beam_size'], flux = self.args['flux'], Lss = self.args['Lss'], beam_profile = self.args['beam_profile'])
        
    def set_sample(self):
        """
        set the sample properties: including the grain map, voxel size, .
        """
        self.sample_input = sample(self.sample_filename, min_misori = self.args['min_misori'], crystal_system = self.args['crystal_system'],
                                      remove_small_grains = self.args['remove_small_grains'], min_vol = self.args['min_vol'])    
        self.sample_input.set_rou(rou = self.args['rou'])
        self.sample_input.set_mass_abs(mass_abs = self.args['mass_abs'])
    
    def set_pars(self):
        """
        set the parameter file.
        """
        self.pars = ImageD11.parameters.read_par_file(self.pars_filename)
        phases = ImageD11.unitcell.Phases(self.pars_filename)
        if self.verbose >= 1:
            logging.info("Got phases: {} and phase name: {}".format(phases, self.phase_name))
        self.ucell = ImageD11.unitcell.unitcell(phases.unitcells[self.phase_name].lattice_parameters, phases.unitcells[self.phase_name].spacegroup)
        
    def set_sample_mask(self):
        """
        set sample mask for computation by grain IDs (gids_mas)
        Args:
            gids (list): a list of grain IDs for making the mask, i.e. only these grains are considered for computation
        """
        if self.gids_mask is not None:
            for i, gid_mask in enumerate(self.gids_mask):
                if i == 0:
                    mask = (self.sample_input.DS['labels'] == gid_mask) & (~np.isnan(self.sample_input.DS['U'][:, :, :, 0, 0]))
                else:
                    mask += (self.sample_input.DS['labels'] == gid_mask) & (~np.isnan(self.sample_input.DS['U'][:, :, :, 0, 0]))
            mask = np.moveaxis(mask, 0, 2)
            self.sample_mask = mask
        else:
            self.sample_mask = None
    
    def set_acquisition(self):
        """
        set acquisition parameters: omega angles and dty range for scanning 3DXRD or box-beam 3DXRD
        """
        # set omega angles, for computing the peaks omega step can be larger, e.g. 0.5 deg for accelerating the computation
        self.omega_angles = np.arange(self.args['rot_start'], self.args['rot_end']+self.args['rot_step']/2.0, self.args['rot_step'])
        self.omega_angles_sparse = np.arange(self.args['rot_start'], self.args['rot_end']+self.args['rot_step']/2.0, np.max([self.args['rot_step'], 0.5]))
        
        if not hasattr(self, 'beam_input'):
            self.set_beam()            
        if self.beam_input.name == 'pencil_beam':
            self.set_dty(scanning_mode = True)   # a list of dty position for scanning
        else:
            self.set_dty(scanning_mode = False)  # only 1 dty position if it is not in scanning mode
            
    def set_dty(self, scanning_mode = True):
        """
        set dty range for data collection
        """
        if scanning_mode:
            halfy = np.abs(self.args['halfy'])
            dty_step = np.abs(self.args['dty_step'])
            self.dtys = np.arange(-halfy, halfy+dty_step/2.0, dty_step)
            if self.verbose >= 1:
                logging.info("dty range: [{} {}] with {} step size (um)".format(self.dtys[0], self.dtys[-1], dty_step))
        else:
            self.dtys = [0.0,]
            if self.verbose >= 1:
                logging.info("Not scanning mode !!! dtys = {} (um)".format(self.dtys))
                
    def set_opts_seg(self, bg = None, cut = 1, pixels_in_spot = 1):
        """
        set segmentation parameters
        """
        self.opts_seg = get_opts_seg(mask = self.detector_mask, bg = bg, cut = cut, pixels_in_spot = pixels_in_spot)
        
    def set_image_size(self):
        """
        set image size for the projection images
        """
        p = pars_conversion.convert_ImageD11pars2p(self.pars)
        self.image_size = (p['detysize'], p['detzsize'])
    
    def run_all(self):
        """
        run all dty positions with forward projection either on cluster or on local machine
        """
        if self.args['use_cluster']:
            bash_script_path = self.generate_slurm_script()
            print("Submitted {} SLURM jobs for dty computation.".format(bash_script_path))
            slurm_submit_and_wait(bash_script_path, wait_time_sec=60)
        else:
            self.run_series_dty()
        
    def run_series_dty(self, dtys = None):
        """
        run multiple dty positions with forward projection
        Args:
            dty (None or list): None means run all the dty positions in the class
        """
        if dtys is None:
            dtys = self.dtys
        for i, dty in enumerate(dtys):
            self.run_single_dty(dty)
            if self.verbose >=1 and (i > 1 and i % 10 == 0):
                logging.info("Done for {} / {} dty steps".format(i, len(dtys)))
                
    def run_single_dty(self, dty = 0.0):
        """
        run single dty forward projection
        computing time: ~12 s for computing fwd_peaks and ~25 s for segmenting and generating peak sparse file for 3620 omega angles for 361*361 grain map
        """
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        if self.args['sparse_omega']:
            omega_angles_to_calc = self.omega_angles_sparse
        else:
            omega_angles_to_calc = self.omega_angles
        # check existence of fwd_peaks_file and reads it from file if exists
        fwd_peaks_exists = False
        fwd_peaks_3d_exists = False
        if self.read_fwd_peaks_from_file:
            fwd_peaks_file = os.path.join(self.output_folder, 'fpks_dty_' + str(round(dty, 2)).replace('.', 'p') + '.h5')
            if os.path.exists(fwd_peaks_file):
                fwd_peaks = io.read_fwd_peaks(fwd_peaks_file, verbose=self.verbose)
                fwd_peaks_exists = True
                
            fwd_peaks_3d_file = os.path.join(self.output_folder, 'fpks_3d_dty_' + str(round(dty, 2)).replace('.', 'p') + '.h5')
            if os.path.exists(fwd_peaks_3d_file):
                fwd_peaks_3d = io.read_fwd_peaks(fwd_peaks_3d_file, verbose=self.verbose)
                fwd_peaks_3d_exists = True
        
        # Compute fwd_peaks
        if not fwd_peaks_exists:
            fwd_peaks = forward_peaks_voxels(self.beam_input, self.sample_input, omega_angles_to_calc, self.ucell, self.pars, dty = dty, mask = self.sample_mask,
                         ds_max = self.args['ds_max'], exp_time = self.args['exp_time'], plot_peaks = self.args['plot_peaks'],
                         verbose = self.verbose, use_cluster = self.args['use_cluster'], y0_offset = self.args['y0_offset'], plot_flag = self.args['plot_flag'], detector = self.args['detector'])
        # write fwd_peaks
        if fwd_peaks.size > 0:
            # merge 2D peaks to 3D peaks
            if not fwd_peaks_3d_exists:
                fwd_peaks_3d, clusters = merge_rows(fwd_peaks, columns_to_cluster=[0, 9, 12, 13, 14, 18, 19],
                                                    epsilons=[0.1, 2, 0.1, 0.1, 0.1, self.args['dety_merge_range'], self.args['detz_merge_range']],
                                                    special_col_index=23, min_samples = 1, weighted_avg=True)
                io.write_fwd_peaks(fwd_peaks_3d, output_folder=self.output_folder, fname_prefix='fpks_3d', verbose=self.verbose)
            if not fwd_peaks_exists:
                io.write_fwd_peaks(fwd_peaks, output_folder=self.output_folder, fname_prefix=None, verbose=self.verbose)
                
            if self.to_sparse:
                # processing frame-by-frmae by making the projection and immediately segmenting sparse peaks
                # ~32 s for 3620 omega angles, 3.6 times faster
                destname = os.path.join(self.output_folder, 'fsparse_dty_' + str(round(dty, 2)).replace('.', 'p') + '.h5')
                make_projs_and_sparse_file(
                    fwd_peaks,
                    destname,
                    self.omega_angles,
                    self.opts_seg,
                    detector_mask = self.detector_mask,
                    detector = self.args['detector'],
                    image_size=self.image_size,
                    int_factors=self.args['int_factors']
                )

                # # an older version:  generate frames over rotation angles and then segmenting, memory consuming and too slower
                # # ~80 s (making projs) + 34 s (segmenting peaks) for 3620 projs
                # frms, _ = make_projections_with_psf(
                #     fwd_peaks,
                #     self.omega_angles,
                #     image_size=self.image_size,
                #     detector=self.args['detector'],
                #     int_factors=self.args['int_factors']
                # )
                # segment_frms(frms, destname, detector=self.args['detector'], opts_seg = self.opts_seg)
        else:
            print("No peaks expected for dty = {}".format(dty))
    

    def generate_slurm_script(self, script_name="fwd_projector_slurm.sh"):
        """
        Generate a SLURM submission script to run `run_series_dty` with grouped dty values.

        Args:
            script_name (str): Name of the SLURM script file.

        Returns:
            str: Path to the generated SLURM script.
        """
        slurm_folder = self.args['slurm_folder']
        if not os.path.exists(slurm_folder):
            os.makedirs(slurm_folder)
        id11_code_path = get_ImageD11_code_path()
        # Split dty values into `num_groups` parts
        num_groups = np.min([int(len(self.dtys) / 12), 40])
        dty_groups = np.array_split(self.dtys, num_groups)
        if self.verbose >= 1:
            logging.info('Got {} dty values'.format(len(self.dtys)))
            logging.info("All dty positions will be divided into {} jobs and each job computes about {} dty values".format(num_groups, dty_groups[0].shape[0]))

        bash_script_path = os.path.join(slurm_folder, script_name)
        outfile_path = os.path.join(slurm_folder, "fp_slurm_%A_%a.out")
        errfile_path = os.path.join(slurm_folder, "fp_slurm_%A_%a.err")
        log_path = os.path.join(slurm_folder, "fp_slurm")

        # Create a temporary file to store dty groups for safer parsing
        dty_file_list = os.path.join(slurm_folder, "dty_groups.txt")
        with open(dty_file_list, "w") as f:
            for group in dty_groups:
                f.write(" ".join(map(str, group)) + "\n")
        
        bash_script_string = """#!/bin/bash
#SBATCH --job-name=fwd_projector
#SBATCH --output={outfile_path}
#SBATCH --error={errfile_path}
#SBATCH --partition=nice
#SBATCH --array=1-{num_groups}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={num_cores}
#SBATCH --mem-per-cpu=8000
#SBATCH --time=01:00:00

date
source /cvmfs/hpc.esrf.fr/software/packages/linux/x86_64/jupyter-slurm/latest/envs/jupyter-slurm/bin/activate
log_path="{log_path}_${{SLURM_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}.log"
dty_file="{dty_file_list}"
dty_values=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" $dty_file)
echo "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH={id11_code_path}" >> $log_path 2>&1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH={id11_code_path} \
python3 -c "from ImageD11.forward_model import forward_projector; \
fp = forward_projector.forward_projector(sample_filename='{sample_filename}', pars_filename='{pars_filename}', \
phase_name='{phase_name}', output_folder='{output_folder}', gids_mask={gids_mask}, \
verbose={verbose}, auto_set={auto_set}, detector_mask={detector_mask}, to_sparse={to_sparse}, \
read_fwd_peaks_from_file={read_fwd_peaks_from_file}, **{args}); \
fp.run_series_dty([float(v) for v in '$dty_values'.split()])" >> $log_path 2>&1
date""".format(
            outfile_path=outfile_path,
            errfile_path=errfile_path,
            num_groups=num_groups,
            num_cores=self.args["num_cores"],
            dty_file_list=dty_file_list,
            log_path=log_path,
            id11_code_path=id11_code_path,
            sample_filename=self.sample_filename,
            pars_filename=self.pars_filename,
            phase_name=self.phase_name,
            output_folder=self.output_folder,
            gids_mask=self.gids_mask,
            verbose=self.verbose,
            auto_set=self.auto_set,
            detector_mask=self.detector_mask,
            to_sparse=self.to_sparse,
            read_fwd_peaks_from_file=self.read_fwd_peaks_from_file,
            args=self.args
        )
        with open(bash_script_path, "w") as f:
            f.writelines(bash_script_string)

        print("Done with writing bash script in", bash_script_path)
        return bash_script_path

    
    def get_cf(self, dtys = None):
        """
        gather peaks from file for a given list of dty values and convert them to cf object using ImageD11.sinograms.dataset.colfile_from_dict
        by default, gather all peaks from all the files existing in the self.output_folder
        """
        self.peaks_2d = []
        self.peaks_3d = []
        if dtys is None:
            # gather all peaks
            dtys = self.dtys
        for dty in dtys:
            # 2D peaks
            fwd_peaks_file = os.path.join(self.output_folder, 'fpks_dty_' + str(round(dty, 2)).replace('.', 'p') + '.h5')
            if os.path.exists(fwd_peaks_file):
                self.peaks_2d.append( io.read_fwd_peaks(fwd_peaks_file, verbose=0) )
            else:
                if self.verbose >= 1:
                    logging.info("{} not found ...".format(fwd_peaks_file))
            # 3D peaks
            fwd_peaks_3d_file = os.path.join(self.output_folder, 'fpks_3d_dty_' + str(round(dty, 2)).replace('.', 'p') + '.h5')
            if os.path.exists(fwd_peaks_3d_file):
                self.peaks_3d.append( io.read_fwd_peaks(fwd_peaks_3d_file, verbose=0) )
            else:
                if self.verbose >= 1:
                    logging.info("{} not found ...".format(fwd_peaks_3d_file))

        if len(self.peaks_2d) > 0:
            self.peaks_2d = np.vstack(self.peaks_2d)
            self.cf_2d = io.convert_fwd_peaks_to_cf(self.peaks_2d)
            if self.verbose >= 1:
                logging.info("Got cf_2d !")
        else:
            print("No peaks_2d is found")
        if len(self.peaks_3d) > 0:
            self.peaks_3d = np.vstack(self.peaks_3d)
            self.cf_3d = io.convert_fwd_peaks_to_cf(self.peaks_3d)
            if self.verbose >= 1:
                logging.info("Got cf_3d !")
        else:
            print("No peaks_3d is found")
            
    def get_cf_4d(self, source = '3d', factor = 50):
        """
        merge 2D/3D peaks to 4D by merging dty
        source = '2d' or '3d'
        """
        ystep = self.args['dty_step']
        if source == '2d':
            if not (hasattr(self, 'peaks_2d') and self.peaks_2d.size > 0):
                self.read_cf()
            # assemble to 4D peaks by merging dty, omega for 2D peaks
            print('Merging peaks_2d to render peaks_4d ...')
            self.peaks_4d, _ = merge_rows(self.peaks_2d, columns_to_cluster=[0, 2, 3, 8, 9, 10, 18, 19],
                                          epsilons=[0.1, factor*ystep, factor*ystep, factor*ystep, 3, 0.2, 20, 20],
                                          special_col_index=23, min_samples = 1, weighted_avg=True)
        else:
            if not (hasattr(self, 'peaks_3d') and self.peaks_3d.size > 0):
                self.read_cf()
            # assemble to 4D peaks by merging dty, omega for 3D peaks
            print('Merging peaks_3d to render peaks_4d ...')
            self.peaks_4d, _ = merge_rows(self.peaks_3d, columns_to_cluster=[0, 2, 3, 8, 9, 10, 18, 19],
                                      epsilons=[0.1, factor*ystep, factor*ystep, factor*ystep, 3, 0.2, 20, 20],
                                      special_col_index=23, min_samples = 1, weighted_avg=True)
        self.cf_4d = io.convert_fwd_peaks_to_cf(self.peaks_4d)
        if not os.path.exists(self.cf_4d_file):
            colfile_to_hdf(self.cf_4d, self.cf_4d_file)
        else:
            print('{} has already existed; not overwriting it ...'.format(self.cf_4d_file))

    def read_cf(self):
        """
        read peaks from h5 file, using ImageD11.columnfile.colfile_from_hdf
        """
        if os.path.exists(self.cf_2d_file):
            self.cf_2d = colfile_from_hdf(self.cf_2d_file)
            self.peaks_2d = io.convert_cf_to_fwd_peaks(self.cf_2d)
            if self.verbose >= 1:
                logging.info('loaded 2D peaks from {}'.format(self.cf_2d_file))
            if os.path.exists(self.cf_2d_file):
                self.cf_3d = colfile_from_hdf(self.cf_3d_file)
                self.peaks_3d = io.convert_cf_to_fwd_peaks(self.cf_3d)
                if self.verbose >= 1:
                    logging.info('loaded 3D peaks from {}'.format(self.cf_3d_file))
                if os.path.exists(self.cf_4d_file):
                    self.cf_4d = colfile_from_hdf(self.cf_4d_file)
                    self.peaks_4d = io.convert_cf_to_fwd_peaks(self.cf_4d)
                    if self.verbose >= 1:
                        logging.info('loaded 4D peaks from {}'.format(self.cf_4d_file))   
        else:
            self.get_cf()

    def write_cf(self):
        """
        write 2D, 3D and 4D peaks to h5 file, using ImageD11.columnfile.colfile_to_hdf
        """
        if hasattr(self, 'cf_2d') and self.cf_2d.dty.size > 0:
            if not os.path.exists(self.cf_2d_file):
                colfile_to_hdf(self.cf_2d, self.cf_2d_file)
                if self.verbose >= 1:
                    logging.info("Done with writing cf_2d !")
            else:
                print('{} already exists; Not overwriting it.'.format(self.cf_2d_file))
        if hasattr(self, 'cf_3d') and self.cf_3d.dty.size > 0:
            if not os.path.exists(self.cf_3d_file):
                colfile_to_hdf(self.cf_3d, self.cf_3d_file)
                if self.verbose >= 1:
                    logging.info("Done with writing cf_3d !")
            else:
                print('{} already exists; Not overwriting it.'.format(self.cf_3d_file))
        if hasattr(self, 'cf_4d') and self.cf_4d.dty.size > 0:
            if not os.path.exists(self.cf_4d_file):
                colfile_to_hdf(self.cf_4d, self.cf_4d_file)
                if self.verbose >= 1:
                    logging.info("Done with writing cf_4d !")
            else:
                print('{} already exists; Not overwriting it.'.format(self.cf_4d_file))
                
    def plot_cf(self, cf_type = '2d', m=None):
        """
        plot the sinogram of cf peaks, i.e. dty vs omega colored by peak intensity
        """
        plt.figure(figsize=(10,8))
        if cf_type == '2d':
            if m is None:
                m = (self.cf_2d.grainID > -1)
            if self.verbose >= 1:
                logging.info('Got {} / {} peaks to plot'.format(np.where(m==True)[0].shape[0], self.cf_2d.nrows))
            plt.hist2d( self.cf_2d.omega[m], self.cf_2d.dty[m], weights= np.log(self.cf_2d.sum_intensity[m]),
                   bins = (self.omega_angles.shape[0], self.dtys.shape[0]), norm='log')
        if cf_type == '3d':
            if m is None:
                m = (self.cf_3d.grainID > -1)
            if self.verbose >= 1:
                logging.info('Got {} / {} peaks to plot'.format(np.where(m==True)[0].shape[0], self.cf_3d.nrows))
            plt.hist2d( self.cf_3d.omega[m], self.cf_3d.dty[m], weights= np.log(self.cf_3d.sum_intensity[m]),
                   bins = (self.omega_angles.shape[0], self.dtys.shape[0]), norm='log')
        if cf_type == '4d':
            if m is None:
                m = (self.cf_4d.grainID > -1)
            if self.verbose >= 1:
                logging.info('Got {} / {} peaks to plot'.format(np.where(m==True)[0].shape[0], self.cf_4d.nrows))
            plt.hist2d( self.cf_4d.omega[m], self.cf_4d.dty[m], weights= np.log(self.cf_4d.sum_intensity[m]),
                   bins = (self.omega_angles.shape[0], self.dtys.shape[0]), norm='log')
        plt.colorbar(label="Log Intensity")
        plt.xlabel("Omega ($^{o}$)")
        plt.ylabel("dty ($\mu$m)")
        plt.show()


def forward_peaks_voxels(beam, sample, omega_angles, ucell, pars, dty=0.0,
                         mask=None, ds_max=1.2, exp_time = 0.002, plot_peaks = False,
                         verbose=1, use_cluster = False, **kwargs):
    """
    Forward calculating the expected peaks on the peaks for a given beam and sample
    For a given dty and beam size, first computed the list of intersected sample voxels with the ray;
    Then, computed the possible peaks on the detector for each voxel
    Suitable for both s3DXRD and box-beam 3DXRD
    
    Benchmark shows that this parallelized version is 11.4 times faster than non-parallelized version:
    360 projs -> 16 s, 24 s with new beam profile
    720 projs -> 40 s, 60 s with new beam profile
    3600 projs -> 180 s, 312 s with new beam profile
    
    Args:
        beam (class): beam class defining the beam properties
        sample (class): sample class defining the sample properties
        omega_angles (array): rotation angles sorted from small to big [deg]
        ucell (object): ImageD11.unitcell.unitcell object, e.g. ucell = ImageD11.unitcell.unitcell([4.04761405231186, 4.04761405231186, 4.04761405231186, 90.0, 90.0, 90.0], 225)
        pars (objected): ImageD11.parameters object, e.g. pars = ImageD11.parameters.read_par_file('/data/visitor/ma6288/id11/20241119/PROCESSED_DATA/nscope_pars/pars.json')
        dty (float): rotation axis position, dty [um]
        mask (3D array): mask of sample effective voxels
        ds_max (float):  maximum 1/d
        exp_time (float): exposure time [s]
        plot_peaks (logic): flag for plotting peaks
        **kwargs:
        y0_offset (float): offset value of the sample when dty = 0, [um]
        weight (array): weights of considering fractional contribution of the voxel for scaled to the distance-to-ray-center, determined by the beam FWHM and profile
        weight_pos (array): distances corresponding to the weights, scaled by beam FHWM
        plot_flag (logic): plot or not for intersected voxels
        detector (string): detector name, 'eiger' by default
        
    Returns:
        results: forward peaks expected on the detector
    """
    beampol = beam.polarization[1]
    DS = sample.DS
    if mask is None:
        mask = (DS['labels'] > -1) & (~np.isnan(DS['U'][:, :, :, 0, 0]))
        mask = np.moveaxis(mask, 0, 2)
    Lsam2det = pars.get_parameters()['distance'] / 1000.0
    rot_step = omega_angles[1] - omega_angles[0]
    if rot_step < 0:
        rot_step = -rot_step

    args = {
        "dty": dty,
        "y0_offset": kwargs.get("y0_offset", 0.0),
        "voxel_size": [DS['voxel_size'][2], DS['voxel_size'][1], DS['voxel_size'][0]],  # change to follow x, y, z order
        "ray_size": np.mean(beam.FWHM) * 1000, # [um]
        "weight": beam.weight,
        "weight_pos": beam.weight_pos,
        "plot_flag": kwargs.get("plot_flag", False),
        "detector": kwargs.get("detector", "eiger"),
    }
    
    # "loky" uses process-based parallelism, requiring data pickling (which is failing).
    # "threading" uses thread-based parallelism, avoiding memory-mapped files, but very SLOW !!!
    t0 = time.time()
    results = Parallel(n_jobs=-1, backend="loky")(delayed(process_omega)(
        omega, mask, DS, beam, sample, pars, ucell, rot_step, ds_max, Lsam2det, args, verbose
    ) for omega in tqdm(omega_angles, desc="Processing omega angles"))
    # # Regarding multiprocessing.Pool, tried chunk, shared memory etc, it does not help to improve the speed
    # set_tmp_dir()
    # with multiprocessing.Pool(processes=min(multiprocessing.cpu_count(), 24)) as pool:
    #     results = list(tqdm(pool.starmap(process_omega, [
    #         (omega, mask, DS, beam, sample, pars, ucell, rot_step, ds_max, Lsam2det, args, verbose)
    #         for omega in omega_angles
    #     ]), total=len(omega_angles), desc="Processing omega angles"))
    
    # Flatten the nested list of results
    # 0-4 grainID, voxel_indices(ZYX), weight(contributing to fraction of voxel and fraction of source intensity)
    # 5-7 pos
    # 8-14 dty, rot, tth, eta, hkl
    # 15-19 gx, gy, gz, dety_pixel(fc), detz_pixel(sc)
    # 20-22 Lorentz, Polarization, transmission factors
    # 23-24 intensity, ds, append to the last column
    fwd_peaks = [item for sublist in results for item in sublist]
    fwd_peaks = np.array(fwd_peaks)
    t1 = time.time()
    if verbose >= 1:
        logging.info("Computing dty = {} took {} seconds ...".format(dty, t1-t0))
    
    # sum of intensity: scalable with an arbitrary factor
    Ahkls = ucell.gethkls(ds_max)  # a list of structure_factor + hkl
    ds_values = np.array(find_ds_for_multiple_targets(Ahkls, fwd_peaks[:, 12:15]), dtype='float') # 1 / d-spacing
    # TODO: make a structure_factor calculation given atomic sites and sint, see structure_factor.m in DCT code
    structure_factor = 1.0/ds_values # there is no function to calculate structure factor, so take it proportional to 1/d* first for the moment
    
    Vcell = cellvolume(ucell.lattice_parameters)    # unit cell volume [angstrom^3]
    Vvoxel = DS['voxel_size'][0] * DS['voxel_size'][1] * DS['voxel_size'][2] # voxel volume [um^3]
    K1 = K1_calc()                                  # square of Thomson scattering lenth r0^2 [mm^2]
    K2 = (pars.parameters['wavelength']**3) * fwd_peaks[:, 4]*Vvoxel * 1e12 / (Vcell**2)  # K2, dimensionless [-]
    
    # K1 * K2 * I0 * Lorentz * Polarization * trans_factor * structure_factor * exp_time
    intensity = K1 * K2 * beam.flux * fwd_peaks[:, 4] * fwd_peaks[:, 20] * fwd_peaks[:, 21] * fwd_peaks[:, 22] * structure_factor * exp_time # [-]
    
    fwd_peaks = np.hstack([fwd_peaks, intensity.reshape(-1, 1)]) # append the 'intensity' to the last column for fwd_peaks
    fwd_peaks = np.hstack([fwd_peaks, ds_values.reshape(-1, 1)]) # append 'ds_values' to the last column for fwd_peaks

    # shift sampos_y back, remember to account for dty + y0_offset shift for reproducing fwd_peaks calculation
    fwd_peaks[:,6] = fwd_peaks[:,6] + (-args['dty'] - args['y0_offset'])/1000.0 # sampos_y [mm]
    
    if plot_peaks:
        plot_fwd_peaks(fwd_peaks)
        
    return fwd_peaks


def process_omega(omega, mask, DS, beam, sample, pars, ucell, rot_step, ds_max, Lsam2det, args, verbose):
    """
    sub-function for processing one omega angle for forward_peaks_voxels
    """
    
    intersected_sampos, intersected_labpos, intersected_voxels = intersected_sample_pos(mask, omega = omega, **args)
    fwd_peaks = []
    
    if intersected_voxels is not None:
        rot_start = omega - rot_step / 2
        rot_end = omega + rot_step
        for i in range(intersected_sampos.shape[0]):
            pos = intersected_sampos[i, :]  # [mm]
            ind_voxel = [int(intersected_voxels[i, 2]), int(intersected_voxels[i, 1]), int(intersected_voxels[i, 0])]
            if DS['labels'][ind_voxel[0], ind_voxel[1], ind_voxel[2]] > -1:
                U = DS['U'][ind_voxel[0], ind_voxel[1], ind_voxel[2], :, :]
                B = DS['B'][ind_voxel[0], ind_voxel[1], ind_voxel[2], :, :]
                fwd, Nr_simu = forward_model.forward_comp(pos, U, B, ucell, pars, ds_max=ds_max, rot_start=rot_start, rot_end=rot_end, rot_step=rot_step, verbose=0)

                if Nr_simu > 0 and fwd[0][8]:
                    Gt = fwd[0][3]
                    eta_rad = np.arccos(np.dot(np.array([0, Gt[1], Gt[2]]) / np.linalg.norm([0, Gt[1], Gt[2]]), np.array([0, 0, 1])))
                    if Gt[1] > 0:
                        eta_rad = 2 * np.pi - eta_rad
                    eta = np.rad2deg(eta_rad)
                    rho_rad = eta_rad + np.pi / 2

                    tth_rad = np.deg2rad(fwd[0][2])
                    Lorentz = 1 / np.sin(tth_rad)
                    Polarization = (1 + np.cos(tth_rad) ** 2 - beam.polarization[1] * np.cos(2 * rho_rad) * np.sin(tth_rad) ** 2) / 2

                    dety_mm = fwd[0][6]
                    detz_mm = fwd[0][7]
                    trans_factor, L_total = forward_model.beam_attenuation(intersected_labpos[i, :], Lsam2det, dety_mm, detz_mm,
                                                                         Rsample=sample.Rsample, beam_name=beam.name, rou=sample.rou, mass_abs=sample.mass_abs)

                    fwd_peaks.append([
                        DS['labels'][ind_voxel[0], ind_voxel[1], ind_voxel[2]], ind_voxel[0], ind_voxel[1], ind_voxel[2], intersected_voxels[i, 4],  # 0-4 grainID, voxel_indices(ZYX), weight
                        intersected_sampos[i, 0], intersected_sampos[i, 1], intersected_sampos[i, 2],  # 5-7 pos
                        args["dty"], fwd[0][0], fwd[0][2], eta, fwd[0][1][0], fwd[0][1][1], fwd[0][1][2],  # 8-14 dty, rot, tth, eta, hkl
                        fwd[0][3][0], fwd[0][3][1], fwd[0][3][2], fwd[0][4], fwd[0][5],  # 15-19 gx, gy, gz, dety_pixel(fc), detz_pixel(sc)
                        Lorentz, Polarization, trans_factor  # 20-22 Lorentz, Polarization, transmission factors
                    ])

    return fwd_peaks


def intersected_sample_pos(mask, dty = 0.0, y0_offset = 0.0, voxel_size = [1.0, 1.0, 1.0], omega = 0.0, ray_size = 1.0,
                           weight = [1.0, 0.8, 0.5, 0.13], weight_pos = [0.5, 1.0, 1.5, 2.3], plot_flag = False, detector = 'eiger'):
    """
    Given a dty and a sample rotation angle omega, compute a list of voxels intersected with the X-ray beam and convert them to both sample and lab positions
    In lab system, the sample rotates by omega, corresponding to rotating the incoming X-ray source by -omega without rotating the sample.
    This has the benefit to keep the sample mask shape constant throughout the whole calculation
    Note that all the calculations are based on the sample coordinate system and do conversions after
    
    Args:
        mask (3D array): sample (or grain) mask Nx*Ny*Nz, note that sample axis normally follow Z, Y, X, so permute the axis may be required
        dty (float): horizontal distance between rotation axis and X-ray beam [um]
        y0_offset (float): offset of y0, if dty = 0 going through the sample center, then y0_offset = 0 [um]
        voxel size (float): voxel size in x, y, z axis [um]
        omega (float): sample rotation angle [deg]
        ray_size (float): FWHM of the incoming X-ray beam [um]       
        weight (array): weights of considering fractional contribution of the voxel for scaled to the distance-to-ray-center, determined by the beam FWHM and profile
        weight_pos (array): distances corresponding to the weights, scaled by beam FHWM
        plot_flag (logic): plotting flag
        detector (string): name of the detector, eiger or frelon related to the flip story of the detector
        
    Returns:
        intersected_sampos: intersected positions of the sample voxels in sample coordinate system after shifting by dty+y0_offset [x y z] [mm]
        intersected_labpos: intersected positions of the sample voxels in lab coordinate system [x y z] [mm]
        intersected_voxels: intersected voxels follow [x y z]
    """
    
    if not isinstance(voxel_size, np.ndarray):
        voxel_size = np.array(voxel_size, dtype='float')
    Omega = forward_model.get_omega_matrix(np.deg2rad(-omega))
    ray_direction = np.dot(Omega, [1.0, 0.0, 0.0])
    grid_size = mask.shape
    assert len(grid_size) == 3, "The input mask must be in 3D"
    ray_origin = [0.0, -dty - y0_offset, 0.0]

    intersected_voxels = intersected_voxels_3d(grid_size, ray_origin, ray_direction, voxel_size = voxel_size,
                                               ray_size=ray_size, weight = weight, weight_pos = weight_pos, mask = mask, plot_flag = plot_flag)
    
    if intersected_voxels.size > 0:
        intersected_sampos = np.zeros((intersected_voxels.shape[0], 3), dtype = 'float')
        if detector in ['Eiger', 'eiger']:
            for i in range(3):
                intersected_sampos[:,i] = (intersected_voxels[:,i] - grid_size[i]/2 + 0.5) * voxel_size[i]
        elif detector in ['Frelon', 'frelon']:
            for i in range(3):
                intersected_sampos[:,i] = (intersected_voxels[:,i] - grid_size[i]/2 + 0.5) * voxel_size[i] # to be tested

        # convert to lab positions
        intersected_sampos[:,1] = intersected_sampos[:,1] + dty + y0_offset
        intersected_sampos /= 1000.0                                           # [mm]
        Omega = forward_model.get_omega_matrix(np.deg2rad(omega))
        intersected_labpos = np.einsum('ij,nj->ni', Omega, intersected_sampos) # [mm]
    else:
        intersected_sampos, intersected_labpos, intersected_voxels = None, None, None
    
    return intersected_sampos, intersected_labpos, intersected_voxels


def intersected_voxels_3d(grid_size, ray_origin, ray_direction, voxel_size = [1.0, 1.0, 1.0], ray_size = 1.0,
                          weight = [1.0, 0.8, 0.5, 0.13], weight_pos = [0.5, 1.0, 1.5, 2.3], mask = None, plot_flag = False):
    """
    Calculate and visualize the intersected voxels for a given ray in 3D.
    
    Args:
        grid_size (tuple): Dimensions of the grid (Nx, Ny, Nz).
        ray_origin (array): Origin of the ray [x, y, z].
        ray_direction (array): Direction of the ray [dx, dy, dz].
        voxel_size (array): voxel size of the grain map [um/pixel]
        ray_size (float): FWHM of the ray [um]
        weight (array): weights of considering fractional contribution of the voxel for scaled to the distance-to-ray-center, determined by the beam FWHM and profile
        weight_pos (array): distances corresponding to the weights, scaled by beam FHWM
        mask (array): 0 and 1 values, 1 corresponds to activated sample voxels
        plot_flag (logic): to plot or not.
    
    Returns:
        intersected (list): List of intersected voxel coordinates [x, y, z].
    """
    
    if not isinstance(voxel_size, np.ndarray):
        voxel_size = np.array(voxel_size, dtype='float')
    ray_direction = ray_direction / np.linalg.norm(ray_direction)  # Normalize the direction
    grid_extent = np.array(grid_size) * voxel_size
    max_extent = np.linalg.norm(grid_extent)

    t_values = np.linspace(-0.707 * max_extent + ray_origin[1], 0.707 * max_extent + ray_origin[1], int(max(grid_size)*1.414))
    ray_shift = np.array([grid_size[1] / 2 + 0.5, grid_size[1] / 2 + 0.5, 0.0]) * voxel_size
    ray_path = np.array(ray_origin) / voxel_size + ray_shift + np.outer(t_values, ray_direction)

    tol_distances = np.array(weight_pos)*ray_size
    weights       = np.array(weight)

    if mask is None:
        x, y, z = np.meshgrid(
            np.arange(grid_size[0]),
            np.arange(grid_size[1]),
            np.arange(grid_size[2]),
            indexing="ij",
        )
        voxel_indices = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    else:
        voxel_indices = np.argwhere(mask > 0)

    chunk_size = np.min([t_values.shape[0], 100])  # Set an appropriate chunk size based on available memory
    intersected = []

    for start in range(0, voxel_indices.shape[0], chunk_size):
        end = min(start + chunk_size, voxel_indices.shape[0])
        chunk = voxel_indices[start:end]

        intersected_chunk = compute_intersections_single_thread(ray_path, chunk, tol_distances, weights, voxel_size)
        intersected.extend(intersected_chunk)

    intersected = np.array(intersected)
    
    if intersected.size == 0:
        plot_flag = False
        
    if plot_flag:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection="3d")
        if mask is not None:
            indices = np.array(voxel_indices)
            ax.scatter(indices[:, 0], indices[:, 1], indices[:, 2], c="blue", s=5, alpha=0.03, label="Masked Voxels")

        ax.scatter(
            intersected[:, 0],
            intersected[:, 1],
            intersected[:, 2],
            c=intersected[:, 4] * 100,
            cmap="viridis",
            s=8,
            alpha=0.9,
            label="Intersected Voxels",
        )

        ax.plot(ray_path[:, 0]/voxel_size[0], ray_path[:, 1]/voxel_size[1], ray_path[:, 2]/voxel_size[2], "r-", label="Ray Path", linewidth=2)
        ax.set_xlim(0, grid_size[0])
        ax.set_ylim(0, grid_size[1])
        ax.set_zlim(0, grid_size[2])
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title("3D Ray-Voxel Intersection with Ray Size {:.2f}".format(ray_size))
        ax.legend()
        plt.show()
        
    return intersected


@njit
def compute_intersections_single_thread(ray_path, voxel_indices, tol_distances, weights, voxel_size):
    """
    sub-function for intersected_voxels_3d
    """
    intersected_chunk = []
    
    for i in range(voxel_indices.shape[0]):  # Process voxel indices sequentially
        x, y, z = voxel_indices[i]
        voxel_center = np.array([x + 0.5, y + 0.5, z + 0.5]) * voxel_size
        distances = np.sqrt(np.sum((ray_path - voxel_center) ** 2, axis=1))
        min_distance = np.min(distances)

        # Assign weights based on distance thresholds
        if min_distance <= tol_distances[-1]:
            for j, tol_distance in enumerate(tol_distances):
                if j == 0 and min_distance <= tol_distance:
                    weight = weights[j]
                elif j > 0:
                    if min_distance <= tol_distance and min_distance > tol_distances[j-1]:
                        weight = weights[j]
                    else:
                        continue
            intersected_chunk.append((x, y, z, min_distance, weight))
    
    return intersected_chunk


def make_projections_with_psf(fwd_peaks, omega_angles, image_size=(2162, 2068), detector='eiger', int_factors=(0.1065, 0.7807, 0.1065), sum_flag = False):
    """
    Make a series of projection images across different rotation angles for one rotation considering point spread function (psf).
    
    Args:
        fwd_peaks (numpy array): Forward computed peaks N*24 array.
        omega_angles (numpy array): List of rotation angles for projections. [deg]
        image_size (tuple): Size of the output image (rows, cols).
        detector (string): Detector name, "eiger" or "frelon".
        int_factors (tuple): Intensity distribution factors for projection left, middle and right across rotation angles.
    
    Returns:
        projs (array): a 3D array containing projection images [num_angles, rows, cols].
        projs_sum: a 2D array summing up all the projections along rotation angle
        
    3620 projections takes about 91 seconds
    """
    
    rot_step = omega_angles[1] - omega_angles[0]
    
    # Initialize 3D array for projections
    projs = np.zeros((len(omega_angles), image_size[0], image_size[1]), dtype='float')
    
    # Iterate over rotation angles with tqdm for progress tracking
    for i, rot_angle in enumerate(tqdm(omega_angles, desc="Generating projections over omega angles", unit="projection")):
        mask = (fwd_peaks[:, 9] >= rot_angle - 3 * rot_step) & (fwd_peaks[:, 9] <= rot_angle + 3 * rot_step)
        fwd_peaks_selected = fwd_peaks[mask]
        if fwd_peaks_selected.size > 0:
            projs[i, :, :] = make_one_projection_with_psf(
                fwd_peaks=fwd_peaks_selected,
                rot_angle=rot_angle,
                rot_step=rot_step,
                image_size=image_size,
                detector=detector,
                int_factors=int_factors
            )
    
    if sum_flag:
        projs_sum = np.sum(projs, axis = 0)
    else:
        projs_sum = None
    return projs, projs_sum


def make_one_projection_with_psf(fwd_peaks, rot_angle, rot_step = 0.05, image_size=(2162, 2068), detector = 'eiger', int_factors = (0.1065, 0.7807, 0.1065)):
    """
    Create one single projection image from position and intensity data with a point spread function (psf).
    Instead of directly applying psf on 3D image, we decompose the psf to 2D-filters across the neiboring projections:
    for psf_sigma = [1.0, 1.0, 0.5],
    psf [omega -] = [[0.008, 0.0132, 0.008], [0.0132, 0.0217, 0.0132],[0.008, 0.0132, 0.008]]     intensity fraction: 0.1065
    psf [omega]   = [[0.0591, 0.0975, 0.0591],[0.0975, 0.1607, 0.0975, 0.0591, 0.0975, 0.0591]]   intensity fraction: 0.7870
    psf [omega +] = [[0.008, 0.0132, 0.008], [0.0132, 0.0217, 0.0132],[0.008, 0.0132, 0.008]]     intensity fraction: 0.1065
    So, we directly distribute the sum_intensity over 3 frames across omega-, omega and omega+, corresponding to [0.1065, 0.7870, 0.1065] of the sum_intensity,
    while the shape of the filter can be remained as the same, i.e. Gaussian filter with a sigma of 1
    
    

    Args:
        fwd_peaks (numpy array): forward computed peaks N*24 array
        rot_angle (float): the ideal rotation angle corresponding to the projection no.
        rot_step (float): step size of the omega rotation [deg]
        image_size (tuple): size of the output image (rows, cols)
        detector (string): detector name, "eiger" or "frelon"
        int_factors (tuple): intensity distribution factors for projection left, middle and right across rotation angles
        
    Returns:
        image: generated projection [x, y]
        
    One projection takes about 0.005 seconds
    
    TODO:
    psf should also depend on the incident angle between the diffracted beam and the detector plane, i.e. more extended along radial direction
    """
    if 'eiger' in detector:
        psf_values = np.array([0.0591, 0.0975, 0.0591,
                               0.0975, 0.1607, 0.0975,
                               0.0591, 0.0975, 0.0591], dtype=np.float64) # eiger detector, psf 3*3, sigma = 1
        neigb_inds = np.array([[-1, -1], [0, -1], [1, -1],
                               [-1, 0], [0, 0], [1, 0],
                               [-1, 1], [0, 1], [1, 1]], dtype=np.int32)
    else:
        psf_values = np.array([0.0144, 0.0281, 0.0351, 0.0281, 0.0144,
                               0.0281, 0.0547, 0.0683, 0.0547, 0.0281,
                               0.0351, 0.0683, 0.0853, 0.0683, 0.0351,
                               0.0281, 0.0547, 0.0683, 0.0547, 0.0281,
                               0.0144, 0.0281, 0.0351, 0.0281, 0.0144], dtype=np.float64) # frelon detector, psf 5*5, sigma = 1.5
        neigb_inds = np.array([[-2, -2], [-1, -2], [0, -2], [1, -2], [2, -2],
                               [-2, -1], [-1, -1], [0, -1], [1, -1], [2, -1],
                               [-2, 0], [-1, 0], [0, 0], [1, 0], [2, 0],
                               [-2, 1], [-1, 1], [0, 1], [1, 1], [2, 1],
                               [-2, 2], [-1, 2], [0, 2], [1, 2], [2, 2]], dtype=np.int32)
    
    # Extract data
    ind_left = np.where((fwd_peaks[:, 9] >= rot_angle - rot_step*1.5) & (fwd_peaks[:, 9] < rot_angle - rot_step*0.5))[0]
    ind_middle = np.where((fwd_peaks[:, 9] >= rot_angle - rot_step*0.5) & (fwd_peaks[:, 9] < rot_angle + rot_step*0.5))[0]
    ind_right = np.where((fwd_peaks[:, 9] >= rot_angle + rot_step*0.5) & (fwd_peaks[:, 9] < rot_angle + rot_step*1.5))[0]
    
    positions = np.vstack([  np.vstack([fwd_peaks[ind_left, 18:20], fwd_peaks[ind_middle, 18:20]]), fwd_peaks[ind_right, 18:20]  ])
    intensities = np.hstack([  np.hstack([fwd_peaks[ind_left, 23]*int_factors[0], fwd_peaks[ind_middle, 23]*int_factors[1]]), fwd_peaks[ind_right, 23]*int_factors[2]  ])

#     # make it compatible with numba, but turned out to be slower --------------- to be investigated
#     positions = np.zeros((np.sum(ind_left) + np.sum(ind_middle) + np.sum(ind_right), 2), dtype=np.float64)
#     intensities = np.zeros(np.sum(ind_left) + np.sum(ind_middle) + np.sum(ind_right), dtype=np.float64)

#     # Fill positions and intensities
#     count = 0
#     for i in range(fwd_peaks.shape[0]):
#         if ind_left[i]:
#             positions[count] = fwd_peaks[i, 18:20]
#             intensities[count] = fwd_peaks[i, 23] * int_factors[0]
#             count += 1
#         elif ind_middle[i]:
#             positions[count] = fwd_peaks[i, 18:20]
#             intensities[count] = fwd_peaks[i, 23] * int_factors[1]
#             count += 1
#         elif ind_right[i]:
#             positions[count] = fwd_peaks[i, 18:20]
#             intensities[count] = fwd_peaks[i, 23] * int_factors[2]
#             count += 1

#     image = np.zeros(image_size, dtype = 'float')
#     # Accumulate intensities in the image, for clarity and compatibility
#     for i in range(positions.shape[0]):
#         x = int(positions[i, 0])
#         y = int(positions[i, 1])
#         intensity = intensities[i]
#         if 0 <= x < image_size[1] and 0 <= y < image_size[0]:
#             # distribute the intensity over neighboring pixels, equivalent to gaussian convolution but much faster
#             for j in range(len(neigb_inds)):
#                 x_neigh = x + neigb_inds[j, 0]
#                 y_neigh = y + neigb_inds[j, 1]
#                 if 0 <= x_neigh < image_size[1] and 0 <= y_neigh < image_size[0]:
#                     image[y_neigh, x_neigh] += intensity * psf_values[j]
            
    # Initialize the image
    image = np.zeros(image_size, dtype = 'float')
    
    # Map positions to pixel indices
    x_indices = (positions[:, 0]).astype(int)
    y_indices = (positions[:, 1]).astype(int)
    
    # Accumulate intensities in the image    
    for x, y, intensity in zip(x_indices, y_indices, intensities):
        if 0 <= x < image_size[1] and 0 <= y < image_size[0]:
            # distribute the intensity over neighboring pixels, equivalent to gaussian convolution but much faster
            for (neigb_ind, psf_value) in zip(neigb_inds, psf_values):
                pos_list = [x + neigb_ind[0], y + neigb_ind[1]]
                if 0 <= pos_list[0] < image_size[1] and 0 <= pos_list[1] < image_size[0]:
                    image[pos_list[1], pos_list[0]] += intensity * psf_value

    # Apply Gaussian filter for PSF, this is 17 times slower
    # image = gaussian_filter(image, sigma=psf_sigma)

    return image


def cellvolume(latticepar):
    """
    calculate unit cell volume [angstrom^3]
    """
    a = latticepar[0]
    b = latticepar[1]
    c = latticepar[2]
    calp = np.cos(np.deg2rad(latticepar[3]))
    cbet = np.cos(np.deg2rad(latticepar[4]))
    cgam = np.cos(np.deg2rad(latticepar[5]))
    
    angular = np.sqrt(1 - calp**2 - cbet**2 - cgam**2 + 2*calp*cbet*cgam)
    
    Vcell = a*b*c*angular
    
    return Vcell


def K1_calc():
    """
    calculate square of Thomson scattering lenth r0^2, [mm^2]
    K1 = 7.940788332083948e-24 [mm^2]
    """
    emass = 9.1093826e-31
    echarge = 1.60217653e-19
    pi4eps0 = 1.11265e-10
    c = 299792458
    K1 = (echarge**2/(pi4eps0 * emass * c * c) *1000)**2  # [mm^2]
    return K1

    
def find_ds_for_multiple_targets(Ahkls, target_hkls):
    """
    Find the ds for a list of target hkls in the list Ahkls.

    Args:
        Ahkls (list): A list of [structure_factor, hkl], e.g. Ahkls = ucell.gethkls(ds_max)
        target_hkls (list): a list of target hkl to search for, e.g. [[0, -2, 2], [0, 2, -2]]

    Returns:
        list: A list of ds values for each target hkl (or 0 if not found).
    """
    # Ensure target_hkls is a 2D list
    target_hkls = ensure_2d_array(target_hkls)

    # Convert Ahkls to a dictionary for fast lookups
    Ahkls_dict = {tuple(hkl): ds for ds, hkl in Ahkls}

    return [Ahkls_dict.get(tuple(target_hkl), 0.0) for target_hkl in target_hkls]


def ensure_2d_array(target_hkls):
    """
    Ensures the input target_hkls is a 2D array.

    Args:
        target_hkls (list or array): A single 1D array or a list of N x 3 arrays.

    Returns:
        list: A 2D list of N x 3 arrays.
    """
    if isinstance(target_hkls[0], (int, float)):
        return [target_hkls]
    return target_hkls


def plot_fwd_peaks(fwd_peaks):
    assert fwd_peaks.shape[1] == 25, "fwd_peaks shapes are not qualified"
    f, ax = plt.subplots(1, 2, figsize=(15, 9))

    sc = ax[0].scatter(fwd_peaks[:, 18], fwd_peaks[:, 19], c=fwd_peaks[:, 23], cmap='viridis', s=8)
    ax[0].set_aspect('equal', 'box')
    cb = f.colorbar(sc, ax=ax[0])
    # cb.set_label('Intensity', fontsize = 20)
    cb.ax.tick_params(labelsize=14)
    ax[0].set_xlabel('fc (pixel)', fontsize = 20)
    ax[0].set_ylabel('sc (pixel)', fontsize = 20)
    ax[0].set_title('(a) Forward peaks on detector', fontsize = 20)
    ax[0].tick_params(width=1.5, length=6, labelsize=14)
    ax[0].invert_yaxis()

    sc = ax[1].scatter(fwd_peaks[:, 5], fwd_peaks[:, 6], c=fwd_peaks[:, 23], cmap='viridis', s=8)
    ax[1].set_aspect('equal', 'box')
    cb = f.colorbar(sc, ax=ax[1])
    cb.set_label('Intensity', fontsize = 20)
    cb.ax.tick_params(labelsize=14)
    ax[1].set_xlabel('X (mm)', fontsize = 20)
    ax[1].set_ylabel('Y (mm)', fontsize = 20)
    ax[1].set_title('(b) Forward peaks from the sample', fontsize = 20)
    ax[1].tick_params(width=1.5, length=6, labelsize=14) 

    f.tight_layout()
    plt.show()
    
    
def make_intensity_map(x_positions, y_positions, intensities, x_range = None, y_range = None, pixel_size=1e-3, plot_flag=True, cmap="viridis"):
    """
    Generate and plot an intensity map based on input x, y positions and the associated intensities.
    
    Args:
        x_positions (array): x positions
        y_positions (array): y positions
        intensities (array): intensities
        x_range (float): min and max x range values
        y_range (float): min and max y range values
        pixel_size (float): pixel size for generating the image
        plot_flag (logic): for making a plot
        cmap (str): color map for plotting
    
    Returns:
        intensity_map
    """
    assert (x_positions.size == y_positions.size) and (x_positions.size == intensities.size), "The sizes of x, y, intensities must be the same"

    # Compute the bounds of the grid
    if x_range is None:
        x_min, x_max = np.min(x_positions), np.max(x_positions)
    else:
        x_min = x_range[0]
        x_max = x_range[1]
    if y_range is None:
        y_min, y_max = np.min(y_positions), np.max(y_positions)
    else:
        y_min = y_range[0]
        y_max = y_range[1]

    # Determine the grid size
    x_bins = int(np.ceil((x_max - x_min) / pixel_size))
    y_bins = int(np.ceil((y_max - y_min) / pixel_size))

    intensity_map = np.zeros((y_bins, x_bins))

    # Map positions to pixel indices
    x_indices = np.floor((x_positions - x_min) / pixel_size).astype(int)
    y_indices = np.floor((y_positions - y_min) / pixel_size).astype(int)

    # Accumulate intensities
    for x_idx, y_idx, intensity in zip(x_indices, y_indices, intensities):
        intensity_map[y_idx, x_idx] += intensity

    if plot_flag:
        plt.figure(figsize=(10, 8))
        plt.imshow(intensity_map, origin='lower', extent=[x_min, x_max, y_min, y_max], cmap=cmap,
                  norm=LogNorm(vmin=max(10, intensity_map.min()), vmax=intensity_map.max()), interpolation="nearest")
        # plt.imshow(intensity_map, origin='lower', extent=[x_min, x_max, y_min, y_max], cmap=cmap)
        plt.colorbar(label="Intensity")
        plt.xlabel("X position")
        plt.ylabel("Y position")
        plt.title("Intensity Map, pixel size = {}".format(pixel_size))
        plt.show()
    
    return intensity_map


"""
The following functions are related to segment the peaks from the projections, which are all adapted from ImageD11.sinograms.lima_segmenter
"""

def get_opts_seg(mask = None, bg = None, cut = 1, pixels_in_spot = 1):
    """
    set segmentation parameters
    Args:
        mask (2D array): mask for the detector
        bg (2D array):   background for the detector
        cut (float):     keep values abuve cut in first look at image
        pixels_in_spot:  minimum number of pixels in a single spot
    """
    thresholds = tuple([cut * pow(2, i) for i in range(6)])
    opts_seg = {"mask": mask,
        "bg": bg,
        "cut": cut,           # keep values abuve cut in first look at image
        "overwrite": False,   # overwrite existing sparse h5 file
        "howmany": 100000,    # max pixels per frame to keep
        "thresholds": thresholds,
        "pixels_in_spot": pixels_in_spot}
    return opts_seg


def make_projs_and_sparse_file(fwd_peaks, destname, omega_angles, opts_seg, detector = 'eiger', detector_mask = None, image_size=(2162, 2068), int_factors=(0.1065, 0.7807, 0.1065)): 
    """
    Make projection image and do segmentation on it to generate sparse peaks h5 file
    Combine make_one_projection_with_psf and adaptation of "segment_lima" in ImageD11.sinograms.lima_segmenter
    The advantage with this is processing frame by frame, instead of producing a series of frames and then doing sementation (slower and more memory consuming) 
    Args:
        fwd_peaks (numpy array): Forward computed peaks N*24 array.
        dsetname (str):  output filename
        omega_angles (numpy array): List of rotation angles for projections. [deg]
        opts_seg (dict): segmenting options, e.g.:
        opts_seg = {"mask": None,
            "bg": None,
            "cut": 1,           # keep values abuve cut in first look at image
            "overwrite": False,
            "howmany": 100000,  # max pixels per frame to keep
            "thresholds": 55,
            "pixels_in_spot": 1}
        detector (str):  detector name, e.g. "eiger" or "frelon3"
        detector_mask (2D array):  detector mask image, 1-active, 0-non-active
        image_size (int array): dimensions of the 2D projection image
        int_factors (tuple): Intensity distribution factors for projection left, middle and right across rotation angles.
    Returns:
        dsetname (str):  output filename with saved sparse peaks
    """
    
    dataset = "/entry_0000/ESRF-ID11/" + detector + "/data"
    # saving compression style:
    opts = {
        "chunks": (10000,),
        "maxshape": (None,),
        "compression": "lzf",
        "shuffle": True,
    }
    assert fwd_peaks.size > 0, "fwd_peaks must not be empty"
    if detector_mask is None:
        detector_mask = np.ones(image_size, dtype = 'float')
    opts_seg["mask"] = detector_mask
    
    if os.path.exists(destname) and opts_seg['overwrite']:
        os.remove(destname)
    elif os.path.exists(destname):
        print('{} exists already; skipping ...'.format(destname))
        return destname
    
    rot_step = omega_angles[1] - omega_angles[0]
    dty_values = d = np.full((omega_angles.shape[0],), fwd_peaks[0, 8])   # save dty as the same shape as omega_angles
    start = time.time()
    with h5py.File(destname, "w") as hout:
        print("# ", fwd_peaks.shape, destname)
        print("# time now", time.ctime(), "\n#")
        g = hout.require_group(dataset)
        g.create_dataset("rot", data = omega_angles, dtype='float')
        g.create_dataset("rot_center", data = omega_angles, dtype='float')
        g.create_dataset("dty", data = dty_values, dtype='float')
        row = g.create_dataset("row", (0,), dtype=np.uint16, **opts)
        col = g.create_dataset("col", (0,), dtype=np.uint16, **opts)
        sig = g.create_dataset("intensity", (0,), dtype=np.float64, **opts)
             
        nframes = omega_angles.shape[0]
        nnz = g.create_dataset("nnz", (nframes,), dtype=np.uint32)
        g.attrs["itype"] = np.dtype(np.uint16).name
        g.attrs["nframes"] = nframes
        g.attrs["shape0"] = image_size[0]
        g.attrs["shape1"] = image_size[1]
        
        npx = 0
        # Iterate over rotation angles with tqdm for progress tracking
        for i, rot_angle in enumerate(tqdm(omega_angles, desc="Making projection and segmenting sparse peaks over omega angles", unit="projection")):
            peaks_mask = (fwd_peaks[:, 9] >= rot_angle - 3 * rot_step) & (fwd_peaks[:, 9] <= rot_angle + 3 * rot_step)
            fwd_peaks_selected = fwd_peaks[peaks_mask]
            if fwd_peaks_selected.size > 0:
                frm = make_one_projection_with_psf(
                    fwd_peaks=fwd_peaks_selected,
                    rot_angle=rot_angle,
                    rot_step=rot_step,
                    image_size=image_size,
                    detector=detector,
                    int_factors=int_factors
                )
                fun = lima_segmenter.frmtosparse(detector_mask, frm.dtype)
                if opts_seg['bg'] is not None:
                    frm = frm.astype(np.float32) - opts_seg['bg']
                frm_npx, frm_row, frm_col, frm_val = fun(frm, opts_seg['cut'])
                spf = clean_frms(frm_npx, frm_row, frm_col, frm_val, opts_seg)
                if spf is None or spf.nnz == 0:
                    nnz[i] = 0
                    continue
                if spf.nnz + npx > len(row):
                    row.resize(spf.nnz + npx, axis=0)
                    col.resize(spf.nnz + npx, axis=0)
                    sig.resize(spf.nnz + npx, axis=0)
            
                row[npx:] = spf.row[:]
                col[npx:] = spf.col[:]
                sig[npx:] = spf.pixels["intensity"]
                nnz[i] = spf.nnz
                npx += spf.nnz
        g.attrs["npx"] = npx
        
    end = time.time()
    print("\n# Done", nframes, "frames", npx, "pixels  fps", nframes / (end - start))
    return destname


def segment_frms(frms, destname, opts_seg, detector = 'eiger'):
    """
    Adapted from "segment_lima" in ImageD11.sinograms.lima_segmenter
    Does segmentation on a series of frames, e.g. generated from fwd_peaks
    Modified based on the function segment_lima in ImageD11/sinograms/lima_segmenter.py
    Args:
        frms (3D array): projections generated from fwd_peaks using forward_projector.make_projections_with_psf
        dsetname (str):  output filename
        detector (str):  detector name, e.g. "eiger" or "frelon3"
        opts_seg (dict): segmenting options, e.g.:
        opts_seg = {"mask": None,
            "bg": None,
            "cut": 1,           # keep values abuve cut in first look at image
            "overwrite": False,
            "howmany": 100000,  # max pixels per frame to keep
            "thresholds": 55,
            "pixels_in_spot": 1}
    Returns:
        dsetname (str):  output filename
    """
    dataset = "/entry_0000/ESRF-ID11/" + detector + "/data"
    # saving compression style:
    opts = {
        "chunks": (10000,),
        "maxshape": (None,),
        "compression": "lzf",
        "shuffle": True,
    }
    assert len(frms.shape) == 3, "frms must have 3D shape"
    if opts_seg["mask"] is None:
        opts_seg["mask"] = np.ones((frms.shape[1], frms.shape[2]), dtype = "float")
    detector_mask = opts_seg['mask']
    
    if os.path.exists(destname) and opts_seg['overwrite']:
        os.remove(destname)
    elif os.path.exists(destname):
        print('{} exists already; skipping ...'.format(destname))
        return destname
    start = time.time()
    with h5py.File(destname, "a") as hout:
        print("# ", frms.shape, destname)
        print("# time now", time.ctime(), "\n#")
        g = hout.require_group(dataset)
        row = g.create_dataset("row", (0,), dtype=np.uint16, **opts)
        col = g.create_dataset("col", (0,), dtype=np.uint16, **opts)
        # can go over 65535 frames in a scan
        # num = g.create_dataset("frame", (1,), dtype=np.uint32, **opts)
        sig = g.create_dataset("intensity", (0,), dtype=frms.dtype, **opts)
        nnz = g.create_dataset("nnz", (frms.shape[0],), dtype=np.uint32)
        g.attrs["itype"] = np.dtype(np.uint16).name
        g.attrs["nframes"] = frms.shape[0]
        g.attrs["shape0"] = frms.shape[1]
        g.attrs["shape1"] = frms.shape[2]
        npx = 0
        nframes = frms.shape[0]
        for i, spf in enumerate(reader_frms(frms, detector_mask, opts_seg)):
            if i % 100 == 0:
                if spf is None:
                    print("%4d 0," % (i))
                else:
                    print("%4d %d," % (i, spf.nnz))
                sys.stdout.flush()
            if spf is None:
                nnz[i] = 0
                continue
            if spf.nnz + npx > len(row):
                row.resize(spf.nnz + npx, axis=0)
                col.resize(spf.nnz + npx, axis=0)
                sig.resize(spf.nnz + npx, axis=0)
            row[npx:] = spf.row[:]
            col[npx:] = spf.col[:]
            sig[npx:] = spf.pixels["intensity"]
            # num[npx:] = i
            nnz[i] = spf.nnz
            npx += spf.nnz
        g.attrs["npx"] = npx
    end = time.time()
    print("\n# Done", nframes, "frames", npx, "pixels  fps", nframes / (end - start))
    return destname

    # the output file should be flushed and closed when this returns
    
    
def reader_frms(frms, mask, opts_seg, start=0):
    """
    Adapted from "reader" in ImageD11.sinograms.lima_segmenter
    iterator to read frames and segment them
    returns sparseframes
    
    Args:
        frms (3D array): projections to be segmented
        mask (2D array): mask image for the detector
        opts_seg (dict): segmenting options, e.g.:
        opts_seg = {"mask": None,
            "bg": None,
            "cut": 1,           # keep values abuve cut in first look at image
            "overwrite": False,
            "howmany": 100000,  # max pixels per frame to keep
            "thresholds": 55,
            "pixels_in_spot": 1}
        start (int): first indice for projections to be segmented
    Returns:
        dsetname (str):  output filename
    """
    
    assert start < len(frms)

    fun = lima_segmenter.frmtosparse(mask, frms.dtype)
    for i in range(start, frms.shape[0]):
        frm = frms[i]
        if opts_seg['bg'] is not None:
            frm = frm.astype(np.float32) - opts_seg['bg']
        npx, row, col, val = fun(frm, opts_seg['cut'])
        spf = clean_frms(npx, row, col, val, opts_seg)
        yield spf


def clean_frms(nnz, row, col, val, opts_seg):
    """
    Adapted from "clean" in ImageD11.sinograms.lima_segmenter
    """
    if nnz == 0:
        return None
    if nnz > opts_seg['howmany']:
        nnz = lima_segmenter.top_pixels(nnz, row, col, val, opts_seg['howmany'], opts_seg['thresholds'])
        # Now get rid of the single pixel 'peaks'
        #   (for the mallocs, data is copied here)
        s = lima_segmenter.sparseframe.sparse_frame(
            row[:nnz].copy(), col[:nnz].copy(), opts_seg['mask'].shape
        )
        s.set_pixels("intensity", val[:nnz].copy())
    else:
        s = lima_segmenter.sparseframe.sparse_frame(row, col, opts_seg['mask'].shape)
        s.set_pixels("intensity", val)
    if opts_seg['pixels_in_spot'] <= 1:
        return s
    # label them according to the connected objects
    s.set_pixels("f32", s.pixels["intensity"].astype(np.float32))
    npk = lima_segmenter.sparseframe.sparse_connected_pixels(
        s, threshold=0, data_name="f32", label_name="cp"
    )
    # only keep spots with more than 3 pixels ...
    #    mom = sparseframe.sparse_moments( s,
    #                                     intensity_name="f32",
    #                                     labels_name="cp" )
    #    npx = mom[:, cImageD11.s2D_1]
    npx = np.bincount(s.pixels["cp"], minlength=npk)
    pxcounts = npx[s.pixels["cp"]]
    pxmsk = pxcounts >= opts_seg['pixels_in_spot']
    if pxmsk.sum() == 0:
        return None
    sf = s.mask(pxmsk)
    return sf


def assemble_sparsefiles(sparsefiles_folder, dtys, outname_sparse = "fwd_sparse.h5", image_size=(2162, 2068)):
    """
    assemble sparsefiles from each dty position to a single sparse file
    adapted from "harvest_masterfile" in ImageD11.sinograms.assemble_label

    Args:
        sparsefiles_folder (str): folder holds individual dty sparse files
        dtys (array):  a list of dty positions [um]
        outname_sparse (str): sparse file to write
        image_size (int array): dimensions of the 2D projection image
    Returns:
        outname_sparse (str): sparse file to write
    """
    if os.path.exists(outname_sparse):
        print('{} already exists; exit ...'.format(outname_sparse))
        return outname_sparse
    scanmotors = ['dty', 'rot', 'rot_center']
    headermotors = ['dty', 'rot']
    with h5py.File(outname_sparse, "a") as hout:
        for i, dty in enumerate(dtys):
            scan = str(i+1) + '.1'
            sparsefile = os.path.join(sparsefiles_folder, 'fsparse_dty_' + str(round(dty, 2)).replace('.', 'p') + '.h5')
            if not os.path.exists(sparsefile):
                print('{} not found ...'.format(sparsefile))
                continue
            else:
                fsparse_pks = io.read_fsparse(sparsefile)
                # write basic motor info
                g = hout.require_group(scan)
                gm = g.require_group("measurement")
                for m in scanmotors:
                    if m in fsparse_pks.keys():
                        data = fsparse_pks[m][:]
                        ds = gm.require_dataset(m, shape=data.shape, dtype=data.dtype)
                        ds[()] = data
                gip = g.require_group("instrument/positioners")
                for m in headermotors:
                    if m in fsparse_pks.keys():
                        data = fsparse_pks[m][:]
                        ds = gip.require_dataset(m, shape=data.shape, dtype=data.dtype)
                        ds[()] = data
                g.attrs["itype"] = fsparse_pks['intensity'].dtype.name
                g.attrs["nframes"] = fsparse_pks['nnz'].shape
                g.attrs["shape0"] = image_size[0]
                g.attrs["shape1"] = image_size[1]

                # write peaks info
                g.create_dataset("row", data = fsparse_pks['row'], dtype=np.uint16)
                g.create_dataset("col", data = fsparse_pks['col'], dtype=np.uint16)
                g.create_dataset("nnz", data = fsparse_pks['nnz'], dtype=np.uint32)
                g.create_dataset("intensity", data = fsparse_pks['intensity'], dtype=g.attrs["itype"])
            print(scan, ", ")
        print()
    return outname_sparse


def merge_rows(data, columns_to_cluster=[0, 9, 12, 13, 14], epsilons=[0.1, 2, 0.1, 0.1, 0.1],
               special_col_index=23, weighted_avg=True, min_samples=1):
    """
    Merge rows based on proximity of specific columns with column-specific epsilon values.
    The default values are suitable for fwd_peaks, where columns_to_cluster=[0, 9, 12, 13, 14, 18, 19] correspond to
    grainID, rot, hkl, dety_pixel(fc), detz_pixel(sc)

    Args:
        data (ndarray): Input array with shape (n_rows, n_cols).
        columns_to_cluster (list): Indices of columns to check for closeness.
        epsilons (list): List of epsilon values for each column in columns_to_cluster.
        special_col_index (int): Index of the special column to sum instead of averaging.
        weighted_avg (bool): Whether to use weighted averaging (weights from a specific column).
        min_samples (int): Minimum samples to form a cluster.

    Returns:
        ndarray: Array with merged rows.
        clusters: cluster ID
    """
    # Input validation
    if not isinstance(data, np.ndarray) or data.ndim != 2:
        raise ValueError("Input data must be a 2D numpy array.")
    if special_col_index >= data.shape[1]:
        raise ValueError("special_col_index is out of bounds.")
    if not all(0 <= col < data.shape[1] for col in columns_to_cluster):
        raise ValueError("Invalid column indices in columns_to_cluster.")
    if len(columns_to_cluster) != len(epsilons):
        raise ValueError("columns_to_cluster and epsilons must have the same length.")

    # Normalize the columns to cluster by their respective epsilon
    subset = data[:, columns_to_cluster]
    normalized_subset = subset / np.array(epsilons)

    # Perform clustering using DBSCAN
    dbscan = DBSCAN(eps=1.0, min_samples=min_samples, metric='euclidean')  # `eps=1.0` because of normalization
    clusters = dbscan.fit_predict(normalized_subset)

    # Merge rows based on clusters
    merged_data = []
    for cluster_id in np.unique(clusters):
        if cluster_id == -1:  # Noise points (optional: handle them separately)
            continue
        cluster_rows = data[clusters == cluster_id]

        if not weighted_avg:
            # Average all columns except the special column
            averaged = np.mean(cluster_rows[:, :], axis=0)
        else:
            # Use weighted average (weights from the special column)
            weights = cluster_rows[:, special_col_index]
            averaged1 = np.average(cluster_rows[:, 0:special_col_index], axis=0, weights=weights)
            averaged2 = np.average(cluster_rows[:, special_col_index+1:], axis=0, weights=weights)

        # Sum the special column
        special_col_sum = np.sum(cluster_rows[:, special_col_index], axis=0)

        # Combine averaged and summed columns
        merged_row = np.append(averaged1, special_col_sum)
        merged_row = np.append(merged_row, averaged2)
        merged_data.append(merged_row)

    # Convert merged data to a numpy array
    merged_data = np.array(merged_data)

    print("Original shape:", data.shape)
    print("Merged shape:", merged_data.shape)
    return merged_data, clusters


def get_ImageD11_code_path():
    """
    get the code path of ImageD11
    """
    return os.path.split(os.path.split(os.path.split(inspect.getfile(forward_projector))[0])[0])[0]


def set_tmp_dir(tmp_dir = "/tmp_14_days/slurm/log"):
    """
    set a temporary directory for multiprocessing to a writable location
    """
    tempdir = tempfile.mkdtemp(dir=tmp_dir)
    multiprocessing.set_start_method("fork", force=True)  # Use 'fork' for better HPC support
    os.environ["TMPDIR"] = tempdir  # Ensure workers use the correct temp dir
