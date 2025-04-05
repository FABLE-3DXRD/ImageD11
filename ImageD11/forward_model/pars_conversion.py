# pars_conversion.py
# convert parameters defined in different codes to a common parameter dictionary (p), which are compatible with fwd-DCT, this includes:
# 1) convert ImageD11 pars to p
# 2) convert DCT parameters.mat to p
# 3) convert pyFAI poni to p
# 4) build_ImageD11par_from_poni
# 5) build_poni_from_ImageD11par
# Haixing Fang, haixing.fang@esrf.fr
# Oct 27th, 2024
# Updates on Jan 19th, 2025: introduce logging

import os
import numpy as np
import ImageD11.transformer
import pyFAI
from datetime import datetime
from ImageD11.forward_model.io import read_matlab_file, read_h5_file
from scipy.optimize import minimize
import logging

econst = 12.3984
logging.basicConfig(level=logging.INFO, force=True)

def convert_ImageD11pars2p(pars, verbose = 1):
    """
    Convert ImageD11 parameters to specific geometry and detector parameters (unit in mm), which are compatible with fwd-DCT
    
    Arguments:
    pars -- ImageD1.parameters object, e.g. pars = ImageD11.parameters.read_par_file(ds.parfile)
    
    Returns:
    a dictionary of p, which has keys of rawfile, Lsam2det, tilt_xyz, RotDet, dety0, detz0, dety00, detz00, wedge, pixelysize, pixelzsize,
                                         detysize, detzsize, BeamStopY, BeamStopZ, S, flip_lr_flag, flip_ud_flag, wavelength, Energy etc.
    """
    p = {}
    try:
        p['rawfile'] = pars.parameters['filename']
    except Exception as e:
        # sometime there is no "filename" key
        if verbose > 1:
            logging.info("An unexpected error occurred: {}".format(e))
    p['Lsam2det'] = pars.parameters['distance']/1000.0
    p['tilt_xyz'] = [np.rad2deg(pars.parameters['tilt_x']), np.rad2deg(pars.parameters['tilt_y']), np.rad2deg(pars.parameters['tilt_z'])]  # [deg]
    p['RotDet'] = get_det_R(pars.parameters['tilt_x'], pars.parameters['tilt_y'], pars.parameters['tilt_z'], verbose)
    p['pixelysize'] = pars.parameters['y_size']/1000.0
    p['pixelzsize'] = pars.parameters['z_size']/1000.0
    if pars.parameters['o11'] == -1 and pars.parameters['o22'] == -1:
        p['detysize'] = 2162
        p['detzsize'] = 2068
        p['flip_lr_flag'] = False
        p['flip_ud_flag'] = False
    elif pars.parameters['o11'] == 1 and pars.parameters['o22'] == -1:
        p['detysize'] = 2048
        p['detzsize'] = 2048
        p['flip_lr_flag'] = False
        p['flip_ud_flag'] = True
    else:
        if verbose >= 1:
            logging.info('Note this detector configuration is not supported yet; Set to not have any flips - equivalent to Eiger !')
        p['detysize'] = 2162
        p['detzsize'] = 2068
        p['flip_lr_flag'] = False
        p['flip_ud_flag'] = False
    p['dety0'] = pars.parameters['y_center']
    p['detz0'] = pars.parameters['z_center']
    p['dety00'] = 0.0
    p['detz00'] = 0.0
    p['wedge'] = pars.parameters['wedge']
    if pars.parameters['omegasign'] > 0:
        p['S'] = np.array([[1.0, 0, 0],
             [0, 1.0, 0],
             [0, 0, 1.0]])
    else:
        p['S'] = np.array([[1.0, 0, 0],
             [0, -1.0, 0],
             [0, 0, 1.0]])
    p['wavelength'] = pars.parameters['wavelength'] # [Angstrom]
    p['Energy'] = econst / pars.parameters['wavelength'] # [keV]
    bs_margin = 2.5/(pars.parameters['y_size']/1000.0)   # assume a 5 mm beamstop [pixel]
    p['BeamStopY'] = [pars.parameters['y_center']-bs_margin, pars.parameters['y_center']+bs_margin]
    p['BeamStopZ'] = [pars.parameters['z_center']-bs_margin, pars.parameters['z_center']+bs_margin]
    
    return p
    

def convert_DCTpars2p(par_path, Binning = None, verbose = 1):
    """
    Convert DCT parameters.mat to specific geometry and detector parameters (unit in mm), which are compatible with fwd-DCT
    
    Arguments:
    par_path -- DCT parameters.mat
    Binning  -- binning of the detector, None by default
    Equivalent to get_para_from_parameters.m in DCT code
    
    Returns:
    a dictionary of p, which has keys of rawfile, Lsam2det, tilt_xyz, RotDet, dety0, detz0, dety00, detz00, wedge, pixelysize, pixelzsize,
                                         detysize, detzsize, BeamStopY, BeamStopZ, S, flip_lr_flag, flip_ud_flag, wavelength, Energy etc.
    """
    if par_path.endswith('.mat'):
        mat_data = read_matlab_file(par_path)
        if 'parameters' in mat_data.keys():
            parameters = mat_data['parameters']  # old parameters file version: parameters saved as a one single struct
        else:
            parameters = mat_data.copy()         # new parameters file version: parameters saved as multiple structs, i.e. acq, labgeo, detgeo, rec, index etc.
    elif par_path.endswith('.h5'):
        parameters = read_h5_file(par_path)
    else:
        raise Exception('At present, the input parameters file must be in .h5 or .mat format')
    
    p = {}
    p['rawfile'] = par_path
    if len(parameters['acq']['energy']) == 1:
        p['Energy'] = parameters['acq']['energy']  # [keV]
        p['wavelength'] = econst / parameters['acq']['energy']  # [Angstrom]
    else:
        p['Energy'] = parameters['acq']['energy'][0]  # [keV]
        p['wavelength'] = econst / parameters['acq']['energy'][0]  # [Angstrom]
    
    # detector information
    p['pixelysize'] = parameters['detgeo']['pixelsizeu'][0]
    p['pixelzsize'] = parameters['detgeo']['pixelsizev'][0]
    p['dety0']      = parameters['detgeo']['detrefu'][0]
    p['detz0']      = parameters['detgeo']['detrefv'][0]
    p['detysize']   = parameters['detgeo']['detsizeu'][0]
    p['detzsize']   = parameters['detgeo']['detsizev'][0]
    if len(parameters['acq']['bbdir']) == 4:
        p['BeamStopY']  = [ parameters['acq']['bbdir'][0], parameters['acq']['bbdir'][0] + parameters['acq']['bbdir'][2]]
        p['BeamStopZ']  = [ parameters['acq']['bbdir'][1], parameters['acq']['bbdir'][1] + parameters['acq']['bbdir'][3]]
    else:
        p['BeamStopY']  = [ parameters['acq']['bbdir'][0][0], parameters['acq']['bbdir'][0][0] + parameters['acq']['bbdir'][0][2]]
        p['BeamStopZ']  = [ parameters['acq']['bbdir'][0][1], parameters['acq']['bbdir'][0][1] + parameters['acq']['bbdir'][0][3]] 
    p['Qdet']       = parameters['detgeo']['Qdet']    
    
    # acquisition / geometry information
    VoxSize         = parameters['recgeo']['voxsize'] # [mm]
    detrefpos       = parameters['detgeo']['detrefpos']
    detorig         = parameters['detgeo']['detorig'] # top left corner of the detector [mm]
    acquisition_center=[0, 0, 0]
    det_center      = detrefpos
    # Note: shift in y is not needed, because the setup for tomo recon has been forced to bring y centered. shift Z needs to be considered
    if len(parameters['acq']['bb']) == 4:
        det_center[2]   = det_center[2] + (parameters['acq']['bb'][1] + parameters['acq']['bb'][3]/2 - p['detz0'])*p['pixelzsize'] # because parameters.acq.bb defines FOV for tomo recon [mm]
    else:
        det_center[2]   = det_center[2] + (parameters['acq']['bb'][0][1] + parameters['acq']['bb'][0][3]/2 - p['detz0'])*p['pixelzsize'] 
    if det_center[2] != detrefpos[2]:
        print('Shift in vertical axis (Z) is found: det_center changes from {:.5f} to {:.5f} mm'.format(detrefpos[2], det_center[2]))
    p['Lsam2det']   = detrefpos[0] # more precise than parameters.acq.dist
    p['dety00']     = det_center[1] - acquisition_center[1] # [mm]
    p['detz00']     = det_center[2] - acquisition_center[2] # [mm]
    p['RotAxisOffset'] = detrefpos[1]                       # not used but save this value [mm]
    
    det_normal      = parameters['detgeo']['detnorm'].ravel()
    if np.abs(det_normal[0]+1) < 0.5:
        det_normal = -det_normal # should be along the beam
        if verbose >= 1:
            logging.info('Note: the original parameters.detgeo.detnorm is opposite to the beam direction.')
            logging.info('Now det_normal is along the beam direction')
    detdiru = parameters['detgeo']['detdiru'].ravel()
    if np.abs(detdiru[1] - 1) < 0.5:
        detdiru = -detdiru
        p['Qdet'][0,:] = -p['Qdet'][0,:]
        flip_lr_flag = True
        if verbose >= 1:
            logging.info('Note: the original parameters.detgeo.detdiru is pointing to the left')
            logging.info('Now detdiru points to the right => flip_lr_flag = True')
    else:
        flip_lr_flag = False
    p['flip_lr_flag'] = flip_lr_flag
    
    detdirv = parameters['detgeo']['detdirv'].ravel()
    if np.abs(detdirv[2] - 1) < 0.5:
        detdirv = -detdirv
        p['Qdet'][1,:] = -p['Qdet'][1,:]
        flip_lr_flag = True
        if verbose >= 1:
            logging.info('Note: the original parameters.detgeo.detdirv is pointing to upwards')
            logging.info('Now detdirv points downwards')
    p['flip_ud_flag'] = False
                                 
    tilt_x, tilt_y, tilt_z = get_det_tilt(det_normal, detdiru, degree = False)
    p['tilt_xyz'] = [np.rad2deg(tilt_x), np.rad2deg(tilt_y), np.rad2deg(tilt_z)] # [deg]
    p['RotDet'] = get_det_R(tilt_x, tilt_y, tilt_z, verbose)           
    p['wedge'] = 0.0
                                 
    # sample information
    if len(parameters['acq']['dir']) == 2:
        p['FileFolder'] = parameters['acq']['dir'][0]
    else:
        p['FileFolder'] = parameters['acq']['dir']
    
    if len(parameters['acq']['refon']) > 1:
        parameters['acq']['refon'] = parameters['acq']['refon'][0]
    if not ('labgeo' in parameters and 'omstep' in parameters['labgeo']):
        if parameters['acq']['type'] == '360degree':
            parameters['labgeo']['omstep'] = 360.0/parameters['acq']['refon']
        else:
            parameters['labgeo']['omstep'] = 180.0/parameters['acq']['refon']
    p['rot_step']  = parameters['labgeo']['omstep']
    if parameters['acq']['type'] == '360degree':
        p['rot_start'] = 0.0
        p['rot_end']   = 360.0
        p['rot_angles'] = np.arange(p['rot_start'], p['rot_end'], p['rot_step'])
    if len(parameters['acq']['name']) == 2:
        p['prefix']   = parameters['acq']['name'][0]
    else:
        p['prefix']   = parameters['acq']['name']
    p['SampleName'] = p['prefix']
    p['ImageNr']  = parameters['acq']['refon']
                                 
    if Binning is not None and Binning != 1:
        p['pixelysize'] = p['pixelysize']*Binning
        p['pixelzsize'] = p['pixelzsize']*Binning
        p['dety0']      = p['dety0']/Binning
        p['detz0']      = p['detz0']/Binning
        p['detysize']   = p['detysize']/Binning
        p['detzsize']   = p['detzsize']/Binning
        p['BeamStopY'][0]   = np.floor(p['BeamStopY'][0]/Binning)
        p['BeamStopY'][1]   = np.ceil(p['BeamStopY'][1]/Binning)
        p['BeamStopZ'][0]   = np.floor(p['BeamStopZ'][0]/Binning)
        p['BeamStopZ'][1]   = np.ceil(p['BeamStopZ'][1]/Binning)
        if verbose >= 1:
            logging.info('Detector images are binned to {:.1f} * {:.1f}'.format(Binning, Binning))
    return p
    
    
def get_det_R(tilt_x, tilt_y, tilt_z, verbose = 1):
    """
    Calculate detector rotation matrix from tilt angles around x, y, z axes in radians
    Arguments:
    tilt_x, tilt_y, tilt_z in radians
    Returns:
    RotDet, 3*3 numpy array
    """
    if verbose >= 2:
        logging.debug('Detector tilt angles = [{:.4f} {:.4f} {:.4f}] deg'.format(np.rad2deg(tilt_x), np.rad2deg(tilt_y), np.rad2deg(tilt_z)))
    # Define rotation matrices
    RotX = np.array([[1, 0, 0], 
                     [0, np.cos(tilt_x), -np.sin(tilt_x)], 
                     [0, np.sin(tilt_x), np.cos(tilt_x)]])
    
    RotY = np.array([[np.cos(tilt_y), 0, np.sin(tilt_y)], 
                     [0, 1, 0], 
                     [-np.sin(tilt_y), 0, np.cos(tilt_y)]])
    
    RotZ = np.array([[np.cos(tilt_z), -np.sin(tilt_z), 0], 
                     [np.sin(tilt_z), np.cos(tilt_z), 0], 
                     [0, 0, 1]])
    RotDet = np.dot(RotZ, np.dot(RotY, RotX))
    
    return RotDet


def get_det_tilt(det_normal, detdiru, degree = False):
    """
    Calculate detector tilt angles around x, y, z axis from detector normal and detdiru
    Arguments:
    det_normal   --  unit vector of detector normal
    detdiru      --  unit vector of detector u
    
    Returns:
    tilt_x, tilt_y, tilt_z in radians if degree = False
    """
    tilt_x = (np.arccos(np.dot(detdiru, [0, 0, 1]))) - np.pi/2
    tilt_y = (np.arccos(np.dot(det_normal * [1, 0, 1], [0, 0, 1]))) - np.pi/2
    tilt_z = np.pi/2 - (np.arccos(np.dot(det_normal * [1, 1, 0], [0, 1, 0])))
    if degree:
        return np.rad2deg(tilt_x), np.rad2deg(tilt_y), np.rad2deg(tilt_z)
    else:
        return tilt_x, tilt_y, tilt_z


def build_ImageD11par_from_poni(lattice_par, sgno, poni_file, par_path = 'demo.pars', detector = 'Eiger', verify = True):
    """
    build a ImageD11 par from a poni file
    
    Arguments:
    lattice_par   -- lattice parameters
    sgno          -- space group number
    poni_file     -- poni file
    par_path      -- the output ImageD11 par filename
    detector      -- detector name, 'Eiger' or 'Frelon3'
    verify        -- BOOL for verify the errors
    
    Returns:
    an ImageD11 parameter object
    
    Example:
    lattice_par = [5.43094, 5.43094, 5.43094, 90, 90, 90] # Si cube
    sgno = 227
    poni_file = 'demo.poni'
    build_par_from_poni(lattice_par, sgno, poni_file, par_path = 'demo.pars', detector = 'Eiger')
    """
    poni_results = pyFAI.load(poni_file).getImageD11()
    write_ImageD11pars_from_poni(lattice_par, sgno, poni_results, par_path = par_path, detector = detector)
    trn = ImageD11.transformer.transformer()
    trn.loadfileparameters(par_path)
    if verify:
        err = calc_err(trn.pars, poni_file)
    return trn.pars


def build_poni_from_ImageD11par(par_file, spatial_file = None, poni_path = 'demo.poni', detector = 'Eiger', verify = True):
    
    """
    convert ImageD11 .par to poni:
    Poin1 = z_center*pixel_size [m]
    Poin2 = y_center*pixel_size [m]
    Rot1 = -tilt_z [radian]
    Rot2 = tilt_y [radian]
    Rot3 = tilt_x [radian]
    
    In ImageD11  --- flip story: eiger has no flip, while frelon3 is flipped up and down
    For Eiger:
    'o11': -1,
    'o12': 0,
    'o21': 0,
    'o22': -1,
    For Frelon:
    'o11': 1,
    'o12': 0,
    'o21': 0,
    'o22': -1,
    """    
    
    if spatial_file is None:
        spatial_file = auto_load_spatial_file(detector)
    
    trn = ImageD11.transformer.transformer()
    trn.loadfileparameters(par_file)

    wavelength = trn.pars['wavelength']
    distance = float(trn.pars['distance'])
    y_center = trn.pars['y_center']
    z_center = trn.pars['z_center']
    pixel_sizey = trn.pars['y_size']
    pixel_sizez = trn.pars['z_size']
    tilt_x = trn.pars['tilt_x']
    tilt_y = trn.pars['tilt_y']
    tilt_z = trn.pars['tilt_z']
    poni1 = z_center*pixel_sizez/1e6
    poni2 = y_center*pixel_sizey/1e6
    print('Initial poni1 = {}'.format(poni1))
    print('Initial poni2 = {}'.format(poni2))
    current_time = datetime.now()
    formatted_time = current_time.strftime("%a %b %d %H:%M:%S %Y")
    # write a first poni file
    write_poni(poni_path, formatted_time, spatial_file, distance, poni1, poni2, tilt_x, tilt_y, tilt_z, wavelength)
    
    pyFAI_pars = pyFAI.load(poni_path)
    print(pyFAI_pars)
    poni_results = pyFAI.load(poni_path).getImageD11()
    print(poni_results)
    
    # minimize the errors of distance, poni1 and poni2
    params0 = [trn.pars['distance'], poni1, poni2]
    print('Initial guess: {}'.format(params0))
    result = minimize(error_function, x0 = params0, args=(trn.pars, poni_path, formatted_time, spatial_file, tilt_x, tilt_y, tilt_z, wavelength), method = 'BFGS')
    

    print('Initial guess: {}'.format(params0))
    print('Final fitted results: {}'.format(result.x))
    
    # check the errors again
    err = calc_err(trn.pars, poni_path)
    print('{} poni has been exported.'.format(poni_path))
    
    return err    
    
    
def auto_load_spatial_file(detector):
    if detector in ['Eiger', 'eiger']:
        spatial_file = '/data/id11/nanoscope/Eiger/newSpatial_E-08-0144_20240205.h5'
    elif detector in ['Frelon36', 'frelon36']:
        spatial_file = "/data/id11/3dxrd/inhouse/Frelon36/frelon36_spline_20240604_full.h5"
    elif detector in ['Frelon4M', 'frelon4M', 'frelon4m']:
        spatial_file = "/data/id11/3dxrd/inhouse/Frelon4M/frelon4m.spline"
    else:
        print('{} has not been identified; currently Eiger, Frelon36 and Frelon4M are supported only.'.format(detector))
        spatial_file = None
    return spatial_file


def write_ImageD11pars_from_poni(lattice_par, sgno, poni_results, par_path = 'demo.pars', detector = 'Eiger'):
    """
    write ImageD11 pars from a pyFAI poni file   

    Example:
    lattice_par = [5.43094, 5.43094, 5.43094, 90, 90, 90] # Si cube
    sgno = 227
    poni_results = pyFAI.load('demo.poni').getImageD11()
    write_pars_from_poni(lattice_par, sgno, poni_results, par_path = 'demo.pars', detector = 'Eiger')
    """
    
    if detector in ['Eiger', 'eiger']:
        poni_results["o11"] = -poni_results["o11"]

    par_string = """cell__a {0}
cell__b {1}
cell__c {2}
cell_alpha {3}
cell_beta {4}
cell_gamma {5}
cell_lattice_[P,A,B,C,I,F,R] {6}
chi 0.0
distance {7}
fit_tolerance 0.05
min_bin_prob 1e-05
no_bins 10000
o11 {8}
o12 {9}
o21 {10}
o22 {11}
omegasign 1.0
t_x 0
t_y 0
t_z 0
tilt_x {12}
tilt_y {13}
tilt_z {14}
wavelength {15}
wedge 0.0
weight_hist_intensities 0
y_center {16}
y_size {17}
z_center {18}
z_size {19}""".format(
    lattice_par[0],
    lattice_par[1],
    lattice_par[2],
    lattice_par[3],
    lattice_par[4],
    lattice_par[5],
    sgno,
    poni_results["distance"],
    poni_results["o11"],
    poni_results["o12"],
    poni_results["o21"],
    poni_results["o22"],
    poni_results["tilt_x"],
    poni_results["tilt_y"],
    poni_results["tilt_z"],
    poni_results["wavelength"] * 10,
    poni_results["y_center"] - 0.5,
    poni_results["y_size"],
    poni_results["z_center"] - 0.5,
    poni_results["z_size"])

    # Save it to file
    with open(par_path, "w") as f:
        f.writelines(par_string)
        print('write par file: {}'.format(par_path))
        

def write_poni(poni_path, formatted_time, spatial_file, distance, poni1, poni2, tilt_x, tilt_y, tilt_z, wavelength):
    
    poni_string = """# Nota: C-Order, 1 refers to the Y axis, 2 to the X axis 
# Calibration done at {0}
poni_version: 2
Detector: NexusDetector
Detector_config: {{"filename": "{1}"}}
Distance: {2}
Poni1: {3}
Poni2: {4}
Rot1: {5}
Rot2: {6}
Rot3: {7}
Wavelength: {8}""".format(
    formatted_time,
    spatial_file,
    distance * 1e-6,
    poni1,
    poni2,
    -tilt_z,
    tilt_y,
    tilt_x,
    wavelength * 1e-10)
    # Save it to file
    with open(poni_path, "w") as f:
        f.writelines(poni_string)
        print('write poni file: {}'.format(poni_path))


def error_function(params, ImageD11_pars, poni_path, formatted_time, spatial_file, tilt_x, tilt_y, tilt_z, wavelength):
    """
    calculate errors of distance, poni1 and poni2 between ImageD11_pars and poni
    Arguments:
    params         -- distance, poni1, poni2
    ImageD11_pars  -- ImageD11 parameter object
    others:        -- poni_path, formatted_time, spatial_file, distance, poni1, poni2, tilt_x, tilt_y, tilt_z, wavelength
    
    Returns:
    err  -- sum of squared errors
    
    Example:
    trn = ImageD11.transformer.transformer()
    trn.loadfileparameters(par_path)
    err = calc_err(trn.pars, poni_path) or err = calc_err(trn.pars, poni_path, parameter = 'distance')
    """ 
    
    distance, poni1, poni2 = params
    keys = ['distance', 'y_center']
    write_poni(poni_path, formatted_time, spatial_file, distance, poni1, poni2, tilt_x, tilt_y, tilt_z, wavelength)
    err = calc_err(ImageD11_pars, poni_path)
    return err['distance']**2 + err['y_center']**2 + err['z_center']**2


def calc_err(ImageD11_pars, poni_path, parameter = None):
    """
    calculate difference between ImageD11_pars and poni for each parameter attribute
    Arguments:
    ImageD11_pars -- ImageD11 parameter object
    poni_path      -- pyFAI poni file name
    
    Returns:
    err  -- difference between each parameter attribute
    
    Example:
    trn = ImageD11.transformer.transformer()
    trn.loadfileparameters(par_path)
    err = calc_err(trn.pars, poni_path) or err = calc_err(trn.pars, poni_path, parameter = 'distance')
    """    
    poni_results = pyFAI.load(poni_path).getImageD11()
    if parameter is None:
        err = {}
        for key in poni_results.keys():
            if key == 'wavelength':
                err[key] = poni_results[key]*10 - ImageD11_pars[key]
            elif key == 'y_center' or key == 'z_center':
                err[key] = poni_results[key] - ImageD11_pars[key] - 0.5  # poni has 0.5 pixels bigger than ImageD11, which is due to difference in counting from the side or the center of the first pixel
            else:
                err[key] = poni_results[key] - ImageD11_pars[key]
            print('{} error: {:.6f}'.format(key, err[key]))
    else:
        err = poni_results[parameter] - ImageD11_pars[parameter]
    return err
