# io.py, a batch of function to read and write files in the following format
# h5 files
# edf images
# tif images
# xml (to be added)
# fwd_peaks
# bridge with cf
# write to dream3d and xdmf files
# read motor positions
# Haixing Fang, haixing.fang@esrf.fr
# March 20, 2024
# updated on January 30, 2025

import os, shutil
import glob

import numpy as np
import fabio
import matplotlib.pyplot as plt
import h5py
import matplotlib
import ImageD11.columnfile
import ImageD11.unitcell
import ImageD11.transformer
from ImageD11.sinograms.dataset import colfile_from_dict

import xml.etree.ElementTree as ET
import scipy.io
from scipy import ndimage
import tifffile as tiff
from PIL import Image

from tqdm import tqdm
import logging

logging.basicConfig(level=logging.INFO, force=True)


def read_images_from_h5(h5name, scan = '1.1', detector = None, StartIndex = None, EndIndex = None):
    # read images from h5 file
    # StartIndex starts from 0
    # number of images = EndIndex - StartIndex
    if detector is None:
        detector = guess_camera_name(h5name, scan)
    with h5py.File(h5name, 'r') as hin:
        dataset = hin[scan]['measurement/'][detector]
        print('data dimensions: {}'.format(dataset.shape))
        if StartIndex is not None and EndIndex is not None:
            if StartIndex < 0:
                StartIndex = 0
            if EndIndex > dataset.shape[0]:
                EndIndex = dataset.shape[0]
            imgs = dataset[StartIndex:EndIndex, :, :]
        else:
             imgs = dataset[()]
    return imgs
    
    
def read_tomo_from_hdf5(hdf5name, ctrl_name = '/entry0000/reconstruction/results/data', StartIndex = None, EndIndex = None, scale = (1.0,), diffrz_offset = None):
    # read tomography volume from nabu reconstructed hdf5 file
    # usually the file 
    with h5py.File(hdf5name, 'r') as hin:
        dataset = hin[ctrl_name]
        print('data dimensions: {}'.format(dataset.shape))
        if StartIndex is not None and EndIndex is not None:
            if StartIndex < 0:
                StartIndex = 0
            if EndIndex > dataset.shape[0]:
                EndIndex = dataset.shape[0]
            pct_vol = dataset[StartIndex:EndIndex, :, :]
        else:
            pct_vol = dataset[()]
    if scale[0] != 1.0:
        print("I am going to make a scaling on pct_vol because scaling = {}.".format(scale))
        pct_vol = im_scaling(pct_vol, scale)
    
    from scipy import ndimage

    if diffrz_offset is not None:
        print("I am going to make a rotation about Z axis for pct_vol because diffrz_offset = {:.4f} degrees.".format(diffrz_offset))
        for i in range(pct_vol.shape[0]):
            pct_vol[i, :, :] = ndimage.rotate(pct_vol[i, :, :], diffrz_offset, reshape=False, mode='nearest', cval=0)
        
    return pct_vol


def im_fabio(data):
    assert isinstance(data, np.ndarray), "Dataset is not a NumPy array."
    if data.ndim < 3:
        obj = fabio.edfimage.EdfImage(data)
        return obj


def write_edf(data, edfname):
    obj = im_fabio(data)
    if not edfname.endswith('.edf'):
        edfname = edfname + '.edf'
    obj.write(edfname)


def im_scaling(im_in, scale = (0.55/1.47,)):
    # example input for scaling: 0.61/1.63 (20x/7.5x) before 14/11/2023 and 0.55/1.47 (20x/7.5x) since 14/11/2023
    if len(scale) != len(im_in.shape):
        scale = [scale[0]] * len(im_in.shape)
        scale = tuple(scale)
        print('Scaling the image data with a factor of {}'.format(scale))
    im_out = ndimage.zoom(im_in, scale, order=0)
    return im_out
   

def read_motor_posisitions_from_h5(h5name, scan = '1.1'):
    # read all the motor positions from h5 file
    with h5py.File(h5name, 'r') as hin:
        MotorPos = {}
        for motor in hin[scan]['instrument/positioners'].keys():
            MotorPos[motor] = hin[scan]['instrument/positioners'][motor][()]
    return MotorPos


def convert_h52_tif(h5name, save_path, prefix_name = 'proj', save_stack = False, ctrl_name = '1.1/measurement/marana3', data_format = 'uint16', clip_range = None):
    '''
    Read images from an h5/hdf5 file and convert to a specified data format (uint16 or uint8) and save to tif images
    Suitable for convert e.g. nabu_recon hdf5 to tif images
    '''
    print('Loading data from {}'.format(h5name))
    with h5py.File(h5name, 'r') as hin:
        dataset = hin[ctrl_name][:]
        print("Dataset shape:", dataset.shape)
        
    if clip_range is None:
        dataset_mean = np.mean(dataset)
        dataset_std  = np.std(dataset)
        clip_range = [dataset_mean - 3*dataset_std, dataset_mean + 3*dataset_std]
    print("clip_range = {}".format(clip_range))
    dataset = np.clip(dataset, clip_range[0], clip_range[1])  # Clip values to the given range

    if data_format == 'uint16':
        print('Converting to uint16 format ...')
        normalized_data = (65535 * (dataset - np.min(dataset)) / (np.max(dataset) - np.min(dataset))).astype(np.uint16)
    elif data_format == 'uint8':
        print('Converting to uint8 format ...')
        normalized_data = (255 * (dataset - np.min(dataset)) / (np.max(dataset) - np.min(dataset))).astype(np.uint8)
    else:
        print('Converting to uint16 format ...')
        normalized_data = (65535 * (dataset - np.min(dataset)) / (np.max(dataset) - np.min(dataset))).astype(np.uint16)
    
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    
    print("Writing tif images to {}".format(save_path))
    if save_stack:
        # save the entire 3D dataset as a single multi-page TIFF image
        images = [Image.fromarray(normalized_data[i]) for i in range(normalized_data.shape[0])]
        filename = os.path.join(save_path, prefix_name + "_stack.tif")
        images[0].save(filename, save_all=True, append_images=images[1:])
    else:
        # Save each slice of the 3D dataset as an individual TIFF image
        for i in range(normalized_data.shape[0]):
            image = Image.fromarray(normalized_data[i])
            filename = os.path.join(save_path, "{}{:04d}.tif".format(prefix_name, i))
            image.save(filename)
    print('Done')
    return normalized_data


def convert_to_list(input_data):
    if isinstance(input_data, str):
        return [input_data]
    return input_data
    

def get_scan_h5(experiment, scan_type = ['dct',]):
    # get all the samples with a particular scan type
    # scan_type = 'dct', 'pct', 'ff', 'topotomo'
    scan_type = convert_to_list(scan_type)
    subfolders = glob.glob(os.path.join('/data/visitor/',experiment,'id11','*'))
    for subfolder in subfolders:
        sub_folders = glob.glob(os.path.join(subfolder,'*'))
        for sub_folder in sub_folders:
            if 'RAW_DATA' not in os.path.basename(sub_folder):
                continue
            else:
                raw_path = sub_folder
    
    print('...................................... Experiment {} ............................'.format(experiment))
    print('Raw data path: ', raw_path)
    
    # find relevant h5 files
    samples = {}
    skips = 'align', 'sample','beam', 'set_up', 'setup', 'test', 'optics', 'slits', 'CeO2', 'ceO2', 'gpfs', 'check_distance', 'calib', 'temp', 'change_energy', 'energy', 'clean'
    for item in sorted( os.listdir( raw_path ) ):
        folder = os.path.join( raw_path, item )
        if not os.path.isdir( folder ) or folder in skips:
            continue
        datasets = os.listdir( folder )
        for dset in datasets:
            for skip in skips:
                if dset.find(skip)>=0:
                    # print("# Skipping", dset, 'because', skip)
                    break
            else:
                scanfolder = os.path.join( folder, dset )
                # print('scanfolder: ', scanfolder)
                if os.path.isdir(scanfolder):
                    if len(scan_type)==1:
                        scans = glob.glob(scanfolder + '/*' + scan_type[0] + '*.h5')
                    elif len(scan_type)>1:
                        scans = []
                        for scan_type_sub in scan_type:
                            scan = glob.glob(scanfolder + '/*' + scan_type_sub + '*.h5')
                            scans.append(scan)
                    if any(sublist for sublist in scans if sublist):
                        for scan_type_sub in scan_type:
                            if dset.find(scan_type_sub)>=0:
                                 # print('# TODO', item, dset, len(scans))
                                if item in samples:
                                    samples[item].append( dset )
                                else:
                                    samples[item] = [dset,]

    # pprint.pprint(samples)
    print('Found {} samples with scan type of {}'.format(len(samples), scan_type))
    print('Samples of {} are as follows:'.format(scan_type))
    for sample in samples:
        print(sample)
    return samples, raw_path

   
def read_matlab_file(file_path):
    # Load the .mat file
    if is_mat_v7p3(file_path):
        mat_data = read_h5_struct(file_path)
    else:
        mat_data = scipy.io.loadmat(file_path)
    return mat_data
    
    
def is_mat_v7p3(file_path):
    with open(file_path, 'rb') as f:
        header = f.read(19)
        print(header)
        # Check if the file starts with HDF5 magic bytes
        if header == b'MATLAB 7.3 MAT-file':
            return True
        else:
            return False

        
def read_h5_file(file_path):
    '''
    Read all data inside an h5/hdf5 file and return the output as a dictionary
    '''
    data_dict = {}

    def visit_group(name, obj):
        keys = name.split('/')  # Split by group level
        current = data_dict

        # Traverse to the correct location in the dictionary for nested groups
        for key in keys[:-1]:
            current = current.setdefault(key, {})

        # Add dataset or group
        if isinstance(obj, h5py.Dataset):
            current[keys[-1]] = obj[()]
        elif isinstance(obj, h5py.Group):
            current[keys[-1]] = {}

    with h5py.File(file_path, 'r') as file:
        file.visititems(visit_group)
    
    return data_dict


def read_h5_struct(file_path):
    '''
    suitable for reading matlab file saved as a struct, e.g. parameters.mat
    '''
    data_dict = {}

    def load_h5_group(obj, data_dict):
        for key, item in obj.items():
            if isinstance(item, h5py.Group):
                # Initialize a dictionary for this group
                data_dict[key] = {}
                load_h5_group(item, data_dict[key])  # Recursively load sub-groups

            elif isinstance(item, h5py.Dataset):
                # Initialize an entry in the dictionary for datasets
                if item.dtype.kind == 'S':  # ASCII string
                    data_dict[key] = item.asstr()[()]
                elif item.dtype.kind == 'O':  # UTF-8 string or other objects
                    if not item.name == '/parameters/cryst/symm':   # not able to read multiple p.cryst.symm yet
                        data_dict[key] = []
                        #print(item.shape[0], item.name)
                        for i in range(item.shape[0]):
                            entry = item[i][0]
                            referenced_item = obj[entry][()]
                            # print(referenced_item.dtype)
                            if isinstance(entry, h5py.Group):
                                ref_dict = {}
                                load_h5_group(entry, ref_dict)
                                data_dict[key].append(ref_dict)
                            else:
                                if referenced_item.dtype == 'u2':  # Unsigned integer                                   
                                    # data_content = h5py.File(file_path, 'r')[entry][()]
                                    data_content = obj[entry][()]
                                    ascii_codes = data_content.flatten().astype(int)
                                    char_list = [chr(code) for code in ascii_codes]
                                    merged_string = ''.join(char_list)
                                    data_dict[key].append(merged_string)
                                else:
                                    # data_content = h5py.File(file_path, 'r')[entry][()]
                                    data_content = obj[entry][()]
                                    data_dict[key].append(data_content)
                    else:
                        print('NOTE: not able to read multiple parameters.cryst.symm yet!')
                elif item.dtype.kind == 'u':  # Unsigned integer (ASCII codes)
                    data_content = item[:]  # Read all data
                    ascii_codes = data_content.flatten().astype(int)
                    char_list = [chr(code) for code in ascii_codes]
                    merged_string = ''.join(char_list)
                    data_dict[key] = merged_string
                else:
                    data_dict[key] = item[()]

    with h5py.File(file_path, 'r') as file:
        load_h5_group(file, data_dict)
    
    return data_dict


def write_pct_vol(h5name, data, ctrl_name = '/entry0000/reconstruction/results/data'):
    if not ctrl_name.startswith('/'):
        ctrl_name = '/' + ctrl_name
    if not h5name.endswith(('.h5', '.hdf5')):
        h5name = h5name + '.h5'
    with h5py.File(h5name, 'w') as f:
        f.create_dataset(ctrl_name, data = data)
        print('Finished writing data with a dimension of: {}'.format(data.shape))


def write_xdmf(h5name, xdmf_filename = None, ctrl = 'recon_mlem', attributes = ['recon_mlem',], voxelsize = [0.0002, 0.0002, 0.0002]):
    # write .xdmf file for visualing 3D/4D data with ParaView
    if xdmf_filename == None:
        xdmf_filename = h5name.replace('.h5','.xdmf')
    print('Writing {}'.format(xdmf_filename))
    with h5py.File(h5name, 'r') as hin:
        recon_3d = hin[ctrl][()]
    dim = recon_3d.shape
    dim4 = list(dim[:3])
    dim4.append(1)
    
    # MeshDimensions
    MeshDimensions = (dim[0] + 1, dim[1] + 1, dim[2] + 1)
    MeshDimensionsStr = 'Dimensions="%d %d %d"' % (MeshDimensions[0], MeshDimensions[1], MeshDimensions[2])

    # ScalarDimensions
    ScalarDimensions = dim4
    ScalarDimensionsStr = 'Dimensions="%d %d %d %d"' % (ScalarDimensions[0], ScalarDimensions[1], ScalarDimensions[2], ScalarDimensions[3])

    # VectorDimensions
    VectorDimensions = list(dim[:3])
    VectorDimensions.append(3)
    VectorDimensionsStr = 'Dimensions="%d %d %d %d"' % (VectorDimensions[0], VectorDimensions[1], VectorDimensions[2], VectorDimensions[3])
    
    # Write XDMF file
    with open(xdmf_filename, 'wt') as fileID:
        fileID.write('<?xml version="1.0"?>\n')
        fileID.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd"[]>\n')
        fileID.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
        fileID.write(' <Domain>\n')
        fileID.write('  <!-- *************** START OF GM3D *************** -->\n')
        fileID.write('  <Grid Name="GM3D" GridType="Uniform">\n')
        fileID.write('   <Topology TopologyType="3DCoRectMesh" Dimensions="{0} {1} {2}"></Topology>\n'.format(MeshDimensions[0], MeshDimensions[1], MeshDimensions[2]))
        fileID.write('    <Geometry Type="ORIGIN_DXDYDZ">\n')
        fileID.write('     <!-- Origin  Z, Y, X -->\n')
        fileID.write('     <DataItem Format="XML" Dimensions="3">0 0 0</DataItem>\n')
        fileID.write('     <!-- DxDyDz (Spacing/Resolution) Z, Y, X -->\n')
        fileID.write('     <DataItem Format="XML" Dimensions="3">{0:.6f} {1:.6f} {2:.6f}</DataItem>\n'.format(voxelsize[0], voxelsize[1], voxelsize[2]))
        fileID.write('    </Geometry>\n')

        for attr in attributes:
            fileID.write('    <Attribute Name="{0}" AttributeType="Scalar" Center="Cell">\n'.format(attr))
            fileID.write('      <DataItem Format="HDF" Dimensions="{0} {1} {2} {3}" NumberType="Float" Precision="6" >{4}:/{5}</DataItem>\n'.format(ScalarDimensions[0], ScalarDimensions[1], ScalarDimensions[2], ScalarDimensions[3], h5name, attr))
            fileID.write('    </Attribute>\n')

        fileID.write('  </Grid>\n')
        fileID.write('  <!-- *************** END OF GM3D *************** -->\n')
        fileID.write(' </Domain>\n')
        fileID.write('</Xdmf>\n')

    print('Done writing the xdmf file')

    
def camera_names_ID11():
    return ['marana3', 'marana1', 'marana2', 'marana', 'frelon1', 'frelon3', 'frelon6', 'eiger', 'basler_eh32', 'frelon16']


def guess_camera_name(h5name, scan = '1.1'):
    scan = convert_to_list(scan)
    camera_list = camera_names_ID11()
    camera_name = None
    with h5py.File(h5name, 'r') as hin:
        for key in hin[scan[0]]['measurement'].keys():
            if key in camera_list:
                print('I got this camera name: {}'.format(key))
                camera_name = key
    return camera_name


def delete_file_if_exists(filepath):
    if os.path.exists(filepath):
        os.remove(filepath)
        print('{} has been deleted.'.format(filepath))
    else:
        print('{} does not exist.'.format(filepath))


def delete_group_from_h5(h5name, group_name):
    with h5py.File(h5name, 'a') as h5f:  # Use 'a' mode to open for read/write access
        if group_name in h5f:
            del h5f[group_name]  # Delete the existing group
            return True
        else:
            return False
        
        
def print_all_keys(d, prefix=""):
    """
    Recursively prints all keys in a dictionary, including sub-keys.

    Parameters:
    d (dict): The dictionary to explore.
    prefix (str): The prefix for nested keys (used internally for recursion).
    """
    if not isinstance(d, dict):
        raise TypeError("Expected a dictionary, but got {} instead.".format(type(d)))
    print('Got the following keys:')
    for key, value in d.items():
        full_key = "{}.{}".format(prefix, key) if prefix else key
        print(full_key)
        if isinstance(value, dict):  # If the value is a dictionary, recurse
            print_all_keys(value, prefix=full_key)


def write_fwd_peaks(fwd_peaks, output_folder = None, fname_prefix = None, verbose = 1):
    """
    Write fwd_peaks to an h5 file with each column as a separate dataset, see how to compute fwd_peaks in forward_projector.py.

    Args:
        fwd_peaks (array): an N*25 array containing forward computed peaks
        output_folder (string): output folder
        fname_prefix (string): filename prefix, "fpks_" by default
        verbose: logging level, 0, 1 or 2
        
    Returns:
        None
    """
    assert fwd_peaks.shape[1] == 25, "fwd_peaks must have 25 columns"
    if output_folder is None:
        output_folder = os.getcwd()
    if fname_prefix is None:
        fname_prefix = 'fpks'
        
    dty = fwd_peaks[0, 8]
    outname = os.path.join(output_folder, fname_prefix + '_dty_' + str(round(dty, 2)).replace('.', 'p') + '.h5')  # e.g. dty = 0.1, outname = 'fpks_dty_0p10.h5'
    
    col_names = ['grainID', 'voxel_zyx', 'weight', 'pos_xyz', 'dty', 'omega', 'tth', 'eta', 'hkl',
                 'g_xyz', 'det_pixel_yz', 'lorentz', 'polarization', 'transmission', 'sum_intensity', 'ds']
    
    with h5py.File(outname, 'w') as hout:
        hout.create_dataset(col_names[0], data = fwd_peaks[:, 0], dtype='int32')
        hout.create_dataset(col_names[1], data = fwd_peaks[:, [1, 2, 3]], dtype='int32')
        hout.create_dataset(col_names[2], data = fwd_peaks[:, 4], dtype='float')
        hout.create_dataset(col_names[3], data = fwd_peaks[:, [5, 6, 7]], dtype='float')
        hout.create_dataset(col_names[4], data = fwd_peaks[:, 8], dtype='float')
        hout.create_dataset(col_names[5], data = fwd_peaks[:, 9], dtype='float')
        hout.create_dataset(col_names[6], data = fwd_peaks[:, 10], dtype='float')
        hout.create_dataset(col_names[7], data = fwd_peaks[:, 11], dtype='float')
        hout.create_dataset(col_names[8], data = fwd_peaks[:, [12, 13, 14]], dtype='int32')
        hout.create_dataset(col_names[9], data = fwd_peaks[:, [15, 16, 17]], dtype='float')
        hout.create_dataset(col_names[10], data = fwd_peaks[:, [18, 19]], dtype='float')
        hout.create_dataset(col_names[11], data = fwd_peaks[:, 20], dtype='float')
        hout.create_dataset(col_names[12], data = fwd_peaks[:, 21], dtype='float')
        hout.create_dataset(col_names[13], data = fwd_peaks[:, 22], dtype='float')
        hout.create_dataset(col_names[14], data = fwd_peaks[:, 23], dtype='float')
        hout.create_dataset(col_names[15], data = fwd_peaks[:, 24], dtype='float')
    if verbose >= 1:
        logging.info('fwd_peaks written to {}'.format(outname))

    return None


def read_fwd_peaks(h5_file, verbose=1):
    """
    Read fwd_peaks data from an h5 file and reconstruct it as an N*24 array.
    
    Args:
        h5_file (string): path to the HDF5 file containing fwd_peaks data
        verbose: logging level, 0, 1 or 2
        
    Returns:
        fwd_peaks: an N*25 array containing forward computed peaks
    """
    
    col_names = ['grainID', 'voxel_zyx', 'weight', 'pos_xyz', 'dty', 'omega', 'tth', 'eta', 'hkl',
                 'g_xyz', 'det_pixel_yz', 'lorentz', 'polarization', 'transmission', 'sum_intensity', 'ds']
    if not os.path.exists(h5_file):
        return None
    try:
        with h5py.File(h5_file, 'r') as hin:
            num_rows = hin[col_names[0]].shape[0]  # Determine the number of rows
            fwd_peaks = np.zeros((num_rows, 25), dtype=np.float64)

            # Read each dataset and fill the corresponding columns in fwd_peaks
            fwd_peaks[:, 0] = hin[col_names[0]][:]    # grainID
            fwd_peaks[:, 1:4] = hin[col_names[1]][:]  # voxel_zyx indices
            fwd_peaks[:, 4] = hin[col_names[2]][:]    # weight
            fwd_peaks[:, 5:8] = hin[col_names[3]][:]  # pos_xyz [mm]
            fwd_peaks[:, 8] = hin[col_names[4]][:]    # dty [um]
            fwd_peaks[:, 9] = hin[col_names[5]][:]    # omega [deg]
            fwd_peaks[:, 10] = hin[col_names[6]][:]   # tth [deg]
            fwd_peaks[:, 11] = hin[col_names[7]][:]   # eta [deg]
            fwd_peaks[:, 12:15] = hin[col_names[8]][:]  # hkl
            fwd_peaks[:, 15:18] = hin[col_names[9]][:]  # g_xyz
            fwd_peaks[:, 18:20] = hin[col_names[10]][:] # det_pixel_yz [pixel]
            fwd_peaks[:, 20] = hin[col_names[11]][:]  # Lorentz
            fwd_peaks[:, 21] = hin[col_names[12]][:]  # Polarization
            fwd_peaks[:, 22] = hin[col_names[13]][:]  # transmission
            fwd_peaks[:, 23] = hin[col_names[14]][:]  # sum_intensity
            fwd_peaks[:, 24] = hin[col_names[15]][:]  # ds
        if verbose >= 1:
            logging.info('fwd_peaks read from {}'.format(h5_file))
        return fwd_peaks
    except Exception as e:
        print("An unexpected error occurred: {}".format(e))
        return None


def convert_fwd_peaks_to_cf(fwd_peaks):
    """
    convert fwd_peaks to ImageD11 column file object.
    
    Args:
        fwd_peaks: an N*25 array containing forward computed peaks        
    Returns:
        cf: an ImageD11 cf object, see ImageD11.sinograms.dataset.colfile_from_dict
    """
    col_names = ['grainID', 'voxel_z', 'voxel_y', 'voxel_x', 'weight', 'pos_x', 'pos_y', 'pos_z', 'dty', 'omega', 'tth', 'eta', 'h', 'k', 'l',
             'gx', 'gy', 'gz', 'det_pixel_y', 'det_pixel_z', 'lorentz', 'polarization', 'transmission', 'sum_intensity', 'ds']
    if isinstance(fwd_peaks, dict):
        raise TypeError("Expected a numpy array, but got a dict")
    if not isinstance(fwd_peaks, np.ndarray):
        fwd_peaks = np.vstack(fwd_peaks)
    assert fwd_peaks.shape[1] == 25, "fwd_peaks must have 25 columns"
    fwd_peaks_dict = {}
    for i, col_name in enumerate(col_names):
        fwd_peaks_dict[col_name] = fwd_peaks[:,i]
        
    cf = colfile_from_dict(fwd_peaks_dict)  # convert to a cf object
    return cf


def convert_cf_to_fwd_peaks(cf):
    """
    convert cf (generated by forward_projector) to fwd_peaks numpy array.
    
    Args:
        cf: ImageD11 column file object (generated by forward_projector), e.g. fp.cf_2d or fp.cf_3d
    Returns:
        fwd_peaks: an N*25 array containing forward computed peaks
    """
    num_rows = cf.nrows  # Determine the number of rows
    fwd_peaks = np.zeros((num_rows, 25), dtype=np.float64)

    # Read each dataset and fill the corresponding columns in fwd_peaks
    fwd_peaks[:, 0] = cf.grainID  # grainID
    fwd_peaks[:, 1] = cf.voxel_z  # voxel_zyx indices
    fwd_peaks[:, 2] = cf.voxel_y  # voxel_zyx indices
    fwd_peaks[:, 3] = cf.voxel_x  # voxel_zyx indices
    fwd_peaks[:, 4] = cf.weight   # weight
    fwd_peaks[:, 5] = cf.pos_x    # pos_xyz [mm]
    fwd_peaks[:, 6] = cf.pos_y    # pos_xyz [mm]
    fwd_peaks[:, 7] = cf.pos_z    # pos_xyz [mm]
    fwd_peaks[:, 8] = cf.dty      # dty [um]
    fwd_peaks[:, 9] = cf.omega    # omega [deg]
    fwd_peaks[:, 10] = cf.tth     # tth [deg]
    fwd_peaks[:, 11] = cf.eta     # eta [deg]
    fwd_peaks[:, 12] = cf.h       # hkl
    fwd_peaks[:, 13] = cf.k       # hkl
    fwd_peaks[:, 14] = cf.l       # hkl
    fwd_peaks[:, 15] = cf.gx      # g_xyz
    fwd_peaks[:, 16] = cf.gy      # g_xyz
    fwd_peaks[:, 17] = cf.gz      # g_xyz
    fwd_peaks[:, 18] = cf.det_pixel_y  # det_pixel_yz [pixel]
    fwd_peaks[:, 19] = cf.det_pixel_z  # det_pixel_yz [pixel]
    fwd_peaks[:, 20] = cf.lorentz      # Lorentz
    fwd_peaks[:, 21] = cf.polarization # Polarization
    fwd_peaks[:, 22] = cf.transmission # transmission
    fwd_peaks[:, 23] = cf.sum_intensity# sum_intensity
    fwd_peaks[:, 24] = cf.ds           # ds
    return fwd_peaks


def read_fsparse(h5_file, group_name = "/entry_0000/ESRF-ID11/eiger/data"):
    """
    Read fsparse_pks from a sparse h5 peaks file.
    
    Args:
        h5_file (string): path to the HDF5 file containing the sparse peaks generated from forward_projector.make_projs_and_sparse_file
        group_name: group name
        
    Returns:
        fsparse_pks: a dictionary containing 'col', 'row', 'nnz', 'intensity', 'dty', 'rot'
    """
    names = ['col', 'row', 'nnz', 'intensity', 'dty', 'rot']
    fsparse_pks = {}
    with h5py.File(h5_file, 'r') as hin:
        g = hin.get(group_name)
        if g is None:
            raise ValueError("Group {} not found in file {}".format(group_name, h5_file))
        for name in names:
            if name in g.keys():
                fsparse_pks[name] = g[name][()]
    return fsparse_pks


def copytree_with_progress(src, dst, skip_keys=None, create_subfolder = True, overwrite = False):
    """
    Copies a folder with all subfolders and files, displaying progress. 
    Skips subfolders containing specified keys.

    Main usage is to backup files from the server to local computer by skipping filefolders containing keywords specified by the input "skip_keys"
    e.g. skip_keys = ['0_rawdata', '2_difblob', '2_difspot', '6_rendering', '7_fed', '8_analysis', 'OAR_log', 'RAW_DATA']
    
    Args:
        src (str): Source directory path.
        dst (str): Destination directory path.
        skip_keys (list): List of strings; skip folders containing these keys in their path.
        create_subfolder (logic): using the last name of the source directory to create a new sub directory under the destination direction
        overwrite (logic): flag for overwriting any existing files in the destination directory, False for skipping
    """
    if not os.path.exists(src):
        raise ValueError("Source folder {} does not exist.".format(src))

    if create_subfolder:
        dst = os.path.join(dst, os.path.split(src)[1])
        
    print("Source directory: {}".format(src))
    print("Destination directory: {}".format(dst))
    print("Skipping {}".format(skip_keys))

    if skip_keys is None:
        skip_keys = []
    
    # Get all files and directories in the source folder
    files_to_copy = []
    for root, dirs, files in os.walk(src):
        # Skip subfolders containing any of the keys
        if any(key in root for key in skip_keys):
            continue
        for file in files:
            src_file = os.path.join(root, file)
            dst_file = os.path.join(dst, os.path.relpath(src_file, src))
            # Add file only if it doesn't already exist or the source is newer
            if overwrite:
                files_to_copy.append((src_file, dst_file))
            else:
                if not os.path.exists(dst_file) or os.path.getmtime(src_file) > os.path.getmtime(dst_file):
                    files_to_copy.append((src_file, dst_file))
                    
    # Create destination directory if it does not exist
    os.makedirs(dst, exist_ok=True)
    
    # Use tqdm for progress visualization
    with tqdm(total=len(files_to_copy), desc="Copying files", unit="file") as pbar:
        for src_file, dst_file in files_to_copy:
            os.makedirs(os.path.dirname(dst_file), exist_ok=True)
            shutil.copy2(src_file, dst_file)
            pbar.update(1)    
    print("Folder '{}' has been successfully copied to '{}', skipping subfolders containing {}.".format(src, dst, skip_keys))
    return None