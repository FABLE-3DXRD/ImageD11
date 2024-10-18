# a batch of function to read and write files in the following format
# h5 files
# edf images
# tif images
# xml
# write to dream3d and xdmf files
# read motor positions
# Haixing Fang, haixing.fang@esrf.fr
# March 20, 2024

import os
import glob

import numpy as np
import fabio
import matplotlib.pyplot as plt
import h5py
import matplotlib
import ImageD11.columnfile
import ImageD11.unitcell
import ImageD11.transformer
import xml.etree.ElementTree as ET
import scipy.io
from scipy import ndimage
import tifffile as tiff
from PIL import Image




def read_images_from_h5(h5name, scan = '1.1', detector = None, StartIndex = None, EndIndex = None):
    # read images from h5 file
    # StartIndex starts from 0
    # number of images = EndIndex - StartIndex
    if detector is None:
        detector = guess_camera_name(h5name, scan)
    with h5py.File(h5name, 'r') as hin:
        dataset = hin[scan]['measurement/'][detector]
        print(f'data dimensions: {dataset.shape}')
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
        print(f'data dimensions: {dataset.shape}')
        if StartIndex is not None and EndIndex is not None:
            if StartIndex < 0:
                StartIndex = 0
            if EndIndex > dataset.shape[0]:
                EndIndex = dataset.shape[0]
            pct_vol = dataset[StartIndex:EndIndex, :, :]
        else:
            pct_vol = dataset[()]
    if scale[0] != 1.0:
        print(f"I am going to make a scaling on pct_vol because scaling = {scale}.")
        pct_vol = im_scaling(pct_vol, scale)
    
    from scipy import ndimage

    if diffrz_offset is not None:
        print(f"I am going to make a rotation about Z axis for pct_vol because diffrz_offset = {diffrz_offset:.4f} degrees.")
        for i in range(pct_vol.shape[0]):
            pct_vol[i, :, :] = ndimage.rotate(pct_vol[i, :, :], diffrz_offset, reshape=False, mode='nearest', cval=0)
        
    return pct_vol

def im_fabio(data):
    assert isinstance(data, np.ndarray), "Dataset is not a NumPy array."
    if data.ndim < 3:
        obj = fabio.edfimage.EdfImage(data)
        return obj

def im_scaling(im_in, scale = (0.55/1.47,)):
    # example input for scaling: 0.61/1.63 (20x/7.5x) before 14/11/2023 and 0.55/1.47 (20x/7.5x) since 14/11/2023
    if len(scale) != len(im_in.shape):
        scale = [scale[0]] * len(im_in.shape)
        scale = tuple(scale)
        print(f'Scaling the image data with a factor of {scale}')
    im_out = ndimage.zoom(im_in, scale, order=0)
    return im_out
   

def read_motor_posisitions_from_h5(h5name, scan = '1.1'):
    # read all the motor positions from h5 file
    with h5py.File(h5name, 'r') as hin:
        MotorPos = {}
        for motor in hin[scan]['instrument/positioners'].keys():
            MotorPos[motor] = hin[scan]['instrument/positioners'][motor][()]
    return MotorPos


def convert_h52_tif(h5name, save_path, prefix_name = 'proj', save_stack = False, ctrl_name = '1.1/measurement/marana3', data_format = 'uint16'):
    with h5py.File(h5name, 'r') as hin:
        dataset = hin[ctrl_name][:]
        print("Dataset shape:", dataset.shape)
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
    
    print(f"Writing tif images to {save_path}")
    if save_stack:
        # save the entire 3D dataset as a single multi-page TIFF image
        images = [Image.fromarray(normalized_data[i]) for i in range(normalized_data.shape[0])]
        filename = os.path.join(save_path, prefix_name + "_stack.tif")
        images[0].save(filename, save_all=True, append_images=images[1:])
    else:
        # Save each slice of the 3D dataset as an individual TIFF image
        for i in range(normalized_data.shape[0]):
            image = Image.fromarray(normalized_data[i])
            filename = os.path.join(save_path, prefix_name + f"{i:04d}" + ".tif")
            image.save(filename)
    print('Done')
    
    
    print(f"Writing TIFF images to {save_path}")
    if save_stack:
        # Save the entire 3D dataset as a single multi-page TIFF image
        filename = os.path.join(save_path, prefix_name + "_stack.tif")
        tiff.imwrite(filename, normalized_data, photometric='minisblack')
    else:
        # Save each slice of the 3D dataset as an individual TIFF image
        for i in range(normalized_data.shape[0]):
            filename = os.path.join(save_path, f"{prefix_name}_{i:04d}.tif")
            tiff.imwrite(filename, normalized_data[i], photometric='minisblack')



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


def read_matlab_file(filename):
    # Load the .mat file
    mat_data = scipy.io.loadmat(filename)
    return mat_data
    
    
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
        fileID.write(f'   <Topology TopologyType="3DCoRectMesh" Dimensions="{MeshDimensions[0]} {MeshDimensions[1]} {MeshDimensions[2]}"></Topology>\n')
        fileID.write('    <Geometry Type="ORIGIN_DXDYDZ">\n')
        fileID.write('     <!-- Origin  Z, Y, X -->\n')
        fileID.write('     <DataItem Format="XML" Dimensions="3">0 0 0</DataItem>\n')
        fileID.write('     <!-- DxDyDz (Spacing/Resolution) Z, Y, X -->\n')
        fileID.write(f'     <DataItem Format="XML" Dimensions="3">{voxelsize[0]:.6f} {voxelsize[1]:.6f} {voxelsize[2]:.6f}</DataItem>\n')
        fileID.write('    </Geometry>\n')
        
        for attr in attributes:
            fileID.write(f'    <Attribute Name="{attr}" AttributeType="Scalar" Center="Cell">\n')
            fileID.write(f'      <DataItem Format="HDF" Dimensions="{ScalarDimensions[0]} {ScalarDimensions[1]} {ScalarDimensions[2]} {ScalarDimensions[3]}" NumberType="Float" Precision="6" >{h5name}:/{attr}</DataItem>\n')
            fileID.write('    </Attribute>\n')
        
        fileID.write('  </Grid>\n')
        fileID.write('  <!-- *************** END OF GM3D *************** -->\n')
        fileID.write(' </Domain>\n')
        fileID.write('</Xdmf>\n')
    print('Done writing the xdmf file')

    
def camera_names_ID11():
    return ['marana3', 'marana1', 'marana2', 'marana', 'frelon3', 'frelon6', 'eiger', 'basler_eh32']


def guess_camera_name(h5name, scan = '1.1'):
    scan = convert_to_list(scan)
    camera_list = camera_names_ID11()
    camera_name = None
    with h5py.File(h5name, 'r') as hin:
        for key in hin[scan[0]]['measurement'].keys():
            if key in camera_list:
                print(f'I got this camera name: {key}')
                camera_name = key
    return camera_name


def delete_file_if_exists(filepath):
    if os.path.exists(filepath):
        os.remove(filepath)
        print(f'{filepath} has been deleted.')
    else:
        print(f'{filepath} does not exist.')


def delete_group_from_h5(h5name, group_name):
    with h5py.File(h5name, 'a') as h5f:  # Use 'a' mode to open for read/write access
        if group_name in h5f:
            del h5f[group_name]  # Delete the existing group
            return True
        else:
            return False