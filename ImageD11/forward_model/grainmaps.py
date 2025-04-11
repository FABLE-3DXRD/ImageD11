# grainmaps.py
# Main functions are as follows:
# 1) grainmap class
#    - read the tensor_map from tomo or pbp reconstructed results
#    - build a DS-type grain map from the ImageD11 results
#    - merge regions to have grain IDs, i.e. 'labels', if there is no effective 'labels' in the tensor_map
#    - extract grain boundaries and compute grain boundary character from the DS-type grain map
#    - TODO: operations on the grain map, including cropping, rotating, re-scaling, comparing with another grain map and overlaying on tomographic slice etc
# 2) indexing_iterative: iterative indexing with left peaks
# 3) merge_grains: identify grains by merging voxels with similar UBIs, 'similar' is defined with a user defined tolerance in misorientation and/or unit cells
# 4) pairing_grains: track the same grains between two different grains object, useful for tracking grains across different datasets
# 5) Crystal, symmetry related functions, e.g. calculate misorientations considering the crystal symmetry, compute mean Rodrigues vector etc.
# Haixing Fang, haixing.fang@esrf.fr
# Dec 9, 2024
# updated on January 30, 2025

import os
import enum
import functools
from math import sin, cos, sqrt, gcd
import numpy as np
from numpy import pi, dot, transpose, radians
import matplotlib
from matplotlib import pyplot as plt

from joblib import Parallel, delayed
import time
from tqdm import tqdm
import re
from sklearn.cluster import KMeans
import warnings
from scipy.ndimage import distance_transform_edt

import h5py
import ImageD11.nbGui.nb_utils as utils
import ImageD11.grain
import ImageD11.indexing
import ImageD11.columnfile
from ImageD11.unitcell import Phases
from ImageD11.unitcell import unitcell
from ImageD11.peakselect import select_ring_peaks_by_intensity

from ImageD11.forward_model import forward_model
from ImageD11.forward_model import io
from ImageD11.forward_model import ori_converter


class grainmap:
    '''
    grain map class to create and do operations on a grain map, akin fwd-DCT DS
    example:
    filename = ds.grains_file
    gm = grainmap(filename)
    gm.merge_and_identify_grains(FirstGrainID = 0, dis_tol = np.sqrt(2))
    DS = gm.DS
    '''
    def __init__(self, filename = None, outname = None, min_misori = 3, crystal_system = 'cubic', remove_small_grains = True, min_vol = 2):
        self.filename = filename
        self.outname = outname
        if self.outname is None and self.filename is not None:
            self.outname = os.path.join(os.path.split(self.filename)[0], 'DS.h5')
        if filename is not None:
            self.DS, self.tensor_map = self.convert2DS()
        else:
            self.DS = None
            self.tensor_map = None
        self.min_misori = min_misori
        self.crystal_system = crystal_system
        self.remove_small_grains = remove_small_grains
        self.min_vol = min_vol
        print('****************************************** Parameters for operating the grain map: ')
        print('Output file name: {}'.format(self.outname))
        print('min_misori = {}'.format(min_misori))
        print('crystal_system: {}'.format(crystal_system))
        print('remove_small_grains = {}'.format(remove_small_grains))
        print('min_vol = {} voxels'.format(min_vol))
        
    
    def read_tensor_map(self):
        if os.path.exists(self.filename):
            tensor_map = io.read_h5_file(self.filename)
            if isinstance(tensor_map, dict):
                io.print_all_keys(tensor_map)
            else:
                print("The file did not produce a dictionary.")
            return tensor_map
        else:
            print('{} is not found'.format(self.filename))
            return None
    
    
    def convert2DS(self):
        """
        convert tensor_map to a dictionary like DS, which is similar to fwd-DCT struct, easier for subsequent operation
        """
        tensor_map = self.read_tensor_map()
        
        for key in tensor_map.keys():
            if 'TensorMap' in key:
                print('Got the key name {}'.format(key))
                keyname = key
        
        print('********************* Converting to DS format *****************************')
        keys_list = ['B', 'U', 'UB', 'UBI', 'eps_sample', 'euler', 'ipf_x', 'ipf_y', 'ipf_z', 'mt', 'nuniq', 'phase_ids', 'unitcell', 'intensity', 'labels']
        DS = {}
        for k in keys_list:
            if k in tensor_map[keyname]['maps'].keys():
                print('Loading {} with a shape of {} to DS ...'.format(k, tensor_map[keyname]['maps'][k].shape))
                DS[k] = tensor_map[keyname]['maps'][k]
        if 'step' in tensor_map[keyname].keys():
            DS['voxel_size'] = tensor_map[keyname]['step'] # [um/pixel]
        
        
        DS['Rod'] = np.empty((DS['UBI'].shape[0], DS['UBI'].shape[1], DS['UBI'].shape[2], 3))
        DS['Rod'].fill(np.nan)
                
        # assign labels ID without any merging, i.e. grain ID
        if 'labels' not in DS.keys():
            DS['labels'] = np.zeros_like(DS['phase_ids']) - 1
            indices = np.argwhere(DS['phase_ids'] > -1)

            for label, (i, j, k) in enumerate(indices):
                DS['labels'][i, j, k] = label  # label starts from 0
            
            print('Found a map (dimensions {}) with {} voxels to be labeled with grain ID'.format(DS['labels'].shape, np.max(DS['labels'])))
            print('Please note: current labels are only randomized without any region merging !!!')

        # get U from UBI if no 'U' existed
        if 'U' not in DS.keys() and 'UBI' in DS.keys():
            try:
                latticepar = np.array(tensor_map[keyname]['phases']['0'].decode('utf-8').split(), dtype = float)
                ucell = unitcell(latticepar[0:6], int(latticepar[-1]))
                B = ucell.B
                print('Got lattice parameters: {} for calculating B matrix to decompose UBI and derive U matrix'.format(latticepar))
                
                DS['U'] = np.empty((DS['UBI'].shape[0], DS['UBI'].shape[1], DS['UBI'].shape[2], 3, 3))
                DS['U'].fill(np.nan)
                
                indices = np.argwhere(~np.isnan(DS['UBI'][:, :, :, 0, 0]))
                for (i, j, k) in indices:
                    DS['U'][i, j, k, :, :] = np.dot(B, DS['UBI'][i, j, k, :, :]).T
                print('Done with deriving U matrix !')
            except KeyError as e:
                print("Key error: {}".format(e))
            except Exception as e:
                print("An unexpected error occurred: {}".format(e))
        
        # assign Rodrigues vector
        try:
            for i in range(DS['U'].shape[0]):
                for j in range(DS['U'].shape[1]):
                    for k in range(DS['U'].shape[2]):
                        if not np.isnan(DS['U'][i, j, k, 0, 0]):
                            DS['Rod'][i, j, k] = ori_converter.quat2rod(ori_converter.u2quat(DS['U'][i, j, k, :, :]))
            print('Got Rodrigues vector for each voxel !')
        except KeyError as e:
            print("Key error: {}".format(e))
        except Exception as e:
            print("An unexpected error occurred: {}".format(e))

        return DS, tensor_map
    
    
    def merge_and_identify_grains(self, FirstGrainID = 0, dis_tol = np.sqrt(2), count_max = 100):
        """
        merge regions and identify grains to update grain IDs, i.e. 'labels' in DS
        default option to also remove small grains, e.g. remove grain IDs with size no bigger than 2 voxels
        """
        self.DS = DS_merge_and_identify_grains(self.DS, FirstGrainID = FirstGrainID, min_misori = self.min_misori, dis_tol = dis_tol, crystal_system = self.crystal_system, count_max = count_max)
        if self.remove_small_grains:
            self.DS = DS_remove_small_grains(self.DS, self.min_vol)
        
    
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
        
        
    def write_labels2tmfile(self):
        """
        write labels as new dataset to the original tensor map file'
        """
        with h5py.File(self.filename, 'a') as hout:
            if 'labels' in self.DS.keys():
                hout.create_dataset('labels', data = self.DS['labels'])
                print('Done with appending labels to {}'.format(self.filename))
            else:
                print('No labels in DS keys; Do merge_and_identify_grains first ...')
    
    
    'TODO'
    """
    method to create grain boundaries and compute boundary misorientations    
    method for cropping
    method for rotating
    method for rescaling
    method for plotting
    method for comparing with another DS
    
    """


def DS_merge_and_identify_grains(DS, FirstGrainID = 0, min_misori = 3.0, dis_tol = np.sqrt(2), crystal_system = 'cubic', count_max = 100):

    """
    Merge regions within which the voxel misorientation is smaller than a pre-defined values
    New grain ID will be assigned to DS dictionary 'labels' key
    Optional to keep the first grain IDs static by defining FirstGrainID (by default = 0 means no first grains to be static)

    Arguments:
    DS                -- grain map dictionary, gm = grainmap(tensor_map_file); DS = gm.DS
    FirstGrainID      -- First grain ID to be considered for merging, by default = 0
    min_misori        -- misorientation for merging regions
    dis_tol           -- maximum distance around the target voxel for identifying its neigboring grains, np.sqrt(2) corresponds to the diagonal voxel
    crystal_system    -- crystal system name one of ['cubic', 'hexagonal', 'orthorhombic', 'tetragonal', 'trigonal', 'monoclinic', 'triclinic']
    count_max         -- maximum iterations for merging regions

    Returns:
    DS_new          -- DS grain map dictionary with updated grain IDs
    """
    
    stop_merge = False
    count = 0
    start_time = time.time()
    while not stop_merge:
        DS_new = DS_merge_and_identify_grains_sub(DS, FirstGrainID = FirstGrainID, min_misori = min_misori, dis_tol = dis_tol, crystal_system = crystal_system)

        N_old = np.max(DS['labels'])
        N_new = np.max(DS_new['labels'])

        count = count+1
        print('*************************** Done iterative No. {} for merging ***************************'.format(count))

        if N_new == N_old or count >= count_max:
            stop_merge = True
        else:
            DS = DS_new.copy()

    end_time = time.time()
    print("The whole merging took {:.4f} seconds".format(end_time - start_time))

    return DS_new


def DS_remove_small_grains(DS, min_vol = 2):
    """
    Remove small grains defined by no bigger than min_vol voxels
    
    Arguments:
    DS                -- grain map dictionary, gm = grainmap(tensor_map_file); DS = gm.DS
    min_vol           -- minimum number of voxels to be kept in DS [voxel]

    Returns:
    DS_out            -- DS grain map dictionary with updated grain IDs
    """
    
    assert 'labels' in DS.keys(), "labels must be in DS keys"
    DS_out = DS.copy()
    count1 = 0
    count2 = 0
    for j in range(np.max(DS['labels'])+1):
        ind = np.array(np.where(DS['labels'] == j)).T
        Vregion = ind.shape[0] # number of voxels for each grain
        if Vregion <= min_vol:
            count1 += 1
            DS_out['labels'][ind[:,0], ind[:,1], ind[:,2]] = -1
        elif Vregion > min_vol:
            count2 += 1
    print('Found and removed {} small grains with size <= {} voxels; {} grains left now.'.format(count1, min_vol, count2))
    return DS_out


def DS_merge_and_identify_grains_sub(DS, FirstGrainID = 0, min_misori = 3.0, dis_tol = np.sqrt(2), crystal_system = 'cubic'):
    
    """
    sub function:
    Merge regions within which the voxel misorientation is smaller than a pre-defined values
    New grain ID will be assigned to DS dictionary 'labels' key
    Optional to keep the first grain IDs static by defining FirstGrainID (by default = 0 means no first grains to be static)

    Arguments:
    DS                -- grain map dictionary, gm = grainmap(tensor_map_file); DS = gm.DS
    FirstGrainID      -- First grain ID to be considered for merging, by default = 0
    min_misori        -- misorientation for merging regions
    dis_tol           -- maximum distance around the target voxel for identifying its neigboring grains, np.sqrt(2) corresponds to distance to the diagonal voxel in 2D
    crystal_system    -- crystal system name one of ['cubic', 'hexagonal', 'orthorhombic', 'tetragonal', 'trigonal', 'monoclinic', 'triclinic']

    Returns:
    DS_merge          -- DS grain map dictionary with updated grain IDs
    """
    
    if crystal_system in ['cubic', 'hexagonal', 'orthorhombic', 'tetragonal', 'trigonal', 'monoclinic', 'triclinic']:
        print('{} is OK'.format(crystal_system))
        crystal_structure = Symmetry[crystal_system]
    else:
        raise ValueError('{} is not supported.'.format(crystal_system))
    
    DS_merge = DS.copy()
    DS_merge['labels'] = np.zeros_like(DS['labels']) - 1

    id = np.unique(DS['labels'].ravel())
    id = id[id >= FirstGrainID]  # Filter IDs greater than or equal to FirstGrainID
    id = np.sort(id)
    id0 = len(id)
    print('Initial number of grain IDs is {}'.format(id0))

    # Initialize Vregion array to store the number of voxels for each region
    Vregion=np.zeros((np.max(DS['labels'])+1, 1), dtype = 'int64')
    print('Get the number of voxels for each of the {} regions ...'.format(len(id)) )       
    for j in range(np.max(DS['labels'])+1):
        Vregion[j] = np.array(np.where(DS['labels'] == j)).T.shape[0] # number of voxels for each region
        if j % 20000 == 0 or j == np.max(DS['labels'])-1:
            print('Done for {} regions ...'.format(j))
            
    # Retrieve GrainIDs for the unchanged ones
    if FirstGrainID > 0:
        print('Retrieve the first {} grain IDs, which are supposed to stay static ...'.format(FirstGrainID - 1))
        for i in range(FirstGrainID):
            ind = np.array(np.where(DS['labels'] == i)).T
            if ind.size > 0:
                DS_merge['labels'][ind[:, 0], ind[:, 1], ind[:, 2]] = i

    # merging loop
    ct = 0
    stop_merge = False
    i = -1
    cutbox = np.zeros((3, 2), dtype=int)
    print('{} regions to be merged ...'.format(id0))
    start_time = time.time()
    while not stop_merge:
        i += 1
        id_vol = DS['labels'] == id[0]
        id_bbox = np.array(np.where(id_vol > 0)).T

        cutbox[:, 0] = np.maximum(np.min(id_bbox, axis=0) - 2, 0)
        cutbox[:, 1] = np.minimum(np.max(id_bbox, axis=0) + 2, np.array(DS['labels'].shape) - 1)

        slices = tuple(slice(cutbox[ii, 0], cutbox[ii, 1] + 1) for ii in range(3))
        id_vol_cut = id_vol[slices]
        vol_cut = DS['labels'][slices]

        # Compute distance map and neighbors
        id_dismap = distance_transform_edt(~id_vol_cut)  # Equivalent to bwdist in MATLAB
        id_neigb = (id_dismap > 0) & (id_dismap < dis_tol)
        id_neigb = np.unique(vol_cut[id_neigb])
        id_neigb = id_neigb[id_neigb > -1]  # Filter out -1 values

        id_indices = np.array(np.where(id_vol)).T

        if len(id_indices) > 1:
            # Multiple indices case
            id_indices_ind = np.ravel_multi_index(
                (id_indices[:, 0], id_indices[:, 1], id_indices[:, 2]), id_vol.shape
            )
            r0 = get_mean_rod(DS['Rod'][id_indices[:, 0], id_indices[:, 1], id_indices[:, 2], :])
        else:
            r0 = np.ravel(DS['Rod'][id_indices[:, 0], id_indices[:, 1], id_indices[:, 2]])

        # get the information for the target region
        U0 = ori_converter.quat2u(ori_converter.rod2quat(r0))
        euler_angles0 = ori_converter.u2euler(U0)
        V0 = Vregion[ id[0] ]
        mergeInd = id_indices    # indices for merging regions    

        # find the potential neighboring region to be merged by checking the misorientation
        indices_neigb_tomerge = []    
        if id_neigb.size > 0:
            r1_list = []
            ang_list = []
            indices_neigb_tomerge = []
            V1_list = []

            for id_neigb_j in id_neigb:
                neigb_vol = DS['labels'] == id_neigb_j
                indices_neigb = np.array(np.where(neigb_vol)).T

                r1 = get_mean_rod(DS['Rod'][indices_neigb[:, 0], indices_neigb[:, 1], indices_neigb[:, 2], :], auto_check = False)
                U1 = ori_converter.quat2u(ori_converter.rod2quat(r1))
                the_angle, _, _ = disorientation(U0, U1, crystal_structure=crystal_structure)

                ang_list.append(np.rad2deg(the_angle))
                V1_list.append(Vregion[id_neigb_j])
                r1_list.append(r1)

                if np.rad2deg(the_angle) <= min_misori:
                    indices_neigb_tomerge.append(indices_neigb)

            # Batch computation for neighbors
            ang = np.array(ang_list)
            V1 = np.array(V1_list)
            r1 = np.array(r1_list)

            ID_for_merge_indices = np.where(ang <= min_misori)[0]
            ID_for_merge = id_neigb[ID_for_merge_indices]

            newGrainID = i + FirstGrainID

            if ID_for_merge.size != 0:
                # Inherit_region_nr = np.hstack([id[0], ID_for_merge])
                weights = np.vstack([V0, V1[ID_for_merge_indices]])/np.sum( np.vstack([V0, V1[ID_for_merge_indices]]) ) # weights for updating the orientation and completeness

                r_new = np.dot(np.vstack([r0, r1[ID_for_merge_indices,:].reshape(len(ID_for_merge_indices), 3)]).T, weights).T # new Rodrigues vector
                U_new = ori_converter.quat2u(ori_converter.rod2quat(np.ravel(r_new)))
                EulerZXZ_new = ori_converter.u2euler(U_new)
                
                indices_neigb_tomerge = np.vstack(indices_neigb_tomerge)
                mergeInd = np.vstack([mergeInd, indices_neigb_tomerge])  # all indices to be merged together
                if np.min(ID_for_merge) < FirstGrainID:
                    # print('Detect smaller grain ID:')
                    # print(newGrainID, ID_for_merge)
                    newGrainID = np.min(ID_for_merge)
                    ct += 1
                    # print(ct)
            else:
                # Inherit_region_nr = id[0]
                EulerZXZ_new = euler_angles0
                r_new = r0

            # update DS_merge
            DS_merge['labels'][mergeInd[:,0], mergeInd[:,1], mergeInd[:,2]] = newGrainID

            # remove the id which have been processed
            id = np.setdiff1d(id, np.hstack([id[0], np.ravel(ID_for_merge)]))
        else:
            Inherit_region_nr = id[0]
            U_new = U0
            DS_merge['labels'][mergeInd[:,0], mergeInd[:,1], mergeInd[:,2]] = i + FirstGrainID
            id = np.setdiff1d(id, id[0])

        # Display progress
        if i > 1 and i % 10000 == 0:
            print('{} grains identified. {} regions have been merged and {} regions waiting for merging ...'.format(i+1, id0 - len(id), len(id)))

        stop_merge = id.size == 0

    end_time = time.time()
    print("The whole merging took {:.4f} seconds".format(end_time - start_time))
    print('{} grains identified out of {} regions.'.format(i+1, id0))
    
    return DS_merge


def DS_to_paraview(DS, h5name = 'DS.h5'):
    """
    Write an .xdmf file that lets you read the DS data with ParaView

    Arguments:
    DS                -- a dictionary contains DS-like map from the output of ImageD11.forward_model.grainmap class
    h5name            -- the corresponding h5 filename, if not exist, I will creat one
    """
    assert 'labels' in DS.keys() and 'voxel_size' in DS.keys(), 'DS keys must contain "labels" and "voxel_size"'
    h5_relpath = os.path.split(h5name)[1]
    xdmf_filename = h5name.replace('.h5', '.xdmf')
    # write an .h5 file if it does not exist
    if not os.path.exists(h5name):
        print('{} is not found; I am creating one ...'.format(h5name))       
        with h5py.File(h5name, 'w') as hout:
            for key, value in DS.items():
                hout.create_dataset(key, data = value)
        print('Done with saving DS to {}'.format(h5name))

    dims = DS['labels'].shape
    scalar_dims = dims
    vector_dims = dims + (3,)
    tensor_dims = dims + (3, 3,)
    MeshDimensions = (dims[0] + 1, dims[1] + 1, dims[2] + 1)

    MeshDimensionsStr = 'Dimensions="%d %d %d"' % MeshDimensions
    ScalarDimensionsStr = 'Dimensions="%d %d %d"' % scalar_dims
    VectorDimensionsStr = 'Dimensions="%d %d %d %d"' % vector_dims
    TensorDimensionsStr = 'Dimensions="%d %d %d %d %d"' % tensor_dims

    steps = tuple(DS['voxel_size'])

    # Write .xdmf file
    with open(xdmf_filename, 'wt') as fileID:
        fileID.write('<?xml version="1.0"?>\n')
        fileID.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd"[]>\n')
        fileID.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
        fileID.write(' <Domain>\n')
        fileID.write('  <Grid Name="GM3D" GridType="Uniform">\n')
        fileID.write('   <Topology TopologyType="3DCoRectMesh" %s></Topology>\n' % MeshDimensionsStr)
        fileID.write('    <Geometry Type="ORIGIN_DXDYDZ">\n')
        fileID.write('     <!-- Origin  Z, Y, X -->\n')
        fileID.write('     <DataItem Format="XML" Dimensions="3">0 0 0</DataItem>\n')
        fileID.write('     <!-- DxDyDz (Spacing/Resolution) Z, Y, X -->\n')
        fileID.write('     <DataItem Format="XML" Dimensions="3">%.6f %.6f %.6f</DataItem>\n' % steps)
        fileID.write('    </Geometry>\n')

        # iterate over all the keys
        for key_name in DS.keys():
            array = DS[key_name]

            # work out what sort of array we have
            map_shape = array.shape
            n_dims = len(map_shape)
            if n_dims == 3:
                # scalar field
                fileID.write('    <Attribute Name="%s" AttributeType="Scalar" Center="Cell">\n' % key_name)
                fileID.write('      <DataItem Format="HDF" %s NumberType="Float" Precision="6" >%s:/%s</DataItem>\n' % (
                        ScalarDimensionsStr, h5_relpath, '/' + key_name))
                fileID.write('    </Attribute>\n')
            elif n_dims == 4:
                # vector field (like IPF)
                fileID.write('    <Attribute Name="%s" AttributeType="Vector" Center="Cell">\n' % key_name)
                fileID.write('      <DataItem Format="HDF" %s NumberType="Float" Precision="6" >%s:/%s</DataItem>\n' % (
                        VectorDimensionsStr, h5_relpath, '/' + key_name))
                fileID.write('    </Attribute>\n')
            elif n_dims == 5:
                assert map_shape == tensor_dims, "Tensor {} shape {} does not match {}".format(key_name, map_shape, tensor_dims)
                # Define the 9 tensor components, e.g. (xx, xy, xz, yx, yy, yz, zx, zy, zz) for eps_sample
                if key_name == 'eps_sample':
                    tensor_components = [
                        ('xx', 0, 0), ('xy', 0, 1), ('xz', 0, 2),
                        ('yx', 1, 0), ('yy', 1, 1), ('yz', 1, 2),
                        ('zx', 2, 0), ('zy', 2, 1), ('zz', 2, 2)
                    ]
                else:
                    tensor_components = [
                        ('11', 0, 0), ('12', 0, 1), ('13', 0, 2),
                        ('21', 1, 0), ('22', 1, 1), ('23', 1, 2),
                        ('31', 2, 0), ('32', 2, 1), ('33', 2, 2)
                    ]
                for comp_name, i, j in tensor_components:
                    attr_name = "{}_{}".format(key_name, comp_name)
                    fileID.write('    <Attribute Name="%s" AttributeType="Scalar" Center="Cell">\n' % attr_name)
                    fileID.write('      <DataItem ItemType="HyperSlab" %s>\n' % ScalarDimensionsStr)
                    fileID.write('        <DataItem Dimensions="3 5" Format="XML">\n')
                    fileID.write('         %d %d %d %d %d\n' % (0, 0, 0, i, j))  # Origin: fix i, j for the component
                    fileID.write('         %d %d %d %d %d\n' % (1, 1, 1, 1, 1))  # Stride: 1 in all dims
                    fileID.write('         %d %d %d %d %d\n' % (dims[0], dims[1], dims[2], 1, 1))  # Count: full 3D, 1x1 in tensor dims
                    fileID.write('        </DataItem>\n')
                    fileID.write('        <DataItem Format="HDF" NumberType="Float" Precision="6" %s >%s:/%s</DataItem>\n' % (
                        TensorDimensionsStr, h5_relpath, '/' + key_name))
                    fileID.write('      </DataItem>\n')
                    fileID.write('    </Attribute>\n')
                continue
            else:
                continue
        fileID.write('  </Grid>\n')
        fileID.write(' </Domain>\n')
        fileID.write('</Xdmf>\n')
    print('Done with writing xdmf file to {}'.format(xdmf_filename))
    

def indexing_iterative(cf_strong, grains, ds, ucell, pars, ds_max = 1.6, tol_angle = 0.25, tol_pixel =3, peak_assign_tol = 0.25, tol_misori = 3, crystal_system='cubic', **kwargs):
    """
    Using the rest of peaks that are not matched with already-indexed grains to index additional grains
    Then, try to merge the grains by removing duplicate grains (if misori < tol_misori) and return grains_new
    
    Arguments:
    cf_strong         -- ImageD11 colume file object file containing strong peaks
    grains            -- ImageD11 grains object, already-indexed grains
    ds                -- dataset object,   e.g. ds = ImageD11.sinograms.dataset.load(dset_file)
    ucell             -- ImageD11.unitcell object, e.g. ucell = ds.phases.unitcells[phase_str]
    pars              -- ImageD1.parameters object, e.g. pars = ImageD11.parameters.read_par_file(ds.parfile)
    Optional arguments:
    ds_max    -- maximum 1/d
    tol_angle -- tolerance of angles for matching peaks
    tol_piexl -- tolerance of pixels for matching peaks
    peak_assign_tol -- tolerance for assigning peaks to grains
    tol_misori -- minimal misorientation to distinguish grains
    crystal_system -- string name of the crystal system, 'cubic', 'hexagonal', 'orthorhombic', 'tetragonal', 'trigonal', 'monoclinic', 'triclinic'
    
    Returns:
    grains_new -- ImageD11 grains object, with update grain ids starting from 0
    """
    if crystal_system in ['cubic', 'hexagonal', 'orthorhombic', 'tetragonal', 'trigonal', 'monoclinic', 'triclinic']:
        print('{} is OK'.format(crystal_system))
    else:
        raise ValueError('{} is not supported.'.format(crystal_system))
    
    cosine_tol = np.cos(np.radians(90 - ds.ostep))
    args = {
        "cf": cf_strong,
        "unitcell": ucell,
        "dstol": kwargs.get("indexer_ds_tol", 0.008),
        "forgen": kwargs.get("rings_for_gen", [0, 3, 5]),
        "foridx": kwargs.get("rings_for_scoring", [0, 1, 3, 4, 5]),
        "hkl_tols": kwargs.get("hkl_tols_seq", [0.01, 0.02, 0.03, 0.04, 0.05]),
        "fracs": kwargs.get("fracs", [0.9, 0.7]),
        "cosine_tol": kwargs.get("cosine_tol", cosine_tol),
        "max_grains": kwargs.get("max_grains", 1000),
        }
    print('I got indexing parameters as follows:')
    print(args)
    
    # find all the peaks that have matched with already-indexed grains
    print('****************** step 1: find all the peaks that are matched with already-indexed grains using forward model ***************************')
    cf_matched_all, Comp_all = forward_model.forward_match_peaks(cf_strong, grains, ds, ucell, pars, ds_max = ds_max, tol_angle = tol_angle, tol_pixel=tol_pixel, thres_int=None)
    
    # find all the peaks that left unmatched
    print('****************** step 2: find all the peaks that are left unindexed ***************************')
    print('****************** Removing un-matched peaks ...')
    cf_clean = forward_model.cf_remove_unmatched_peaks(cf_strong, cf_matched_all)
    print('****************** Getting the rest of the peaks by comparing the difference between cf_strong and cf_clean')
    cf_strong_rest, cf_clean_diff = forward_model.cf_set_difference(cf_strong, cf_clean, tol = 0.001)
    
    print('Sinograme of all matched peaks')
    forward_model.cf_plot_sino(cf_clean)
    print('Sinograme of the rest peaks')
    forward_model.cf_plot_sino(cf_strong_rest)
    print('{} / {} = {} matched peaks'.format(cf_clean.nrows, cf_strong.nrows, cf_clean.nrows/cf_strong.nrows))
    print('{} / {} = {} peaks left unindexed'.format(cf_strong_rest.nrows, cf_strong.nrows, cf_strong_rest.nrows/cf_strong.nrows))
    
    # specify our ImageD11 indexer with these peaks
    print('****************** step 3: index additional grains with all the left peaks ***************************')
    if 'cell__a' in pars.parameters.keys():
        indexer_rest = ImageD11.indexing.indexer_from_colfile(cf_strong_rest)
    else:
        # new pars format
        indexer_rest = ImageD11.indexing.indexer_from_colfile_and_ucell(cf_strong_rest, ucell=ucell)
    print("Indexing {} peaks".format(cf_strong_rest.nrows))
    
    # USER: set a tolerance in d-space (for assigning peaks to powder rings)
    indexer_rest.ds_tol = args['dstol']

    # change the log level so we can see what the ring assigments look like
    ImageD11.indexing.loglevel = 1

    # assign peaks to powder rings
    indexer_rest.assigntorings()

    # change log level back again
    ImageD11.indexing.loglevel = 3
    
    # let's plot the assigned peaks
    fig, ax = plt.subplots(layout='constrained', figsize=(10,5))

    # indexer.ra is the ring assignments
    ax.scatter(cf_strong_rest.ds, cf_strong_rest.eta, c=indexer_rest.ra, cmap='tab20', s=1)
    ax.plot( ucell.ringds, [0,]*len(ucell.ringds), '|', ms=90, c="red")
    ax.set_xlim(cf_strong_rest.ds.min()-0.05, cf_strong_rest.ds.max()+0.05)
    ax.set_xlabel("d-star")
    ax.set_ylabel("eta")

    plt.show()
    
    grains_rest, indexer_rest = utils.do_index(**args)
    print('Found {} grains!'.format(len(grains_rest)))

    # add temporary grain IDs to the grains
    for ginc, g in enumerate(grains_rest):
        g.gid = ginc + len(grains)
    
    # plotting unit cell lenth
    mean_unit_cell_lengths_rest = [np.cbrt(np.linalg.det(g.ubi)) for g in grains_rest]
    fig, ax = plt.subplots()
    ax.plot(mean_unit_cell_lengths_rest)
    ax.set_xlabel("Grain ID")
    ax.set_ylabel("Unit cell length")
    plt.show()
    a0_rest = np.median(mean_unit_cell_lengths_rest)
    print('Average lattice parameter of newly indexed additional grains: {} angstrom'.format(a0_rest))
          
    # assign peaks to grains
    utils.assign_peaks_to_grains(grains_rest, cf_strong_rest, tol=peak_assign_tol)
    utils.plot_index_results(indexer_rest, cf_strong_rest, 'First attempt')
    
    # merge duplicated grains
    print('****************** step 4: merge grains by removing duplicated grains based on their misorientations ***************************')
    grains_new = merge_grains(grains, grains_rest, crystal_system = crystal_system, tol_misori = tol_misori)
    print('****************** Number of previously indexed, currently indexed and final merged grains:')
    print(len(grains), len(grains_rest), len(grains_new))
          
    # update grain IDs to the grains
    for ginc, g in enumerate(grains_new):
        g.gid = ginc

    # plotting unit cell lenth for grains_new
    mean_unit_cell_lengths_new = [np.cbrt(np.linalg.det(g.ubi)) for g in grains_new]
    fig, ax = plt.subplots()
    ax.plot(mean_unit_cell_lengths_new)
    ax.set_xlabel("Grain ID")
    ax.set_ylabel("Unit cell length")
    plt.show()

    a0_new = np.median(mean_unit_cell_lengths_new)

    print('Average lattice parameter of all grains: {} angstrom'.format(a0_new))
    print('Done !')
    # if you want to assign peaks to grains_new, do this
    # utils.assign_peaks_to_grains(grains_new, cf_strong, tol=peak_assign_tol)
    # utils.plot_grain_sinograms(grains_new, cf_strong, min(len(grains_new), 25))
    return grains_new


def remove_gid_from_grains(grains, gid_to_remove = None):
    """
    Remove specific grains from the existing grains
    
    Arguments:
    gid_to_remove -- a list of grain IDs to be removed
    
    Returns:
    grains_new -- ImageD11 grains object, with update grain ids starting from 0
    """
    if gid_to_remove is not None:
        grains_new = []
        print('The following grain IDs will be removed and new grains object will be produced (note that thhe grain ID may be altered):')
        print(gid_to_remove)
        for g in grains:
            if g.gid in gid_to_remove:
                continue
            else:
                grains_new.append(g)
        
        # update grain IDs to the grains
        for ginc, g in enumerate(grains_new):
            g.gid = ginc
    else:
        print('Nothing needs to be done')
        grains_new = grains
    return grains_new


def merge_grains(grains, grains_rest, crystal_system = 'cubic', tol_misori=3):
    """
    Identify duplicate grains and then make a new grains_new object
    
    Arguments:
    grains         -- ImageD11 grains object, assuming from the first indexing
    grains_rest    -- ImageD11 grains object, assuming from the second indexing
    crystal_structure -- symmetry class with crystal symmetry operators
    tol_ang            -- tolerance for misorientation to distinguish grains [deg]

    Returns:
    grains_new -- ImageD11 grains object, with update grain ids starting from 0
    """
    
    if crystal_system in ['cubic', 'hexagonal', 'orthorhombic', 'tetragonal', 'trigonal', 'monoclinic', 'triclinic']:
        crystal_structure = Symmetry[crystal_system]
    else:
        raise ValueError('{} is not supported.'.format(crystal_system))
    
    grains_new = []
    oris_duplicate = set()  # Use a set for faster lookups
    tol_rad = np.deg2rad(tol_misori)  # Convert tolerance to radians once
    print('tol_misori = {} or {} rad'.format(tol_misori, tol_rad))
    
    # Precompute inverses for grains and grains_rest
    inv_grains = {g.gid: np.linalg.inv(g.U) for g in grains}
    inv_grains_rest = {g_rest.gid: np.linalg.inv(g_rest.U) for g_rest in grains_rest}
    
    # Check for duplicates and merge
    for g in grains:
        ori1 = inv_grains[g.gid]
        for g_rest in grains_rest:
            if g_rest.gid not in oris_duplicate:
                ori2 = inv_grains_rest[g_rest.gid]
                the_angle, the_axis, the_axis_xyz = disorientation(ori1, ori2, crystal_structure=crystal_structure)
                if the_angle < tol_rad:  # Compare directly in radians
                    print("Angle: {}, Axis: {}, Axis XYZ: {}".format(np.rad2deg(the_angle), the_axis, the_axis_xyz))
                    oris_duplicate.add(g_rest.gid)
                    print('grain no. {} from grains matches with grain no. {} from grains_rest'.format(g.gid, g_rest.gid))
        grains_new.append(g)
    
    # get unique grains from grains_rest
    grains_new.extend(g for g in grains_rest if g.gid not in oris_duplicate)
    
    # update grain IDs to the grains
    for ginc, g in enumerate(grains_new):
        g.gid = ginc
        
    return grains_new


def pairing_grains(grains1, grains2, crystal_system='cubic', tol_misori=3, tol_dis=50, tol_int=[0.0, 0.0]):
    """
    pairing grains from two different grains objects, useful for tracking the same grain across different datasets
    about 3 times faster with parallel, it may be speed up further
    
    Arguments:
    grains1         -- ImageD11 grains object
    grains2         -- ImageD11 grains object assuming from the second indexing
    crystal_structure -- symmetry class with crystal symmetry operators
    tol_ang            -- tolerance for misorientation to distinguish grains [deg]
    tol_dis            -- tolerance for distance-to-volume-center to distinguish grains [um]
    tol_int            -- tolerance for intensity (relative difference in mean intensity) to distinguish grains [-]

    Returns:
    PairIndex -- numpy array with a shape of N*18, column means [i, j, rod_i, rod_j, int_sum1, int_sum2, mis_ori, dis, COM_i, COM_j]
    """   

    if crystal_system in ['cubic', 'hexagonal', 'orthorhombic', 'tetragonal', 'trigonal', 'monoclinic', 'triclinic']:
        crystal_structure = Symmetry[crystal_system]
    else:
        raise ValueError('{} is not supported.'.format(crystal_system))

    print('Got {} grains from grain1 and {} grains from grains2'.format(len(grains1), len(grains2)))
    print('I will use the following tolerances for pairing grains:')
    print('tol_misori = {} deg'.format(tol_misori))
    print('tol_dis = {} um'.format(tol_dis))
    if not tol_int == [0.0, 0.0]:
        print('tol_int = {}'.format(tol_int))
    
    tol_rad = np.deg2rad(tol_misori)
    inv_grains1 = np.array([np.linalg.inv(g.U) for g in grains1])
    inv_grains2 = np.array([np.linalg.inv(g.U) for g in grains2])
    
    translations1 = np.array([g.translation if g.translation is not None else [0.0, 0.0, 0.0] for g in grains1])
    translations2 = np.array([g.translation if g.translation is not None else [0.0, 0.0, 0.0] for g in grains2])
    
    intensities1 = [read_intensity_info(g.intensity_info)['sum_of_all'] for g in grains1]
    intensities2 = [read_intensity_info(g.intensity_info)['sum_of_all'] for g in grains2]

    def process_pair(i, j):
        ori1, ori2 = inv_grains1[i], inv_grains2[j]
        pos1, pos2 = translations1[i], translations2[j]
        intensity1, intensity2 = intensities1[i], intensities2[j]

        # Calculate misorientation
        the_angle, _, _ = disorientation(ori1, ori2, crystal_structure=crystal_structure)

        # Calculate distance
        dis = np.linalg.norm(pos1 - pos2)

        # Check tolerances
        if the_angle < tol_rad and dis < tol_dis:
            if tol_int == [0.0, 0.0]:
                # return [i, j, grains1[i].Rod[0], grains1[i].Rod[1], grains1[i].Rod[2],
                #         grains2[j].Rod[0], grains2[j].Rod[1], grains2[j].Rod[2],
                #         intensity1, intensity2, np.rad2deg(the_angle), dis,
                #         *pos1, *pos2]
                return [i, j, grains1[i].Rod[0], grains1[i].Rod[1], grains1[i].Rod[2],
                        grains2[j].Rod[0], grains2[j].Rod[1], grains2[j].Rod[2],
                        intensity1, intensity2, np.rad2deg(the_angle), dis] + list(pos1) + list(pos2)
            else:
                intensity_diff = abs(intensity1 - intensity2) / intensity1
                if tol_int[0] <= intensity_diff <= tol_int[1]:
                    # return [i, j, grains1[i].Rod[0], grains1[i].Rod[1], grains1[i].Rod[2],
                    #         grains2[j].Rod[0], grains2[j].Rod[1], grains2[j].Rod[2],
                    #         intensity1, intensity2, np.rad2deg(the_angle), dis,
                    #         *pos1, *pos2]
                    return [i, j, grains1[i].Rod[0], grains1[i].Rod[1], grains1[i].Rod[2],
                        grains2[j].Rod[0], grains2[j].Rod[1], grains2[j].Rod[2],
                        intensity1, intensity2, np.rad2deg(the_angle), dis] + list(pos1) + list(pos2)
        return None
    
    total_pairs = len(grains1) * len(grains2)
    print('{} comparisons will be performed.'.format(total_pairs))
    
    start_time = time.time()
    pairs = Parallel(n_jobs=-1, backend='loky')(
        delayed(process_pair)(i, j) for i, j in tqdm([(i, j) for i in range(len(grains1)) for j in range(len(grains2))], total=total_pairs)
    )
    elapsed_time = end_time = time.time() - start_time
    print("Time taken: {} seconds".format(elapsed_time))
    print('PairIndex column info: [i, j, rod_i, rod_j, int_sum1, int_sum2, mis_ori, dis, COM_i, COM_j]')
    
    # Filter and convert to NumPy array
    PairIndex = np.array([pair for pair in pairs if pair is not None])
    print('{} pairs found.'.format(len(PairIndex)))

    return PairIndex


def read_intensity_info(text):
    # Input string example
    # text = 'sum_of_all = 4989843.970096 , middle 343 from 0.000000 to 180.000000 in tth: median = 5725.219460 , min = 133.669141 , max = 90410.322693 , mean = 14010.814810 , std = 19810.559831 , n = 343\n'

    # Regular expression to match key-value pairs
    pattern = r'(\w+(?:_\w+)?) = ([\d\.]+)'

    # Extract matches into a dictionary
    data = {key: float(value) for key, value in re.findall(pattern, text)}

    return data


def get_mean_rod(rod, kmeans_flag=False, auto_check = True):
    """
    Compute the mean Rodrigues vector with optional k-means clustering.
    
    Parameters:
    - rod: numpy array of shape (n, 3) or (3, n) representing Rodrigues vectors
    - kmeans_flag: boolean, whether to perform k-means clustering

    Returns:
    - rod_mean: numpy array of shape (3,) representing the mean Rodrigues vector
    """
    rod = np.array(rod)
    
    # Ensure rod has shape (n, 3)
    if auto_check and (rod.shape[1] != 3 and rod.shape[0] == 3):
        rod = rod.T
    # Remove rows with any NaN values
    rod = rod[~np.isnan(rod).any(axis=1)]
    
    if kmeans_flag and len(rod) > 1:
        # Perform k-means clustering
        kmeans = KMeans(n_clusters=2, random_state=0).fit(rod)
        C = kmeans.cluster_centers_  # Centroids
        idx = kmeans.labels_         # Cluster indices
        
        # Compute the distance between the two cluster centroids
        C_diff = np.linalg.norm(C[0] - C[1])
        
        if C_diff > 0.5:
            # Separate into two clusters
            rod1 = rod[idx == 0]
            rod2 = rod[idx == 1]
            
            # Choose the larger cluster
            if len(rod1) > len(rod2):
                larger_cluster = rod1
            else:
                larger_cluster = rod2
            
            # Compute the median for each coordinate
            mid_ind = len(larger_cluster) // 2
            rod_mean = np.zeros(3)
            for i in range(3):
                rod_temp = np.sort(larger_cluster[:, i])
                rod_mean[i] = rod_temp[mid_ind]
        else:
            # If centroids are close, use the first centroid
            rod_mean = C[0]
    else:
        # Compute the median for the entire set
        mid_ind = len(rod) // 2
        rod_mean = np.zeros(3)
        for i in range(3):
            rod_temp = np.sort(rod[:, i])
            rod_mean[i] = rod_temp[mid_ind]
    
    return rod_mean


# this class comes from pymicro/crystal/lattice.py
class Symmetry(enum.Enum):
    """
    Class to describe crystal symmetry defined by its Laue class symbol.
    """
    cubic = 'm3m'
    hexagonal = '6/mmm'
    orthorhombic = 'mmm'
    tetragonal = '4/mmm'
    trigonal = 'bar3m'
    monoclinic = '2/m'
    triclinic = 'bar1'

    @staticmethod
    def from_string(s):
        if s == 'cubic':
            return Symmetry.cubic
        elif s == 'hexagonal':
            return Symmetry.hexagonal
        elif s == 'orthorhombic':
            return Symmetry.orthorhombic
        elif s == 'tetragonal':
            return Symmetry.tetragonal
        elif s == 'trigonal':
            return Symmetry.trigonal
        elif s == 'monoclinic':
            return Symmetry.monoclinic
        elif s == 'triclinic':
            return Symmetry.triclinic
        else:
            return None

    def to_string(self):
        if self is Symmetry.cubic:
            return 'cubic'
        elif self is Symmetry.hexagonal:
            return 'hexagonal'
        elif self is Symmetry.orthorhombic:
            return 'orthorhombic'
        elif self is Symmetry.tetragonal:
            return 'tetragonal'
        elif self is Symmetry.trigonal:
            return 'trigonal'
        elif self is Symmetry.monoclinic:
            return 'monoclinic'
        elif self is Symmetry.triclinic:
            return 'triclinic'
        else:
            return None

    @staticmethod
    def from_dream3d(n):
        if n in [1, 3]:
            return Symmetry.cubic
        elif n in[0, 2]:
            return Symmetry.hexagonal
        elif n == 6:
            return Symmetry.orthorhombic
        elif n in [7, 8]:
            return Symmetry.tetragonal
        elif n in [9, 10]:
            return Symmetry.trigonal
        elif n == 5:
            return Symmetry.monoclinic
        elif n == 4:
            return Symmetry.triclinic
        else:
            return None

    @staticmethod
    def from_space_group(space_group_number):
        """Create an instance of the `Symmetry` class from a TSL symmetry
        number.

        :raise ValueError: if the space_group_number is not between 1 and 230.
        :param int space_group_number: the number asociated with the
        space group (between 1 and 230).
        :return: an instance of the `Symmetry` class
        """
        if space_group_number < 1 or space_group_number > 230:
          raise ValueError('space_group_number must be between 1 and 230')
          return None
        if space_group_number <= 2:
            return Symmetry.triclinic
        elif space_group_number <= 15:
            return Symmetry.monoclinic
        elif space_group_number <= 74:
            return Symmetry.orthorhombic
        elif space_group_number <= 142:
            return Symmetry.tetragonal
        elif space_group_number <= 167:
            return Symmetry.trigonal
        elif space_group_number <= 194:
            return Symmetry.hexagonal
        else:
            return Symmetry.cubic

    @staticmethod
    def from_tsl(tsl_number):
        """Create an instance of the `Symmetry` class from a TSL symmetry
        number.

        :return: an instance of the `Symmetry` class
        """
        if tsl_number == 43:
            return Symmetry.cubic
        elif tsl_number == 62:
            return Symmetry.hexagonal
        elif tsl_number == 22:
            return Symmetry.orthorhombic
        elif tsl_number == 42:
            return Symmetry.tetragonal
        elif tsl_number == 32:
            return Symmetry.trigonal
        elif tsl_number == 2:
            return Symmetry.monoclinic
        elif tsl_number == 1:
            return Symmetry.triclinic
        else:
            return None

    @staticmethod
    def to_tsl(symmetry):
        """Convert a type of crystal `Symmetry` in the corresponding TSL number.

        :return: the TSL number for this symmetry.
        """
        if symmetry is Symmetry.cubic:
            return 43
        elif symmetry is Symmetry.hexagonal:
            return 62
        elif symmetry is Symmetry.orthorhombic:
            return 22
        elif symmetry is Symmetry.tetragonal:
            return 42
        elif symmetry is Symmetry.trigonal:
            return 32
        elif symmetry is Symmetry.monoclinic:
            return 2
        elif symmetry is Symmetry.triclinic:
            return 1
        else:
            return None

    def symmetry_operators(self, use_miller_bravais=False):
        """Define the equivalent crystal symmetries.
        Note that it takes the highest symmetry in each of the 7 crystal systems
        This may not be quite right if a specific point group with less symmetry is desired

        Those come from Randle & Engler, 2000. For instance in the cubic
        crystal struture, for instance there are 24 equivalent cube orientations.

        :return array: A numpy array of shape (n, 3, 3) where n is the \
        number of symmetries of the given crystal structure.
        """
        if self is Symmetry.cubic:
            sym = np.zeros((24, 3, 3), dtype=float)
            sym[0] = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
            sym[1] = np.array([[0., 0., -1.], [0., -1., 0.], [-1., 0., 0.]])
            sym[2] = np.array([[0., 0., -1.], [0., 1., 0.], [1., 0., 0.]])
            sym[3] = np.array([[-1., 0., 0.], [0., 1., 0.], [0., 0., -1.]])
            sym[4] = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
            sym[5] = np.array([[1., 0., 0.], [0., 0., -1.], [0., 1., 0.]])
            sym[6] = np.array([[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]])
            sym[7] = np.array([[1., 0., 0.], [0., 0., 1.], [0., -1., 0.]])
            sym[8] = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
            sym[9] = np.array([[-1., 0., 0.], [0., -1., 0.], [0., 0., 1.]])
            sym[10] = np.array([[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]])
            sym[11] = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
            sym[12] = np.array([[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]])
            sym[13] = np.array([[0., 0., -1.], [-1., 0., 0.], [0., 1., 0.]])
            sym[14] = np.array([[0., -1., 0.], [0., 0., 1.], [-1., 0., 0.]])
            sym[15] = np.array([[0., 1., 0.], [0., 0., -1.], [-1., 0., 0.]])
            sym[16] = np.array([[0., 0., -1.], [1., 0., 0.], [0., -1., 0.]])
            sym[17] = np.array([[0., 0., 1.], [-1., 0., 0.], [0., -1., 0.]])
            sym[18] = np.array([[0., -1., 0.], [0., 0., -1.], [1., 0., 0.]])
            sym[19] = np.array([[0., 1., 0.], [1., 0., 0.], [0., 0., -1.]])
            sym[20] = np.array([[-1., 0., 0.], [0., 0., 1.], [0., 1., 0.]])
            sym[21] = np.array([[0., 0., 1.], [0., -1., 0.], [1., 0., 0.]])
            sym[22] = np.array([[0., -1., 0.], [-1., 0., 0.], [0., 0., -1.]])
            sym[23] = np.array([[-1., 0., 0.], [0., 0., -1.], [0., -1., 0.]])
        elif self is Symmetry.hexagonal:
            if use_miller_bravais:
              # using the Miller-Bravais representation here
              sym = np.zeros((12, 4, 4), dtype=float)
              sym[0] = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
              sym[1] = np.array([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
              sym[2] = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1]])
              sym[3] = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
              sym[4] = np.array([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, -1]])
              sym[5] = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, -1]])
              sym[6] = np.array([[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]])
              sym[7] = np.array([[0, 0, -1, 0], [-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 0, 1]])
              sym[8] = np.array([[0, -1, 0, 0], [0, 0, -1, 0], [-1, 0, 0, 0], [0, 0, 0, 1]])
              sym[9] = np.array([[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
              sym[10] = np.array([[0, 0, -1, 0], [-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 0, -1]])
              sym[11] = np.array([[0, -1, 0, 0], [0, 0, -1, 0], [-1, 0, 0, 0], [0, 0, 0, -1]])
            else:
              sym = np.zeros((12, 3, 3), dtype=float)
              s60 = np.sin(60 * np.pi / 180)
              sym[0] = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
              sym[1] = np.array([[0.5, s60, 0.], [-s60, 0.5, 0.], [0., 0., 1.]])
              sym[2] = np.array([[-0.5, s60, 0.], [-s60, -0.5, 0.], [0., 0., 1.]])
              sym[3] = np.array([[-1., 0., 0.], [0., -1., 0.], [0., 0., 1.]])
              sym[4] = np.array([[-0.5, -s60, 0.], [s60, -0.5, 0.], [0., 0., 1.]])
              sym[5] = np.array([[0.5, -s60, 0.], [s60, 0.5, 0.], [0., 0., 1.]])
              sym[6] = np.array([[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]])
              sym[7] = np.array([[0.5, s60, 0.], [s60, -0.5, 0.], [0., 0., -1.]])
              sym[8] = np.array([[-0.5, s60, 0.], [s60, 0.5, 0.], [0., 0., -1.]])
              sym[9] = np.array([[-1., 0., 0.], [0., 1., 0.], [0., 0., -1.]])
              sym[10] = np.array([[-0.5, -s60, 0.], [-s60, 0.5, 0.], [0., 0., -1.]])
              sym[11] = np.array([[0.5, -s60, 0.], [-s60, -0.5, 0.], [0., 0., -1.]])               
        elif self is Symmetry.trigonal:
              sym = np.zeros((6, 3, 3), dtype=float)
              s60 = np.sin(60 * np.pi / 180)
              sym[0] = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
              sym[1] = np.array([[-0.5, -s60, 0.], [s60, -0.5, 0.], [0., 0., 1.]])
              sym[2] = np.array([[-0.5, s60, 0.], [-s60, -0.5, 0.], [0., 0., 1.]])
              sym[3] = np.array([[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]])    
              sym[4] = np.array([[-0.5, s60, 0.], [s60, 0.5, 0.], [0., 0., -1.]])
              sym[5] = np.array([[-0.5, -s60, 0.], [-s60, 0.5, 0.], [0., 0., -1.]])
        elif self is Symmetry.orthorhombic:
            sym = np.zeros((4, 3, 3), dtype=float)
            sym[0] = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
            sym[1] = np.array([[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]])
            sym[2] = np.array([[-1., 0., 0.], [0., 1., 0.], [0., 0., -1.]])
            sym[3] = np.array([[-1., 0., 0.], [0., -1., 0.], [0., 0., 1.]])
        elif self is Symmetry.tetragonal:
            sym = np.zeros((8, 3, 3), dtype=float)
            sym[0] = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
            sym[1] = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
            sym[2] = np.array([[-1., 0., 0.], [0., -1., 0.], [0., 0., 1.]])
            sym[3] = np.array([[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]])
            sym[4] = np.array([[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]])
            sym[5] = np.array([[-1., 0., 0.], [0., 1., 0.], [0., 0., -1.]])
            sym[6] = np.array([[0., 1., 0.], [1., 0., 0.], [0., 0., -1.]])
            sym[7] = np.array([[0., -1., 0.], [-1., 0., 0.], [0., 0., -1.]])
        elif self is Symmetry.monoclinic:
            sym = np.zeros((2, 3, 3), dtype=float)
            sym[0] = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
            sym[1] = np.array([[-1., 0., 0.], [0., -1., 0.], [0., 0., 1.]])
        elif self is Symmetry.triclinic:
            sym = np.zeros((1, 3, 3), dtype=float)
            sym[0] = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        else:
            raise ValueError('warning, symmetry not supported: %s' % self)
        return sym

    def move_vector_to_FZ(self, v):
        """
        Move the vector to the Fundamental Zone of a given `Symmetry` instance.

        :param v: a 3 components vector.
        :return: a new 3 components vector in the fundamental zone.
        """
        omegas = []  # list to store all the rotation angles
        syms = self.symmetry_operators()
        for sym in syms:
            # apply symmetry to the vector and compute the corresponding angle
            v_sym = np.dot(sym, v)
            omega = 2 * np.arctan(np.linalg.norm(v_sym)) * 180 / np.pi
            omegas.append(omega)
        # the fundamental zone corresponds to the minimum angle
        index = np.argmin(omegas)
        return np.dot(syms[index], v)

    def move_rotation_to_FZ(self, g, verbose=False):
        """Compute the rotation matrix in the Fundamental Zone of a given
        `Symmetry` instance.

        The principle is to apply all symmetry operators to the rotation matrix
        and identify which one yield the smallest rotation angle. The
        corresponding rotation is then returned. This computation is vectorized
        to save time.

        :param g: a 3x3 matrix representing the rotation.
        :param bool verbose: flag for verbose mode.
        :return: a new 3x3 matrix for the rotation in the fundamental zone.
        """
        syms = self.symmetry_operators()
        g_syms = np.dot(syms, g)
        traces = np.trace(g_syms, axis1=1, axis2=2)
        omegas = np.arccos(0.5 * (traces - 1))
        index = np.argmin(omegas)
        if verbose:
            print(traces)
            print(omegas)
            print('moving to FZ, index = %d' % index)
        return g_syms[index]

    def lattice_parameters_number(self):
        """Return the number of parameter associated with a lattice of this
        symmetry.

        :return: the number of parameters.
        """
        if self is Symmetry.cubic:
            return 1
        elif self in [Symmetry.hexagonal, Symmetry.trigonal, Symmetry.tetragonal]:
            return 2
        elif self is Symmetry.orthorhombic:
            return 3
        elif self is Symmetry.monoclinic:
            return 4
        else:  # triclinic case
            return 6

    def elastic_constants_number(self):
        """Return the number of independent elastic constants for this symmetry.

        :return: the number of elastic constants.
        """
        if self is Symmetry.cubic:
            return 3
        elif self is Symmetry.hexagonal:
            return 5
        elif self is Symmetry.tetragonal:
            return 6
        elif self is Symmetry.orthorhombic:
            return 9
        elif self is Symmetry.monoclinic:
            return 13
        else:  # triclinic case
            return 21

    def stiffness_matrix(self, elastic_constants):
        """Build the stiffness matrix for this symmetry using Voigt convention.

        The Voigt notation contracts 2 tensor indices into a single index:
        11 -> 1, 22 -> 2, 33 -> 3, 23 -> 4, 31 -> 5, 12 -> 6

        :param list elastic_constants: the elastic constants (the number must
            correspond to the type of symmetry, eg 3 for cubic).
        :return ndarray: a numpy array of shape (6, 6) representing
            the stiffness matrix.
        """
        if self is Symmetry.cubic:
            if len(elastic_constants) != 3:
                raise ValueError('Error: need 3 elastic constants for cubic '
                                 'symmetry, got %d' % len(elastic_constants))
            C11, C12, C44 = elastic_constants
            C = np.array([[C11, C12, C12,   0,   0,   0],
                          [C12, C11, C12,   0,   0,   0],
                          [C12, C12, C11,   0,   0,   0],
                          [  0,   0,   0, C44,   0,   0],
                          [  0,   0,   0,   0, C44,   0],
                          [  0,   0,   0,   0,   0, C44]])
            return C
        elif self is Symmetry.hexagonal:
            if len(elastic_constants) != 5:
                raise ValueError('Error: need 5 elastic constants for hexagonal '
                                 'symmetry, got %d' % len(elastic_constants))
            C11, C12, C13, C33, C44 = elastic_constants
            C66 = (C11 - C12) / 2
            C = np.array([[C11, C12, C13,   0,   0,   0],
                          [C12, C11, C13,   0,   0,   0],
                          [C13, C13, C33,   0,   0,   0],
                          [  0,   0,   0, C44,   0,   0],
                          [  0,   0,   0,   0, C44,   0],
                          [  0,   0,   0,   0,   0, C66]])
            return C
        elif self is Symmetry.tetragonal:
            if len(elastic_constants) != 6:
                raise ValueError('Error: need 6 elastic constants for tetragonal '
                                 'symmetry, got %d' % len(elastic_constants))
            C11, C12, C13, C33, C44, C66 = elastic_constants
            C = np.array([[C11, C12, C13,   0,   0,   0],
                          [C12, C11, C13,   0,   0,   0],
                          [C13, C13, C33,   0,   0,   0],
                          [  0,   0,   0, C44,   0,   0],
                          [  0,   0,   0,   0, C44,   0],
                          [  0,   0,   0,   0,   0, C66]])
            return C
        elif self is Symmetry.orthorhombic:
            if len(elastic_constants) != 9:
                raise ValueError('Error: need 9 elastic constants for tetragonal '
                                 'symmetry, got %d' % len(elastic_constants))
            C11, C12, C13, C22, C23, C33, C44, C55, C66 = elastic_constants
            C = np.array([[C11, C12, C13,   0,   0,   0],
                          [C12, C22, C23,   0,   0,   0],
                          [C13, C23, C33,   0,   0,   0],
                          [  0,   0,   0, C44,   0,   0],
                          [  0,   0,   0,   0, C55,   0],
                          [  0,   0,   0,   0,   0, C66]])
            return C
        elif self is Symmetry.monoclinic:
            if len(elastic_constants) != 13:
                raise ValueError('Error: need 13 elastic constants for monoclinic '
                                 'symmetry, got %d' % len(elastic_constants))
            C11, C12, C13, C16, C22, C23, C26, C33, C36, C44, C45, \
            C55, C66 = elastic_constants
            C = np.array([[C11, C12, C13,   0,   0, C16],
                          [C12, C22, C23,   0,   0, C26],
                          [C13, C23, C33,   0,   0, C36],
                          [  0,   0,   0, C44, C45,   0],
                          [  0,   0,   0, C45, C55,   0],
                          [C16, C26, C36,   0,   0, C66]])
            return C
        elif self is Symmetry.triclinic:
            if len(elastic_constants) != 21:
                raise ValueError('Error: need 21 elastic constants for triclinic '
                                 'symmetry, got %d' % len(elastic_constants))
            C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, \
            C34, C35, C36, C44, C45, C46, C55, C56, C66 = elastic_constants
            C = np.array([[C11, C12, C13, C14, C15, C16],
                          [C12, C22, C23, C24, C25, C26],
                          [C13, C23, C33, C34, C35, C36],
                          [C14, C24, C34, C44, C45, C46],
                          [C15, C25, C35, C45, C55, C56],
                          [C16, C26, C36, C46, C56, C66]])
            return C
        else:
            raise ValueError('warning, symmetry not supported: %s' % self)

    @staticmethod
    def orthotropic_constants_from_stiffness(C):
        """Return orthotropic elastic constants from stiffness matrix.

        :param ndarray C: a numpy array of shape (6, 6) representing
            the stiffness matrix.
        :return dict ortho_elas: a dictionary of orthotropic elastic constants
            corresponding to the input stiffness matrix. Keys are
            'E1','E2','E3','nu12','nu13','nu23','G12','G13','G23'
        """
        # compute the compliance matrix
        S = np.linalg.inv(C)
        # compute the orthotropic elastic constants
        ortho_elas = dict()
        ortho_elas['E1'] = 1 / S[0, 0]
        ortho_elas['E2'] = 1 / S[1, 1]
        ortho_elas['E3'] = 1 / S[2, 2]
        ortho_elas['Nu12'] = -ortho_elas['E1'] * S[1, 0]
        ortho_elas['Nu13'] = -ortho_elas['E1'] * S[2, 0]
        ortho_elas['Nu23'] = -ortho_elas['E2'] * S[2, 1]
        ortho_elas['G12'] = 1 / S[5, 5]
        ortho_elas['G13'] = 1 / S[4, 4]
        ortho_elas['G23'] = 1 / S[3, 3]
        # return a dictionary populated with the relevant values
        return ortho_elas


# these functions are adapted from pymicro/crystal/microstructure.py
def disorientation_list(ori1_list, ori2_list, crystal_structure=Symmetry.triclinic):
    """
    Compute the misorientation for lists of crystal orientations.
    
    Arguments:
    ori1_list: A list or array of 3x3 rotation matrices (shape: [N, 3, 3])
        describing the first set of orientations.
    ori2_list: A list or array of 3x3 rotation matrices (shape: [N, 3, 3])
        describing the second set of orientations.
    crystal_structure: An instance of the `Symmetry` class describing the crystal symmetry.
    
    Returns:
    A list of numpy arrays:
    - angles: Misorientation angles in radians (shape: [N]).
    - axes: Misorientation axes in crystal coordinates (shape: [N, 3]).
    - axes_xyz: Misorientation axes in sample coordinates (shape: [N, 3]).
    """
    # Ensure inputs are numpy arrays
    ori1_list = np.asarray(ori1_list)
    ori2_list = np.asarray(ori2_list)
    assert len(ori1_list) == len(ori2_list), "ori1_list must be the same length as ori2_list"
    symmetries = crystal_structure.symmetry_operators()  # Shape: [num_sym_ops, 3, 3]
    if len(ori1_list.shape) == 2:
        num_orientations = 1
    elif len(ori1_list.shape) == 3:
        num_orientations = len(ori1_list)
    else:
        print('Only supports 2 or 3 dimensional inputs.')

    # Initialize outputs
    angles = []
    axes = []
    axes_xyz = []

    # Iterate through all orientations
    if num_orientations > 1:
        for idx in range(num_orientations):
            ori1 = ori1_list[idx]
            ori2 = ori2_list[idx]
            print(ori1)

            the_angle, the_axis, the_axis_xyz = disorientation(ori1, ori2, crystal_structure=crystal_structure)
            angles.append(the_angle)
            axes.append(the_axis)
            axes_xyz.append(the_axis_xyz)
    else:
        the_angle, the_axis, the_axis_xyz = disorientation(ori1_list, ori2_list, crystal_structure=crystal_structure)
        angles.append(the_angle)
        axes.append(the_axis)
        axes_xyz.append(the_axis_xyz)

    return angles, axes, axes_xyz


def disorientation(ori1, ori2, crystal_structure=Symmetry.triclinic):
    """
    Compute the disorientation between two orientations using vectorized symmetry operations.

    Arguments:
    ori1: A 3x3 rotation matrix representing the first orientation.
    ori2: A 3x3 rotation matrix representing the second orientation.
    crystal_structure: A `Symmetry` class instance describing the crystal symmetry.
    
    Returns: Tuple (misorientation angle in radians, misorientation axis in crystal coordinates,
              misorientation axis in sample coordinates).
                  
    Replace the function disorientation_deprecated, speed up by 10x faster
    """
    if (ori1 == ori2).all():
        the_angle = 0.0
        the_axis = np.array([0.0, 0.0, 1.0])
        the_axis_xyz = np.array([0.0, 0.0, 1.0])
    else:
        the_angle = np.pi
        the_axis = np.array([0., 0., 1.])
        the_axis_xyz = np.array([0., 0., 1.])

        # Get symmetry operators as an array of shape (N, 3, 3)
        symmetries = crystal_structure.symmetry_operators()  # Shape: [N, 3, 3]
        num_sym_ops = symmetries.shape[0]

        # Apply symmetries to ori1 and ori2
        oj_list = np.einsum('nij,jk->nik', symmetries, ori1)
        oi_list = np.einsum('nij,jk->nik', symmetries, ori2)

        # Compute all combinations of symmetrized ori1 and ori2
        delta_list_AB = np.einsum('aij,bkj->abik', oi_list, oj_list)
        delta_list_BA = np.einsum('aij,bkj->abik', oj_list, oi_list)
        
        # Flatten (2*Nsymmetry*Nsymmetry, 3, 3)
        delta_list = np.concatenate([
            delta_list_AB.reshape(-1, 3, 3),
            delta_list_BA.reshape(-1, 3, 3)
        ])

        # Calculate misorientation angles for all delta matrices
        mis_angles = np.array([misorientation_angle_from_delta(delta) for delta in delta_list])

        # Find the minimum non-zero misorientation angle
        nonzero_mask = mis_angles > 0
        try:
            if np.any(nonzero_mask):
                min_idx = np.argmin(mis_angles[nonzero_mask])
                min_idx = np.where(nonzero_mask)[0][min_idx]
        except ValueError:
            min_idx = 0
            print("No non-zero misorientation angles found.")
        except Exception as e:
            min_idx = 0
            print("An unexpected error occurred: {}".format(e))            

        the_angle = mis_angles[min_idx]

        # Compute misorientation axis for the minimum angle
        delta_min = delta_list[min_idx]
        the_axis = misorientation_axis_from_delta(delta_min)

        # Calculate the axis in sample coordinates
        oi_min = oi_list[min_idx // num_sym_ops]
        the_axis_xyz = np.dot(oi_min.T, the_axis)

    return the_angle, the_axis, the_axis_xyz


def disorientation_deprecated(ori1, ori2, crystal_structure=Symmetry.triclinic):
    """This function is deprecated. Use `new_function` instead."""
    
    """Compute the disorientation another crystal orientation.

    Considering all the possible crystal symmetries, the disorientation
    is defined as the combination of the minimum misorientation angle
    and the misorientation axis lying in the fundamental zone, which
    can be used to bring the two lattices into coincidence.

    .. note::

     Both orientations are supposed to have the same symmetry. This is not
     necessarily the case in multi-phase materials.

    :param orientation: an instance of
        :py:class:`~pymicro.crystal.microstructure.Orientation` class
        describing the other crystal orientation from which to compute the
        angle.
    :param crystal_structure: an instance of the `Symmetry` class
        describing the crystal symmetry, triclinic (no symmetry) by
        default.
    :returns tuple: the misorientation angle in radians, the axis as a
        numpy vector (crystal coordinates), the axis as a numpy vector
        (sample coordinates).
    """
    the_angle = np.pi
    the_axis = np.array([0., 0., 1.])
    the_axis_xyz = np.array([0., 0., 1.])
    symmetries = crystal_structure.symmetry_operators()
    (gA, gB) = (ori1, ori2)  # nicknames
    for (g1, g2) in [(gA, gB), (gB, gA)]:
        for j in range(symmetries.shape[0]):
            sym_j = symmetries[j]
            oj = np.dot(sym_j, g1)  # the crystal symmetry operator is left applied
            for i in range(symmetries.shape[0]):
                sym_i = symmetries[i]
                oi = np.dot(sym_i, g2)
                delta = np.dot(oi, oj.T)
                mis_angle = misorientation_angle_from_delta(delta)
                if mis_angle < the_angle:
                    # now compute the misorientation axis, should check if it lies in the fundamental zone
                    mis_axis = misorientation_axis_from_delta(delta)
                    # here we have np.dot(oi.T, mis_axis) = np.dot(oj.T, mis_axis)
                    # print(mis_axis, mis_angle*180/np.pi, np.dot(oj.T, mis_axis))
                    the_angle = mis_angle
                    the_axis = mis_axis
                    the_axis_xyz = np.dot(oi.T, the_axis)
    return the_angle, the_axis, the_axis_xyz


def misorientation_angle_from_delta(delta):
    """Compute the misorientation angle from the misorientation matrix.

    Compute the angle associated with this misorientation matrix :math:`\\Delta g`.
    It is defined as :math:`\\omega = \\arccos(\\text{trace}(\\Delta g)/2-1)`.
    To avoid float rounding point error, the value of :math:`\\cos\\omega`
    is clipped to [-1.0, 1.0].

    .. note::

      This does not account for the crystal symmetries. If you want to
      find the disorientation between two orientations, use the
      :py:meth:`~pymicro.crystal.microstructure.Orientation.disorientation`
      method.

    :param delta: The 3x3 misorientation matrix.
    :returns float: the misorientation angle in radians.
    """
    cw = np.clip(0.5 * (delta.trace() - 1), -1., 1.)
    omega = np.arccos(cw)
    return omega


def misorientation_axis_from_delta(delta):
    """Compute the misorientation axis from the misorientation matrix.

    :param delta: The 3x3 misorientation matrix.
    :returns: the misorientation axis (normalised vector).
    """
    n = np.array([delta[1, 2] - delta[2, 1], delta[2, 0] -
                  delta[0, 2], delta[0, 1] - delta[1, 0]])
    n /= np.sqrt((delta[1, 2] - delta[2, 1]) ** 2 +
                 (delta[2, 0] - delta[0, 2]) ** 2 +
                 (delta[0, 1] - delta[1, 0]) ** 2)
    return n