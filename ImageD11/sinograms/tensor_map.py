import os

import h5py
import numpy as np
import numba

from ImageD11.sinograms import dataset
from ImageD11.sinograms.sinogram import save_array
from ImageD11 import unitcell
from ImageD11.sinograms.point_by_point import nb_inv_3d


class TensorMap:
    """This is a class to store a contiguous voxel-based representation of a sample.
    At its core is the self.maps attribute, which is a dictionary of Numpy arrays.
    Each Numpy array should represent a 3D voxel grid of the sample.
    The dimensions of the array are aligned with the laboratory reference frame in the order (Z, Y, X, ...)
    The shape of the first three axes should therefore be (NZ, NY, NX, ...)
    The total number of dimensions can vary
    E.g a phase ID map might be (1, 15, 20) but a UBI map might be (1, 15, 20, 3, 3)
    """
    
    # TODO:
    # MT, RMT, unitcell, B, U methods

    def __init__(self, maps=None, phases=None, steps=None):
        """maps: dict of Numpy arrays, each with shape (NZ, NY, NX, ...), with string keys for the map names
           phases: dict of ImageD11.unitcell.unitcell objects with integer keys for the phase IDs
           steps: [zstep, ystep, xtep] step sizes in um of voxels"""
        
        # Initialise an empty dictionary of voxel maps
        if maps is None:
            maps = dict()
        self.maps = maps
        
        # dict to store the ImageD11.unitcell.unitcell objects for each phase ID
        if phases is None:
            phases = dict()
        self.phases = phases
        
        if steps is None:
            steps = [1.0, 1.0, 1.0]
        self.steps = steps

        # dict to store the meta orix orientations for each phase ID
        self._meta_orix_oriens = dict()
        
        # dict to store grain merges
        # e.g when we merge together multiple TensorMap layers
        # if we want to detect and merge duplicate grains in multiple layers
        # then we will have new merged grain ids in self.labels
        # we need to store the mapping from new merged grain ids to original grain ids
        self.merged_mapping = dict()
        
        self._shape = None
        
        # Ensure that all the maps have the same shape (will also fill in self._shape)
        self.check_shape()

    def __getattribute__(self, item):
        """lets you call self.UBI for example
        called whenever getattr(item, 'attr') or item.attr is called
        for speed, look in self.maps FIST before giving up"""
        if item == 'maps':
            return object.__getattribute__(self, item)
        if item in self.maps.keys():
            return super(TensorMap, self).__getattribute__('maps')[item]
        return super(TensorMap, self).__getattribute__(item)

    def __getitem__(self, map_name):
        """allows you to use dictionary syntax like self["UBI"]"""
        return self.maps[map_name]

    def __setitem__(self, name, array):
        """allows you to use dictionary syntax like self["UBI"]"""
        self.add_map(name, array)

    def keys(self):
        """allows you to use dictionary syntax"""
        return self.maps.keys()

    def clear_cache(self):
        # Clear calculated maps
        for name in ("U", "UB"):
            if name in self.keys():
                del self.maps[name]
    
    def check_shape(self):
        """Checks that all the maps in self.maps have equal shape for their first 3 dimensions"""
        
        if len(self.maps) > 0:
            map_dims = [array.shape[:3] for array in self.maps.values()]
            map_dims_set = set(map_dims)
            if len(map_dims_set) > 1:
                raise ValueError("Not all the maps in self.maps have the right shape!")
            self._shape = list(map_dims_set)[0]
    
    @property
    def shape(self):
        """The shape of the map (NZ, NY, NX)"""
        if self._shape is None:
            self.check_shape()
        return self._shape

    def add_map(self, name, array):
        """Required to clear caches if UBI is set"""
        if name == "UBI":
            self.clear_cache()
        self.maps[name] = array
    
    def plot(self, map_name, z_layer=0):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.imshow(self.maps[map_name][z_layer, ...], origin="lower")
        ax.set_xlabel('Lab X axis --->')
        ax.set_ylabel('Lab Y axis --->')
        ax.set_title(map_name)
        plt.show()

    @property
    def UB(self):
        """The UB matrix"""
        if "UB" in self.keys():
            return self.maps['UB']
        else:
            if 'UBI' not in self.keys():
                raise KeyError('No UBI to calculate from!')
            # calculate it
            UB_map = np.zeros_like(self.maps['UBI'])
            UBI_map = self.maps['UBI']
            nb_inv_3d(UBI_map, UB_map)
            self.add_map('UB', UB_map)
            return UB_map

    def get_meta_orix_orien(self, phase_id=0):
        """Get a meta orix orientation for all voxels of a given phase ID"""
        if phase_id in self._meta_orix_oriens.keys():
            return self._meta_orix_oriens[phase_id]
        else:
            # calculate it
            if 'phase_ids' not in self.keys():
                raise ValueError('No phase map to select from!')

            if phase_id not in self.phases:
                raise KeyError('Phase ID ' + phase_id + ' not in self.phases!')

            # Get a mask for UBI from this phase ID
            phase_mask = self.phase_ids == phase_id

            # Get a mask for non-nan UBs
            non_nan_mask = ~np.isnan(self.UB[:, :, :, 0, 0])

            # Combine masks
            total_mask = phase_mask & non_nan_mask
            UBcalc = self.UB[total_mask]

            # Get the reference orientation for this phase
            ref_unitcell = self.phases[phase_id]

            # Calculate meta orien
            meta_orien = ref_unitcell.get_orix_orien(UBcalc)

            self._meta_orix_oriens[phase_id] = meta_orien
            return self._meta_orix_oriens[phase_id]

    def get_ipf_maps(self):
        """Calculate all IPF maps and add them to self.maps"""
        if 'phase_ids' not in self.keys():
            raise ValueError('No phase map to select from!')

        try:
            from orix.vector.vector3d import Vector3d
        except ImportError:
            raise ImportError("Missing diffpy and/or orix, can't compute orix phase!")

        shape = self.phase_ids.shape

        # iterate over IPF directions and xyz strings
        for axis, letter in zip(np.eye(3), ["x", "y", "z"]):
            rgb_map = np.zeros(shape + (3,))

            if len(self.phases) == 0:
                raise KeyError("No phases in self.phases to compute IPF colours for!")

            rgb_map.shape = -1, 3

            # iterate over phases
            for phase_id in self.phases.keys():
                phase_mask = self.phase_ids == phase_id
                inds = np.mgrid[0:shape[0] * shape[1] * shape[2]][phase_mask.ravel()]
                ipf_direction = Vector3d(axis)
                rgb_flat = self.phases[phase_id].get_ipf_colour_from_orix_orien(self.get_meta_orix_orien(phase_id),
                                                                                axis=ipf_direction)
                rgb_map[inds] = rgb_flat

            rgb_map.shape = shape + (3,)

            self.add_map('ipf_' + letter, rgb_map)

    def to_h5(self, h5file, h5group='TensorMap'):
        """Write all maps to an HDF5 file (h5file) with a parent group h5group. Creates h5group if doesn't already exist"""
        with h5py.File(h5file, "a") as hout:
            parent_group = hout.require_group(h5group)
            
            maps_group = parent_group.require_group('maps')
            
            for map_name in self.keys():
                array = self.maps[map_name]
                
                # save array to H5 file in reverse order (Z Y X rather than X Y Z) because XDMF importer requires it)
                
                # array = np.moveaxis(array, (0,1,2), (2,1,0))
                saved_array = save_array(maps_group, map_name, array)
                if "ipf" in map_name:
                    # set 'IMAGE" attribute
                    saved_array.attrs['CLASS'] = 'IMAGE'
            
            # store the phases
            if len(self.phases) > 0:
                phase_group = parent_group.require_group('phases')
                for phase_id in self.phases.keys():
                    phase_group[str(phase_id)] = self.phases[phase_id].tostring()
            
            # store the step sizes
            parent_group.create_dataset("step", data=np.array(self.steps))

    
    def to_paraview(self, h5name, h5group='TensorMap'):
        """Exports to H5, then writes an XDMF file that lets you read the data with ParaView"""
        # Write H5 first
        if not os.path.exists(h5name):
            self.to_h5(h5name, h5group=h5group)
        
        h5_relpath = os.path.split(h5name)[1]
        
        xdmf_filename = h5name.replace('.h5','.xdmf')
        
        dims = self.shape
        scalar_dims = dims
        vector_dims = dims + (3,)
        tensor_dims = dims + (3,3,)

        MeshDimensions = (dims[0] + 1, dims[1] + 1, dims[2] + 1)
        
        MeshDimensionsStr = 'Dimensions="%d %d %d"' % MeshDimensions
        ScalarDimensionsStr = 'Dimensions="%d %d %d"' % scalar_dims
        VectorDimensionsStr = 'Dimensions="%d %d %d %d"' % vector_dims
        TensorDimensionsStr = 'Dimensions="%d %d %d %d %d"' % tensor_dims
        
        steps = self.steps

        # Write XDMF file
        with open(xdmf_filename, 'wt') as fileID:
            fileID.write('<?xml version="1.0"?>\n')
            fileID.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd"[]>\n')
            fileID.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
            fileID.write(' <Domain>\n')
            fileID.write('  <Grid Name="GM3D" GridType="Uniform">\n')
            fileID.write('   <Topology TopologyType="3DCoRectMesh" %s></Topology>\n' % MeshDimensionsStr)
            fileID.write('    <Geometry Type="ORIGIN_DXDYDZ">\n')
            fileID.write('     <DataItem Format="XML" Dimensions="3">0 0 0</DataItem>\n')
            fileID.write('     <DataItem Format="XML" Dimensions="3">%.6f %.6f %.6f</DataItem>\n' % steps)
            fileID.write('    </Geometry>\n')
            
            # iterate over all the maps
            for map_name in self.keys():
                array = self.maps[map_name]
                
                # work out what sort of array we have
                map_shape = array.shape
                n_dims = len(map_shape)
                if n_dims == 3:
                    # scalar field
                    fileID.write('    <Attribute Name="%s" AttributeType="Scalar" Center="Cell">\n' % map_name)
                    fileID.write('      <DataItem Format="HDF" %s NumberType="Float" Precision="6" >%s:/%s</DataItem>\n' % (ScalarDimensionsStr, h5_relpath, h5group + '/maps/' + map_name))
                    fileID.write('    </Attribute>\n')
                elif n_dims == 4:
                    # vector field (like IPF)
                    fileID.write('    <Attribute Name="%s" AttributeType="Vector" Center="Cell">\n' % map_name)
                    fileID.write('      <DataItem Format="HDF" %s NumberType="Float" Precision="6" >%s:/%s</DataItem>\n' % (VectorDimensionsStr, h5_relpath, h5group + '/maps/' + map_name))
                    fileID.write('    </Attribute>\n')
                elif n_dims == 5:
                    # tensor field (like UBI)
                    # fileID.write('    <Attribute Name="%s" AttributeType="Tensor" Center="Cell">\n' % map_name)
                    # fileID.write('      <DataItem Format="HDF" %s NumberType="Float" Precision="6" >%s:/%s</DataItem>\n' % (TensorDimensionsStr, h5_relpath, h5group + '/maps/' + map_name))
                    # fileID.write('    </Attribute>\n')
                    continue
                else:
                    continue

            fileID.write('  </Grid>\n')
            fileID.write(' </Domain>\n')
            fileID.write('</Xdmf>\n')
        
    @staticmethod
    def recon_order_to_map_order(recon_arr):
        """Transform a 2D array from reconstruction space (first axis X, second axis is -Y)
           to TensorMap space (Z, Y ,X)
           The input array must have the first two dimensions (X, -Y) as is the case
           for reconstructed grain shapes from sinogram or PBP methods"""
        # we are currently in (X, -Y, ...)
        nx, ny = recon_arr.shape[:2]
        
        # flip the second axis
        micro_arr = np.flip(recon_arr, 1)

        # now we are (X, Y, ...)
        # add the Z axis in the front
        micro_arr = np.expand_dims(micro_arr, 0)

        # now we are (Z, X, Y, ...)
        # now swap X and Y
        micro_arr = np.swapaxes(micro_arr, 1, 2)

        # now we are (Z, Y, X)
        
        assert micro_arr.shape[:3] == (1, ny, nx)
        
        return micro_arr
    
    @staticmethod
    def map_index_to_recon(mj, mk, yshape):
        """From a 3D array index (mi, mj, mk) and the Y shape of the map,
        determine the corresponding position in reconstruction space"""
        return (mk, yshape-mj-1)
    
    @classmethod
    def from_ubis(cls, ubi_array):
        """Make simplest possible TensorMap object from a UBI array in reconstuction space (X, -Y)"""
        
        ubi_array = cls.recon_order_to_map_order(ubi_array)
        
        # just make a simple maps container with the UBI array
        maps = {'UBI': ubi_array}
        return TensorMap(maps=maps)

    @classmethod
    def from_h5(cls, h5file, h5group='TensorMap'):
        """Load TensorMap object from an HDF5 file"""
        
        maps = dict()
        phases = dict()
        
        with h5py.File(h5file, 'r') as hin:
            parent_group = hin[h5group]
            
            maps_group = parent_group['maps']

            for map_name in maps_group.keys():
                array = maps_group[map_name][:]  # load the array from disk

                maps[map_name] = array
            
            if 'phases' in parent_group.keys():
                phase_group = parent_group['phases']
                
                for phase_id in phase_group.keys():
                    phases[int(phase_id)] = unitcell.cellfromstring(phase_group[phase_id][()].decode('utf-8'))
            
            steps = parent_group['step'][:]
        
        tensor_map = cls(maps=maps, phases=phases, steps=steps)
        
        return tensor_map

    @classmethod
    def from_pbpmap(cls, pbpmap, steps=None):
        """Create TensorMap from a pbpmap object"""
        
        maps = dict()
        
        # see if we have a ubibest to take
        if hasattr(pbpmap, 'ubibest'):
            ubi_map = pbpmap.ubibest
        else:
            ubi_map = pbpmap.ubi
        
        # reshape ubi map and add it to the dict
        maps['UBI'] = cls.recon_order_to_map_order(ubi_map)
        
        # add npks to the dict
        if hasattr(pbpmap, 'npks'):
            maps['npks'] = cls.recon_order_to_map_order(pbpmap.npks)
        
        # add nuniq to the dict
        if hasattr(pbpmap, 'nuniq'):
            maps['nuniq'] = cls.recon_order_to_map_order(pbpmap.nuniq)
                
        tensor_map = cls(maps=maps, steps=steps)

        return tensor_map

    @classmethod
    def from_grainsinos(cls, grainsinos, method="iradon", use_gids=True, cutoff_level=0.1, steps=None):
        """Build a TensorMap object from a list of GrainSinos.
        method is the recon that we look for inside each grainsino
        use_gids will look for grainsino.grain.gid inside each grainsino to use as the label
        if it can't find it, it will use the increment"""
        
        # make empty maps container
        maps = dict()
        
        # make empty phases container
        phases = dict()

        # get the labels for each grain
        grain_labels = [inc for inc, _ in enumerate(grainsinos)]
        if use_gids:
            try:
                grain_labels = [gs.grain.gid for gs in grainsinos]
            except AttributeError:
                print("Some/all GIDs are missing! Using increments instead:")

        # check that each grain has a reference unitcell
        # if it doesn't, 

        # work out the phase for each grain to decide how many phases we have
        # use a set to remove duplicate phases
        # if this fails, we can't continue
        try:
            phase_set = {gs.grain.ref_unitcell for gs in grainsinos}
            phases_list = list(phase_set)
            
            # add the phases we found to the phases dictionary
            for phase_inc, phase in enumerate(phases_list):
                phases[phase_inc] = phase
        except NameError:
            raise AttributeError("Some/all grains are missing reference unit cells! Can't continue")

        # work out the phase ID for each grain
        phase_ids = [phases_list.index(gs.grain.ref_unitcell) for gs in grainsinos]
        
        # construct the maps in reconstruction space
        # we will convert them all at the end before adding
        
        map_shape = grainsinos[0].recons[method].shape

        # make an empty grain label map
        grain_labels_map = np.full(map_shape, -1)

        # make an empty intensity map
        raw_intensity_map = np.zeros(map_shape)

        # make an empty phase ID map
        phase_id_map = np.full(map_shape, -1)

        # make an empty UBI map
        ubi_map = np.full((map_shape[:2] + (3, 3)), 0.0)

        # check if we have IPF colours

        have_ipfs = False
        if all([all([hasattr(gs.grain, attr) for attr in ("rgb_x", "rgb_y", "rgb_z")]) for gs in grainsinos]):
            have_ipfs = True

        # IPF maps
        if have_ipfs:
            redx = np.zeros(map_shape)
            grnx = np.zeros(map_shape)
            blux = np.zeros(map_shape)

            redy = np.zeros(map_shape)
            grny = np.zeros(map_shape)
            bluy = np.zeros(map_shape)

            redz = np.zeros(map_shape)
            grnz = np.zeros(map_shape)
            bluz = np.zeros(map_shape)

        # normalisation function
        def norm(r):
            m = r > r.max() * 0.2
            return (r / r[m].mean()).clip(0, 1)

        for label, gs, phase_id in zip(grain_labels, grainsinos, phase_ids):
            g_raw_intensity = norm(gs.recons[method])
            g_raw_intensity_mask = g_raw_intensity > raw_intensity_map
            g_raw_intensity_map = g_raw_intensity[g_raw_intensity_mask]
            raw_intensity_map[g_raw_intensity_mask] = g_raw_intensity_map
            grain_labels_map[g_raw_intensity_mask] = label
            phase_id_map[g_raw_intensity_mask] = phase_id
            ubi_map[g_raw_intensity_mask] = gs.grain.ubi

            if have_ipfs:
                redx[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_x[0]
                grnx[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_x[1]
                blux[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_x[2]

                redy[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_y[0]
                grny[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_y[1]
                bluy[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_y[2]

                redz[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_z[0]
                grnz[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_z[1]
                bluz[g_raw_intensity_mask] = g_raw_intensity_map * gs.grain.rgb_z[2]

        raw_intensity_map[raw_intensity_map <= cutoff_level] = 0.0
        
        maps["intensity"] = cls.recon_order_to_map_order(raw_intensity_map)
        maps["labels"] = cls.recon_order_to_map_order(grain_labels_map)
        maps["phase_ids"] = cls.recon_order_to_map_order(phase_id_map)
        maps["UBI"] = cls.recon_order_to_map_order(ubi_map)

        if have_ipfs:
            rgb_x_map = np.transpose((redx, grnx, blux), axes=(1, 2, 0))
            rgb_y_map = np.transpose((redy, grny, bluy), axes=(1, 2, 0))
            rgb_z_map = np.transpose((redz, grnz, bluz), axes=(1, 2, 0))
            
            maps["ipf_x"] = cls.recon_order_to_map_order(rgb_x_map)
            maps["ipf_y"] = cls.recon_order_to_map_order(rgb_y_map)
            maps["ipf_z"] = cls.recon_order_to_map_order(rgb_z_map)
        
        # get the step size from the dataset of one of the grainsino objects
        if steps is None:
            ystep = grainsinos[0].ds.ystep
            steps = (1.0, ystep, ystep)
        
        tensor_map = cls(maps=maps, phases=phases, steps=steps)
        
        return tensor_map