import h5py
import numpy as np
import numba

from ImageD11.sinograms import dataset
from ImageD11.sinograms.sinogram import save_array
from ImageD11 import unitcell


@numba.njit
def nb_choose_best(i, j, u, n, NY, ubiar,
                   minpeaks=6):
    # map of the unique scores
    uniq = np.ones((NY, NY), dtype='q')
    uniq.fill(minpeaks)  # peak cutorr
    npk = np.zeros((NY, NY), dtype='q')
    ubi = np.zeros((NY, NY, 3, 3), dtype='d')
    ubi.fill(np.nan)
    for k in range(i.size):
        ip = i[k]
        jp = j[k]
        #        if ip == 96 and jp == 510:
        #            print(ip,jp,k, ubiar[:,:,k])
        if u[k] > uniq[ip, jp]:
            uniq[ip, jp] = u[k]
            npk[ip, jp] = n[k]
            for ii in range(3):
                for jj in range(3):
                    ubi[ip, jp, ii, jj] = ubiar[ii, jj, k]
    return uniq, npk, ubi


@numba.njit
def nb_inv(mats, imats):
    for i in range(mats.shape[0]):
        for j in range(mats.shape[1]):
            for k in range(mats.shape[2]):
                if np.isnan(mats[i, j, k, 0, 0]):
                    imats[i, j, k] = np.nan
                else:
                    try:
                        imats[i, j, k] = np.linalg.inv(mats[i, j, k])
                    except:
                        imats[i, j, k] = np.nan


class pbpmap:
    """Class to load and manipulate point-by-point indexing results"""

    def __init__(self, fname):
        pbp_array = np.loadtxt(fname).T
        n = len(pbp_array[0])
        self.i = pbp_array[0].astype(int)
        self.i -= self.i.min()
        self.j = pbp_array[1].astype(int)
        self.j -= self.j.min()
        self.n = pbp_array[2].astype(int)  # total peaks indexing with hkl==int with 0.03
        self.u = pbp_array[3].astype(int)  # unique (h,k,l) labels on indexed peaks
        self.NI = int(self.i.max() - self.i.min()) + 1
        self.NJ = int(self.j.max() - self.j.min()) + 1
        self.NY = max(self.NI, self.NJ)
        self.ubi = pbp_array[4:].astype(float)
        self.ubi.shape = 3, 3, -1

    def choose_best(self, minpeaks=6):
        self.muniq, self.npks, self.ubibest = nb_choose_best(
            self.i, self.j,
            self.u, self.n,
            self.NY, self.ubi, minpeaks)


class GrainMap:
    # TODO:
    # MT, RMT, unitcell, B, U methods

    def __init__(self, ds):
        """Class to store 3D grain map data produced from Scanning 3DXRD datasets"""
        self.maps = {}
        self.ds = ds

        # dict to store the ImageD11.unitcell.unitcell objects for each phase ID
        self.phases = dict()

        # dict to store the meta orix orientations for each phase
        self._meta_orix_oriens = dict()

    def __getattribute__(self, item):
        """lets you call gmap.UBI for example
        called whenever getattr(item, 'attr') or item.attr is called
        for speed, look in self.maps FIST before giving up"""
        if item == 'maps':
            return object.__getattribute__(self, item)
        if item in self.maps.keys():
            return super(GrainMap, self).__getattribute__('maps')[item]
        return super(GrainMap, self).__getattribute__(item)

    def __getitem__(self, map_name):
        """allows you to use dictionary syntax like gmap["UBI"]"""
        return self.maps[map_name]

    def __setitem__(self, name, array):
        """allows you to use dictionary syntax like gmap["UBI"]"""
        self.add_map(name, array)

    def keys(self):
        """allows you to use dictionary syntax"""
        return self.maps.keys()

    def clear_cache(self):
        # Clear calculated maps
        for name in ("U", "UB"):
            if name in self.keys():
                del self.maps[name]

    def add_map(self, name, array, suppress_clear=False):
        """Required to clear caches if UBI is set"""
        if name == "UBI" and not suppress_clear:
            self.clear_cache()
        self.maps[name] = array

    def plot(self, map_name, z_layer=0):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        map_to_plot = self.maps[map_name]
        if len(map_to_plot.shape) == 4:
            # 4D so it's an IPF
            map_to_plot = map_to_plot[..., z_layer, :]
        else:
            map_to_plot = map_to_plot[..., z_layer]
        ax.imshow(map_to_plot, origin="lower", interpolation="nearest")
        ax.set_title(map_name)
        ax.set_xlabel('<--- Lab Y')
        ax.set_ylabel('Lab X')
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
            nb_inv(UBI_map, UB_map)
            self.add_map('UB', UB_map)
            return UB_map

    def get_meta_orix_orien(self, phase_id=0):
        """Get a meta orix orientation for all voxels of a given phase ID"""
        if phase_id in self._meta_orix_oriens.keys():
            return self._meta_orix_oriens[phase_id]
        else:
            # calculate it
            if 'phase' not in self.keys():
                raise ValueError('No phase map to select from!')

            if phase_id not in self.phases:
                raise KeyError('Phase ID ' + phase_id + ' not in self.phases!')

            # Get a mask for UBI from this phase ID
            phase_mask = self.phase == phase_id

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
        if 'phase' not in self.keys():
            raise ValueError('No phase map to select from!')

        try:
            from orix.vector.vector3d import Vector3d
        except ImportError:
            raise ImportError("Missing diffpy and/or orix, can't compute orix phase!")

        shape = self.phase.shape

        # iterate over IPF directions and xyz strings
        for axis, letter in zip(np.eye(3), ["x", "y", "z"]):
            rgb_map = np.zeros(shape + (3,))

            if len(self.phases) == 0:
                raise KeyError("No phases in self.phases to compute IPF colours for!")

            rgb_map.shape = -1, 3

            # iterate over phases
            for phase_id in self.phases.keys():
                phase_mask = self.phase == phase_id
                inds = np.mgrid[0:shape[0] * shape[1] * shape[2]][phase_mask.ravel()]
                ipf_direction = Vector3d(axis)
                rgb_flat = self.phases[phase_id].get_ipf_colour_from_orix_orien(self.get_meta_orix_orien(phase_id),
                                                                                axis=ipf_direction)
                rgb_map[inds] = rgb_flat

            rgb_map.shape = shape + (3,)

            self.add_map('ipf_' + letter, rgb_map)

    def to_h5(self, h5file, h5group="GrainMap"):
        """Write all maps to an HDF5 file (h5file) with a parent group h5group. Creates h5group if doesn't already exist"""
        with h5py.File(h5file, "a") as hout:
            parent_group = hout.require_group(h5group)
            
            maps_group = parent_group.require_group('maps')
            
            for map_name in self.keys():
                saved_array = save_array(maps_group, map_name, self.maps[map_name])
                if "ipf" in map_name:
                    # set 'IMAGE" attribute
                    saved_array.attrs['CLASS'] = 'IMAGE'
            
            # store the dataset path
            maps_group.attrs['dsetfile'] = self.ds.dsfile
            
            # store the phases
            if len(self.phases) > 0:
                phase_group = parent_group.require_group('phases')
                for phase_id in self.phases.keys():
                    phase_group[str(phase_id)] = self.phases[phase_id].tostring()

            

    @classmethod
    def from_h5(cls, h5file, h5group="GrainMap", ds=None):
        """Load map object from an HDF5 file"""
        with h5py.File(h5file, 'r') as hin:
            parent_group = hin[h5group]
            
            maps_group = parent_group['maps']

            if ds is None:  # load the DataSet from disk if none provided
                dsfile = maps_group.attrs['dsetfile']
                ds = dataset.load(dsfile)

            gmap = GrainMap(ds=ds)
            for map_name in maps_group.keys():
                loaded_array = maps_group[map_name][:]  # load the array from disk
                # add array to the class, suppressing cache clearing (so we don't have to recompute UB if it's saved to disk)
                gmap.add_map(map_name, loaded_array, suppress_clear=True)
            
            if 'phases' in parent_group.keys():
                phase_group = parent_group['phases']
                
                for phase_id in phase_group.keys():
                    gmap.phases[int(phase_id)] = unitcell.cellfromstring(phase_group[phase_id][()].decode('utf-8'))

        return gmap

    @classmethod
    def from_pbpmap(cls, pbpmap, ds):
        """From a pbpmap object and a dataset"""

        # make an empty GrainMap object
        gmap = GrainMap(ds=ds)

        # see if we have a ubibest to take
        if hasattr(pbpmap, 'ubibest'):
            ubi_map = pbpmap.ubibest
        else:
            ubi_map = pbpmap.ubi

        # reshape
        ubi_map = ubi_map.reshape(ubi_map.shape[:2] + (1, 3, 3))
        gmap.add_map("UBI", ubi_map)

        return gmap

    @classmethod
    def from_grainsinos(cls, grainsinos, method="iradon", use_gids=True, cutoff_level=0.1):
        """Build a GrainMap from a list of GrainSinos.
        method is the recon that we look for inside each grainsino
        use_gids will look for grainsino.grain.gid inside each grainsino to use as the label
        if it can't find it, it will use the increment"""
        # get the dataset
        ds = grainsinos[0].ds
        # make an empty GrainMap object
        gmap = GrainMap(ds=ds)

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
            phases = {gs.grain.ref_unitcell for gs in grainsinos}
            phases_list = list(phases)
            gmap.phases = {inc: phase for inc, phase in enumerate(phases_list)}
        except NameError:
            raise AttributeError("Some/all grains are missing reference unit cells! Can't continue")

        # work out the phase ID for each grain
        phase_ids = [phases_list.index(gs.grain.ref_unitcell) for gs in grainsinos]

        # get the shape of the maps
        map_shape = grainsinos[0].recons[method].shape + (1,)

        # make an empty grain label map
        grain_labels_map = np.full(map_shape, -1)

        # make an empty intensity map
        raw_intensity_map = np.zeros(map_shape)

        # make an empty phase ID map
        phase_id_map = np.full(map_shape, -1)

        # make an empty UBI map
        ubi_map = np.full((map_shape[:2] + (1, 3, 3)), 0.0)

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
            g_raw_intensity = norm(gs.recons[method])[..., np.newaxis]
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

        gmap.add_map("intensity", raw_intensity_map)
        gmap.add_map("labels", grain_labels_map)
        gmap.add_map("phase", phase_id_map)
        gmap.add_map("UBI", ubi_map)

        if have_ipfs:
            rgb_x_map = np.transpose((redx, grnx, blux), axes=(1, 2, 3, 0))
            rgb_y_map = np.transpose((redy, grny, bluy), axes=(1, 2, 3, 0))
            rgb_z_map = np.transpose((redz, grnz, bluz), axes=(1, 2, 3, 0))

            gmap.add_map("ipf_x", rgb_x_map)
            gmap.add_map("ipf_y", rgb_y_map)
            gmap.add_map("ipf_z", rgb_z_map)

        return gmap
