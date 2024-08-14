import numpy
import numba 


class GrainMap:
    # TODO:
    # Read/write to/from HDF5
    # MT, RMT, unitcell, B, U methods
    # method to calculate IPF colours if we don't have them already
    #     use the phase ID map
    #     we'll have the reference unit cell for each voxel that way
    #     then make a meta orix orien - cached property?
    #     compute the IPF colours for the meta orix orien
    #     fill in the maps
    def __init__(self, ds):
        """Class to store 3D grain map data produced from Scanning 3DXRD datasets"""
        self.maps = {}
        self.ds = ds

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
        return self.maps[map_name]
    
    def __setitem__(self, name, array):
        self.add_map(name, array)
    
    def keys(self):
        return self.maps.keys()
    
    def clear_cache(self):
        # Clear calculated maps
        for name in ("U", "UB"):
            if name in self.keys():
                del self.maps[name]
    
    def add_map(self, name, array):
        """Required to clear caches if UBI is set"""
        if name == "UBI":
            self.clear_cache()
        self.maps[name] = array
    
    def plot(self, map_name, z_layer=0):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.imshow(self.maps[map_name][..., z_layer], origin="lower", interpolation="nearest")
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
            gmap.phases = list(phases)
        except NameError:
            raise AttributeError("Some/all grains are missing reference unit cells! Can't continue")
        
        # work out the phase ID for each grain
        phase_ids = [gmap.phases.index(gs.grain.ref_unitcell) for gs in grainsinos]
        
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

            # ubi_map[g_raw_intensity_mask] = gs.grain.ubi
            
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
            rgb_x_map = np.transpose((redx, grnx, blux), axes=(1, 2, 0, 3))
            rgb_y_map = np.transpose((redy, grny, bluy), axes=(1, 2, 0, 3))
            rgb_z_map = np.transpose((redz, grnz, bluz), axes=(1, 2, 0, 3))
            
            gmap.add_map("ipf_x", rgb_x_map)
            gmap.add_map("ipf_y", rgb_y_map)
            gmap.add_map("ipf_z", rgb_z_map)
        
        return gmap

    
@numba.njit
def nb_inv(mats, imats):
    for i in range(mats.shape[0]):
        for j in range(mats.shape[1]):
            for k in range(mats.shape[2]):
                if np.isnan( mats[i,j,k,0,0] ):
                    imats[i,j,k] = np.nan
                else:
                    try:
                        imats[i,j,k] = np.linalg.inv( mats[i,j,k] )
                    except:
                        imats[i,j,k] = np.nan
