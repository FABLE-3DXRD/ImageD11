import os

import numba
import h5py
import numpy as np

from ImageD11.sinograms.geometry import recon_to_step
from ImageD11.sinograms.sinogram import save_array
from ImageD11 import unitcell


# from ImageD11.sinograms.point_by_point import nb_inv_3d


# Numba vectorized functions to work on maps of any dimensionality, provided the last dimensions are the ones of the array
# e.g these can be called and broadcasted for a UBI map of (1,200,500,3,3)

# take in a nxn float64 array, return a nxn float64 array
@numba.guvectorize([(numba.float64[:, :], numba.float64[:, :])], '(n,n)->(n,n)', nopython=True)
def fast_invert(mat, res):
    """Fast broadcastable matrix invert.
    To use (no need to pre-populate res):
    imat = fast_invert(mat)"""
    if np.isnan(mat[0, 0]):
        res[...] = np.nan
    else:
        res[...] = np.linalg.inv(mat)


# take in a nxn float64 array, return a nxn float64 array
@numba.guvectorize([(numba.float64[:, :], numba.float64[:, :])], '(n,n)->(n,n)', nopython=True)
def ubi_to_mt(ubi, res):
    """Get metric tensor from UBI, the same way as in ImageD11.grain
    To use (no need to pre-populate res):
    mt = ubi_to_mt(ubi)"""
    if np.isnan(ubi[0, 0]):
        res[...] = np.nan
    else:
        res[...] = np.dot(ubi, ubi.T)


# take in a nxn float64 array and a dummy variable p-length float64 vector (required by numba), return a p-length float64 vector
@numba.guvectorize([(numba.float64[:, :], numba.float64[:], numba.float64[:])], '(n,n),(p)->(p)', nopython=True)
def mt_to_unitcell(mt, dum, res):
    """Get unitcell from mt, the same way as in ImageD11.grain
    To use (no need to pre-populate res):
    dummy_var = np.arange(6)  # ensure dummy variable has the right length
    unitcell = mt_to_unitcell(mt, dummy_var)"""
    if np.isnan(mt[0, 0]):
        res[...] = np.nan
    else:
        G = mt
        a, b, c = np.sqrt(np.diag(G))
        al = np.degrees(np.arccos(G[1, 2] / b / c))
        be = np.degrees(np.arccos(G[0, 2] / a / c))
        ga = np.degrees(np.arccos(G[0, 1] / a / b))
        res[..., 0] = a
        res[..., 1] = b
        res[..., 2] = c
        res[..., 3] = al
        res[..., 4] = be
        res[..., 5] = ga


# take in a p-length float64 vector and a dummy variable nxn float64 array (required by numba), return a nxn float64 array
@numba.guvectorize([(numba.float64[:], numba.float64[:, :], numba.float64[:, :])], '(p),(n,n)->(n,n)', nopython=True)
def unitcell_to_b(unitcell, dum, res):
    """Get B matrix from unitcell, the same way as in ImageD11.grain (via ImageD11.unitcell)
    To use (no need to pre-populate res):
    dummy_var = np.eye(3)  # ensure dummy variable has the right shape (3x3)
    B = unitcell_to_b(unitcell, dummy_var)"""
    if np.isnan(unitcell[0]):
        res[...] = np.nan
    else:
        a, b, c = unitcell[:3]
        ralpha, rbeta, rgamma = np.radians(unitcell[3:])  # radians
        ca = np.cos(ralpha)
        cb = np.cos(rbeta)
        cg = np.cos(rgamma)
        g = np.full((3, 3), np.nan, float)
        g[0, 0] = a * a
        g[0, 1] = a * b * cg
        g[0, 2] = a * c * cb
        g[1, 0] = a * b * cg
        g[1, 1] = b * b
        g[1, 2] = b * c * ca
        g[2, 0] = a * c * cb
        g[2, 1] = b * c * ca
        g[2, 2] = c * c
        gi = np.linalg.inv(g)
        astar, bstar, cstar = np.sqrt(np.diag(gi))
        betas = np.degrees(np.arccos(gi[0, 2] / astar / cstar))
        gammas = np.degrees(np.arccos(gi[0, 1] / astar / bstar))

        res[..., 0, 0] = astar
        res[..., 0, 1] = bstar * np.cos(np.radians(gammas))
        res[..., 0, 2] = cstar * np.cos(np.radians(betas))
        res[..., 1, 0] = 0.0
        res[..., 1, 1] = bstar * np.sin(np.radians(gammas))
        res[..., 1, 2] = -cstar * np.sin(np.radians(betas)) * ca
        res[..., 2, 0] = 0.0
        res[..., 2, 1] = 0.0
        res[..., 2, 2] = 1.0 / c


# take in a nxn float64 array, return a nxn float64 array
@numba.guvectorize([(numba.float64[:, :], numba.float64[:, :], numba.float64[:, :])], '(n,n),(n,n)->(n,n)',
                   nopython=True)
def ubi_and_b_to_u(ubi, b, res):
    """Get U matrix from UBI and B matrices, the same way as in ImageD11.grain
    To use (no need to pre-populate res):
    U = ubi_and_b_to_u(ubi, b)"""
    if np.isnan(ubi[0, 0]):
        res[...] = np.nan
    elif np.isnan(b[0, 0]):
        res[...] = np.nan
    else:
        res[...] = np.dot(b, ubi).T


@numba.guvectorize([(numba.float64[:, :], numba.float64[:], numba.float64[:, :])], '(n,n),(p)->(n,n)', nopython=True)
def ubi_and_unitcell_to_eps_sample(ubi, dzero_cell, res):
    if np.isnan(ubi[0, 0]):
        res[...] = np.nan
    elif np.isnan(dzero_cell[0]):
        res[...] = np.nan
    else:
        a, b, c = dzero_cell[:3]
        ralpha, rbeta, rgamma = np.radians(dzero_cell[3:])  # radians
        ca = np.cos(ralpha)
        cb = np.cos(rbeta)
        cg = np.cos(rgamma)
        g = np.full((3, 3), np.nan, float)
        g[0, 0] = a * a
        g[0, 1] = a * b * cg
        g[0, 2] = a * c * cb
        g[1, 0] = a * b * cg
        g[1, 1] = b * b
        g[1, 2] = b * c * ca
        g[2, 0] = a * c * cb
        g[2, 1] = b * c * ca
        g[2, 2] = c * c
        gi = np.linalg.inv(g)
        astar, bstar, cstar = np.sqrt(np.diag(gi))
        betas = np.degrees(np.arccos(gi[0, 2] / astar / cstar))
        gammas = np.degrees(np.arccos(gi[0, 1] / astar / bstar))

        B = np.full((3, 3), np.nan, float)
        B[0, 0] = astar
        B[0, 1] = bstar * np.cos(np.radians(gammas))
        B[0, 2] = cstar * np.cos(np.radians(betas))
        B[1, 0] = 0.0
        B[1, 1] = bstar * np.sin(np.radians(gammas))
        B[1, 2] = -cstar * np.sin(np.radians(betas)) * ca
        B[2, 0] = 0.0
        B[2, 1] = 0.0
        B[2, 2] = 1.0 / c

        F = np.dot(ubi.T, B.T)
        w, sing, vh = np.linalg.svd(F)
        V = np.dot(w, np.dot(np.diag(sing), w.T))
        em = V - np.eye(3)
        res[...] = em


@numba.guvectorize([(numba.float64[:, :], numba.float64[:], numba.float64[:, :])], '(n,n),(p)->(n,n)', nopython=True)
def ubi_and_unitcell_to_eps_crystal(ubi, dzero_cell, res):
    if np.isnan(ubi[0, 0]):
        res[...] = np.nan
    elif np.isnan(dzero_cell[0]):
        res[...] = np.nan
    else:
        a, b, c = dzero_cell[:3]
        ralpha, rbeta, rgamma = np.radians(dzero_cell[3:])  # radians
        ca = np.cos(ralpha)
        cb = np.cos(rbeta)
        cg = np.cos(rgamma)
        g = np.full((3, 3), np.nan, float)
        g[0, 0] = a * a
        g[0, 1] = a * b * cg
        g[0, 2] = a * c * cb
        g[1, 0] = a * b * cg
        g[1, 1] = b * b
        g[1, 2] = b * c * ca
        g[2, 0] = a * c * cb
        g[2, 1] = b * c * ca
        g[2, 2] = c * c
        gi = np.linalg.inv(g)
        astar, bstar, cstar = np.sqrt(np.diag(gi))
        betas = np.degrees(np.arccos(gi[0, 2] / astar / cstar))
        gammas = np.degrees(np.arccos(gi[0, 1] / astar / bstar))

        B = np.full((3, 3), np.nan, float)
        B[0, 0] = astar
        B[0, 1] = bstar * np.cos(np.radians(gammas))
        B[0, 2] = cstar * np.cos(np.radians(betas))
        B[1, 0] = 0.0
        B[1, 1] = bstar * np.sin(np.radians(gammas))
        B[1, 2] = -cstar * np.sin(np.radians(betas)) * ca
        B[2, 0] = 0.0
        B[2, 1] = 0.0
        B[2, 2] = 1.0 / c

        F = np.dot(ubi.T, B.T)
        w, sing, vh = np.linalg.svd(F)
        S = np.dot(vh.T, np.dot(np.diag(sing), vh))
        em = S - np.eye(3)
        res[...] = em


@numba.guvectorize([(numba.float64[:, :], numba.float64[:, :], numba.float64[:], numba.float64[:, :])],
                   '(n,n),(k,k),()->(n,n)', nopython=True)
def strain_to_stress(strain, C, phase_mask, res):
    mvvec = np.zeros(6, dtype=np.float64)
    mvvec[0] = strain[0, 0]
    mvvec[1] = strain[1, 1]
    mvvec[2] = strain[2, 2]
    mvvec[3] = np.sqrt(2.) * strain[1, 2]
    mvvec[4] = np.sqrt(2.) * strain[0, 2]
    mvvec[5] = np.sqrt(2.) * strain[0, 1]

    stress_MV = np.dot(C, mvvec)

    symm = np.full((3, 3), np.nan, dtype=np.float64)
    if phase_mask:
        symm[0, 0] = stress_MV[0]
        symm[1, 1] = stress_MV[1]
        symm[2, 2] = stress_MV[2]
        symm[1, 2] = stress_MV[3] * (1. / np.sqrt(2.))
        symm[0, 2] = stress_MV[4] * (1. / np.sqrt(2.))
        symm[0, 1] = stress_MV[5] * (1. / np.sqrt(2.))
        symm[2, 1] = stress_MV[3] * (1. / np.sqrt(2.))
        symm[2, 0] = stress_MV[4] * (1. / np.sqrt(2.))
        symm[1, 0] = stress_MV[5] * (1. / np.sqrt(2.))

    res[...] = symm


@numba.guvectorize([(numba.float64[:, :], numba.float64[:])], '(n,n)->()', nopython=True)
def sig_to_hydro(sig, res):
    """Get hydrostatic stress scalar from stress tensor (frame invariant)"""
    res[...] = np.sum(np.diag(sig)) / 3


@numba.guvectorize([(numba.float64[:, :], numba.float64[:])], '(n,n)->()', nopython=True)
def sig_to_vm(sig, res):
    """Get von-Mises stress scalar from stress tensor"""
    sig11 = sig[0, 0]
    sig22 = sig[1, 1]
    sig33 = sig[2, 2]
    sig12 = sig[0, 1]
    sig23 = sig[1, 2]
    sig31 = sig[2, 0]
    vm = np.sqrt(((sig11 - sig22) ** 2 + (sig22 - sig33) ** 2 + (sig33 - sig11) ** 2 + 6 * (
                sig12 ** 2 + sig23 ** 2 + sig31 ** 2)) / 2.)
    res[...] = vm


# @numba.njit
# def UB_map_to_U_B_map(UB_map):
#     """QR decomp from xfab.tools"""
#     U_map = np.zeros_like(UB_map)
#     U_map.fill(np.nan)
#     B_map = np.zeros_like(UB_map)
#     B_map.fill(np.nan)
#     for i in range(UB_map.shape[0]):
#         for j in range(UB_map.shape[1]):
#             for k in range(UB_map.shape[2]):
#                 UB = UB_map[i, j, k]
#                 if np.isnan(UB[0, 0]):
#                     U_map[i, j, k] = np.nan
#                     B_map[i, j, k] = np.nan
#                 else:
#                     U, B = np.linalg.qr(UB)
#
#                     if B[0, 0] < 0:
#                         B[0, 0] = -B[0, 0]
#                         B[0, 1] = -B[0, 1]
#                         B[0, 2] = -B[0, 2]
#                         U[0, 0] = -U[0, 0]
#                         U[1, 0] = -U[1, 0]
#                         U[2, 0] = -U[2, 0]
#                     if B[1, 1] < 0:
#                         B[1, 1] = -B[1, 1]
#                         B[1, 2] = -B[1, 2]
#                         U[0, 1] = -U[0, 1]
#                         U[1, 1] = -U[1, 1]
#                         U[2, 1] = -U[2, 1]
#                     if B[2, 2] < 0:
#                         B[2, 2] = -B[2, 2]
#                         U[0, 2] = -U[0, 2]
#                         U[1, 2] = -U[1, 2]
#                         U[2, 2] = -U[2, 2]
#
#                     U_map[i, j, k] = U
#                     B_map[i, j, k] = B
#
#     return U_map, B_map


@numba.njit
def _arctan2(y, x):
    """From xfab.tools"""
    tol = 1e-8
    if np.abs(x) < tol: x = 0
    if np.abs(y) < tol: y = 0

    if x > 0:
        return np.arctan(y / x)
    elif x < 0 <= y:
        return np.arctan(y / x) + np.pi
    elif x < 0 and y < 0:
        return np.arctan(y / x) - np.pi
    elif x == 0 and y > 0:
        return np.pi / 2
    elif x == 0 and y < 0:
        return -np.pi / 2
    elif x == 0 and y == 0:
        raise ValueError('Local function _arctan2() does not accept arguments (0,0)')


# take in a nxn float64 array and a dummy variable p-length float64 vector (required by numba), return a p-length float64 vector
@numba.guvectorize([(numba.float64[:, :], numba.float64[:], numba.float64[:])], '(n,n),(p)->(p)', nopython=True)
def u_to_euler(u, dum, res):
    """Get Euler angles (radians) from U matrix, the same way as in xfab.tools
    To use (no need to pre-populate res):
    dummy = np.arange(3)
    euler = u_to_euler(u, dummy)"""
    if np.isnan(u[0, 0]):
        res[...] = np.nan
    else:
        tol = 1e-8
        PHI = np.arccos(u[2, 2])
        if np.abs(PHI) < tol:
            phi1 = _arctan2(-u[0, 1], u[0, 0])
            phi2 = 0
        elif np.abs(PHI - np.pi) < tol:
            phi1 = _arctan2(u[0, 1], u[0, 0])
            phi2 = 0
        else:
            phi1 = _arctan2(u[0, 2], -u[1, 2])
            phi2 = _arctan2(u[2, 0], u[2, 1])

        if phi1 < 0:
            phi1 = phi1 + 2 * np.pi
        if phi2 < 0:
            phi2 = phi2 + 2 * np.pi

        res[..., 0] = phi1
        res[..., 1] = PHI
        res[..., 2] = phi2


# @numba.njit
# def U_map_to_euler_map(U_map):
#     """From xfab.tools"""
#     tol = 1e-8
#     euler_map = np.zeros((U_map.shape[0], U_map.shape[1], U_map.shape[2], 3))
#     euler_map.fill(np.nan)
#     for i in range(U_map.shape[0]):
#         for j in range(U_map.shape[1]):
#             for k in range(U_map.shape[2]):
#                 U = U_map[i, j, k]
#                 if np.isnan(U[0, 0]):
#                     euler_map[i, j, k] = np.nan
#                 else:
#                     PHI = np.arccos(U[2, 2])
#                     if np.abs(PHI) < tol:
#                         phi1 = _arctan2(-U[0, 1], U[0, 0])
#                         phi2 = 0
#                     elif np.abs(PHI - np.pi) < tol:
#                         phi1 = _arctan2(U[0, 1], U[0, 0])
#                         phi2 = 0
#                     else:
#                         phi1 = _arctan2(U[0, 2], -U[1, 2])
#                         phi2 = _arctan2(U[2, 0], U[2, 1])
#
#                     if phi1 < 0:
#                         phi1 = phi1 + 2 * np.pi
#                     if phi2 < 0:
#                         phi2 = phi2 + 2 * np.pi
#
#                     euler_map[i, j, k, 0] = phi1
#                     euler_map[i, j, k, 1] = PHI
#                     euler_map[i, j, k, 2] = phi2
#
#                     # euler_map[i, j, k] = np.array([phi1, PHI, phi2])
#
#     return euler_map


class TensorMap:
    """This is a class to store a contiguous voxel-based representation of a sample.
    At its core is the self.maps attribute, which is a dictionary of Numpy arrays.
    Each Numpy array should represent a 3D voxel grid of the sample.
    The dimensions of the array are aligned with the laboratory reference frame in the order (Z, Y, X, ...)
    The shape of the first three axes should therefore be (NZ, NY, NX, ...)
    The total number of dimensions can vary
    E.g a phase ID map might be (1, 15, 20) but a UBI map might be (1, 15, 20, 3, 3)
    """

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
        for name in ("U", "B", "UB", "mt", "unitcell", "euler"):
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
    def UBI(self):
        """The UBI matrix"""
        return self.maps['UBI']

    @UBI.setter
    def UBI(self, value):
        self.clear_cache()
        self.maps['UBI'] = value

    @property
    def UB(self):
        """The UB matrix"""
        if "UB" in self.keys():
            return self.maps['UB']
        else:
            if 'UBI' not in self.keys():
                raise KeyError('No UBI to calculate from!')
            # calculate it
            UBI_map = self.maps['UBI']
            # old way (non-broadcast, fixed dimensionality):
            # UB_map = np.zeros_like(self.maps['UBI'])
            # nb_inv_3d(UBI_map, UB_map)
            # new way (broadcastable)
            UB_map = fast_invert(UBI_map)
            self.add_map('UB', UB_map)
            return UB_map

    @property
    def mt(self):
        """The metric tensor"""
        if "mt" in self.keys():
            return self.maps['mt']
        else:
            # calculate it
            UBI_map = self.maps["UBI"]
            mt_map = ubi_to_mt(UBI_map)
            self.add_map('mt', mt_map)
            return mt_map

    @property
    def unitcell(self):
        """The unitcell - from metric tensor"""
        if "unitcell" in self.keys():
            return self.maps['unitcell']
        else:
            # calculate it
            mt_map = self.mt
            dummy = np.arange(6)  # dummy variable, not used
            unitcell_map = mt_to_unitcell(mt_map, dummy)
            self.add_map('unitcell', unitcell_map)
            return unitcell_map

    @property
    def B(self):
        """The B matrix via the unitcell"""
        if "B" in self.keys():
            return self.maps['B']
        else:
            # calculate it
            unitcell_map = self.unitcell
            dummy = np.eye(3)  # dummy 3x3 variable, not used
            B_map = unitcell_to_b(unitcell_map, dummy)
            self.add_map('B', B_map)
            return B_map

    # @property
    # def B(self):
    #     # """The B matrix - using QR decomposition from UB"""
    #     """The B matrix - via unitcell"""
    #     if "B" in self.keys():
    #         return self.maps['B']
    #     else:
    #         # calculate it
    #         U_map, B_map = UB_map_to_U_B_map(self.UB)
    #         self.add_map('U', U_map)
    #         self.add_map('B', B_map)
    #         return B_map

    # @property
    # def U(self):
    #     """The U matrix - using QR decomposition from UB"""
    #     if "U" in self.keys():
    #         return self.maps['U']
    #     else:
    #         # calculate it
    #         U_map, B_map = UB_map_to_U_B_map(self.UB)
    #         self.add_map('U', U_map)
    #         self.add_map('B', B_map)
    #         return U_map

    @property
    def U(self):
        """The U matrix - via UBI and B"""
        if "U" in self.keys():
            return self.maps['U']
        else:
            # calculate it
            UBI_map = self.UBI
            B_map = self.B
            U_map = ubi_and_b_to_u(UBI_map, B_map)
            self.add_map('U', U_map)
            return U_map

    @property
    def euler(self):
        """The euler angles - using Numba-accelerated vectorised conversion"""
        if "euler" in self.keys():
            return self.maps['euler']
        else:
            # calculate it
            # euler_map = U_map_to_euler_map(self.U)
            U_map = self.U
            dummy = np.arange(3)
            euler_map = u_to_euler(U_map, dummy)
            self.add_map('euler', euler_map)
            return euler_map

    @property
    def dzero_unitcell(self):
        """The dzero unitcell for each voxel from the reference phase for that voxel"""
        if 'dzero_unitcell' in self.keys():
            return self.maps['dzero_unitcell']
        else:
            # calculate it
            dzero_unitcell_map = np.full(self.shape + (6,), np.nan, dtype=float)
            for phase_id, ucell in self.phases.items():
                dzero_unitcell_map[self.phase_ids == phase_id] = ucell.lattice_parameters
            self.add_map('dzero_unitcell', dzero_unitcell_map)
            return dzero_unitcell_map

    @property
    def eps_sample(self):
        """The per-voxel strain in sample frame, relative to the B0 of the unitcell of the reference phase for that voxel."""
        if 'eps_sample' in self.keys():
            return self.maps['eps_sample']
        else:
            # calculate it
            eps_sample_map = ubi_and_unitcell_to_eps_sample(self.UBI, self.dzero_unitcell)
            self.add_map('eps_sample', eps_sample_map)
            return eps_sample_map

    @property
    def eps_crystal(self):
        """The per-voxel strain in crystal frame, relative to the B0 of the unitcell of the reference phase for that
        voxel."""
        if 'eps_crystal' in self.keys():
            return self.maps['eps_crystal']
        else:
            # calculate it
            eps_crystal_map = ubi_and_unitcell_to_eps_crystal(self.UBI, self.dzero_unitcell)
            self.add_map('eps_crystal', eps_crystal_map)
            return eps_crystal_map

    # TODO - make multiphase - store C map for each voxel
    def get_stress(self, C, phase_id):
        """Fill in sig_crystal and sig_sample using 6x6 stiffness tensor C"""
        print('Warning! This is currently single-phase only - calling this more than once with different phases will '
              'overwrite the previous result')
        phase_mask = self.phase_ids == phase_id
        sig_crystal_map = strain_to_stress(self.eps_crystal, C, phase_mask)
        sig_sample_map = strain_to_stress(self.eps_sample, C, phase_mask)
        self.add_map('sig_crystal', sig_crystal_map)
        self.add_map('sig_sample', sig_sample_map)

    @property
    def sig_hydro(self):
        """The per-voxel hydrostatic stress (frame invariant)"""
        if 'sig_hydro' in self.keys():
            return self.maps['sig_hydro']
        else:
            sig_hydro_map = sig_to_hydro(self.sig_crystal)
            self.add_map('sig_hydro', sig_hydro_map)
            return sig_hydro_map

    @property
    def sig_mises(self):
        """The per-voxel von-Mises stress"""
        if 'sig_mises' in self.keys():
            return self.maps['sig_mises']
        else:
            sig_mises_map = sig_to_vm(self.sig_sample)
            self.add_map('sig_mises', sig_mises_map)
            return sig_mises_map

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
                    # store the phase name as an attribute
                    phase_group[str(phase_id)].attrs['phase_name'] = self.phases[phase_id].name

                    # store the step sizes
            parent_group.create_dataset("step", data=np.array(self.steps))

    def to_paraview(self, h5name, h5group='TensorMap'):
        """Exports to H5, then writes an XDMF file that lets you read the data with ParaView"""
        # Write H5 first
        if not os.path.exists(h5name):
            self.to_h5(h5name, h5group=h5group)

        h5_relpath = os.path.split(h5name)[1]

        xdmf_filename = h5name.replace('.h5', '.xdmf')

        dims = self.shape
        scalar_dims = dims
        vector_dims = dims + (3,)
        tensor_dims = dims + (3, 3,)

        MeshDimensions = (dims[0] + 1, dims[1] + 1, dims[2] + 1)

        MeshDimensionsStr = 'Dimensions="%d %d %d"' % MeshDimensions
        ScalarDimensionsStr = 'Dimensions="%d %d %d"' % scalar_dims
        VectorDimensionsStr = 'Dimensions="%d %d %d %d"' % vector_dims
        TensorDimensionsStr = 'Dimensions="%d %d %d %d %d"' % tensor_dims

        steps = tuple(self.steps)

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
                    fileID.write(
                        '      <DataItem Format="HDF" %s NumberType="Float" Precision="6" >%s:/%s</DataItem>\n' % (
                            ScalarDimensionsStr, h5_relpath, h5group + '/maps/' + map_name))
                    fileID.write('    </Attribute>\n')
                elif n_dims == 4:
                    # vector field (like IPF)
                    fileID.write('    <Attribute Name="%s" AttributeType="Vector" Center="Cell">\n' % map_name)
                    fileID.write(
                        '      <DataItem Format="HDF" %s NumberType="Float" Precision="6" >%s:/%s</DataItem>\n' % (
                            VectorDimensionsStr, h5_relpath, h5group + '/maps/' + map_name))
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

    def to_ctf_mtex(self, ctf_path, z_index=0):
        """Export a Z slice to a CTF file for MTEX processing.
        The resulting ctf file can be loaded in MTEX with the command:
        ebsd = EBSD.load(ctf_path)
        Note that no Euler or spatial conversions are needed"""
        euler_slice = np.degrees(self.euler[z_index, :, :])  # covert to degrees
        phase_slice = self.phase_ids[z_index, :, :] + 1  # unindexed should be 0

        euler_slice[np.isnan(euler_slice)] = 0.0  # nans to 0

        # get XY placements
        NY, NX = self.shape[1:]
        ystep, xstep = self.steps[1:]
        X, Y = np.meshgrid(np.arange(NY) * ystep,
                           np.arange(NX) * xstep)  # XY flip in MTEX - map only (orientations are fine)

        # flatten arrays
        npx = NY * NX
        euler_flat = euler_slice.reshape((npx, -1))
        phases_flat = phase_slice.reshape((npx, -1))
        Y_flat = Y.reshape((npx, -1))
        X_flat = X.reshape((npx, -1))
        XY_flat = np.hstack((X_flat, Y_flat))

        bands = np.zeros_like(Y_flat)
        error = np.zeros_like(Y_flat)
        mad = np.zeros_like(Y_flat)
        bc = np.zeros_like(Y_flat)
        bs = np.zeros_like(Y_flat)

        combined_array = np.hstack((phases_flat, XY_flat, bands, error, euler_flat, mad, bc, bs))

        header_lines = [
            "Channel Text File",
            "Prj unnamed",
            "Author\t[Unknown]",
            "JobMode\tGrid",
            "XCells\t%d" % NX,
            "YCells\t%d" % NY,
            "XStep\t%f" % xstep,
            "YStep\t%f" % ystep,
            "AcqE1\t0",
            "AcqE2\t0",
            "AcqE3\t0",
            "Euler angles refer to Sample Coordinate system (CS0)!	Mag	2E3	Coverage	100	Device	0	KV	1.5E1	TiltAngle	70	TiltAxis	0",
            "Phases\t%d" % len(self.phases)
        ]

        for phase in self.phases.values():
            phase_string = "%f;%f;%f\t%f;%f;%f\t%s\t11\t%s" % (
                phase.lattice_parameters[0], phase.lattice_parameters[1], phase.lattice_parameters[2],
                phase.lattice_parameters[3], phase.lattice_parameters[4], phase.lattice_parameters[5], phase.name,
                phase.symmetry)
            header_lines.extend([phase_string])

        header_lines.extend(["Phase\tX\tY\tBands\tError\tEuler1\tEuler2\tEuler3\tMAD\tBC\tBS"])

        header_text = "\n".join(header_lines)

        np.savetxt(ctf_path, combined_array, delimiter='\t', header=header_text, comments='',
                   fmt=['%d', '%06.6s', '%06.6s', '%d', '%d', '%06.6s', '%06.6s', '%06.6s', '%06.6s', '%d', '%d'])

        print('CTF exported!')
        print('In MTEX, run the command:')
        print("import_wizard('EBSD')")
        print("Click the '+', choose file %s and click 'Open'" % ctf_path)
        print("Click 'Next >>'")
        print("Click 'Next >>' though the phases, changing if necessary (should be right though)")
        print("Choose 'apply rotation to Euler angles and spatial coordinates' with an angle of [0,0,0]")
        print("Choose the default MTEX plotting convention (X east, Y north)")
        print("Click 'Next >>' then 'script (m-file)' for the 'Import to', then 'Finish'")
        print("Click 'Run' at the top, save the import script somewhere, then it should run")
        print("Your EBSD should now be imported into MTEX. You can try to plot it with 'plot(ebsd)'")

    @staticmethod
    def map_order_to_recon_order(map_arr, z_layer=0):
        """Transform a 2D slice of a TensorMap to reconstruction space. The opposite of recon_order_to_map_order"""
        # we are currently in (Z, Y, X, ...)
        nz, ny, nx = map_arr.shape[:3]
        # get the 2D slice
        slice_2D = map_arr[z_layer, ...].copy().squeeze()

        # now we are (Y, X, ...)
        # swap X and Y
        slice_2D = np.swapaxes(slice_2D, 0, 1)

        # now we are (X, Y, ...)
        # flip the second axis
        slice_2D = np.flip(slice_2D, 1)

        assert slice_2D.shape[:2] == (nx, ny)
        # now we are (X, -Y, ...)
        return slice_2D

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
        return (mk, yshape - mj - 1)

    @staticmethod
    def recon_index_to_map(ri, rj, yshape):
        """From a 2D reconstruction space index (ri, rj) and the Y shape of the map, determine
        the corresponding position in TensorMap space"""

        return (0, yshape - 1 - rj, ri)

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
                    phases[int(phase_id)].name = phase_group[phase_id].attrs['phase_name']

            steps = parent_group['step'][:]

        tensor_map = cls(maps=maps, phases=phases, steps=steps)

        return tensor_map

    def to_pbpmap(self, z_layer=0, default_npks=20, default_nuniq=20):
        """Get a PBPMap object from this TensorMap object (good for refining strains)"""
        ubi_pbpmap_order = self.map_order_to_recon_order(self.UBI, z_layer=z_layer)
        # this is in (ri, rj, 3, 3)

        # make the columns to hold the (si, sj) indices for the pbpmap
        si_col = np.zeros(ubi_pbpmap_order.shape[0] * ubi_pbpmap_order.shape[1])
        sj_col = np.zeros(ubi_pbpmap_order.shape[0] * ubi_pbpmap_order.shape[1])
        # make some columns to hold npks and nuniq
        npks_col = np.zeros(ubi_pbpmap_order.shape[0] * ubi_pbpmap_order.shape[1])
        nuniq_col = np.zeros(ubi_pbpmap_order.shape[0] * ubi_pbpmap_order.shape[1])
        # we need to make some UBI columns accordingly
        # first we transpose to (3, 3, ri, rj)
        # then when we reshape, we are doing (:2.ravel(), 2:4.ravel()) and structure is preserved
        ubi_cols = ubi_pbpmap_order.transpose(2, 3, 0, 1).reshape(9, len(si_col))
        row_idx = 0
        for ri in np.arange(ubi_pbpmap_order.shape[0]):
            for rj in np.arange(ubi_pbpmap_order.shape[1]):
                # get the si, sj values
                si, sj = recon_to_step(ri, rj, ubi_pbpmap_order.shape[:2])
                # get the mi, mj, mk values
                mi, mj, mk = self.recon_index_to_map(ri, rj, ubi_pbpmap_order.shape[1])
                si_col[row_idx] = si
                sj_col[row_idx] = sj
                # do we have anything at this reconstruction point?
                # otherwise we have nothing
                if self.phase_ids[mi, mj, mk] > -1:
                    npks_col[row_idx] = default_npks
                    nuniq_col[row_idx] = default_nuniq
                row_idx += 1

        # unpack ubi columns
        ubi00, ubi01, ubi02, ubi10, ubi11, ubi12, ubi20, ubi21, ubi22 = ubi_cols

        from ImageD11.sinograms.point_by_point import PBPMap
        output_map = PBPMap(new=True)
        output_map.nrows = ubi_pbpmap_order.shape[0] * ubi_pbpmap_order.shape[1]
        output_map.addcolumn(si_col, 'i')
        output_map.addcolumn(sj_col, 'j')
        output_map.addcolumn(npks_col, 'ntotal')
        output_map.addcolumn(nuniq_col, 'nuniq')
        output_map.addcolumn(ubi00, 'ubi00')
        output_map.addcolumn(ubi01, 'ubi01')
        output_map.addcolumn(ubi02, 'ubi02')
        output_map.addcolumn(ubi10, 'ubi10')
        output_map.addcolumn(ubi11, 'ubi11')
        output_map.addcolumn(ubi12, 'ubi12')
        output_map.addcolumn(ubi20, 'ubi20')
        output_map.addcolumn(ubi21, 'ubi21')
        output_map.addcolumn(ubi22, 'ubi22')

        return output_map

    @classmethod
    def from_pbpmap(cls, pbpmap, steps=None, phases=None):
        """Create TensorMap from a pbpmap object"""

        maps = dict()

        # see if we have a best UBI map to take
        if not hasattr(pbpmap, 'best_ubi'):
            raise ValueError('PBPMap has no best UBI selection to take from!')

        ubi_map = pbpmap.best_ubi

        # create a mask from ubi_map
        ubi_mask = np.where(np.isnan(ubi_map[:, :, 0, 0]), 0, 1).astype(bool)

        # reshape ubi map and add it to the dict
        maps['UBI'] = cls.recon_order_to_map_order(ubi_map)

        # add npks to the dict
        if hasattr(pbpmap, 'npks'):
            maps['npks'] = cls.recon_order_to_map_order(np.where(ubi_mask, pbpmap.best_npks, 0))

        # add nuniq to the dict
        if hasattr(pbpmap, 'nuniq'):
            maps['nuniq'] = cls.recon_order_to_map_order(np.where(ubi_mask, pbpmap.best_nuniq, 0))

        tensor_map = cls(maps=maps, steps=steps, phases=phases)

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
        raw_intensity_map = np.full(map_shape, cutoff_level)

        # make an empty phase ID map
        phase_id_map = np.full(map_shape, -1)

        # make an empty UBI map full of nan (invalid value)
        ubi_map = np.full((map_shape[:2] + (3, 3)), np.nan, float)

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

    @classmethod
    def from_stack(cls, tensormaps, zstep=1.0):
        """Stack multiple TensorMaps together along Z. The order of tensormaps will determine the Z-order, lowest-index first.
        All Tensor Maps in tensormaps must have the same phase mapping (i.e phase 0 is always the same unitcell)"""
        combined_maps = {}
        tm0 = tensormaps[0]

        # work out what map names we have in all of our tensormaps
        common_map_names = [map_name for map_name in tm0.keys() if all([map_name in tm.keys() for tm in tensormaps])]

        for map_name in common_map_names:
            # just concatenate the maps together (defaults to first axis which is Z)
            combined_maps[map_name] = np.concatenate([tm.maps[map_name] for tm in tensormaps])

        # make the tensormap object
        combined_tensormap = cls(maps=combined_maps, phases=tm0.phases, steps=(zstep, tm0.steps[1], tm0.steps[2]))
        return combined_tensormap

    @classmethod
    def from_combine_phases(cls, tensormaps):
        """Combine multiple mono-phase TensorMaps with different phases into one TensorMap.
           For now, this handles grain label collisions by offsetting the grain labels of
           subsequent tensormaps before adding"""
        combined_maps = {}
        tm0 = tensormaps[0]

        # make sure each input TensorMap only has one phase
        for tm in tensormaps:
            if len(tm.phases) > 1:
                raise ValueError("Each input TensorMap should only have one phase!")

        # combine phases
        combined_phases = {inc: tm.phases[0] for inc, tm in enumerate(tensormaps)}

        # work out what map names we have in all of our tensormaps
        common_map_names = [map_name for map_name in tm0.keys() if all([map_name in tm.keys() for tm in tensormaps])]

        for map_name in common_map_names:
            # get the base array (the first tensormap)
            base_arr = tm0[map_name].copy()

            # iterate over the other tensormaps
            for tm_inc, tm in enumerate(tensormaps[1:]):
                tm_inc += 1  # because we start at 1:
                # update the base array where the other tensormaps have nonzero phases
                if map_name == 'phase_ids':
                    # we are updating the phase id map for the combined TensorMap
                    # so the new array we put in is just the tm_inc

                    new_arr = tm_inc
                elif map_name == 'labels':
                    # for now, we need to make sure there's no collisions between grain IDs
                    # for now we are shifting grain labels of subsequent maps
                    # get the highest current grain ID
                    max_prev_gid = np.max(base_arr)

                    # make a version of labels for this tensormap that's shifted by max_prev_gid + 1
                    shifted_labels = tm['labels'].copy()
                    shifted_labels = np.where(shifted_labels > -1, shifted_labels + (max_prev_gid + 1), shifted_labels)

                    # make sure that aside from -1, no grain labels intersect between base_arr and the shifted labels we will introduce
                    assert len(np.intersect1d(np.unique(base_arr)[1:], np.unique(shifted_labels)[1:])) == 0

                    new_arr = shifted_labels
                else:
                    new_arr = tm[map_name]

                # selectively overwrite base_arr with new_arr
                # the array slicing stuff below is to allow arbitrary rightward broadcasting
                # phase_id is (NZ, NY, NX) but base_arr might be UBI for example (NZ, NY, NX, 3, 3)
                # Numpy can auto-broadcast leftwards (e.g (NZ, NY, NX) to (3, 3, NZ, NY, NX))
                # but not rightwards!
                # so we need to slice like this (NZ, NY, NX)[..., np.newaxis, np.newaxis]
                # In Python 2, slicing grammar is different, so we can't invoke ... directly inside a tuple
                base_arr = np.where((tm['phase_ids'] > -1)[(Ellipsis,) + (np.newaxis,) * (base_arr.ndim - 3)], new_arr,
                                    base_arr)

            combined_maps[map_name] = base_arr

        combined_tm = cls(maps=combined_maps, phases=combined_phases, steps=tm0.steps)

        return combined_tm
