## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py
# ImageD11_v0.4 Software for beamline ID11
# Copyright (C) 2005  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

from __future__ import print_function

import h5py
import numpy as np
import xfab.tools

import ImageD11.finite_strain
import ImageD11.unitcell


# helpers : put these into xfab.tools at some point?
def e6_to_symm(e):
    """ Follows the eps convention from xfab 
    eps = [e11, e12, e13, e22, e23, e33]
    """
    e11, e12, e13, e22, e23, e33 = e
    return np.array(((e11, e12, e13),
                     (e12, e22, e23),
                     (e13, e23, e33)))


def symm_to_e6(m):
    """ Follows the eps convention from xfab 
    eps = [e11, e12, e13, e22, e23, e33]
    """
    return np.array((m[0, 0], m[0, 1], m[0, 2],
                     m[1, 1], m[1, 2],
                     m[2, 2]))


class grain:
    def __init__(self, ubi, translation=None, **kwds):
        if translation is None:
            # If translation has not been read from ubi file make it 
            # be None to avoid confusion with a grain which is known
            # to be at [0,0,0]
            self.translation = None
        else:
            self.translation = np.array(translation, float)
        self.set_ubi(ubi)

    def set_ubi(self, ubi):
        """ Update the orientation and clear cached values """
        self.ubi = np.array(ubi, float)
        assert np.linalg.det(self.ubi) >= 0, 'Left handed axis system!'
        self.clear_cache()

    def clear_cache(self):
        # We will cache a bunch of things and access them
        # via properties. Usually we return copies of these
        # so we avoid people overwriting cached values by accident.
        # e.g. cell = grain.unitcell
        #      cell[4:] = 90 # degree angle for orthorhombic
        #      oups : grain._unitcell was corrupted!
        self._UB = None
        self._U = None
        self._B = None
        self._rod = None
        self._mt = None
        self._rmt = None
        self._unitcell = None
        self._ref_unitcell = None
        self._orix_orien = None

    @property
    def UB(self):
        """ The UB matrix from Busing and Levy
        columns are the reciprocal space lattice vectors
        """
        if self._UB is None:
            self._UB = np.linalg.inv(self.ubi)
        return self._UB.copy()

    @property
    def ub(self):
        return self.UB

    @property
    def B(self):
        if self._B is None:
            self._B = ImageD11.unitcell.unitcell(self.unitcell).B.copy()
        return self._B.copy()

    @property
    def U(self):
        """ The orientation matrix (U) from Busing and Levy """
        if self._U is None:
            # ubi = inv(UB) = inv(B)inv(U)
            self._U = np.dot(self.B, self.ubi).T
        return self._U.copy()

    @property
    def u(self):
        return self.U

    @property
    def Rod(self):
        """ A Rodriguez vector. 
        Length proportional to angle, direction is axis"""
        if self._rod is None:
            self._rod = xfab.tools.u_to_rod(self.U)
        return self._rod.copy()

    @property
    def mt(self):
        """Metric tensor """
        if self._mt is None:
            self._mt = np.dot(self.ubi, self.ubi.T)
        return self._mt.copy()

    @property
    def rmt(self):
        """ Reciprocal metric tensor """
        if self._rmt is None:
            self._rmt = np.linalg.inv(self.mt)
        return self._rmt.copy()

    @property
    def unitcell(self):
        """ a,b,c,alpha,beta,gamma """
        if self._unitcell is None:
            G = self.mt
            a, b, c = np.sqrt(np.diag(G))
            al = np.degrees(np.arccos(G[1, 2] / b / c))
            be = np.degrees(np.arccos(G[0, 2] / a / c))
            ga = np.degrees(np.arccos(G[0, 1] / a / b))
            self._unitcell = np.array((a, b, c, al, be, ga))
        return self._unitcell.copy()

    def eps_grain_matrix(self, dzero_cell, m=0.5):
        """ dzero_cell can be another grain or cell parameters
                [a,b,c,alpha,beta,gamma]
            m is the exponent for the Seth-Hill finite strain tensors
            E = 1/2m (U^2m - I)
            (Negative m reverses the reference and lab systems).
        Returns eps as a symmetric matrix
        ... in the grain reference system of dzero_cell
        """
        if hasattr(dzero_cell, "UB"):
            B = dzero_cell.UB
        else:
            B = ImageD11.unitcell.unitcell(dzero_cell).B
        F = ImageD11.finite_strain.DeformationGradientTensor(self.ubi, B)
        eps = F.finite_strain_ref(m)
        return eps

    def eps_grain(self, dzero_cell, m=0.5):
        """ dzero_cell can be another grain or cell parameters
                [a,b,c,alpha,beta,gamma]
            m is the exponent for the Seth-Hill finite strain tensors
            E = 1/2m (U^2m - I)
            (Negative m reverses the reference and lab systems).
        Returns eps as a the 6 independent entries in the symmetric matrix
         e11 e12 e13 e22 e23 e33
        ... in the grain reference system of dzero_cell
        """
        E = self.eps_grain_matrix(dzero_cell, m)
        return symm_to_e6(E)

    def eps_sample_matrix(self, dzero_cell, m=0.5):
        """ dzero_cell can be another grain or cell parameters:
                [a,b,c,alpha,beta,gamma]
            m is the exponent for the Seth-Hill finite strain tensors
            E = 1/2m (V^2m - I)
            (Negative m reverses the reference and lab systems).
        Returns eps as a symmetric matrix
        ... in the sample system (z up, x along the beam at omega=0)
        """
        if hasattr(dzero_cell, "UB"):
            B = dzero_cell.UB
        else:
            B = ImageD11.unitcell.unitcell(dzero_cell).B
        F = ImageD11.finite_strain.DeformationGradientTensor(self.ubi, B)
        eps = F.finite_strain_lab(m)
        return eps

    def eps_sample(self, dzero_cell, m=0.5):
        """ dzero_cell can be another grain or cell parameters:
                [a,b,c,alpha,beta,gamma]
            m is the exponent for the Seth-Hill finite strain tensors
            E = 1/2m (V^2m - I)
            (Negative m reverses the reference and lab systems).
        Returns eps as a the 6 independent entries in the symmetric matrix
         e11 e12 e13 e22 e23 e33
        ... in the sample system (z up, x along the beam at omega=0)
        """
        E = self.eps_sample_matrix(dzero_cell, m)
        return symm_to_e6(E)

    # TODO: we need I/O for this
    @property
    def ref_unitcell(self):
        """The reference (strain-free) unitcell object for this grain"""
        if self._ref_unitcell is None:
            raise NameError("Reference unitcell not set for this grain!")
        else:
            return self._ref_unitcell

    @ref_unitcell.setter
    def ref_unitcell(self, value):
        if not isinstance(value, ImageD11.unitcell.unitcell):
            raise TypeError("ref_unitcell must be an ImageD11.unitcell.unitcell instance!")
        self._ref_unitcell = value
        self._orix_orien = None  # must be recomputed because we changed the reference unitcell

    @property
    def orix_phase(self):
        """The orix phase for the grain, taken straight from the reference unitcell"""
        return self.ref_unitcell.orix_phase

    @property
    def orix_orien(self):
        if self._orix_orien is None:
            # get the orix orientation from the reference unit cell
            self._orix_orien = self.ref_unitcell.get_orix_orien_fast(self.U)
        return self._orix_orien

    # TODO: classmethod for from_orix_orien to make a grain from an orientation

    def get_ipf_colour(self, axis=np.array([0., 0., 1.])):
        """Calls ImageD11.unitcell.unitcell.get_ipf_colour with the grain's UB"""

        return self.ref_unitcell.get_ipf_colour_fast(self.U, axis)

    def to_h5py_group(self, parent_group, group_name):
        """Creates a H5Py group for this grain.
           parent_group is the parent H5py Group
           group_name is the name of the H5py Group for this name
           Very useful for saving lists of grains to an H5 file.
           Uses require_group to modify existing data if present"""

        grain_group = parent_group.require_group(group_name)

        # essential attributes:
        save_array(grain_group, 'ubi', self.ubi)

        # optional attributes:
        for attr in STRINGATTRS + NUMATTRS:
            try:
                value = getattr(self, attr)
                if value is not None:
                    grain_group[attr] = value
            except (NameError, AttributeError, TypeError, OSError):
                continue

        # write array attributes
        for attr in ARRATTRS:
            if hasattr(self, attr):
                if getattr(self, attr) is not None:
                    save_array(grain_group, attr, getattr(self, attr))

        return grain_group

    @classmethod
    def from_h5py_group(cls, grain_group):
        """Creates a grain object from an h5py group"""
        # read essential attributes"
        ubi = grain_group["ubi"][:]

        # make grain object:
        g = grain(ubi=ubi)

        # read optional attributes:
        # strings:
        for attr in STRINGATTRS:
            if attr in grain_group.keys():
                setattr(g, attr, grain_group.get(attr)[()].decode())

        # numbers:
        for attr in NUMATTRS:
            if attr in grain_group.keys():
                setattr(g, attr, grain_group.get(attr)[()])

        # arrays:
        for attr in ARRATTRS:
            if attr in grain_group.keys():
                setattr(g, attr, grain_group.get(attr)[:])

        return g


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

def write_grain_file(filename, list_of_grains):
    f = open(filename, "w")
    for g in list_of_grains:
        t = g.translation
        if t is not None:
            f.write("#translation: %g %g %g\n" % (t[0], t[1], t[2]))
        if hasattr(g, "name"):
            f.write("#name %s\n" % (g.name.rstrip()))
        if hasattr(g, "intensity_info"):
            f.write("#intensity_info %s\n" % (g.intensity_info.rstrip()))
        if hasattr(g, "npks"):
            f.write("#npks %d\n" % (int(g.npks)))
        if hasattr(g, "nuniq"):
            f.write("#nuniq %d\n" % (int(g.nuniq)))
        if hasattr(g, "Rod"):
            try:
                f.write("#Rod %f %f %f\n" % tuple([float(r) for r in g.Rod]))
            except:
                f.write("#Rod %s" % (g.Rod))
        f.write("#UBI:\n")
        u = g.ubi
        # More than float32 precision
        f.write("%.9g %.9g %.9g\n" % (u[0, 0], u[0, 1], u[0, 2]))
        f.write("%.9g %.9g %.9g\n" % (u[1, 0], u[1, 1], u[1, 2]))
        f.write("%.9g %.9g %.9g\n\n" % (u[2, 0], u[2, 1], u[2, 2]))
    f.close()


def read_grain_file(filename):
    """read ubifile and return a list of ubi arrays """
    f = open(filename, "r")
    grainsread = []
    u = []
    t = None
    p = {}
    for line in f:
        if line.find("#translation:") == 0:
            t = [float(x) for x in line.split()[1:]]
            continue
        if line[0] == "#" and line.find("UBI") < 0:
            k, v = line[1:].split(" ", 1)
            p[k] = v
            continue
        if line[0] == "#" and line.find("intensity_info") > -1:
            p["intensity_info"] = line.split("intensity_info")[1].rstrip()
        if line.find("#") == 0: continue
        vals = [float(x) for x in line.split()]
        if len(vals) == 3:
            u = u + [vals]
        if len(u) == 3:
            grainsread.append(grain(u, t))
            for k in ["name", "npks", "nuniq", "intensity_info"]:  # Rod - is recomputed when needed
                if k in p:
                    setattr(grainsread[-1], k, p[k])
            p = {}
            u = []
            t = None
    f.close()
    return grainsread


STRINGATTRS = ["intensity_info", "name"]
NUMATTRS = ["npks", "nuniq"]
ARRATTRS = ["translation"]
MATATTRS = ["UB", "ubi", "B", "U", "mt", "rmt"]

def write_grain_file_h5(filename, list_of_grains, group_name='grains'):
    """Write list of grains to H5py file."""
    # TODO: Use Silx Nexus IO instead?

    with h5py.File(filename, 'a') as hout:
        grains_group = hout.create_group(group_name)

        # the H5 grain file should be order-preserving
        # i.e it should always be possible to read the grains in the same order as they are written
        # so we will use the list positions as the names for the H5 groups
        # but we will try to read the names if we have them

        for ginc, g in enumerate(list_of_grains):
            group_name = str(ginc)
            g.to_h5py_group(parent_group=grains_group, group_name=group_name)


def read_grain_file_h5(filename, group_name='grains'):
    """Read list of grains from H5py file"""
    with h5py.File(filename, 'r') as hin:
        grains_group = hin[group_name]
        grains = []

        # take all the keys in the grains group, sort them by integer value, iterate
        for gid_string in sorted(grains_group.keys(), key=lambda x: int(x)):
            grain_group = grains_group[gid_string]

            g = grain.from_h5py_group(grain_group)

            grains.append(g)

    return grains


def export_to_vtp_paraview(grains, output_path, colour='npks'):
    """
    Export list of grains to vtp file for Paraview visualisation, with helpful arrays included.
    Thanks to Marilyn Sarkis for the idea.
    """
    import vtk
    from vtk.util import numpy_support as ns
    from scipy.spatial.transform import Rotation as R 

    def add_scalar_attribute(pd, attr, np_array):
        """Add scalar attribute array to polydata"""
        attr_array = ns.numpy_to_vtk(np_array, deep=True)
        attr_array.SetName(attr)
        pd.GetPointData().AddArray(attr_array)
    
    def add_vector_attribute(pd, attr, np_array):
        """Add vector attribute array to polydata"""
        attr_array = ns.numpy_to_vtk(np_array, deep=True)
        attr_array.SetNumberOfComponents(np_array.shape[1])  # required for vectors
        attr_array.SetName(attr)
        pd.GetPointData().AddArray(attr_array)
    
    # ensure we have ipf colours
    rgbattr = 'rgb_z'
    try:
        col = [getattr(grain, rgbattr) for grain in grains]  # IPF colour
    except AttributeError:
        # couldn't get the IPF attributes
        # try to compute it first
        # will still fail if we don't have reference unitcells
        from ImageD11.nbGui.nb_utils import get_rgbs_for_grains
        get_rgbs_for_grains(grains)
        col = [getattr(grain, rgbattr) for grain in grains]  # IPF colour

    translations = np.array([grain.translation for grain in grains])
    
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()
    
    for i, t in enumerate(translations):
        pid = points.InsertNextPoint(t[0], t[1], t[2])
        verts.InsertNextCell(1)
        verts.InsertCellPoint(pid)
    
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetVerts(verts)

    # add scalar attributes to polydata
    for numattr in NUMATTRS:
        np_array = np.array([float(getattr(grain, numattr)) for grain in grains])
        add_scalar_attribute(polydata, numattr, np_array)
    
    # compute volumes
    prop_volumes = np.array([float(grain.intensity_info.split("mean = ")[1].split(" , ")[0].replace("'", "")) for grain in grains])
    prop_radii = np.cbrt(3*prop_volumes/(4*np.pi))
    
    add_scalar_attribute(polydata, 'prop_volume', prop_volumes)
    add_scalar_attribute(polydata, 'prop_radius', prop_radii)
    
    # add vector attributes (e.g. colours)
    for vecattr in ["rgb_x", "rgb_y", "rgb_z", "Rod", "unitcell"]:
        np_array = np.array([getattr(grain, vecattr) for grain in grains])
        add_vector_attribute(polydata, vecattr, np_array)

    # compute quaternions for glyph orientation in Paraview

    quats = R.from_matrix([g.U for g in grains]).as_quat(scalar_first=True)  # WXYZ order to match Paraview
    
    add_vector_attribute(polydata, 'quat', quats)

    # it's also useful to have lab-frame vector fields for (x_c, y_c, z_c) and (a, b, c) and (astar, bstar, cstar) for viz
    # x_c, y_c, z_c are columns of U
    add_vector_attribute(polydata, "x_c", np.array([g.U[:, 0] for g in grains]))
    add_vector_attribute(polydata, "y_c", np.array([g.U[:, 1] for g in grains]))
    add_vector_attribute(polydata, "z_c", np.array([g.U[:, 2] for g in grains]))

    # a, b, c are rows of UBI (because we invert)
    add_vector_attribute(polydata, "a", np.array([g.ubi[0, :] for g in grains]))
    add_vector_attribute(polydata, "b", np.array([g.ubi[1, :] for g in grains]))
    add_vector_attribute(polydata, "c", np.array([g.ubi[2, :] for g in grains]))

    # a*, b*, c* are columns of UB
    add_vector_attribute(polydata, "astar", np.array([g.UB[:, 0] for g in grains]))
    add_vector_attribute(polydata, "bstar", np.array([g.UB[:, 1] for g in grains]))
    add_vector_attribute(polydata, "cstar", np.array([g.UB[:, 2] for g in grains]))
    
    # add matrix attributes
    for matattr in MATATTRS:
        tensor_components = [
            ('11', 0, 0), ('12', 0, 1), ('13', 0, 2),
            ('21', 1, 0), ('22', 1, 1), ('23', 1, 2),
            ('31', 2, 0), ('32', 2, 1), ('33', 2, 2)
        ]
        for comp_name, i, j in tensor_components:
            attr_name = "{}_{}".format(matattr, comp_name)
            np_array = np.array([getattr(grain, matattr)[i, j] for grain in grains])
            add_scalar_attribute(polydata, attr_name, np_array)

    eps_s = np.array([grain.eps_sample_matrix(dzero_cell=grain.ref_unitcell.lattice_parameters) for grain in grains])
    eps_c = np.array([grain.eps_grain_matrix(dzero_cell=grain.ref_unitcell.lattice_parameters) for grain in grains])

    # it's a strain tensor, so we can name them xx, yy, zz etc
    tensor_components = [
        ('xx', 0, 0), ('xy', 0, 1), ('xz', 0, 2),
        ('yx', 1, 0), ('yy', 1, 1), ('yz', 1, 2),
        ('zx', 2, 0), ('zy', 2, 1), ('zz', 2, 2)
    ]
    for comp_name, i, j in tensor_components:
        add_scalar_attribute(polydata, "eps_s_{}".format(comp_name), eps_s[:,i,j])
        add_scalar_attribute(polydata, "eps_c_{}".format(comp_name), eps_c[:,i,j])
    
    # set the active scalars
    polydata.GetPointData().SetActiveScalars(colour)
    
    # export to file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_path)
    writer.SetInputData(polydata)
    writer.Write()