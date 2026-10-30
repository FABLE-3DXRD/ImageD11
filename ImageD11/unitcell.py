
from __future__ import print_function, division

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

#
#
#
import logging, math
from math import fabs

import numpy as np
from numpy.linalg import inv
from ImageD11 import cImageD11
from xfab import tools
from ImageD11.parameters import parameters, AnalysisSchema


def _inv3(m):
    """Closed-form 3x3 inverse. Raises LinAlgError on singular, like np.linalg.inv."""
    a, b, c = m[0]
    d, e, f = m[1]
    g, h, i = m[2]
    A =  e*i - f*h
    B = -(d*i - f*g)
    C =  d*h - e*g
    det = a*A + b*B + c*C
    if det == 0.0:
        raise np.linalg.LinAlgError("Singular matrix")
    r = 1.0 / det
    return np.array((
        (A*r,            -(b*i - c*h)*r,   (b*f - c*e)*r),
        (B*r,             (a*i - c*g)*r,  -(a*f - c*d)*r),
        (C*r,            -(a*h - b*g)*r,   (a*e - b*d)*r)))

def radians(x):
    return x* math.pi / 180.


def degrees(x):
    return x * 180. / math.pi


def cross(a, b):
    """
    a x b has length |a||b|sin(theta)
    """
    return np.array([a[1] * b[2] - a[2] * b[1],
                     a[2] * b[0] - b[2] * a[0],
                     a[0] * b[1] - b[0] * a[1]], float)


def norm2(a):
    """
    Compute the unit 2 norm
    """
    return np.sqrt(np.dot(a, a))


def unit(a):
    """
    Normalise vector a to unit length
    """
    try:
        return a / norm2(a)
    except:
        logging.error("cannot normalise to unit length a=%s" % (str(a)))
        raise


# Systematic absences

def P(h, k, l):
    return False


def A(h, k, l):
    return (k + l) % 2 != 0


def B(h, k, l):
    return (h + l) % 2 != 0


def C(h, k, l):
    return (h + k) % 2 != 0


def I(h, k, l):
    return (h + k + l) % 2 != 0


def F(h, k, l):
    return (h + k) % 2 != 0 or (h + l) % 2 != 0 or (k + l) % 2 != 0


def R(h, k, l):
    return (-h + k + l) % 3 != 0


outif = {
    "P": P,
    "A": I,
    "B": B,
    "C": C,
    "I": I,
    "F": F,
    "R": R}


def orient_BL(B, h1, h2, g1, g2):
    """Algorithm was intended to follow this one using 2 indexed
    reflections.
    W. R. Busing and H. A. Levy. 457. Acta Cryst. (1967). 22, 457
    """
    h1c = np.dot(B, h1)  # cartesian H1
    h2c = np.dot(B, h2)  # cartesian H2
    t1c = unit(h1c)  # unit vector along H1
    t3c = unit(np.cross(h1c, h2c))
    t2c = unit(np.cross(h1c, t3c))
    t1g = unit(g1)
    t3g = unit(np.cross(g1, g2))
    t2g = unit(np.cross(g1, t3g))
    T_g = np.transpose(np.array([t1g, t2g, t3g]))  # Array are stored by rows and
    T_c = np.transpose(np.array([t1c, t2c, t3c]))  # these are columns
    U = np.dot(T_g, np.linalg.inv(T_c))
    UB = np.dot(U, B)
    UBI = np.linalg.inv(UB)
    return UBI, UB


def cosangles_many(ha, hb, gi):
    """ Finds the cosines of angles between two lists of hkls in
    reciprocal metric gi """
    assert len(ha[0]) == 3 and len(hb[0]) == 3
    na = len(ha)
    nb = len(hb)
    hag = np.dot(ha, gi)
    hbg = np.dot(hb, gi)
    hagha = np.sqrt((hag * ha).sum(axis=1))
    hbghb = np.sqrt((hbg * hb).sum(axis=1))
    haghb = np.dot(ha, hbg.T)
    ca = haghb / np.outer(hagha, hbghb)
    assert ca.shape == (na, nb)
    return ca


def cellfromstring(s):
    items = s.split()
    # print items
    latt = [float(x) for x in items[0:6]]
    try:
        symm = items[6]
    except IndexError:
        symm = 'P'
    return unitcell(latt, symm)


class Phases(AnalysisSchema):
    """
    Phases class - extends AnalysisSchema from xfab.parameters
    Reads analysis parameters from file
    Contains self.unitcells which is a dict of unitcell objects
    """
    def __init__(self, filename=None):
        super(Phases, self).__init__(filename=filename)
        self.unitcells = {}
        if filename is not None:
            self.get_unitcells()
    
    def add_phase_from_unitcell(self, phase_name, unitcell):
        """
        Add phase from unitcell object
        
        phase_name: the name of the phase
        unitcell: the unitcell object
        """
        super(Phases, self).add_phase_from_unitcell(phase_name, unitcell)
        self.get_unitcells()
    
    def get_unitcells(self):
        # dict of parameter objects for each phase
        for phase_name, phase_pars_obj in self.phase_pars_obj_dict.items():
            self.unitcells[phase_name] = unitcell_from_parameters(phase_pars_obj)
            # set the name of the unitcell from the phase name
            self.unitcells[phase_name].name = phase_name

# R space groups can be on hexagonal (standard) or rhombohedral axes
RHOMBOHEDRAL_SPACEGROUPS = (146, 148, 155, 160, 161, 166, 167)
# relative change of the metric tensor accepted as symmetric
METRIC_TOL = 1e-5


def laue_misfit(g, spacegroup, cell_choice="standard"):
    """
    Largest change of the real space metric tensor g under the Laue class
    rotations of spacegroup, relative to g. Zero when the cell has the symmetry.
    """
    from ImageD11 import sym_u
    ops = sym_u.laue_rotations(spacegroup, cell_choice).group
    return max(np.abs(np.dot(np.dot(o, g), o.T) - g).max() for o in ops) / np.abs(g).max()


class unitcell:
    # Unit cell stuff
    # Generate a list of peaks from a unit cell
    def __init__(self, lattice_parameters, symmetry="P", verbose=0, name=None, sgname=None):
        """
        Unit cell class
        supply a list (tuple etc) of a,b,c,alpha,beta,gamma
        optionally a lattice symmetry, one of "P","A","B","C","I","F","R"
        or : a space group number (integer)
        optionally sgname : an xfab space group name (cell_sg in a .par file)
          e.g. "Fm-3m" or "R-3mr" (rhombohedral axes). It overrides symmetry
          for the space group and the hkls.
        With a space group, a warning is logged if the cell does not have its symmetry.
        """
        self.lattice_parameters = np.array(lattice_parameters)
        if self.lattice_parameters.shape[0] != 6:
            raise Exception("You must supply 6 lattice parameters\n" + \
                            "      a,b,c,alpha,beta,gamma")

        if isinstance(symmetry, str):
            if symmetry.isdigit():
                symmetry = int(symmetry)

        if isinstance(symmetry, int):
            # this is an integer spacegroup
            self.symmetry = symmetry
            self.spacegroup = symmetry
        else:

            if symmetry not in ["P", "A", "B", "C", "I", "F", "R"]:
                raise Exception("Your symmetry " + symmetry + \
                                " was not recognised")
            else:
                self.symmetry = symmetry
            # assigning a function here!
            self.absent = outif[self.symmetry]

        a = self.lattice_parameters[0]
        b = self.lattice_parameters[1]
        c = self.lattice_parameters[2]
        self.alpha = radians(self.lattice_parameters[3])
        ca = math.cos(radians(self.lattice_parameters[3]))
        cb = math.cos(radians(self.lattice_parameters[4]))
        cg = math.cos(radians(self.lattice_parameters[5]))
        if verbose == 1: print("Unit cell", self.lattice_parameters)
        self.g = np.array([[a * a, a * b * cg, a * c * cb],
                           [a * b * cg, b * b, b * c * ca],
                           [a * c * cb, b * c * ca, c * c]], float)
        if verbose == 1: print("Metric tensor\n", self.g)
        try:
            self.gi = inv(self.g)
        except:
            raise Exception("Unit cell was degenerate, could not determine" + \
                            "reciprocal metric tensor")
        if verbose == 1: print("Reciprocal Metric tensor\n", self.gi)
        self.astar = np.sqrt(self.gi[0, 0])
        self.bstar = np.sqrt(self.gi[1, 1])
        self.cstar = np.sqrt(self.gi[2, 2])

        self.alphas = degrees(math.acos(self.gi[1, 2] / self.bstar / self.cstar))
        self.betas = degrees(math.acos(self.gi[0, 2] / self.astar / self.cstar))
        self.gammas = degrees(math.acos(self.gi[0, 1] / self.astar / self.bstar))
        if verbose == 1: print("Reciprocal cell")
        if verbose == 1:
            print(self.astar, self.bstar, self.cstar, \
                  self.alphas, self.betas, self.gammas)
        # Equation 3 from Busing and Levy
        self.B = np.array(
            [[self.astar,
              self.bstar * math.cos(radians(self.gammas)),
              self.cstar * math.cos(radians(self.betas))],
             [0,
              self.bstar * math.sin(radians(self.gammas)),
              -self.cstar * math.sin(radians(self.betas)) * ca],
             [0, 0,
              1. / c]], float)
        if verbose == 1: print(self.B)
        if verbose == 1:
            print(np.dot(np.transpose(self.B),
                         self.B) - self.gi)  # this should be zero
        self.hkls = None
        self.peaks = None
        self.limit = 0

        self.ringtol = 0.001
        self.ringds = []  # a list of floats
        self.ringhkls = {}  # a dict of lists of integer hkl
        self._rings_for = None
        # used for caching
        self.anglehkl_cache = {"ringtol": self.ringtol,
                               "B": self.B,
                               "BI": np.linalg.inv(self.B)}
        if name is not None:
            self.name = name
        
        # orix stuff
        self._orix_phase = None
        self._orix_alignment = None

        self.sgname = sgname
        self.cell_choice = "standard"
        if sgname is not None:
            self._set_sgname(sgname)
        elif isinstance(self.symmetry, int):
            if (self.spacegroup in RHOMBOHEDRAL_SPACEGROUPS and
                    laue_misfit(self.g, self.spacegroup) > METRIC_TOL and
                    laue_misfit(self.g, self.spacegroup, "rhombohedral") <= METRIC_TOL):
                self.cell_choice = "rhombohedral"
        if hasattr(self, "spacegroup"):
            misfit = laue_misfit(self.g, self.spacegroup, self.cell_choice)
            if misfit > METRIC_TOL:
                logging.warning("Unit cell %s does not have the symmetry of space group %s (%s axes): "
                                "the metric changes by %.2g under its rotations. Check the lattice "
                                "parameters, or the setting (monoclinic: unique axis b)." % (
                                    str(self.lattice_parameters), str(self.spacegroup),
                                    self.cell_choice, misfit))

    def _set_sgname(self, sgname):
        """ Space group from an xfab name, checked against self.symmetry """
        from xfab import sg
        spg = sg.sg(sgname=sgname)
        self.cell_choice = spg.cell_choice
        if isinstance(self.symmetry, int):
            if self.symmetry != spg.no:
                logging.warning("Space group %d does not match %s (number %d): using %s" % (
                    self.symmetry, sgname, spg.no, sgname))
        else:
            centring = "P" if spg.cell_choice == "rhombohedral" else spg.name[0]
            if self.symmetry != centring:
                logging.warning("Lattice %s does not match %s (lattice %s): using %s" % (
                    self.symmetry, sgname, centring, sgname))
        self.spacegroup = spg.no

    def __repr__(self):
        if (not hasattr(self, 'name')) or (self.name is None):
            return "Unitcell" + " | " + str(self.lattice_parameters) + " | " + str(self.symmetry)
        else:
            return str(self.name) + " | " + str(self.lattice_parameters) + " | " + str(self.symmetry)
    
    @property
    def orix_phase(self):
        """Get this unitcell as an Orix Phase object. Useful for pole figures / IPF colours"""
        if self._orix_phase is None:
            # compute the orix phase

            try:
                from orix.quaternion.symmetry import get_point_group
                from orix.crystal_map import Phase
                from diffpy.structure import Lattice, Structure
            except ImportError:
                raise ImportError("Missing orix, can't compute orix phase!")

            try:
                spacegroup = self.spacegroup
            except AttributeError:
                raise AttributeError("You must initialise this unitcell with an integer spacegroup symmetry for Orix to work!")

            lattice = Lattice(*tuple(self.lattice_parameters))
            structure = Structure(lattice=lattice)
            pg = get_point_group(space_group_number=spacegroup)
            phase = Phase(point_group=pg, structure=structure)
            self._orix_phase = phase
        return self._orix_phase

    @property
    def orix_alignment(self):
        """Rotation C of the crystal frame (as in orix_phase) to where orix has the
           point group of this cell: orientations go to orix as U.C^T.
           orix fixes the lattice frame (x along a, z along c*) and each point group
           to its standard axes. These disagree for b-unique 2/m (the 2-fold is along y),
           -31m, and R groups on rhombohedral axes. C is the identity otherwise.
           Several C do that and give different colours (e.g. any turn about the 2-fold
           of 2/m, or +-30 degrees for -31m): the symmetry does not choose. Here:
           - R groups on rhombohedral axes are coloured as the same lattice on its
             hexagonal (obverse) axes, where orix needs no turn.
           - otherwise C takes a* as close as it can to orix's x (a* is x in ImageD11's
             frame, and MTEX's default), then b as close as it can to +z. For -31m
             this gives MTEX's colours. For 2/m, orix's colour key is not MTEX's."""
        if self._orix_alignment is None:
            from ImageD11 import sym_u
            phase = self.orix_phase
            M = self.B.dot(np.linalg.inv(phase.structure.lattice.recbase))
            Mi = np.linalg.inv(M)
            if self.cell_choice == "rhombohedral":
                # a_h = a - b, b_h = b - c, c_h = a + b + c : rows of ubi_h = P.ubi
                P = np.array([[1., -1., 0.], [0., 1., -1.], [1., 1., 1.]])
                a, alpha = self.lattice_parameters[0], radians(self.lattice_parameters[3])
                hexcell = unitcell([2 * a * math.sin(alpha / 2), 2 * a * math.sin(alpha / 2),
                                    a * math.sqrt(3 * (1 + 2 * math.cos(alpha))), 90., 90., 120.],
                                   self.spacegroup)
                # U_h = U.N, then as get_orix_orien_fast does for hexcell
                N = hexcell.B.dot(P).dot(np.linalg.inv(self.B)).T
                Mh = hexcell.B.dot(np.linalg.inv(hexcell.orix_phase.structure.lattice.recbase))
                self._orix_alignment = Mi.dot(N).dot(Mh).dot(hexcell.orix_alignment.T).T
                return self._orix_alignment
            # the cell's own rotations, as they act on orix orientations
            ours = [Mi.dot(S).dot(M).T for S in
                    sym_u.cartesian(sym_u.laue_rotations(self.spacegroup, self.cell_choice), self.B)]
            theirs = []
            for O in phase.point_group.laue.proper_subgroup.to_matrix():
                if not any(np.abs(O - T).max() < 1e-6 for T in theirs):
                    theirs.append(np.asarray(O, float))
            # a* and b in orix's crystal frame: orix vector = M^T . ImageD11 vector
            astar = M[0]
            b = M.T.dot(np.linalg.inv(self.B)[1])
            try:
                self._orix_alignment = sym_u.align(ours, theirs,
                                                   refs=[(astar, (1., 0., 0.)), (b, (0., 0., 1.))])
            except ValueError as e:
                logging.warning("Cell %s, space group %s: %s. IPF colours of equivalent "
                                "orientations may differ" % (str(self.lattice_parameters),
                                                             str(self.spacegroup), str(e)))
                self._orix_alignment = np.eye(3)
        return self._orix_alignment

    def get_orix_orien(self, UBs):
        """Get an Orix orientation (containing one or more orientations)
           using self.orix_phase from an array of ImageD11 UB matrices.
           U is taken as grain.U does: B from the cell of each UB"""
        UBs = np.squeeze(UBs)
        # B is upper triangular (Busing and Levy) with B^T.B = UB^T.UB
        UBTs = np.ascontiguousarray(np.swapaxes(UBs, -1, -2))  # matmul is slow on the transposed view
        B = np.swapaxes(np.linalg.cholesky(np.matmul(UBTs, UBs)), -1, -2)
        return self.get_orix_orien_fast(np.matmul(UBs, np.linalg.inv(B)))

    def get_orix_orien_fast(self, Umats):
        """NOTE - Supply U not UB!
           Get an Orix orientation (containing one or more orientations)
           using self.orix_phase from an array of ImageD11 U matrices.
           """
        try:
            from orix.quaternion import Orientation
        except ImportError:
            raise ImportError("Missing orix, can't compute orix orien!")

        # U.B is invariant
        # U_id11 B_id11 == U_orix B_orix
        # U_orix = U_id11 B_id11 (B_orix)^-1
        # then in orix, orientations are lab -> crystal
        # in ImageD11, orientations are crystal -> lab
        # given U_orix with shape (N, 3, 3) we transpose the last 2 axes accordingly
        # .C^T puts the crystal frame where orix has the point group (orix_alignment)
        
        phase = self.orix_phase
        B0_id11 = self.B  # get the ImageD11 B0 matrix
        B0_orix = phase.structure.lattice.recbase  # get the orix B0 matrix
        U_orix = Umats.dot(B0_id11.dot(np.linalg.inv(B0_orix)).dot(self.orix_alignment.T))

        # swapaxes performs the transposition, works for (3,3) or (N,3,3)
        final_ori = Orientation.from_matrix(np.swapaxes(U_orix, -1, -2), symmetry=phase.point_group)

        return final_ori

    def get_ipf_colour_from_orix_orien(self, orix_orien, axis=np.array([0., 0., 1.])):
        """Gets the TSL IPF colour(s) (Nolze, 2016) in the laboratory frame from an orix orientation
           10.1107/S1600576716012942
           axis is an axis vector in the laboratory frame
           The colours are orix's (IPFColorKeyTSL). They agree with MTEX's, except for
           Laue class 2/m, where the two keys differ."""
        try:
            from orix.vector.vector3d import Vector3d
            from orix.plot import IPFColorKeyTSL
        except ImportError:
            raise ImportError("Missing diffpy and/or orix, can't get IPF colour!")

        ipf_direction = Vector3d(axis)

        ipfkey = IPFColorKeyTSL(self.orix_phase.point_group, direction=ipf_direction)
        rgb = ipfkey.orientation2color(orix_orien)

        return rgb

    def get_ipf_colour_fast(self, Us, axis=np.array([0., 0., 1.])):
        """Gets the TSL IPF colour (Nolze, 2016) in the laboratory frame for a list of U matrices
           10.1107/S1600576716012942
           axis is an axis vector in the laboratory frame
           optionally return the orien which is useful for IPF plotting"""

        orien = self.get_orix_orien_fast(Us)

        rgb = self.get_ipf_colour_from_orix_orien(orien, axis=axis)

        return np.squeeze(rgb)
    
    def get_ipf_colour(self, UBs, axis=np.array([0., 0., 1.])):
        """Gets the TSL IPF colour (Nolze, 2016) in the laboratory frame for a list of U matrices
           10.1107/S1600576716012942
           axis is an axis vector in the laboratory frame
           optionally return the orien which is useful for IPF plotting"""

        orien = self.get_orix_orien(UBs)

        rgb = self.get_ipf_colour_from_orix_orien(orien, axis=axis)

        return np.squeeze(rgb)

    def tostring(self):
        """
        Write out a line containing unit cell information
        """
        return "%f %f %f %f %f %f %s" % (self.lattice_parameters[0],
                                         self.lattice_parameters[1],
                                         self.lattice_parameters[2],
                                         self.lattice_parameters[3],
                                         self.lattice_parameters[4],
                                         self.lattice_parameters[5],
                                         self.symmetry)

    def gethkls_xfab(self, dsmax, spg=None):
        """
        Generate hkl list
        Argument dsmax is the d* limit (eg 1/d)
        Argument spg is the space group name, e.g. 'R3-c'
        """
        stl_max = dsmax / 2.
        if isinstance(self.symmetry, int):
            sgnum = self.symmetry
        else:
            sgnum = None
        # a name carries its own axes (e.g. 'R-3mr'), a number uses those found for the cell
        cell_choice = self.cell_choice if spg is None else "standard"
        raw_peaks = tools.genhkl_all(self.lattice_parameters,
                                     0, stl_max,
                                     sgname=spg,
                                     sgno=sgnum,
                                     cell_choice=cell_choice,
                                     output_stl=True)
        peaks = []
        for i in range(len(raw_peaks)):
            peaks.append([raw_peaks[i, 3] * 2,
                          (int(raw_peaks[i, 0]), int(raw_peaks[i, 1]), int(raw_peaks[i, 2]))])
        # Sort the list: l (index 2) is primary, k (index 1) secondary, h (index 0) tertiary
        peaks.sort(key=lambda x: (x[0], x[1][2], x[1][1], x[1][0]))
        self.peaks = peaks
        self.limit = dsmax
        return peaks

    def gethkls(self, dsmax):
        """
        Generate hkl list
        Argument dsmax is the d* limit (eg 1/d)
        Default of zero gives only the (000) reflection

        assumes [h|k|l] < 200
        """
        if dsmax == self.limit and self.peaks is not None:
            return self.peaks
        if self.sgname is not None:
            return self.gethkls_xfab(dsmax, spg=self.sgname)
        if isinstance(self.symmetry, int):
            return self.gethkls_xfab(dsmax, spg=None)
        h = k = 0
        l = 1  # skip 0,0,0
        hs = ks = ls = 1
        b = 0
        peaks = []
        while abs(h) < 200:  # H
            while abs(k) < 200:  # K
                while abs(l) < 200:  # L
                    ds = self.ds([h, k, l])
                    if ds < dsmax:
                        if not self.absent(h, k, l):
                            peaks.append([ds, (h, k, l)])
                        else:
                            pass
                        b = 0
                    else:
                        if ls == 1:
                            ls = -1
                            l = 0
                        else:
                            ls = 1
                            l = 0
                            b = b + 1
                            break
                    l = l + ls
                k = k + ks
                # l is always zero here
                if b > 1:
                    if ks == 1:
                        ks = -1
                        k = -1
                    else:
                        ks = 1
                        k = 0
                        b = b + 1
                        break
            h = h + hs
            if b > 3:
                if hs == 1:
                    hs = -1
                    h = -1
                else:
                    hs = 1
                    h = 0
                    break

        peaks.sort()

        self.peaks = peaks
        self.limit = dsmax
        return peaks
    
    def ds(self, h):
        """ computes 1/d for this hkl = hgh """
        return math.sqrt(np.dot(h, np.dot(self.gi, h)))  # 1/d or d*

    def makerings(self, limit, tol=0.001):
        """
        Makes a list of computed powder rings
        The tolerance is the difference in d* to decide
        if two peaks overlap
        """
        # Cache if limit and tol are the same
        # this is called each time by assigntorings
        # slow if we are doing pbp indexing
        if (getattr(self, "_rings_for", None) == (limit, tol)
                and getattr(self, "ringds", None) is not None):
            return
        self.peaks = self.gethkls(limit + tol)  # [ ds, [hkl] ]
        if len(self.peaks) == 0:
            raise ValueError(
                "No peaks found for limit=%s; try increasing limit" % (limit,)
            )
        self.ringds = []  # a list of floats
        self.ringhkls = {}  # a dict of lists of integer hkl
        # Append first peak
        peak = self.peaks[0]
        self.ringds.append(peak[0])
        self.ringhkls[peak[0]] = [peak[1]]
        for peak in self.peaks[1:]:
            if abs(peak[0] - self.ringds[-1]) < tol:
                self.ringhkls[self.ringds[-1]].append(peak[1])
            else:
                self.ringds.append(peak[0])
                self.ringhkls[self.ringds[-1]] = [peak[1]]
        self.ringtol = tol
        self._rings_for = (limit, tol)

    def ringtths(self, wavelength):
        """
        Computes the two theta values of powder rings, given a wavelength.
        `makerings` must be called beforehand.
        """
        ring_ds = self.ringds
        if len(ring_ds) == 0:
            raise ValueError('No rings found. Call `makerings` before calling this function.')
        tth = []
        for ds in ring_ds:
            tth.append(2 * np.degrees(np.arcsin(ds * wavelength / 2)))
        return tth

    def anglehkls(self, h1, h2):
        """
        Compute the angle between reciprocal lattice vectors h1, h2
        """
        g1 = np.dot(h1, np.dot(self.gi, h1))
        g2 = np.dot(h2, np.dot(self.gi, h2))
        g12 = np.dot(h1, np.dot(self.gi, h2))
        costheta = g12 / math.sqrt(g1 * g2)
        try:
            return degrees(math.acos(costheta)), costheta
        except:
            if abs(costheta - 1) < 1e-6:
                return 0., 1.0
            if abs(costheta + 1) < 1e-6:
                return 180., -1.0
            print("Error in unit cell class determining angle")
            print("h1", h1, "h2", h2, "Costheta=", costheta)
            raise

    def getanglehkls(self, ring1, ring2):
        """
        Cache the previous pairs called for
        sorted by cos2angle
        """
        if self.ringtol != self.anglehkl_cache['ringtol'] or \
                (self.B != self.anglehkl_cache['B']).any():
            self.anglehkl_cache = {'ringtol': self.ringtol,
                                   'B': self.B,
                                   'BI': np.linalg.inv(self.B)}
        key = (ring1, ring2)
        B = self.anglehkl_cache['B']
        BI = self.anglehkl_cache['BI']
        if key not in self.anglehkl_cache:
            h1 = self.ringhkls[self.ringds[ring1]]
            h2 = self.ringhkls[self.ringds[ring2]]
            cangs = cosangles_many(h1, h2, self.gi)
            val = filter_pairs(h1, h2, cangs, B, BI)
            self.anglehkl_cache[key] = val
        else:
            val = self.anglehkl_cache[key]
        return val
    
    def getcoses(self, ring1, ring2, tol=1e-5):
        """
        Distinct cosines between hkls of ring1 and ring2, descending,
        with 180 degree pairs removed and near-duplicates collapsed.

        Cached per ring pair -- indexer.find rebuilt this on every call, and
        it depends only on the rings and the cell.
        """
        key = ('coses', ring1, ring2, tol)
        if key in self.anglehkl_cache:
            return self.anglehkl_cache[key]
        hkls1 = self.ringhkls[self.ringds[int(ring1)]]
        hkls2 = self.ringhkls[self.ringds[int(ring2)]]
        cosangles = []
        for h1 in hkls1:
            for h2 in hkls2:
                cosangles.append(self.anglehkls(h1, h2)[1])
        cosangles.sort()
        coses = []
        while len(cosangles) > 0:
            a = cosangles.pop()
            if abs(a - 1.0) < tol or abs(a + 1.0) < tol:
                continue
            if len(coses) == 0:
                coses.append(a)
                continue
            if abs(coses[-1] - a) > tol:
                coses.append(a)
        self.anglehkl_cache[key] = coses
        return coses
    
    def orient(self, ring1, g1, ring2, g2, verbose=0, crange=-1.):
        """
        Compute an orientation matrix using cell parameters and the indexing
        of two reflections

        Orientation matrix:
        A matrix such that h = A^1 g
        define 3 vectors t1,t2,t3
        t1 is parallel to first peak (unit vector along g1)
        t2 is in the plane of both   (unit vector along g1x(g1xg2))
        t3 is perpendicular to both  (unit vector along g1xg2)
        """
        # normalized dot product
        costheta = np.dot(g1, g2) / np.sqrt((g1 * g1).sum() * (g2 * g2).sum())
        hab, c2ab, matrs = self.getanglehkls(ring1, ring2)
        if crange > 0:
            best = np.arange(len(c2ab), dtype=int)[abs(c2ab - costheta) < crange]
            if verbose == 1:
                print("possible best", best, len(c2ab))
        else:
            i = np.searchsorted(c2ab, costheta, side='left')
            if i > 0 and (i == len(c2ab) or
                          (fabs(costheta - c2ab[i - 1]) < fabs(costheta - c2ab[i]))):
                best = [i - 1, ]
            else:
                best = [i, ]
        if verbose == 1:
            print("g1, g2", g1, g2)
            print("observed cos2theta", costheta)
            print("hab, c2ab", hab, c2ab)
            print("best", best)
        self.UBIlist = []
        UBlist = []
        for b in best:
            h1, h2 = hab[b]
            if verbose == 1:
                print("Assigning h1", h1, g1, self.ds(h1), \
                      math.sqrt(np.dot(g1, g1)), \
                      self.ds(h1) - math.sqrt(np.dot(g1, g1)))
                print("Assigning h2", h2, g2, self.ds(h2), \
                      math.sqrt(np.dot(g2, g2)), \
                      self.ds(h1) - math.sqrt(np.dot(g1, g1)))
                print("Cos angle calc", self.anglehkls(h1, h2),
                      "obs", costheta, "c2ab", c2ab[b])
            BT = matrs[b]
            UBI = np.empty((3, 3), float)
            UBI[0] = g1
            UBI[1] = g2
            cImageD11.quickorient(UBI, BT)
            if verbose == 1:
                print("UBI")
                print(UBI)
                h = np.dot(UBI, g1)
                print("(%9.3f, %9.3f, %9.3f)" % (h[0], h[1], h[2]))
                h = np.dot(UBI, g2)
                print("(%9.3f, %9.3f, %9.3f)" % (h[0], h[1], h[2]))
            self.UBI = UBI
            self.UB = _inv3(UBI)
            self.UBIlist.append(UBI)
            UBlist.append(self.UB)
        # trim to uniq list? What about small distortions...
        self.UBIlist = ubi_equiv(self.UBIlist, UBlist)
    
    @classmethod
    def from_cif(cls, filename, name=None):
        """Import a unitcell object from a CIF file"""
        import xfab.structure, xfab.sg
        o = xfab.structure.build_atomlist()
        o.CIFread(filename)
        atomlist = o.atomlist
        lattice_parameters = atomlist.cell
        spacegroup = xfab.sg.sg(sgname=atomlist.sgname).no
        return cls(lattice_parameters, spacegroup, name=name)
    
    @classmethod
    def from_pars(cls, pars):
        """
        Produce a unitcell from a parameter dict/object
        """
        return unitcell_from_parameters(pars)
    
    @classmethod
    def from_par_file(cls, filename):
        """
        Produce a unitcell from a .par file
        """
        par_obj = parameters.from_file(filename)
        return cls.from_pars(par_obj)
    
    def to_par_dict(self):
        """
        Return a parameter dict
        """
        parnames = ["cell_" + name for name in "_a _b _c alpha beta gamma lattice_[P,A,B,C,I,F,R]".split()]
        pars_dict = {}
        for inc, parname in enumerate(parnames):
            if inc < 6:
                pars_dict[parname] = self.lattice_parameters[inc]
            else:
                pars_dict[parname] = self.symmetry
        if self.sgname is not None:
            pars_dict['cell_sg'] = self.sgname

        if hasattr(self, 'name') and self.name is not None:
            pars_dict['phase_name'] = self.name
        return pars_dict
    
    def to_par_obj(self):
        """
        Return an ImageD11.parameters.parameters object
        """
        pars_dict = self.to_par_dict()
        pars_obj = parameters.from_dict(pars_dict)
        return pars_obj
    
    def to_par_file(self, filename):
        """
        Write lattice parameters to file as an ImageD11 .par
        """
        pars_obj = self.to_par_obj()
        pars_obj.saveparameters(filename)

    
def BTmat(h1, h2, B, BI):
    """ used for computing orientations
    """
    g1 = np.dot(B, h1)  # gvectors for these hkl
    g2 = np.dot(B, h2)
    g3 = np.cross(g1, g2)
    u1 = unit(g1)  # normalised
    u3 = unit(g3)
    u2 = np.cross(u1, u3)
    BT = np.dot(BI, np.transpose((u1, u2, u3)))
    return BT


HKL0 = np.array([[0, 0, 1, 1, -1, 1, -1, 0, 0, 1, -1, 1, 1, 3, 11],
                 [0, 1, 0, 1, 1, 0, 0, 1, -1, 1, 1, -1, 1, 2, 12],
                 [1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, -1, 1, 13]], float)  # first unit cell


def filter_pairs(h1, h2, c2a, B, BI, tol=1e-5):
    """ remove duplicate pairs for orientation searches
    h1 = reflections of ring1, N1 peaks
    h2 = reflections of ring2, N2 peaks
    c2a  = cos angle between them, N1xN2
    B = B matrix in reciprocal space
    BI = inverse in real space
    """
    assert c2a.shape == (len(h1), len(h2))
    order = np.argsort(c2a.ravel())  # increasing in cosine of angle
    c2as = c2a.flat[order]
    hi, hj = np.mgrid[0:len(h1), 0:len(h2)]
    hi = hi.ravel()[order]  # to get back the peaks
    hj = hj.ravel()[order]
    # Results holders:
    pairs = []
    cangs = []
    matrs = []
    # cluster the cangs assuming a sensible threshold
    dc = (c2as[1:] - c2as[:-1]) > 1e-8  # differences
    inds = list(np.arange(1, len(dc) + 1, dtype=int)[dc]) + [len(c2as) - 1, ]
    p = 0  # previous
    for i in inds:
        c = c2as[p:i]  # block is p:i
        if abs(c2as[p]) < 0.98:  # always keep the first one
            ha = h1[hi[p]]
            hb = h2[hj[p]]
            pairs.append((ha, hb))
            cangs.append(c2as[p])
            BT = BTmat(ha, hb, B, BI)
            matrs.append(BT)
        else:
            p = i
            continue
        if len(c) == 1:
            p = i
            continue
        assert (c.max() - c.min()) < 2.1e-8, "Angles blocking error in filter_pairs"
        # here we have a bunch of hkl pairs which all give the same angle
        # between them. They are not all the same. We generate a pair of peaks
        # from the first one and see which other pairs index differently
        ga = np.dot(B, ha)
        gb = np.dot(B, hb)
        assert abs(np.dot(ga, gb) / np.sqrt(np.dot(ga, ga) * np.dot(gb, gb)) - c2as[p]) < 2e-8, "mixup in filter_pairs"
        gobs = np.array((ga, gb, (0, 0, 0)), float)
        UBI = gobs.copy()
        cImageD11.quickorient(UBI, BT)
        gtest = [np.dot(np.linalg.inv(UBI), HKL0).T.copy(), ]
        for j in range(p + 1, i):
            ha = h1[hi[j]]
            hb = h2[hj[j]]
            BT = BTmat(ha, hb, B, BI)
            newpair = True
            for gt in gtest:
                UBI = gobs.copy()
                cImageD11.quickorient(UBI, BT)
                npk = cImageD11.score(UBI, gt, 1e-6)
                if npk == len(HKL0[0]):
                    newpair = False
                    break
            if newpair:
                pairs.append((ha, hb))
                cangs.append(c2as[j])
                matrs.append(BT)
                gtest.append(np.dot(np.linalg.inv(UBI), HKL0).T.copy())
        p = i
    return pairs, cangs, matrs


def ubi_equiv(ubilist, ublist, tol=1e-8):
    """ Two ubi are considered equivalent if they both index the peaks
    in the HKL0 array exactly"""
    if len(ubilist) < 2:
        return ubilist
    order = np.argsort([np.trace(ubi) for ubi in ubilist])  # low to high
    uniq = [ubilist[order[-1]], ]
    refgv = [np.dot(ublist[order[-1]], HKL0), ]
    for i in order[:-1][::-1]:
        ubi = ubilist[i]
        score = 1
        for pks in refgv:  # pks is (n,3) for an orientation
            hcalc = np.dot(ubi, pks)
            score = min(score, np.abs(np.round(hcalc) - hcalc).sum())
            if score <= tol:  # matches a previous grain
                break
        if score > tol:  # is a new orientation
            uniq.append(ubi)
            refgv.append(np.dot(np.linalg.inv(ubi), HKL0))
    return uniq


def unitcell_from_parameters(pars):
    parnames = "_a _b _c alpha beta gamma".split()
    try:
        sgname = pars.get("cell_sg")  # optional
    except KeyError:
        sgname = None
    cell = unitcell([pars.get("cell_%s" % (s)) for s in parnames],
                    pars.get("cell_lattice_[P,A,B,C,I,F,R]"),
                    sgname=sgname)

    return cell


if __name__ == "__main__":
    import sys, time

    start = time.time()
    cell = unitcell([float(x) for x in sys.argv[1:7]], sys.argv[7])
    cell.makerings(2)
    cell.getanglehkls(11, 12)
