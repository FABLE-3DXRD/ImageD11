from __future__ import print_function, division


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


import numpy as np
from . import cImageD11, unitcell
from ImageD11.grain import grain
from xfab.tools import ubi_to_u, u_to_rod, ubi_to_rod

import math, time, sys

loglevel = 1


class clogging(
    object
):  # because multiprocessing. FIXME. object level logging rather than module ?
    def log(self, *args):
        print(" ".join(str(a) for a in args))

    def debug(self, *args):
        if loglevel <= 0:
            logging.log("debug:", *args)

    def info(self, *args):
        if loglevel <= 1:
            logging.log("info:", *args)

    def warning(self, *args):
        if loglevel <= 2:
            logging.log("warning:", *args)

    def error(self, *args):
        if loglevel <= 3:
            logging.log("error:", *args)


logging = clogging()


def ubi_fit_2pks(ubi, g1, g2):
    """
    Refine a ubi matrix so it matches the pair of g-vectors supplied
    (almost) accounts for cell parameters not being quite right
    """
    ub = np.linalg.inv(ubi)
    h1 = np.round(np.dot(ubi, g1))
    h2 = np.round(np.dot(ubi, g2))
    g3 = np.cross(g1, g2)
    h3 = np.dot(ubi, g3)  # do not round to integer
    g1c = np.dot(ub, h1)
    g2c = np.dot(ub, h2)
    g3c = np.dot(ub, h3)
    R = np.outer(g1, h1) + np.outer(g2, h2) + np.outer(g3, h3)
    H = np.outer(h1, h1) + np.outer(h2, h2) + np.outer(h3, h3)
    ubfit = np.dot(R, np.linalg.inv(H))
    ubifit = np.linalg.inv(ubfit)
    return ubifit


def myhistogram(data, bins):
    """
    The numpy histogram api was changed
    So here is an api that will not change
    It is based on that from the old Numeric manual
    """
    n = np.searchsorted(np.sort(data), bins)
    n = np.concatenate([n, [len(data)]])
    return n[1:] - n[:-1]


def readubis(ubifile):
    """read ubifile and return a list of ubi arrays"""
    f = open(ubifile, "r")
    ubisread = []
    u = []
    for line in f:
        if line[0] == "#":
            continue
        vals = [float(x) for x in line.split()]
        if len(vals) == 3:
            u = u + [vals]
        if len(u) == 3:
            ubisread.append(np.array(u))
            u = []
    f.close()
    return ubisread


def write_ubi_file(filename, ubilist):
    """save 3x3 matrices into file"""
    f = open(filename, "w")
    for u in ubilist:
        f.write("%f %f %f\n" % (u[0][0], u[0][1], u[0][2]))
        f.write("%f %f %f\n" % (u[1][0], u[1][1], u[1][2]))
        f.write("%f %f %f\n\n" % (u[2][0], u[2][1], u[2][2]))
    f.close()


def ubitocellpars(ubi):
    """convert ubi matrix to unit cell"""
    g = np.dot(ubi, np.transpose(ubi))
    from math import acos, degrees, sqrt

    a = sqrt(g[0, 0])
    b = sqrt(g[1, 1])
    c = sqrt(g[2, 2])
    alpha = degrees(acos(g[1, 2] / b / c))
    beta = degrees(acos(g[0, 2] / a / c))
    gamma = degrees(acos(g[0, 1] / a / b))
    return a, b, c, alpha, beta, gamma


def ubitoU(ubi):
    """
    convert ubi to Grainspotter style U
    The convention is B as being triangular, hopefully as Busing and Levy
    TODO - make some testcases please!!
    """
    # return np.transpose(np.dot(ubitoB(ubi),ubi))
    return ubi_to_u(ubi)


def ubitoRod(ubi):
    """
    TODO Testcases!!!
    """
    #     u = ubitoU(ubi)
    #     w, v = np.linalg.eig(u)
    #     print 'Eigenvalues'
    #     print w
    #     print 'Eigen vectors'
    #     print v
    #     #ehat = v[:,0]
    #     #angle = -1*math.acos(np.clip(w[order[1]].real,-1,1))
    #     order = np.argsort(w.real)
    #     print order
    #     ehat = v[:, order[-1]]
    #     if order.tolist() != range(3):
    #         print 'HHFH'
    #         angle = -1*np.arccos(w[order[1]].real)
    #     else:
    #         angle = np.arccos(w[order[1]].real)
    #     Rod = ehat * math.tan(angle/2)
    #     return Rod.real
    return ubi_to_rod(ubi)


def ubitoB(ubi):
    """give the B matrix from ubi"""
    g = np.dot(ubi, np.transpose(ubi))
    return np.transpose(np.linalg.inv(np.linalg.cholesky(g)))


def mod_360(theta, target):
    """
    Find multiple of 360 to add to theta to be closest to target
    """
    diff = theta - target
    while diff < -180:
        theta = theta + 360
        diff = theta - target
    while diff > 180:
        theta = theta - 360
        diff = theta - target
    return theta


def calc_drlv2(UBI, gv):
    """
    Get the difference from being integer hkls
    UBI: ubi matrix
    gv: list of g-vectors
    returns drlv2 = (h_calc - h_int)^2
    """
    h = np.dot(UBI, np.transpose(gv))
    hint = np.floor(h + 0.5).astype(int)  # rounds down
    diff = h - hint
    drlv2 = np.sum(diff * diff, 0)
    return drlv2


def refine(UBI, gv, tol, quiet=True):
    """
    Refine an orientation matrix and rescore it.

    From Paciorek et al Acta A55 543 (1999)
       UB = R H-1
    where:
       R = sum_n r_n h_n^t
       H = sum_n h_n h_n^t
       r = g-vectors
       h = hkl indices
    """
    #      print "Orientation and unit cell refinement of"
    #      print "UBI\n",UBI
    #      print "Scores before",self.score(UBI)
    # Need to find hkl indices for all of the peaks which are indexed
    h = np.dot(UBI, np.transpose(gv))
    hint = np.floor(h + 0.5).astype(int)  # rounds down
    diff = h - hint
    drlv2 = np.sum(diff * diff, 0)
    tol = float(tol)
    tol = tol * tol
    # Only use peaks which are assigned to rings for refinement
    ind = np.compress(np.less(drlv2, tol), np.arange(gv.shape[0]))
    # scoreb4=ind.shape[0]
    contribs = drlv2[ind]
    try:
        fitb4 = math.sqrt(np.sum(contribs) / contribs.shape[0])
        if not quiet:
            logging.debug("Fit before refinement %.8f %5d" % (fitb4, contribs.shape[0]))
    except:
        logging.error("No contributing reflections for \n%s" % (str(UBI)))
        raise
    # drlv2_old=drlv2
    R = np.zeros((3, 3), float)
    H = np.zeros((3, 3), float)
    for i in ind:
        r = gv[i, :]
        k = hint[:, i].astype(float)
        #           print r,k
        R = R + np.outer(r, k)
        H = H + np.outer(k, k)
    try:
        HI = np.linalg.inv(H)
        UBoptimal = np.dot(R, HI)
        UBIo = np.linalg.inv(UBoptimal)
    except:
        # A singular matrix - this sucks.
        UBIo = UBI
    h = np.dot(UBIo, np.transpose(gv))
    hint = np.floor(h + 0.5).astype(int)  # rounds down
    diff = h - hint
    drlv2 = np.sum(diff * diff, 0)
    ind = np.compress(np.less(drlv2, tol), np.arange(gv.shape[0]))
    # scorelastrefined=ind.shape[0]
    contribs = drlv2[ind]
    if contribs.size != 0:
        fitlastrefined = math.sqrt(np.sum(contribs) / contribs.shape[0])
        if not quiet:
            logging.debug("after %.8f %5d" % (fitlastrefined, contribs.shape[0]))
    else:
        logging.error("\n\n\n")
        logging.error("No contributing reflections for %s\n" % (str(UBI)))
        logging.error("After refinement, it was OK before ???")
        logging.error("\n\n\n")
        return UBI
        # raise
    #      for i in ind:
    #         print "( %-6.4f %-6.4f %-6.4f ) %12.8f %12.8f"%(\
    #            h[0,i],h[1,i],h[2,i],sqrt(drlv2[i]),sqrt(drlv2_old[i]))
    #      print UBIo
    #      print "Scores after", self.score(UBIo,self.hkl_tol)
    #      print "diff\n",UBI-UBIo
    #      print "Mean drlv now",sum(sqrt(drlv2))/drlv2.shape[0],
    #      print "Mean drlv old",sum(sqrt(drlv2_old))/drlv2_old.shape[0]
    return UBIo


def indexer_from_colfile(colfile, **kwds):
    uc = unitcell.unitcell_from_parameters(colfile.parameters)
    w = float(colfile.parameters.get("wavelength"))
    gv = np.array((colfile.gx, colfile.gy, colfile.gz), float)
    kwds.update({"unitcell": uc, "wavelength": w, "gv": gv.T})
    ind = indexer(**kwds)
    if "omega" in colfile.titles:
        ind.omega_fullrange = find_omega_ranges(colfile.omega)
    ind.colfile = colfile
    return ind



def indexer_from_colfile_and_ucell(colfile, ucell, **kwds):
    """Force a specific unitcell to be used, 
    ensuring a wavelength is provided."""

    # Try to get the wavelength from kwds, then colfile.parameters, otherwise raise an error
    w = kwds.get("wavelength")
    if w is None:
        try:
            w = colfile.parameters.get("wavelength")
        except KeyError:
            raise KeyError("Wavelength must be provided either in kwds or colfile.parameters")

    gv = np.array((colfile.gx, colfile.gy, colfile.gz), float)
    kwds.update({"unitcell": ucell, "wavelength": float(w), "gv": gv.T})
    ind = indexer(**kwds)
    if "omega" in colfile.titles:
        ind.omega_fullrange = find_omega_ranges(colfile.omega)
    ind.colfile = colfile
    return ind

class indexer:
    """
    A class for searching for orientation matrices
    """

    def __init__(
        self,
        unitcell=None,
        gv=None,
        cosine_tol=0.002,
        minpks=10,
        hkl_tol=0.01,
        ring_1=1,
        ring_2=2,
        ds_tol=0.005,
        wavelength=-1,
        uniqueness=0.5,
        eta_range=0.0,
        max_grains=100,
    ):
        """
        Unitcell would be a unitcell object for generating hkls peaks
        gv would be a 3*n array of points in reciprocal space
        The rest of the arguments are parameters.
        """
        # This stop variable allows computation to be run in a thread...
        self.stop = False
        self.unitcell = unitcell
        self.gv = gv
        self.ra = None
        if gv is not None:  # do init
            logging.info("gv: %s %s %s" % (str(gv), str(gv.shape), str(gv.dtype)))
            assert gv.shape[1] == 3
            self.gv = gv.astype(float)
            self.ds = np.sqrt((gv * gv).sum(axis=1))
            self.ga = np.zeros(len(self.ds), np.int32) - 1  # Grain assignments
            self.gvflat = np.ascontiguousarray(gv, float)
            self.wedge = 0.0  # Default

        self.cosine_tol = cosine_tol
        self.wavelength = wavelength
        self.hkl_tol = hkl_tol
        self.ring_1 = ring_1
        self.ring_2 = ring_2
        self.uniqueness = uniqueness
        self.minpks = minpks
        self.ds_tol = ds_tol
        self.max_grains = max_grains
        self.eta_range = eta_range
        self.ubis = []
        self.scores = []
        self.hits = []
        self.omega = None
        self.colfile = None
        self.index_needs_debug = 0  # track problems quietly...
        self.omega_fullrange = 0  # flags having no idea
        # it would make more sense to inherit the parameter object - will
        # have to think about this some more - how general is it?
        from ImageD11 import parameters

        self.parameterobj = parameters.parameters(
            cosine_tol=self.cosine_tol,
            hkl_tol=self.hkl_tol,
            ring_1=self.ring_1,
            ring_2=self.ring_2,
            minpks=self.minpks,
            uniqueness=self.uniqueness,
            ds_tol=self.ds_tol,
            wavelength=self.wavelength,
            eta_range=self.eta_range,
            max_grains=self.max_grains,
        )

        # Add a resetting functionality, adapted from
        # stackoverflow.com/questions/4866587/pythonic-way-to-reset-an-objects-variables
        import copy

        self.__pristine_dict = copy.deepcopy(self.__dict__)

    def __getattr__(self, name):
        """got some time lost setting tol which does not exist

        this will never be clean :-(
        """
        if name == "tol":
            raise KeyError("tol not in indexer")
        print("WARNING: creating indexer.%s" % (name))
        setattr(self, name, None)

    def reset(self):
        """
        To get a really clean indexer just create a new one (e.g. via __init__)
        This was added for the gui to help it forget what happened before
        but keep parameters as they were set
        """
        import copy

        self.__dict__ = copy.deepcopy(self.__pristine_dict)
        self.__pristine_dict = copy.deepcopy(self.__dict__)

    def loadpars(self, filename=None):
        if filename is not None:
            self.parameterobj.loadparameters(filename)
        # self.parameterobj.update_other(self) # busted CI for logging in windows + py2.7
        for parname in self.parameterobj.parameters:
            if hasattr(self, parname):
                setattr(self, parname, self.parameterobj.get(parname))

    def updateparameters(self):
        self.savepars()
        self.pars = self.parameterobj.parameters

    def savepars(self, filename=None):
        self.parameterobj.update_yourself(self)
        if filename is not None:
            self.parameterobj.saveparameters(filename)

    def out_of_eta_range(self, eta):
        """decide if an eta is going to be kept"""
        e = mod_360(float(eta), 0)
        if e < abs(self.eta_range) and e > -abs(self.eta_range):
            return True
        if e < -180.0 + abs(self.eta_range) or e > 180.0 - abs(self.eta_range):
            return True
        return False

    def assigntorings(self):
        """
        Assign the g-vectors to hkl rings
        """
        # rings are in self.unitcell
        limit = np.amax(self.ds)
        logging.info("Assign to rings, maximum d-spacing considered: %f" % (limit))
        self.unitcell.makerings(limit, tol=self.ds_tol)
        dsr = self.unitcell.ringds
        # npks
        npks = len(self.ds)
        self.ra = np.zeros(npks, np.int32) - 1
        self.na = np.zeros(len(dsr), np.int32)
        logging.info("Ring assignment array shape", self.ra.shape)
        tol = float(self.ds_tol)
        best = np.zeros(npks, float) + tol
        for j, dscalc in enumerate(dsr):
            dserr = abs(self.ds - dscalc)
            sel = dserr < best
            self.ra[sel] = j
            best[sel] = dserr[sel]
        # Report on assignments
        ds = np.array(self.ds)
        logging.info(
            "Ring     (  h,  k,  l) Mult  total indexed to_index  ubis  peaks_per_ubi   tth"
        )
        minpks = 0
        # try reverse order instead
        for j in range(len(dsr))[::-1]:
            ind = np.compress(np.equal(self.ra, j), np.arange(self.ra.shape[0]))
            self.na[j] = ind.shape[0]
            n_indexed = np.sum(np.where(self.ga[ind] > -1, 1, 0))
            n_to_index = np.sum(np.where(self.ga[ind] == -1, 1, 0))
            # diffs = abs(take(ds,ind) - dsr[j])
            h = self.unitcell.ringhkls[dsr[j]][0]
            Mult = len(self.unitcell.ringhkls[dsr[j]])
            if self.omega_fullrange > 0:
                expected_orients = int(
                    180.0 / self.omega_fullrange * self.na[j] / float(Mult)
                )
                expected_npks = int(self.omega_fullrange / 180.0 * Mult)
                minpks += expected_npks
            else:
                expected_orients = "N/A"
                expected_npks = "N/A"
            tth = 2 * np.degrees(np.arcsin(dsr[j] * self.wavelength / 2))
            logging.info(
                "Ring %-3d (%3d,%3d,%3d)  %3d  %5d   %5d    %5d %5s     %2s  %.2f"
                % (
                    j,
                    h[0],
                    h[1],
                    h[2],
                    Mult,
                    self.na[j],
                    n_indexed,
                    n_to_index,
                    expected_orients,
                    expected_npks,
                    tth,
                )
            )
        if minpks > 0:
            logging.info("\nmin_pks:  - Current  --> %3d" % (self.minpks))
            logging.info("          - Expected --> %3d\n" % (minpks))

        # We will only attempt to index g-vectors which have been assigned
        # to hkl rings (this gives a speedup if there
        # are a lot of spare peaks
        ind = np.compress(np.greater(self.ra, -1), np.arange(self.ra.shape[0]))
        self.gvr = self.gv[ind]
        logging.info(
            "Using only those peaks which are assigned to rings for scoring trial matrices"
        )
        logging.info("Shape of scoring matrix", self.gvr.shape)
        self.gvflat = np.ascontiguousarray(self.gvr, float)  # Makes it contiguous
        self.gv = np.ascontiguousarray(self.gv, float)
        # in memory, hkl fast index

    def friedelpairs(self, filename):
        """
        Attempt to identify Freidel pairs

        Peaks must be assigned to the same powder ring
        Peaks will be the closest thing to being 180 degrees apart
        """
        out = open(filename, "w")
        dsr = self.unitcell.ringds
        nring = len(dsr)
        for j in range(nring):
            ind = np.compress(np.equal(self.ra, j), np.arange(self.ra.shape[0]))
            # ind is the indices of the ring assigment array - eg which hkl is this gv
            #
            if len(ind) == 0:
                continue
            thesepeaks = self.gv[ind]
            #
            h = self.unitcell.ringhkls[dsr[j]][0]
            #
            out.write("\n\n\n# h = %d \n" % (h[0]))
            out.write("# k = %d \n" % (h[1]))
            out.write("# l = %d \n" % (h[2]))
            out.write("# npks = %d \n" % (thesepeaks.shape[0]))
            out.write(
                "# score eta1 omega1 tth1 gv1_x gv1_y gv1_z eta2 omega2 tth2 gv2_x gv2_y gv2_z\n"
            )
            for k in range(thesepeaks.shape[0]):
                nearlyzero = thesepeaks + thesepeaks[k]
                mag = np.sum(nearlyzero * nearlyzero, 1)
                best = np.argmin(mag)
                if best > k:
                    a = ind[k]
                    out.write("%f " % (np.sqrt(mag[best])))
                    out.write(
                        "%f %f %f %f %f %f    "
                        % (
                            self.eta[a],
                            self.omega[a],
                            self.tth[a],
                            self.gv[a][0],
                            self.gv[a][1],
                            self.gv[a][2],
                        )
                    )
                    a = ind[best]
                    out.write(
                        "%f %f %f %f %f %f    "
                        % (
                            self.eta[a],
                            self.omega[a],
                            self.tth[a],
                            self.gv[a][0],
                            self.gv[a][1],
                            self.gv[a][2],
                        )
                    )
                    out.write("\n")

    def score_all_pairs(self, n=None, rmulmax=None, rings_to_use=None):
        """
        Generate all the potential pairs of rings and go score them too

        n = maximum number of pairs to try
        maxmult = max multiplicity of rings to use for generating pairs
        rings_to_use = the rings to use for generating ubis
        """
        self.assigntorings()
        if rings_to_use is not None:
            rings = [r for r in rings_to_use if (self.ra == r).sum() > 0]
        else:
            # Which rings have peaks assigned to them?
            rings = [r for r in set(self.ra) if r >= 0]
        # What are the multiplicities of these rings? We will use low multiplicity first
        mults = {r: len(self.unitcell.ringhkls[self.unitcell.ringds[r]]) for r in rings}
        if rmulmax is not None:
            rings = [r for r in rings if mults[r] <= rmulmax]

        # How many peaks per ring? We will use the sparse rings first...
        # why? We assume these are the strongest peaks on a weak high angle ring
        # occupation = {r:self.na[r] for r in rings}
        pairs = [
            (int(mults[r1] * mults[r2]), int(self.na[r1] * self.na[r2]), r1, r2)
            for r1 in rings
            for r2 in rings
        ]
        pairs.sort()
        self.tried = 0
        self.npairs = len(pairs)
        self.stop = False
        k = 0
        for mu, oc, r1, r2 in pairs:
            k += 1
            try:
                self.ring_1 = r1
                self.ring_2 = r2
                self.find()
                if len(self.hits) == 0:  # skip when nothing is found
                    continue
                self.scorethem()
                self.tried += 1
            except KeyboardInterrupt:
                break
            if self.stop:
                break
            if n is not None and k > n:
                break
            logging.info(
                "Tried r1=%d r2=%d attempt %d of %d, got %d grains"
                % (r1, r2, self.tried, len(pairs), len(self.ubis))
            )
        logging.info(
            "\nTested", self.tried, "pairs and found", len(self.ubis), "grains so far"
        )

    def find(self):
        """
        Dig out the potential hits
        """
        # Optionally only used unindexed peaks here. Make this obligatory
        # Need indices of gvectors to test.
        # Bug out early when there are none
        if self.ra is None:
            self.assigntorings()
        iall = np.arange(self.gv.shape[0])
        i1 = np.compress(
            np.logical_and(np.equal(self.ra, self.ring_1), self.ga == -1), iall
        ).tolist()
        i2 = np.compress(
            np.logical_and(np.equal(self.ra, self.ring_2), self.ga == -1), iall
        ).tolist()
        if len(i1) == 0 or len(i2) == 0:
            logging.info("no peaks left for those rings")
            return
        # Which are the rings being used for indexing
        hkls1 = self.unitcell.ringhkls[self.unitcell.ringds[int(self.ring_1)]]
        hkls2 = self.unitcell.ringhkls[self.unitcell.ringds[int(self.ring_2)]]
        logging.info("hkls of rings being used for indexing")
        logging.info("Ring 1: %s" % (str(hkls1)))
        logging.info("Ring 2: %s" % (str(hkls2)))
        cosangles = []
        for h1 in hkls1:
            for h2 in hkls2:
                ca = self.unitcell.anglehkls(h1, h2)
                cosangles.append(ca[1])
        cosangles.sort()
        coses = []
        while len(cosangles) > 0:
            a = cosangles.pop()
            if (
                abs(a - 1.0) < 1e-5 or abs(a + 1.0) < 1e-5
            ):  # Throw out 180 degree angles
                continue
            if len(coses) == 0:
                coses.append(a)
                continue
            if abs(coses[-1] - a) > 1e-5:
                coses.append(a)
        logging.info("Possible angles and cosines between peaks in rings:")
        for c in coses:
            logging.info("%.6f %.6f" % (math.acos(c) * 180 / math.pi, c))
        #
        #
        logging.info("Number of peaks in ring 1: %d" % (len(i1)))
        logging.info("Number of peaks in ring 2: %d" % (len(i2)))
        logging.info("Minimum number of peaks to identify a grain %d" % (self.minpks))
        # print self.gv.shape
        # ntry=0
        # nhits=0
        self.hits = []
        if len(i1) == 0 or len(i2) == 0:
            # return without crashing please
            return
        tol = float(self.cosine_tol)
        # ng=0
        mp = np.sqrt(np.sum(self.gv * self.gv, 1))
        # print mp.shape
        ps1 = np.take(self.gv, i1, 0)
        mp1 = np.take(mp, i1, 0)
        n1 = ps1.copy()
        ps2 = np.take(self.gv, i2, 0)
        mp2 = np.take(mp, i2, 0)
        n2 = ps2.copy()
        # print "mp1.shape",mp1.shape
        # print "n1[:,1].shape",n1[:,1].shape
        for i in range(3):
            n1[:, i] = n1[:, i] / mp1
            n2[:, i] = n2[:, i] / mp2
        cs = np.array(coses, "d")
        # found=0
        hits = []
        start = time.time()
        self.cosangles = cs
        mtol = -tol  # Ugly interface - set cosine tolerance negative for all
        # instead of best
        for i in range(len(i1)):
            costheta = np.dot(n2, n1[i])
            if tol > 0:  # This is the original algorithm - the closest angle
                best, diff = cImageD11.closest(costheta, cs)
                if diff < tol:
                    hits.append([diff, i1[i], i2[best]])
            else:
                for cval in cs:
                    # 1d   scalar  1d
                    diff = cval - costheta
                    candidates = np.compress(abs(diff) < mtol, i2)
                    for c in candidates:
                        hits.append([0.0, i1[i], c])
        logging.info("Number of trial orientations generated %d" % (len(hits)))
        logging.info("Time taken %.6f /s" % (time.time() - start))
        self.hits = hits

    def histogram_drlv_fit(self, UBI=None, bins=None):
        """
        Generate a histogram of |drlv| for a ubi matrix
        For use in validation of grains
        """
        if UBI is None:
            ubilist = self.ubis
        else:
            ubilist = [UBI]
        if bins is None:
            start = 0.25
            fac = 2
            bins = [start]
            while start > 1e-5:
                start = start / fac
                bins.append(start)
            bins.append(-start)
            bins.reverse()
            bins = np.array(bins)
        hist = np.zeros((len(ubilist), bins.shape[0] - 1), int)
        j = 0
        for UBI in ubilist:
            drlv2 = calc_drlv2(UBI, self.gv)
            drlv = np.sort(np.sqrt(drlv2))  # always +ve
            if drlv[-1] > 0.866:
                print("drlv of greater than 0.866!!!", drlv[-1])
            positions = np.searchsorted(drlv, bins)
            hist[j, :] = positions[1:] - positions[:-1]
            j = j + 1
        self.bins = bins
        self.histogram = hist

    def scorethem(self, fitb4=False):
        """decide which trials listed in hits to keep"""
        if self.hits is None or len(self.hits) == 0:  # no idea how this can be None?
            logging.info("No hits to score")
            return
        start = time.time()
        ng = 0
        tol = float(self.hkl_tol)
        gv = self.gvflat
        all = len(self.hits)
        logging.info("Scoring %d potential orientations" % (all))
        progress = 0
        nuniq = 0
        # for getind mallocs
        drlv2tmp = np.empty(len(self.gv), float)
        labelstmp = np.empty(len(self.gv), np.int32)
        while len(self.hits) > 0 and ng < self.max_grains:
            diff, i, j = self.hits.pop()
            if self.ga[i] > -1 or self.ga[j] > -1 or i == j:
                # skip things which are already assigned or errors
                continue
            try:
                self.unitcell.orient(
                    self.ring_1, self.gv[i, :], self.ring_2, self.gv[j, :], verbose=0
                )
            except:
                logging.error(
                    " ".join([str(x) for x in (i, j, self.ring_1, self.ring_2)])
                )
                logging.error(str(self.gv[i]))
                logging.error(str(self.gv[j]))
                logging.error("Failed to find orientation in unitcell.orient")
                raise
            if fitb4:  # FIXME : this does not work
                self.unitcell.UBI = ubi_fit_2pks(
                    self.unitell.UBI, self.gv[i, :], self.gv[j, :]
                )
            # npk = cImageD11.score(self.unitcell.UBI,gv,tol)
            npk = self.score(self.unitcell.UBI, tol)
            UBI = self.unitcell.UBI.copy()
            if npk > self.minpks:
                # Try to get a better orientation if we can...:
                self.unitcell.orient(
                    self.ring_1,
                    self.gv[i, :],
                    self.ring_2,
                    self.gv[j, :],
                    verbose=0,
                    crange=abs(self.cosine_tol),
                )
                if fitb4:
                    for k in range(len(self.unitell.UBIlist)):
                        self.unitell.UBIlist[k] = ubi_fit_2pks(
                            self.unitell.UBIlist[k], self.gv[i, :], self.gv[j, :]
                        )
                if len(self.unitcell.UBIlist) > 1:
                    npks = [
                        self.score(UBItest, tol) for UBItest in self.unitcell.UBIlist
                    ]
                    choice = np.argmax(npks)
                    if npks[choice] >= npk:
                        UBI = self.unitcell.UBIlist[choice].copy()
                        npk = npks[choice]
                _ = cImageD11.score_and_refine(UBI, gv, tol)
                # See if we already have this grain...
                try:
                    ind = self.getind(
                        UBI,
                        drlv2tmp=drlv2tmp,
                        labelstmp=labelstmp,
                    )  # indices of peaks indexed
                    ga = self.ga[ind]  # previous grain assignments
                    uniqueness = np.sum(np.where(ga == -1, 1, 0)) * 1.0 / ga.shape[0]
                    if uniqueness > self.uniqueness:
                        self.ga[ind] = len(self.scores) + 1
                        self.ubis.append(UBI)
                        self.scores.append(npk)
                        ubistr = (" %.6f" * 9) % tuple(UBI.ravel())
                        logging.info(
                            "new grain %d pks, i %d j %d UBI %s" % (npk, i, j, ubistr)
                        )
                        ng = ng + 1
                    else:
                        nuniq = nuniq + 1
                except:
                    raise

        logging.info(
            "Number of orientations with more than %d peaks is %d"
            % (self.minpks, len(self.ubis))
        )
        logging.info("Time taken %.3f/s" % (time.time() - start))
        if len(self.ubis) > 0:
            bestfitting = np.argmax(self.scores)
            logging.info("UBI for best fitting\n%s" % (str(self.ubis[bestfitting])))
            logging.info(
                "Unit cell: %s\n" % (str(ubitocellpars(self.ubis[bestfitting])))
            )
            self.refine(self.ubis[bestfitting])
            logging.info(
                "Indexes %d peaks, with <drlv2>=%f"
                % (self.scorelastrefined, self.fitlastrefined)
            )
            logging.info("That was the best thing I found so far")
            notaccountedfor = ((self.ga < 0) & (self.ra >= 0)).sum()
            logging.info(
                "Number of peaks assigned to rings but not indexed = %d"
                % (notaccountedfor)
            )
        else:
            logging.info(
                "Try again, either with larger tolerance or fewer minimum peaks"
            )

    def fight_over_peaks(self):
        """
        Get the best ubis from those proposed
        Use all peaks (ring assigned or not)
        """
        self.drlv2 = np.zeros(self.gv.shape[0], float) + 2
        labels = np.ones(self.gv.shape[0], np.int32)
        np.subtract(labels, 2, labels)
        i = -1
        for ubi in self.ubis:
            i += 1
            try:
                npk = cImageD11.score_and_assign(
                    ubi, self.gv, self.hkl_tol, self.drlv2, labels, i
                )
            except:
                print(ubi.shape)
                print(self.gv.shape)
                print(self.hkl_tol)
                print(self.drlv2.shape)
                print(labels.shape)
                print("Error in fight_over_peaks", __file__)
                raise
        self.ga = labels
        # For each grain we want to know how many peaks it indexes
        # This is a histogram of labels
        bins = np.arange(-0.5, len(self.ubis) - 0.99)
        hst = myhistogram(labels, bins)
        self.gas = hst
        assert len(self.gas) == len(self.ubis)

    def saveindexing(self, filename, tol=None):
        """
        Save orientation matrices

        FIXME : refactor this into something to do
            grain by grain }
            peak by peak   }
        """
        f = open(filename, "w")
        i = 0
        from ImageD11 import transform

        self.gv = np.ascontiguousarray(self.gv)

        # grain assignment
        self.fight_over_peaks()

        # Printing per grain
        uinverses = []
        allind = np.array(list(range(len(self.ra))))
        tthcalc = np.zeros(len(self.ra), float)
        etacalc = np.zeros(len(self.ra), float)
        omegacalc = np.zeros(len(self.ra), float)
        i = -1

        for ubi in self.ubis:
            i += 1
            # Each ubi has peaks in self.ga
            uinverses.append(np.linalg.inv(ubi))
            # self.ga was filled in during fight_over_peaks
            npk, mdrlv = cImageD11.refine_assigned(ubi.copy(), self.gv, self.ga, i)
            assert npk == self.gas[i]
            f.write(
                "Grain: %d   Npeaks=%d   <drlv>=%f\n" % (i, self.gas[i], np.sqrt(mdrlv))
            )
            f.write("UBI:\n" + str(ubi) + "\n")
            cellpars = ubitocellpars(ubi)
            f.write("Cell pars: ")
            for abc in cellpars[:3]:
                f.write("%10.6f " % (abc))
            for abc in cellpars[3:]:
                f.write("%10.3f " % (abc))
            f.write("\n")
            # Grainspotter U
            f.write("U:\n" + str(ubitoU(ubi)) + "\n")
            f.write("B:\n" + str(ubitoB(ubi)) + "\n")

            # Compute hkls
            h = np.dot(ubi, self.gv.T)
            hint = np.floor(h + 0.5)
            gint = np.dot(uinverses[-1], hint)
            dr = h - hint

            f.write("Peak   (  h       k       l      )   drlv             x       y ")
            if self.wavelength < 0:
                f.write("\n")
            else:
                f.write(
                    "   Omega_obs Omega_calc   Eta_obs Eta_calc   tth_obs tth_calc\n"
                )

                tc, ec, oc = transform.uncompute_g_vectors(
                    gint, self.wavelength, wedge=self.wedge
                )
                ind = np.compress(self.ga == i, allind)

            for j in ind:
                f.write(
                    "%-6d ( % 6.4f % 6.4f % 6.4f ) % 12.8f "
                    % (j, h[0, j], h[1, j], h[2, j], np.sqrt(self.drlv2[j]))
                )
                f.write(" % 7.1f % 7.1f " % (self.xp[j], self.yp[j]))
                if self.wavelength < 0:
                    f.write("\n")
                else:
                    # # # These should be equal to
                    to = math.asin(self.wavelength * self.ds[j] / 2) * 360 / math.pi
                    # tth observed
                    eo = mod_360(self.eta[j], 0)
                    oo = self.omega[j]
                    tc1 = tc[j]
                    # Choose which is closest in eta/omega,
                    # there are two choices, {eta,omega}, {-eta,omega+180}
                    w = np.argmin([abs(ec[0][j] - eo), abs(ec[1][j] - eo)])
                    ec1 = ec[w][j]
                    oc1 = oc[w][j]
                    # Now find best omega within 360 degree intervals
                    oc1 = mod_360(oc1, oo)
                    f.write(
                        "  % 9.4f % 9.4f     % 9.4f % 9.4f   % 9.4f % 9.4f"
                        % (oo, oc1, eo, ec1, to, tc1)
                    )
                    etacalc[j] = ec1
                    omegacalc[j] = oc1
                    tthcalc[j] = tc1
                if self.ra[j] == -1:
                    f.write(" *** was not assigned to ring\n")
                else:
                    f.write("\n")
            f.write("\n\n")
        # peaks assigned to rings
        in_rings = np.compress(np.greater(self.ra, -1), np.arange(self.gv.shape[0]))
        f.write("\n\nAnd now listing via peaks which were assigned to rings\n")
        nleft = 0
        nfitted = 0
        npk = 0
        for peak in in_rings:
            # Compute hkl for each grain
            h = self.gv[peak, :]
            f.write(
                "\nPeak= %-5d Ring= %-5d gv=[ % -6.4f % -6.4f % -6.4f ]   omega= % 9.4f   eta= % 9.4f   tth= % 9.4f\n"
                % (
                    peak,
                    self.ra[peak],
                    h[0],
                    h[1],
                    h[2],
                    self.omega[peak],
                    self.eta[peak],
                    self.tth[peak],
                )
            )
            if self.ga[peak] != -1:
                m = self.ga[peak]
                hi = np.dot(self.ubis[m], h)
                hint = np.floor(hi + 0.5).astype(int)
                gint = np.dot(uinverses[m], hint)
                f.write("Grain %-5d (%3d,%3d,%3d)" % (m, hint[0], hint[1], hint[2]))
                f.write("  ( % -6.4f % -6.4f % -6.4f )  " % (hi[0], hi[1], hi[2]))
                # Now find best omega within 360 degree intervals
                f.write(
                    " omega= % 9.4f   eta= %9.4f   tth= %9.4f\n"
                    % (omegacalc[peak], etacalc[peak], tthcalc[peak])
                )
                npk = npk + 1
            else:
                if len(self.ubis) > 0:
                    f.write("Peak not assigned\n")
                    # , closest=[ % -6.4f % -6.4f % -6.4f ] for grain %d\n"%(hi[0],hi[1],hi[2],m))
                else:
                    f.write("Peak not assigned, no grains found\n")
                nleft = nleft + 1

        f.write("\n\nTotal number of peaks was %d\n" % (self.gv.shape[0]))
        f.write("Peaks assigned to grains %d\n" % (npk))
        f.write("Peaks assigned to rings but remaining unindexed %d\n" % (nleft))

        f.write(
            "Peaks not assigned to rings at all %d\n"
            % (np.sum(np.where(self.ra == -1, 1, 0)))
        )
        f.close()

    def getind(self, UBI, tol=None, drlv2tmp=None, labelstmp=None):
        """
        Returns the indices of peaks in self.gv indexed by matrix UBI
        """
        if tol == None:
            tol = self.hkl_tol
        ng = len(self.gvflat)
        if drlv2tmp is None:
            drlv2 = np.ones(ng, float)
        else:
            drlv2 = drlv2tmp
            drlv2[:] = 1
        if labelstmp is None:
            labels = np.zeros(ng, np.int32)
        else:
            labels = labelstmp
            labels[:] = 0
        # we only use peaks assigned to rings for scoring here
        # already done in making gvflat in assigntorings
        try:
            npk = cImageD11.score_and_assign(UBI, self.gv, tol, drlv2, labels, 1)
        except:
            print(self.gvflat.shape)
            print("ra", self.ra.shape)
            print(drlv2.shape)
            print(labels.shape)
            print("logic error in getind")
            raise
        return labels == 1

    def score(self, UBI, tol=None):
        """
        Decide which are the best orientation matrices
        """
        if tol is None:
            return cImageD11.score(UBI, self.gv, self.hkl_tol)
        else:
            return cImageD11.score(UBI, self.gv, tol)

    def refine(self, UBI):
        """
        Refine an orientation matrix and rescore it.

        From Paciorek et al Acta A55 543 (1999)
        UB = R H-1
           where:
           R = sum_n r_n h_n^t
           H = sum_n h_n h_n^t
           r = g-vectors
           h = hkl indices
        """
        #      print "Orientation and unit cell refinement of"
        #      print "UBI\n",UBI
        #      print "Scores before",self.score(UBI)
        # Need to find hkl indices for all of the peaks which are indexed
        drlv2 = calc_drlv2(UBI, self.gv)
        h = np.dot(UBI, np.transpose(self.gv))
        hint = np.floor(h + 0.5).astype(int)  # rounds down
        tol = float(self.hkl_tol)
        tol = tol * tol
        # Only use peaks which are assigned to rings for refinement
        ind = np.compress(
            np.logical_and(np.less(drlv2, tol), np.greater(self.ra, -1)),
            np.arange(self.gv.shape[0]),
        )
        # scoreb4=ind.shape[0]
        contribs = drlv2[ind]
        if len(contribs) == 0:
            raise Exception("No contributing reflections for" + str(UBI))
        # try:
        #    fitb4=sum(contribs)/contribs.shape[0]
        # except:
        #    print "No contributing reflections for\n",UBI
        #    raise
        # drlv2_old=drlv2
        R = np.zeros((3, 3), float)
        H = np.zeros((3, 3), float)
        for i in ind:
            r = self.gv[i, :]
            k = hint[:, i].astype(float)
            #           print r,k
            R = R + np.outer(r, k)
            H = H + np.outer(k, k)
        try:
            HI = np.linalg.inv(H)
            UBoptimal = np.dot(R, HI)
            UBIo = np.linalg.inv(UBoptimal)
        except:
            # A singular matrix - this sucks.
            UBIo = UBI
        drlv2 = calc_drlv2(UBIo, self.gv)
        ind = np.compress(
            np.logical_and(np.less(drlv2, tol), np.greater(self.ra, -1)),
            np.arange(self.gv.shape[0]),
        )
        self.scorelastrefined = ind.shape[0]
        contribs = drlv2[ind]
        if contribs.size != 0:
            self.fitlastrefined = math.sqrt(np.sum(contribs) / contribs.shape[0])
        else:
            raise ValueError("\n\n\nNo contributing reflections for ubi:\n{}\nAfter refinement, it was OK before ???\n\n\n".format(UBI))
        #      for i in ind:
        #         print "( %-6.4f %-6.4f %-6.4f ) %12.8f %12.8f"%(h[0,i],h[1,i],h[2,i],sqrt(drlv2[i]),sqrt(drlv2_old[i]))
        #      print UBIo
        #      print "Scores after", self.score(UBIo,self.hkl_tol)
        #      print "diff\n",UBI-UBIo
        #      print "Mean drlv now",sum(sqrt(drlv2))/drlv2.shape[0],
        #      print "Mean drlv old",sum(sqrt(drlv2_old))/drlv2_old.shape[0]
        return UBIo

    def saveubis(self, filename):
        """
        Save the generated ubi matrices into a text file
        """
        write_ubi_file(filename, self.ubis)

    def coverage(self):
        """
        Compute the expected coverage of reciprocal space
        use the min/max obs values of xp/yp/omega to work out what was measured in the scan?
        No lambda or
        """
        pass

    def readgvfile(self, filename, quiet=False):
        f = open(filename, "r")
        # Lattice!!!
        self.unitcell = unitcell.cellfromstring(f.readline())
        while 1:
            line = f.readline()
            if line[0] == "#":
                if line.find("wavelength") > -1:
                    self.wavelength = float(line.split()[-1])
                    if not quiet:
                        print("Got wavelength from gv file of ", self.wavelength)
                    continue
                if line.find("wedge") > -1:
                    self.wedge = float(line.split()[-1])
                    if not quiet:
                        print("Got wedge from gv file of ", self.wedge)
                    continue
                if line.find("ds h k l") > -1:
                    continue  # reads up to comment line
                if line.find("omega") > -1 and line.find("xc") > -1:
                    break
        self.eta = []  # Raw peak information
        self.omega = []
        self.ds = []
        self.xr = []
        self.yr = []
        self.zr = []
        self.xp = []
        self.yp = []
        try:
            for line in f.readlines():
                v = [float(x) for x in line.split()]
                if len(v) == 0:
                    # Skip the blank lines
                    continue
                if self.out_of_eta_range(v[6]):
                    continue
                self.xr.append(v[0])
                self.yr.append(v[1])
                self.zr.append(v[2])
                self.xp.append(v[3])
                self.yp.append(v[4])
                self.ds.append(v[5])
                self.eta.append(v[6])
                self.omega.append(v[7])
        except:
            print("LINE:", line)
            raise
        #            raise "Problem interpreting the last thing I printed"
        f.close()
        self.ds = np.array(self.ds)
        self.omega = np.array(self.omega)
        self.omega_fullrange = find_omega_ranges(self.omega)
        if self.wavelength > 0:
            self.tth = (
                np.arcsin(np.array(self.ds) * self.wavelength / 2) * 360 / math.pi
            )
        else:
            self.tth = np.zeros(len(self.ds))
        self.gv = np.transpose(np.array([self.xr, self.yr, self.zr], float))
        self.allgv = self.gv.copy()
        self.ga = np.zeros(len(self.ds), np.int32) - 1  # Grain assignments

        self.gvflat = np.ascontiguousarray(self.gv, float)
        self.gv = self.gvflat.copy()
        # Makes it contiguous in memory, hkl fast index
        if not quiet:
            print("Read your gv file containing", self.gv.shape)


def find_omega_ranges(omega):
    """
    This looks at self.omega and attempts to figure out coverage
    in degrees.

    We round to 1 degree as this is only used for an approximate
    guess to numbers of peaks anyway
    """
    iomega = np.round(omega % 360).astype(int)
    counts = np.bincount(iomega)
    return (counts > 0).sum()


def index(
    colf,
    npk_tol=[(400, 0.01), (200, 0.02)],
    cosine_tol=np.cos(np.radians(90 - 0.1)),
    ds_tol=0.004,
    max_grains=1000,
    rmulmax=10,
    rings_to_use=None,
    maxpairs=None,
    log_level=3,
):
    """
    Creates an indexer from a colfile
    Uses the unit cell from the colfile.parameters

    Loops over pairs of rings with multiplicity less than rmulmax

    npk_tol : list of (minpks, hkl_tol) to test
    cosine_tol : for finding pairs of peaks to make an orientation
    ds_tol : for assigning peaks to hkl rings in d*

    rmulmax : max multiplicity of rings to use in ubi generation (indexer.find)
    rings_to_use : which rings to consider to ubi generation (default = all)
    maxpairs : how many pairs of rings to use for ubi generation

    loglevel : like the logging module, but not

    returns the indexer object
    """
    global loglevel
    loglevel = log_level
    ind = indexer_from_colfile(
        colf, ds_tol=ds_tol, max_grains=max_grains, cosine_tol=cosine_tol
    )
    ind.assigntorings()
    for minpks, tol in npk_tol:
        ind.minpks = minpks
        ind.hkl_tol = tol
        ind.score_all_pairs(n=maxpairs, rmulmax=rmulmax, rings_to_use=rings_to_use)
    return ind


def do_index(
    cf,
    dstol=0.005,
    hkl_tols=(
        0.01,
        0.025,
        0.05,
    ),
    fracs=(
        0.9,
        0.7,
    ),
    cosine_tol=np.cos(np.radians(90 - 0.25)),
    max_grains=1000,
    forgen=(),
    foridx=(),
    unitcell=None,
    **kwds
):
    """
    Does indexing from a columnfile (cf)

    Returns a list of grains and the indexer object

    The peaks used for indexing (on rings in foridx) are selected and put into
    indexer.colfile

    The orientations in forgen are used to search orientations
    dstol, cosine_tol and max_grains are as usual (max grains is for each search).

    There is a loop over "frac" then "hkl_tols" repeating indexing many times.
    """
    for ringid in forgen:
        if ringid not in foridx:
            raise ValueError("All rings in forgen must be in foridx!")

    if cosine_tol < 0:
        import warnings

        warnings.warn("cosine_tol given as a negative number, are you sure about that?")

    logging.info("Indexing {} peaks".format(cf.nrows))

    global loglevel
    loglevel = 3

    # Figure out the peaks to use from foridx:
    if unitcell is not None:
        indexer = indexer_from_colfile_and_ucell(cf, ucell=unitcell,**kwds)
    else:
        indexer = indexer_from_colfile(cf)
    indexer.ds_tol = dstol
    indexer.assigntorings()
    # Only use the peaks in foridx:
    pkmask = np.zeros(cf.nrows, bool)
    for i in foridx:  # select the peaks in foridx
        pkmask |= indexer.ra == i
    cf_for_indexing = cf.copyrows(pkmask)
    cf_for_indexing.parameters = cf.parameters

    if unitcell:
        indexer = indexer_from_colfile_and_ucell(cf_for_indexing, ucell=unitcell,**kwds)
    else:
        indexer = indexer_from_colfile(cf_for_indexing)
    indexer.ds_tol = dstol
    indexer.assigntorings()

    omega_range = indexer.omega_fullrange
    if omega_range < 0:
        omega_range = 180  # guess it as 180

    n_peaks_expected = 0
    rings = []
    for i, dstar in enumerate(indexer.unitcell.ringds):
        # counts_on_this_ring = (indexer.ra == i).sum() is indexer.na above
        if indexer.na[i] > 0:  # useful peak
            if i in foridx:
                multiplicity = len(
                    indexer.unitcell.ringhkls[indexer.unitcell.ringds[i]]
                )
                n_peaks_expected += int(multiplicity * omega_range / 180.0)
                if i in forgen:  # we are generating orientations from this ring
                    rings.append((indexer.na[i], multiplicity, i))

    rings.sort()

    print("{} peaks expected".format(n_peaks_expected))
    print("Trying these rings (counts, multiplicity, ring number): {}".format(rings))

    indexer.cosine_tol = cosine_tol
    indexer.max_grains = max_grains

    try:
        threadb4 = cImageD11.cimaged11_omp_get_max_threads()
        cImageD11.cimaged11_omp_set_num_threads(1)  # ?

        for frac in fracs:
            indexer.minpks = n_peaks_expected * frac
            for indexer.hkl_tol in hkl_tols:
                for i in range(len(rings)):
                    for j in range(i, len(rings)):
                        indexer.ring_1 = rings[i][2]
                        indexer.ring_2 = rings[j][2]
                        indexer.find()
                        if len(indexer.hits) > 0:
                            indexer.scorethem()
                print(frac, indexer.hkl_tol, len(indexer.ubis))
    finally:
        cImageD11.cimaged11_omp_set_num_threads(threadb4)

    grains = [grain(ubi) for ubi in indexer.ubis]
    # if we supplied a unitcell, we probably want that as the reference for the grain
    if unitcell is not None:
        for g in grains:
            g.ref_unitcell = unitcell

    # set names for the grains
    for ginc, g in enumerate(grains):
        # try to make a name that includes the phase name
        try:
            g.name = g.ref_unitcell.name + ':' + str(ginc)
        except (NameError, KeyError, AttributeError) as e:
            g.name = str(ginc)

    return grains, indexer
