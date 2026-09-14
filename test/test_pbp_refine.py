"""unittest suite for the fast g-vector path and the origin finder in
pbp_refine.

    python3 -m unittest test_pbp_refine -v
    python3 -m unittest test_pbp_refine.TestGVectors -v
    python3 test_pbp_refine.py                        # same thing

The reference is a real ImageD11.columnfile with colf.updateGeometry() run on
it -- so xl/yl/zl and gx/gy/gz are exactly the columns a real analysis
produced, and both updateGeometry paths are covered: fast=True is the C
Ctransform, fast=False the python transform.

The refined origin enters as distance - xpos, so ImageD11's g-vectors at a
shifted origin are just updateGeometry() at a shifted distance.
xpos_refined therefore takes a small number of distinct values in the
fixture, and the reference is one updateGeometry call per value.

The partition comes from pbp_refine.build_partition, which goes through
VoxelSinoMasker -- so the CSR arrays compute_origins is fed are the ones
production builds, and build_partition itself gets exercised. PartitionInput
is the handful of attributes it reads off a PBPRefine.

Nothing in the forward model is re-derived either. The reflection list comes
from ImageD11.unitcell.gethkls, the Ewald condition from
ImageD11.transform.uncompute_g_vectors (which is gv_general.g_to_k), and the
projection onto the detector from ImageD11.transform.compute_xyz_from_tth_eta.
TestOrigins.test_forward_model_indexes checks the loop closes: peaks pushed
back through compute_tth_eta and compute_g_vectors at their true origin index
their own UBI to integers.

The sample-stage geometry is ImageD11.sinograms.geometry, not re-derived here:
step_to_sample for the voxel grid (so sy runs negative, as PBPRefine.setmap
builds it), dty_values_grain_in_beam for the dty that illuminates a voxel,
dty_to_dtyi for the discretisation, and sample_to_lab for the lab-frame x that
compute_origins reports.

compute_origins is checked against a forward model where the answer is known
exactly: give every voxel its own orientation and a peak is indexed by exactly
one voxel, so the weighted mean collapses and

    lx_modified[q] == G.sample_to_lab(sx, sy, y0, dty, omega)[0]

for the voxel the peak actually came from.

Note on omegasign: pbp_refine takes sin/cos of the raw omega column while
ImageD11 applies omegasign first, so the module is only correct for
omegasign = +1. TestOmegaSign pins that down in both directions -- if the
convention is ever fixed in setpeaks, test_raw_omega_is_wrong_at_minus_one
is the test that will start failing, and it should be deleted then.
"""

import unittest

import numpy as np


import sys
if int(sys.version_info.major) == 2:
    raise unittest.SkipTest('Skipping pbp_refine tests on Python 2')
else:
    import ImageD11.columnfile

    import ImageD11.parameters
    import ImageD11.transform as T
    import ImageD11.unitcell
    from ImageD11.sinograms import geometry as G, pbp_refine as M

# --------------------------------------------------------------------------
#  Geometry fixtures
# --------------------------------------------------------------------------
# the keys ImageD11.transform.compute_* accept
TKEYS = ("distance", "y_center", "y_size", "tilt_y", "z_center", "z_size",
         "tilt_z", "tilt_x", "o11", "o12", "o21", "o22", "t_x", "t_y", "t_z",
         "wedge", "chi")

BASE = dict(distance=150000.0, y_center=1024.3, y_size=75.0, tilt_y=0.0021,
            z_center=1019.7, z_size=75.0, tilt_z=-0.0043, tilt_x=0.0007,
            o11=1.0, o12=0.0, o21=0.0, o22=-1.0,
            t_x=0.0, t_y=0.0, t_z=0.0, wedge=0.0, chi=0.0,
            wavelength=0.2846, omegasign=1.0)

CASES = {
    "plain":                dict(),
    "wedge":                dict(wedge=1.7),
    "chi":                  dict(chi=-0.9),
    "wedge_chi":            dict(wedge=1.7, chi=-0.9),
    "t_xyz":                dict(t_x=13.0, t_y=-7.0, t_z=4.0),
    "everything":           dict(wedge=1.7, chi=-0.9, t_x=13.0, t_y=-7.0,
                                 t_z=4.0),
    "flipped_detector":     dict(o11=0.0, o12=-1.0, o21=-1.0, o22=0.0,
                                 wedge=0.4, chi=0.3, t_x=5.0, t_y=2.0),
    "transposed_detector":  dict(o11=0.0, o12=1.0, o21=1.0, o22=0.0,
                                 wedge=-0.6, chi=0.5, t_z=-3.0),
    "big_tilts":            dict(tilt_x=0.05, tilt_y=-0.04, tilt_z=0.03,
                                 wedge=2.0, chi=1.0, t_x=9.0, t_y=-3.0,
                                 t_z=1.0),
    "negative_wedge_chi":   dict(wedge=-1.2, chi=0.7, t_x=-8.0, t_y=5.0),
}

SAMPLE_R = 250.0        # sample radius, same units as distance and dty
YSTEP = 2.0
Y0 = 0.35
ACELL = 4.05

# machine precision is ~1e-16 here; ImageD11's own C and python paths differ
# from each other by ~7e-16 on the same columnfile, so this is the floor
TOL = 1e-12


def pars_for(case):
    return dict(BASE, **CASES[case])


# --------------------------------------------------------------------------
#  A columnfile that looks like the one the refinement is handed
# --------------------------------------------------------------------------
_ICOLF_CACHE = {}


def _f(a):
    return np.ascontiguousarray(a, dtype=np.float64)


def make_colf(sc, fc, omega, pars, fast=True, extra=()):
    """A columnfile with updateGeometry() run on it, so xl/yl/zl, tth, eta, ds
    and gx/gy/gz are ImageD11's."""
    titles = ["sc", "fc", "omega"]
    cols = [_f(sc), _f(fc), _f(omega)]
    for name, arr in extra:
        titles.append(name)
        cols.append(_f(arr))
    c = ImageD11.columnfile.newcolumnfile(titles)
    # a list of 1-D arrays, not a 2-D array: set_bigarray keeps the object it
    # is handed as its column store, and updateGeometry appends to it
    c.set_bigarray(cols)
    c.setparameters(ImageD11.parameters.parameters(**pars))
    c.updateGeometry(fast=fast)
    return c


def fake_icolf(pars, n=5000, seed=0, fast=True, nxpos=8):
    """The columnfile the refinement is handed.

    On top of what updateGeometry writes: sinomega/cosomega as
    PBPRefine.setpeaks adds them, and xpos_refined as get_origins adds it --
    the lab-frame x of the diffraction origin. It is quantised to nxpos values
    so id11_gvectors can build its reference with that many updateGeometry
    calls; nothing in the transform cares where the number came from.

    Cached, because several tests want the same one.
    """
    key = (tuple(sorted(pars.items())), n, seed, fast, nxpos)
    if key in _ICOLF_CACHE:
        return _ICOLF_CACHE[key]

    rng = np.random.default_rng(seed)
    sc = rng.uniform(5.0, 2043.0, n)
    fc = rng.uniform(5.0, 2043.0, n)
    omega = rng.uniform(-180.0, 180.0, n)
    # a voxel somewhere in the sample, and the dty that puts it in the beam
    r = SAMPLE_R * np.sqrt(rng.uniform(0.0, 1.0, n))
    th = rng.uniform(0.0, 2.0 * np.pi, n)
    sx, sy = r * np.cos(th), r * np.sin(th)
    dty = G.dty_values_grain_in_beam(sx, sy, Y0, omega)

    c = make_colf(sc, fc, omega, pars, fast=fast, extra=(
        ("dty", dty),
        ("dtyi", G.dty_to_dtyi(dty, YSTEP, dty.min())),
        ("sum_intensity", rng.uniform(1e3, 1e6, n)),
        ("Number_of_pixels", rng.uniform(4.0, 400.0, n))))

    c.addcolumn(np.sin(np.radians(c.omega)), "sinomega")
    c.addcolumn(np.cos(np.radians(c.omega)), "cosomega")
    # xpos_refined is the lab x of the diffraction origin, i.e.
    # G.sample_to_lab(sx, sy, y0, dty, omega)[0]. Quantised to nxpos values so
    # that id11_gvectors needs that many updateGeometry calls; the transform
    # does not care where the number came from.
    lx = G.sample_to_lab(sx, sy, Y0, dty, omega)[0]
    edges = np.linspace(lx.min(), lx.max(), nxpos + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    c.addcolumn(centres[np.clip(np.searchsorted(edges, lx) - 1,
                                0, nxpos - 1)], "xpos_refined")
    _ICOLF_CACHE[key] = c
    return c


def fast_gvectors(icolf, xpos=None, sinomega=None, cosomega=None):
    """What build_refine_inputs computes, minus the partition permutation."""
    pars = icolf.parameters.get_parameters()
    A, oc, invw, aeye = M.gvec_pars(pars)
    xp = _f(icolf.xpos_refined) if xpos is None else _f(xpos)
    so = _f(icolf.sinomega) if sinomega is None else _f(sinomega)
    co = _f(icolf.cosomega) if cosomega is None else _f(cosomega)
    return M.gvectors_from_lab(_f(icolf.xl), _f(icolf.yl), _f(icolf.zl),
                               xp, so, co, A, oc, invw, aeye, 4)


def id11_gvectors(icolf, fast=True):
    """ImageD11's g-vectors with the origin of diffraction at +xpos along the
    beam: colf.updateGeometry() at distance - xpos, one call per distinct
    xpos. TestOriginShift is what establishes that equivalence."""
    pars = dict(icolf.parameters.get_parameters())
    d0 = float(pars["distance"])
    xp = _f(icolf.xpos_refined)
    out = np.empty((icolf.nrows, 3))
    for v in np.unique(xp):
        c = make_colf(icolf.sc, icolf.fc, icolf.omega,
                      dict(pars, distance=d0 - float(v)), fast=fast)
        sel = xp == v
        for k, name in enumerate(("gx", "gy", "gz")):
            out[sel, k] = _f(getattr(c, name))[sel]
    return out


def relmax(ref, got):
    """max |ref - got|, relative to the largest element of ref."""
    ref = np.asarray(ref, dtype=np.float64)
    got = np.asarray(got, dtype=np.float64)
    scale = np.abs(ref).max()
    return float(np.abs(ref - got).max() / (scale if scale > 0 else 1.0))


# --------------------------------------------------------------------------
#  The detector transform
# --------------------------------------------------------------------------
class TestLabFrame(unittest.TestCase):
    """the linear part of compute_xyz_lab"""

    def test_affine_in_sc_fc(self):
        """compute_xyz_lab is affine in (sc, fc), which is what lets
        refine_map average xl/yl/zl over a merge group instead of averaging
        sc/fc and transforming afterwards."""
        for case in CASES:
            with self.subTest(case):
                pars = pars_for(case)
                kw = {k: float(pars[k]) for k in TKEYS}
                rng = np.random.default_rng(3)
                sc = rng.uniform(5.0, 2043.0, 5000)
                fc = rng.uniform(5.0, 2043.0, 5000)
                w = rng.uniform(0.1, 1e5, 5000)          # sum_intensity
                xyz = np.array(T.compute_xyz_lab((sc, fc), **kw))
                mean_then_transform = np.array(T.compute_xyz_lab(
                    (np.array([(sc * w).sum() / w.sum()]),
                     np.array([(fc * w).sum() / w.sum()])), **kw)).ravel()
                transform_then_mean = (xyz * w).sum(axis=1) / w.sum()
                self.assertLess(
                    relmax(mean_then_transform, transform_then_mean), TOL)


class TestOriginShift(unittest.TestCase):
    """Inside ImageD11: reducing `distance` by d shifts lab x by d and leaves
    y and z alone. This is what lets a per-peak origin along the beam become a
    per-peak subtraction from the xl column."""

    def test_distance_is_a_lab_x_shift(self):
        for case in CASES:
            with self.subTest(case):
                pars = pars_for(case)
                rng = np.random.default_rng(2)
                sc = rng.uniform(5.0, 2043.0, 4000)
                fc = rng.uniform(5.0, 2043.0, 4000)
                om = rng.uniform(-180.0, 180.0, 4000)
                base = make_colf(sc, fc, om, pars)
                for d in (-SAMPLE_R, -1.0, 1.0, SAMPLE_R):
                    with self.subTest(case, offset=d):
                        moved = make_colf(
                            sc, fc, om,
                            dict(pars, distance=float(pars["distance"]) - d))
                        self.assertLess(
                            relmax(_f(base.xl) - d, _f(moved.xl)), 1e-14)
                        # y and z are untouched, so bit-identical
                        self.assertTrue(np.array_equal(_f(base.yl),
                                                       _f(moved.yl)))
                        self.assertTrue(np.array_equal(_f(base.zl),
                                                       _f(moved.zl)))


# --------------------------------------------------------------------------
#  g-vectors
# --------------------------------------------------------------------------
class TestGVectors(unittest.TestCase):

    def test_nominal_origin_vs_columnfile(self):
        """xpos = 0 must reproduce the gx/gy/gz columns updateGeometry wrote,
        on both the C and the python path."""
        for case in CASES:
            for fast in (True, False):
                with self.subTest(case, fast=fast):
                    c = fake_icolf(pars_for(case), fast=fast)
                    ref = np.ascontiguousarray(
                        np.array((c.gx, c.gy, c.gz)).T)
                    got = fast_gvectors(c, xpos=np.zeros(c.nrows))
                    self.assertLess(relmax(ref, got), TOL)

    def test_refined_origin_vs_imaged11(self):
        """The real case: the origin of diffraction moved per peak."""
        for case in CASES:
            for fast in (True, False):
                with self.subTest(case, fast=fast):
                    c = fake_icolf(pars_for(case), fast=fast)
                    self.assertGreater(
                        np.abs(c.xpos_refined).max(), 0.5 * SAMPLE_R,
                        "fixture is not exercising xpos")
                    ref = id11_gvectors(c, fast=fast)
                    self.assertLess(relmax(ref, fast_gvectors(c)), TOL)

    def test_xpos_actually_matters(self):
        """Guard against a test that would pass with xpos ignored: dropping it
        has to move the answer by far more than TOL."""
        c = fake_icolf(pars_for("plain"))
        with_xpos = fast_gvectors(c)
        without = fast_gvectors(c, xpos=np.zeros(c.nrows))
        self.assertGreater(relmax(with_xpos, without), 1e-4)


class TestOmegaSign(unittest.TestCase):
    """pbp_refine feeds the kernel sin/cos of the raw omega column, so it is
    only correct for omegasign = +1. gvec_one itself is sign-agnostic."""

    def test_plus_one_agrees(self):
        c = fake_icolf(dict(BASE, omegasign=1.0, wedge=1.1, chi=-0.5,
                            t_x=6.0, t_y=-2.0))
        self.assertLess(relmax(id11_gvectors(c), fast_gvectors(c)), TOL)

    def test_minus_one_agrees_when_sign_applied(self):
        pars = dict(BASE, omegasign=-1.0, wedge=1.1, chi=-0.5,
                    t_x=6.0, t_y=-2.0)
        c = fake_icolf(pars)
        om = _f(c.omega) * float(pars["omegasign"])
        got = fast_gvectors(c, sinomega=np.sin(np.radians(om)),
                            cosomega=np.cos(np.radians(om)))
        self.assertLess(relmax(id11_gvectors(c), got), TOL)

    def test_raw_omega_is_wrong_at_minus_one(self):
        """Known, pre-existing: setpeaks fills sinomega/cosomega from the raw
        omega column. Delete this test when that is fixed."""
        pars = dict(BASE, omegasign=-1.0, wedge=1.1, chi=-0.5,
                    t_x=6.0, t_y=-2.0)
        c = fake_icolf(pars)
        self.assertGreater(relmax(id11_gvectors(c), fast_gvectors(c)), 1e-3)


# --------------------------------------------------------------------------
#  The conditioning test that replaced np.linalg.svd
# --------------------------------------------------------------------------
class TestConditioning(unittest.TestCase):

    EPS = 2.220446049250313e-16

    @staticmethod
    def _matrices(n, seed=17):
        """A mix from realistic UBIs to rank deficient, including the 1e14
        threshold band itself."""
        rng = np.random.default_rng(seed)
        for k in range(n):
            mode = k % 4
            if mode == 0:                              # a realistic UBI
                e = 10.0 ** rng.uniform(-0.3, 0.3, 3)
            elif mode == 1:                            # moderately ill
                e = 10.0 ** rng.uniform(-6.0, 0.0, 3)
            elif mode == 2:                            # at the threshold
                e = np.array([1.0, 10.0 ** rng.uniform(-16.0, -12.0), 1.0])
            else:                                      # rank deficient
                e = np.array([1.0, rng.uniform(0.1, 1.0), 0.0])
            Q1, _ = np.linalg.qr(rng.normal(size=(3, 3)))
            Q2, _ = np.linalg.qr(rng.normal(size=(3, 3)))
            yield (np.ascontiguousarray(Q1.dot(np.diag(e)).dot(Q2.T))
                   * 10.0 ** rng.uniform(-3.0, 3.0))

    def _lapack_verdict(self, U):
        sv = np.linalg.svd(U)[1]
        smax, smin = sv[0], sv[2]
        rank3 = smin > smax * 3.0 * self.EPS
        cnd = 1e300 if smin <= 0.0 else smax / smin
        return bool(rank3 and cnd < 1e14)

    def test_same_verdict_as_lapack(self):
        """Agree with LAPACK wherever the answer is actually determined.
 
        Right at cond_max it is not. Rounding one of these matrices into
        float64 moves its smallest singular value by about eps*smax, which is
        eps*cond in relative terms -- 2% at cond 1e14, and worse above it. So
        for a matrix sitting within a couple of percent of the cut, "is
        cond < 1e14" has no float64 answer: numpy's LAPACK and the one numba
        links against can each return a different smin, both correct to the
        accuracy available, and land on opposite sides. Asserting exact
        agreement there makes the test pass or fail depending on the runner's
        BLAS build, which is not something the code can control.
 
        Only the threshold band of _matrices comes near the cut at all; the
        realistic, moderately ill and rank deficient modes are orders of
        magnitude clear of it, and the rank deficient ones are rejected on
        condition number regardless of which side of the rank test they fall.
        So skip the matrices float64 cannot place, and require exact agreement
        on all the rest.
        """
        n = 40000
        bad = shortcut = ambiguous = 0
        for U in self._matrices(n):
            hi, lo = M._sv_extremes(U)
            if lo > 0.0 and hi < lo * 1e6:
                shortcut += 1
            sv = np.linalg.svd(U)[1]
            smax, smin = sv[0], sv[2]
            if smin > 0.0:
                cond = smax / smin
                # how well is smin defined at this conditioning, relatively
                slack = 8.0 * self.EPS * cond
                if abs(cond / 1e14 - 1.0) < slack:
                    ambiguous += 1
                    continue
            if bool(M._well_conditioned(U)) != self._lapack_verdict(U):
                bad += 1
        self.assertEqual(bad, 0, "%d/%d disagreements with LAPACK" % (bad, n))
        self.assertLess(ambiguous, n // 20,
                        "%d/%d matrices skipped as undecidable -- the guard "
                        "band has swallowed the test" % (ambiguous, n))
        self.assertGreater(shortcut, n // 10,
                           "the analytic shortcut almost never fires, so the "
                           "cheap path is not being exercised")

    def test_sv_extremes_accurate_when_well_conditioned(self):
        """Inside the guard band the analytic singular values are right, not
        merely conservative -- typically to ~1e-14.

        The tail is heavier than the median: Cardano goes through
        arccos(d)/3, and d approaches +-1 when two eigenvalues nearly
        coincide, where the derivative blows up. Over 20000 well-conditioned
        matrices the median relative error is ~9e-15 and p99 ~8e-11, but the
        worst (relative eigen-gap ~2e-7) reaches a few 1e-6. That does not
        matter for what _well_conditioned uses it for -- deciding whether the
        matrix is comfortably under a 1e6 guard, itself eight orders under the
        1e14 cut -- so assert the bulk tightly and just bound the tail.
        """
        rng = np.random.default_rng(23)
        errs = np.empty(20000)
        for i in range(errs.size):
            e = 10.0 ** rng.uniform(-2.0, 0.0, 3)      # cond <= 100
            Q1, _ = np.linalg.qr(rng.normal(size=(3, 3)))
            Q2, _ = np.linalg.qr(rng.normal(size=(3, 3)))
            U = np.ascontiguousarray(Q1.dot(np.diag(e)).dot(Q2.T))
            sv = np.linalg.svd(U)[1]
            hi, lo = M._sv_extremes(U)
            errs[i] = max(abs(hi - sv[0]) / sv[0], abs(lo - sv[2]) / sv[2])
        self.assertLess(np.median(errs), 1e-12)
        self.assertLess(np.percentile(errs, 99), 1e-9)
        self.assertLess(errs.max(), 1e-4)


# --------------------------------------------------------------------------
#  compute_origins, against a forward model
# --------------------------------------------------------------------------
def reflections(acell=ACELL, lattice="F", s2max=44):
    """(3, n) allowed hkl, from ImageD11.unitcell rather than a parity rule."""
    uc = ImageD11.unitcell.unitcell([acell] * 3 + [90.0] * 3, lattice)
    hkls = uc.gethkls(np.sqrt(s2max) / acell)
    return np.array([hkl for _ds, hkl in hkls], dtype=np.float64).T


HKL = reflections()


class PartitionInput(object):
    """The attributes pbp_refine.build_partition reads off a PBPRefine.

    Not a PBPRefine: that wants a DataSet and files on disk. build_partition
    only touches icolf.omega / .dty / .nrows / .titles, the sample grid and
    mask (for rmax), ystep and ymin -- and choose_omega_bins falls through to
    VoxelSinoMasker's heuristic when there is no frm column and no dset.
    """

    def __init__(self, case, ymin=None):
        c = ImageD11.columnfile.newcolumnfile(["omega", "dty"])
        c.set_bigarray([_f(case["omega"]), _f(case["dty"])])
        self.icolf = c
        self.sx_grid = case["sx_grid"]
        self.sy_grid = case["sy_grid"]
        self.mask = case["mask"]
        self.ystep = YSTEP
        self.ymin = float(case["dty"].min()) if ymin is None else float(ymin)
        self.dset = None          # no frm, no obinedges -> the heuristic


def compute_origins_ref(case, hkl_tol=0.05, weight_reg=YSTEP / 3.0):
    """An oracle for compute_origins: every (peak, voxel) pair, in numpy.

    Deliberately NOT built on the partition or on fill_voxel_idx. Those are
    the accelerations compute_origins is meant to be equivalent to, so an
    oracle that used them could not catch a strip walk that clips a voxel. It
    is the O(npeaks x nvoxels) version the whole scheme exists to avoid.

    Same predicates in the same order, so the two agree bit for bit:
      * beam through the voxel, |dty - dty_in_beam| <= ystep, from
        G.dty_values_grain_in_beam_sincos
      * hkl residual, |UBI g - round(UBI g)|^2 < hkl_tol^2
      * weight w = r / (ydist + r), eq (6)
      * report G.sample_to_lab's lx of the weighted mean voxel

    Returns (lx, nclaim): the origin per peak, and how many voxels claimed it.
    """
    sx_ax, sy_ax = case["sx_ax"], case["sy_ax"]
    singlemap, mask, gve = case["singlemap"], case["mask"], case["gve"]
    dty = case["dty"]
    so = np.sin(np.radians(case["omega"]))
    co = np.cos(np.radians(case["omega"]))
    n = dty.size
    accx = np.zeros(n)
    accy = np.zeros(n)
    accw = np.zeros(n)
    nclaim = np.zeros(n, dtype=np.int64)
    tolsq = hkl_tol * hkl_tol
    for i in range(sx_ax.size):
        for j in range(sy_ax.size):
            if not mask[i, j]:
                continue
            ubi = singlemap[i, j]
            if np.isnan(ubi[0, 0]):
                continue
            yd = np.abs(G.dty_values_grain_in_beam_sincos(
                sx_ax[i], sy_ax[j], Y0, so, co) - dty)
            sel = yd <= YSTEP
            if not sel.any():
                continue
            hf = gve[sel].dot(ubi.T)                     # (m, 3) == ubi @ g
            dh = hf - np.round(hf)
            hit = (dh * dh).sum(axis=1) < tolsq
            idx = np.nonzero(sel)[0][hit]
            w = weight_reg / (yd[idx] + weight_reg)
            accx[idx] += sx_ax[i] * w
            accy[idx] += sy_ax[j] * w
            accw[idx] += w
            nclaim[idx] += 1
    lx = np.zeros(n)
    got = accw > 0.0
    lx[got] = G.sample_to_lab_sincos(accx[got] / accw[got],
                                     accy[got] / accw[got],
                                     Y0, dty[got], so[got], co[got])[0]
    return lx, nclaim


class TestOrigins(unittest.TestCase):
    """compute_origins recovers the voxel a peak came from.

    Peaks are generated with the origin of diffraction AT the voxel, which is
    the physics, while the g-vectors handed to compute_origins use the NOMINAL
    origin, as icolf.gx/gy/gz do. That mismatch is the ~xpos/distance error
    the hkl_tol_origins window exists to absorb, so it is part of the test.
    """

    NI = NJ = 9
    PARS = dict(BASE)

    @classmethod
    def setUpClass(cls):
        cls.case = cls._make_case()
        cls.lx, _ = cls._run(cls.case)
        cls.ref, cls.nclaim = compute_origins_ref(cls.case)

    # ---------------------------------------------------------- fixture
    @classmethod
    def _make_case(cls, one_grain_per_voxel=True, seed=5):
        rng = np.random.default_rng(seed)
        pars = cls.PARS
        ni, nj = cls.NI, cls.NJ
        # the grid PBPRefine.setmap builds: step space through step_to_sample,
        # so sy descends (sy = -sj * ystep) and dsy is negative, which is the
        # sign convention compute_origins' strip walk actually meets
        si = np.arange(ni) - ni // 2
        sj = np.arange(nj) - nj // 2
        sx_ax, sy_ax = G.step_to_sample(si, sj, YSTEP)
        sx_grid, sy_grid = np.meshgrid(sx_ax, sy_ax, indexing="ij")

        def random_U():
            Q, R = np.linalg.qr(rng.normal(size=(3, 3)))
            Q = Q * np.sign(np.diag(R))
            if np.linalg.det(Q) < 0:
                Q[:, 0] *= -1
            return Q

        B = np.eye(3) / ACELL
        if one_grain_per_voxel:
            Us = [random_U() for _ in range(ni * nj)]
            which = np.arange(ni * nj).reshape(ni, nj)
        else:
            Us = [random_U() for _ in range(4)]
            which = rng.integers(0, 4, size=(ni, nj))

        singlemap = np.empty((ni, nj, 3, 3))
        vox = []
        for i in range(ni):
            for j in range(nj):
                UB = Us[which[i, j]].dot(B)
                singlemap[i, j] = np.linalg.inv(UB)
                vox.append((sx_grid[i, j], sy_grid[i, j], UB))

        # ---- forward model, all ImageD11 ---------------------------------
        # uncompute_g_vectors is gv_general.g_to_k: it solves the Ewald
        # condition and returns both omega branches. compute_xyz_from_tth_eta
        # ray-traces onto the detector with the origin of diffraction at the
        # voxel's rotated sample position (compute_grain_origins of t_x, t_y).
        wvln = float(pars["wavelength"])
        detkw = {k: float(pars[k]) for k in TKEYS}
        rows = []
        for iv, (sx, sy, UB) in enumerate(vox):
            tth, etas, omegas = T.uncompute_g_vectors(
                UB.dot(HKL), wvln, wedge=detkw["wedge"], chi=detkw["chi"])
            ok = tth > 0.0        # invalid solutions come back zeroed
            if not ok.any():
                continue
            for eta, om in zip(etas, omegas):
                fc, sc = T.compute_xyz_from_tth_eta(
                    tth[ok], eta[ok], om[ok],
                    **dict(detkw, t_x=sx, t_y=sy, t_z=0.0))
                on = ((sc > 2.0) & (sc < 2046.0)
                      & (fc > 2.0) & (fc < 2046.0))
                npk = int(on.sum())
                if npk == 0:
                    continue
                omd = om[ok][on]
                # dty that illuminates this voxel, jittered inside the beam.
                # Uniform in +-ystep/2, not gaussian: then
                # |dty - dty_in_beam| <= ystep/2 < ystep for every peak, so
                # the voxel it came from always passes compute_origins' window
                # and every peak must be claimed. A gaussian tail would leave a
                # few peaks unclaimed and there would be no exact number to
                # assert.
                dty = (G.dty_values_grain_in_beam(sx, sy, Y0, omd)
                       + rng.uniform(-0.5 * YSTEP, 0.5 * YSTEP, npk))
                rows.append(np.column_stack((
                    sc[on] + rng.normal(0.0, 0.1, npk),
                    fc[on] + rng.normal(0.0, 0.1, npk),
                    omd, dty,
                    # the truth: the lab x of that voxel at that omega
                    G.sample_to_lab(sx, sy, Y0, dty, omd)[0],
                    np.full(npk, sx), np.full(npk, sy),
                    np.full(npk, iv))))
        a = np.concatenate(rows)
        sc, fc, omega, dty, lx_true = (a[:, 0], a[:, 1], a[:, 2], a[:, 3],
                                       a[:, 4])
        vox_of_peak = (a[:, 5], a[:, 6])       # the voxel each peak came from
        vox_idx = a[:, 7].astype(np.int64)

        # g-vectors at the nominal origin: the icolf gx/gy/gz columns, which
        # is what get_origins feeds compute_origins
        gc = make_colf(sc, fc, omega, pars)
        gve = np.ascontiguousarray(
            np.array((_f(gc.gx), _f(gc.gy), _f(gc.gz))).T)

        return dict(sc=sc, fc=fc, omega=omega, dty=dty, lx_true=lx_true,
                    sx_grid=np.ascontiguousarray(sx_grid),
                    sy_grid=np.ascontiguousarray(sy_grid),
                    vox_of_peak=vox_of_peak, vox_idx=vox_idx,
                    vox=vox, pars=pars,
                    gve=gve, singlemap=np.ascontiguousarray(singlemap),
                    mask=np.ones((ni, nj), dtype=bool),
                    sx_ax=np.ascontiguousarray(sx_ax),
                    sy_ax=np.ascontiguousarray(sy_ax),
                    rmax=float(np.hypot(sx_grid, sy_grid).max()),
                    npk=sc.shape[0])

    # ---------------------------------------------------------- driver
    @staticmethod
    def _run(case, hkl_tol=0.05, nchunks=3, omega_binsize=1.0, ymin=None):
        # the module's own partition, through VoxelSinoMasker
        part = M.build_partition(PartitionInput(case, ymin),
                                 omega_binsize=omega_binsize, verbose=False)
        order = part["order"]
        lx_perm = M.compute_origins(
            case["singlemap"], case["mask"],
            np.ascontiguousarray(case["gve"][order]),
            np.ascontiguousarray(np.sin(np.radians(case["omega"]))[order]),
            np.ascontiguousarray(np.cos(np.radians(case["omega"]))[order]),
            np.ascontiguousarray(case["dty"][order]),
            case["sx_ax"], case["sy_ax"],
            float(case["sx_ax"][1] - case["sx_ax"][0]),
            float(case["sy_ax"][1] - case["sy_ax"][0]),
            Y0, YSTEP, hkl_tol, YSTEP / 3.0,
            part["omega_partitions"], part["dty_partitions"],
            part["nom"], part["ndty"], part["usin"], part["ucos"],
            part["ymin"], part["ray_margin"],
            int(np.diff(part["dty_partitions"], axis=1).max()), nchunks)
        lx = np.empty(case["npk"])
        lx[order] = lx_perm
        return lx, part

    # ---------------------------------------------------------- tests
    def test_forward_model_indexes(self):
        """The loop closes: take the generated sc/fc back through ImageD11 at
        the origin they were made with, and the g-vectors index that voxel's
        own UBI to integers. If this fails, the fixture is wrong, not
        compute_origins."""
        case = self.case
        detkw = {k: float(case["pars"][k]) for k in TKEYS}
        wvln = float(case["pars"]["wavelength"])
        worst = 0.0
        for iv in np.unique(case["vox_idx"])[::7]:
            sel = case["vox_idx"] == iv
            sx, sy, UB = case["vox"][iv]
            om = case["omega"][sel]
            tth, eta = T.compute_tth_eta(
                (case["sc"][sel], case["fc"][sel]), omega=om,
                **dict(detkw, t_x=sx, t_y=sy, t_z=0.0))
            g = np.asarray(T.compute_g_vectors(
                tth, eta, om, wvln, wedge=detkw["wedge"],
                chi=detkw["chi"]))
            hkl = np.linalg.inv(UB).dot(g)
            worst = max(worst, float(np.abs(hkl - np.round(hkl)).max()))
        self.assertLess(worst, 0.02,
                        "worst hkl residual %.3g -- the forward model does "
                        "not reproduce its own reflections" % worst)

    # ------------------------------------------------------- comparison
    def _assert_matches_reference(self, lx, ref, nclaim, rmax, label=""):
        """compute_origins against the brute-force reference.

        Up to two claiming voxels the two agree bit for bit: IEEE addition is
        commutative, so the order the voxels are visited in cannot matter.
        From three it is not associative, and compute_origins walks the strip
        i-major or j-major depending on the omega bin while the reference is
        always i-major -- so those sums can land one rounding apart. The bound
        is in ulp of rmax, the scale of the terms being summed, not a
        tolerance picked to make the test pass. Measured worst: 1.2 ulp.
        """
        ulp = float(np.spacing(rmax))
        d = np.abs(lx - ref)
        few = nclaim <= 2
        if few.any():
            self.assertEqual(
                float(d[few].max()), 0.0,
                "%s%d/%d peaks with <=2 claiming voxels are not bit-exact "
                "(max %.3e)" % (label, int((d[few] != 0).sum()),
                                int(few.sum()), d[few].max()))
        many = nclaim >= 3
        if many.any():
            self.assertLessEqual(
                float(d[many].max()), 2.0 * ulp,
                "%s%d peaks with >=3 claiming voxels differ by up to %.1f ulp "
                "of rmax" % (label, int((d[many] != 0).sum()),
                             d[many].max() / ulp))

    # ------------------------------------------------------------ tests
    def test_every_peak_is_claimed(self):
        """Every peak's own voxel passes the beam window by construction (the
        dty jitter is uniform inside +-ystep/2) and its hkl residual against
        its own UBI is ~4e-4, far inside hkl_tol_origins. So no peak may go
        unclaimed -- not almost none."""
        self.assertEqual(int((self.nclaim == 0).sum()), 0)
        self.assertEqual(int(self.nclaim.min()), 1)

    def test_matches_brute_force_reference(self):
        """compute_origins equals the brute-force reference on every peak: bit
        for bit where at most two voxels claimed it, and within 2 ulp of rmax
        where three or more did."""
        self._assert_matches_reference(self.lx, self.ref, self.nclaim,
                                       self.case["rmax"])

    def test_single_claim_peaks_give_the_true_voxel(self):
        """Where exactly one voxel claimed a peak, that voxel is the one the
        peak came from -- every time, since the true voxel always claims. The
        value is then G.sample_to_lab's lx for that voxel, to within the one
        rounding of accx/accw: assert 2 ulp, not a tolerance.

        lx is not tested for being nonzero anywhere: the voxel on the rotation
        axis has lab x of exactly 0, so lx == 0 is a legitimate answer and
        cannot distinguish "unclaimed" from "claimed by the on-axis voxel".
        """
        one = self.nclaim == 1
        self.assertGreater(one.mean(), 0.9)
        sx, sy = self.case["vox_of_peak"]
        want = G.sample_to_lab(sx[one], sy[one], Y0,
                               self.case["dty"][one],
                               self.case["omega"][one])[0]
        got = self.lx[one]
        off = np.abs(got - want)
        # the rounding happens at the scale of the terms, not of the answer:
        # accx/accw round-trips sx through a multiply and a divide, and
        # lx = sxc cos - syc sin can then cancel to near zero. So the bound is
        # ulp of rmax. Measured worst: 2.0 ulp.
        ulp = float(np.spacing(self.case["rmax"]))
        bad = int((off > 2.0 * ulp).sum())
        self.assertEqual(bad, 0,
                         "%d/%d single-claim peaks off by more than 2 ulp of "
                         "rmax (worst %.1f ulp)"
                         % (bad, int(one.sum()), off.max() / ulp))
        # and the on-axis voxel really is in play, so the zero case is covered
        self.assertGreater(int((want == 0.0).sum()), 0)

    def test_bounded_by_sample_radius(self):
        """lx is a convex combination of voxel centroids rotated into the lab,
        so |lx| <= rmax always. A hard bound, not a tolerance."""
        self.assertLessEqual(np.abs(self.lx).max(),
                             self.case["rmax"] * (1.0 + 1e-9))

    def test_independent_of_nchunks(self):
        """Each omega bin is owned by one chunk and each peak by one cell, so
        the result must not depend on how the work is split."""
        for nch in (1, 2, 3, 7):
            with self.subTest(nchunks=nch):
                lx, _ = self._run(self.case, nchunks=nch)
                self.assertTrue(
                    np.array_equal(lx, self.ref),
                    "max |dx| %.3g vs the brute-force reference"
                    % np.abs(lx - self.ref).max())

    def test_independent_of_omega_binsize(self):
        """The omega bins are only a pre-filter: the exact per-peak dty
        predicate is applied inside, and ray_margin widens with the binsize,
        so a coarser binning can only be a superset. Identical results mean
        the margin was already sufficient at the finest binning. "auto" is
        VoxelSinoMasker's own heuristic, which is what production uses."""
        for bs in (0.5, 1.0, 2.0, 5.0, 360.0, "auto"):
            with self.subTest(omega_binsize=bs):
                lx, part = self._run(self.case, omega_binsize=bs)
                self._assert_matches_reference(
                    lx, self.ref, self.nclaim, self.case["rmax"],
                    label="binsize %s (ray_margin %.2f): "
                          % (bs, part["ray_margin"]))

    def test_independent_of_ymin(self):
        """ymin is ds.ybincens[0], and correct_bins_for_half_scan pads
        ybincens about y0 -- so on a half scan it sits several ystep below the
        smallest measured dty and the low dty bins are empty. It only sets
        where the bin origin is, so the answer must not move. The boundary
        VoxelSinoMasker allows, dty.min() - ystep/2, is included."""
        d0 = float(self.case["dty"].min())
        for label, ymin in (("dty.min()", d0),
                            ("dty.min() - ystep/2", d0 - 0.5 * YSTEP),
                            ("dty.min() - 1 ystep", d0 - YSTEP),
                            ("dty.min() - 17.5 ystep", d0 - 17.5 * YSTEP),
                            ("half-scan pad, 200 ystep", d0 - 200.0 * YSTEP)):
            with self.subTest(ymin=label):
                lx, part = self._run(self.case, ymin=ymin)
                self.assertTrue(
                    np.array_equal(lx, self.ref),
                    "%d peaks differ with ymin = %s (%d dty bins)"
                    % (int((lx != self.ref).sum()), label, part["ndty"]))

    def test_ymin_above_the_data_is_rejected(self):
        """The one ymin that is not allowed: above dty.min() by more than half
        a step, which would bin peaks to negative dty rows. VoxelSinoMasker
        raises on it and so does the reference builder."""
        d0 = float(self.case["dty"].min())
        with self.assertRaises(ValueError):
            self._run(self.case, ymin=d0 + 0.6 * YSTEP)

    def test_tolerance_sweep_stays_exact(self):
        """At every tolerance the output equals the reference exactly, and the
        number of peaks co-indexed by a second voxel grows monotonically with
        the tolerance -- which is why the value matters: the pre-fix code
        compared a squared residual to an unsquared tol, so a nominal 0.05 was
        really sqrt(0.05) = 0.224."""
        multi = {}
        for tol in (0.02, 0.05, 0.10, 0.20):
            lx, _ = self._run(self.case, hkl_tol=tol)
            ref, nclaim = compute_origins_ref(self.case, hkl_tol=tol)
            multi[tol] = float((nclaim > 1).mean())
            with self.subTest(hkl_tol=tol):
                self._assert_matches_reference(
                    lx, ref, nclaim, self.case["rmax"],
                    label="hkl_tol %g: " % tol)
                self.assertEqual(int((nclaim == 0).sum()), 0)
        tols = sorted(multi)
        for a, b in zip(tols, tols[1:]):
            self.assertLessEqual(multi[a], multi[b],
                                 "co-indexing fell from tol %g to %g" % (a, b))
        self.assertLess(multi[0.02], 0.01)
        self.assertGreater(multi[0.20], 0.10)

    def test_many_voxels_per_grain(self):
        """The realistic case: a grain spans many voxels, so several claim each
        peak and the answer is a weighted mean rather than the single true
        voxel. Still exact against the reference, still every peak claimed,
        and the mean should sit on the truth rather than beside it."""
        case = self._make_case(one_grain_per_voxel=False)
        lx, _ = self._run(case)
        ref, nclaim = compute_origins_ref(case)
        self._assert_matches_reference(lx, ref, nclaim, case["rmax"])
        self.assertEqual(int((nclaim == 0).sum()), 0)
        self.assertGreater(float((nclaim > 1).mean()), 0.5,
                           "this case is meant to exercise multi-voxel means")
        err = lx - case["lx_true"]
        self.assertLess(abs(err.mean()), 0.25 * YSTEP,
                        "bias %.4g on a %g grid" % (err.mean(), YSTEP))
        self.assertLessEqual(np.abs(lx).max(), case["rmax"] * (1.0 + 1e-9))


if __name__ == "__main__":
    unittest.main(verbosity=2)