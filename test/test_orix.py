import numpy as np
import unittest

from ImageD11.unitcell import unitcell

import sys

if int(sys.version_info.major) == 2:
    # we're on Python 2, so skip these tests
    raise unittest.SkipTest("Skipping Orix tests on Python 2")
else:
    from orix.vector import Miller
    from orix.crystal_map import Phase
    from scipy.spatial.transform import Rotation as R


class TestOrixPhase(unittest.TestCase):
    """Tests unitcell.orix_phase"""

    def test_no_spacegroup(self):
        # set up hexagonal B matrix without a spacegroup
        hex_ucell_array = np.array([3, 3, 4, 90.0, 90.0, 120.0])
        hex_ucell = unitcell(hex_ucell_array, symmetry="P")

        # should fail because we didn't pass a spacegroup for the symmetry
        with self.assertRaises(AttributeError):
            phase = hex_ucell.orix_phase

    def test_with_spacegroup(self):
        # set up hexagonal B matrix with a spacegroup
        hex_ucell_array = np.array([3, 3, 4, 90.0, 90.0, 120.0])
        hex_ucell = unitcell(hex_ucell_array, symmetry=194)

        # we should have an orix phase object

        self.assertIsInstance(hex_ucell.orix_phase, Phase)


class TestGetOrixOrien(unittest.TestCase):
    """Tests unitcell.get_orix_orien"""

    def test_fast(self):
        # make sure fast matches original

        # print("")
        # set up hexagonal B matrix
        hex_ucell_array = np.array([3, 3, 4.0, 90.0, 90.0, 120.0])
        hex_ucell = unitcell(hex_ucell_array, symmetry=194)

        B = hex_ucell.B

        # set up reference unit cell
        ref_hex_ucell_array = np.array([3, 3, 4.0, 90.0, 90.0, 120.0])
        ref_hex_ucell = unitcell(ref_hex_ucell_array, symmetry=194)

        # many random orientations
        Us = R.random(1000).as_matrix()

        UBs = np.dot(Us, B)

        orix_phase = ref_hex_ucell.orix_phase
        orix_orien_slow = ref_hex_ucell.get_orix_orien(UBs)
        orix_orien_fast = ref_hex_ucell.get_orix_orien_fast(Us)
        self.assertTrue(np.allclose(orix_orien_slow.to_matrix(), orix_orien_fast.to_matrix()))

    def test_no_spacegroup(self):
        # set up hexagonal B matrix without a spacegroup
        hex_ucell_array = np.array([3, 3, 4.1, 90.0, 90.0, 120.0])
        hex_ucell = unitcell(hex_ucell_array, symmetry=194)

        B = hex_ucell.B

        # set up reference unit cell
        ref_hex_ucell_array = np.array([3, 3, 4.0, 90.0, 90.0, 120.0])
        ref_hex_ucell = unitcell(ref_hex_ucell_array, symmetry="P")

        U = np.squeeze(R.from_euler("zyx", [-90, 0, -90], degrees=True).as_matrix())

        UB = np.dot(U, B)

        # should fail because we didn't pass a spacegroup for the symmetry
        with self.assertRaises(AttributeError):
            ori = ref_hex_ucell.get_orix_orien(UB)

    def test_hexagonal(self):
        # set up hexagonal B matrix
        hex_ucell_array = np.array([3, 3, 4.1, 90.0, 90.0, 120.0])
        hex_ucell = unitcell(hex_ucell_array, symmetry=194)

        B = hex_ucell.B

        # set up reference unit cell
        ref_hex_ucell_array = np.array([3, 3, 4.0, 90.0, 90.0, 120.0])
        ref_hex_ucell = unitcell(ref_hex_ucell_array, symmetry=194)

        # set up U matrix

        # we want:
        # a* to point in +Z lab
        # b* to point in +X lab
        # c* to point in +Y lab

        U = np.squeeze(R.from_euler("zyx", [-90, 0, -90], degrees=True).as_matrix())

        UB = np.dot(U, B)

        orix_phase = ref_hex_ucell.orix_phase
        orix_orien = ref_hex_ucell.get_orix_orien(UB)

        # a* should be pointing in lab Z

        # get direction of (100) normal in cartesian lab frame as Numpy array
        ra_direc_lab = np.array(
            (~orix_orien * Miller(hkl=(1, 0, 0), phase=orix_phase)).xyz
        ).T
        lab_z = np.array([0, 0, 1])

        ra_direc_lab_unit = ra_direc_lab / np.linalg.norm(ra_direc_lab)

        self.assertTrue(np.allclose(ra_direc_lab_unit, lab_z))

    def test_hexagonal_many(self):
        # print("")
        # set up hexagonal B matrix
        hex_ucell_array = np.array([3, 3, 4.0, 90.0, 90.0, 120.0])
        hex_ucell = unitcell(hex_ucell_array, symmetry=194)

        B = hex_ucell.B

        # set up reference unit cell
        ref_hex_ucell_array = np.array([3, 3, 4.0, 90.0, 90.0, 120.0])
        ref_hex_ucell = unitcell(ref_hex_ucell_array, symmetry=194)

        # many random orientations
        Us = R.random(1000).as_matrix()

        UBs = np.dot(Us, B)

        orix_phase = ref_hex_ucell.orix_phase
        orix_orien = ref_hex_ucell.get_orix_orien(UBs)

        # work out a* direction in the lab frame for each orientation

        # get a random vector of a* in reciprocal space
        astar = np.array([1, 2, 3])
        # get Miller vector for this
        astar_miller = Miller(hkl=astar, phase=orix_phase)
        # work out a* directions in lab frame using ImageD11
        astar_lab_id11 = np.dot(UBs, astar)
        # work out a* directions in lab frame using Orix
        astar_lab_orix = np.array((~orix_orien).outer(astar_miller).xyz).T[0]

        # are they essentially the same:
        self.assertTrue(np.allclose(astar_lab_id11, astar_lab_orix))

    def test_one_equals_many(self):
        """
        Check that the orientation we get back with one orientation is the same
        as for an array (the code path is different)
        """
        Umats = [m.as_matrix() for m in R.random(10)]
        for uc in (
            [3, 3, 3, 90, 90, 90],  # cubic
            [3, 3, 4, 90, 90, 120],  # hexagonal
            [3, 3, 3, 65, 65, 65],  # trigonal   PROBLEM CASE!!! which number/setting?
            [3, 3, 4, 90, 90, 90],  # tetragonal
            [3, 4, 5, 90, 90, 90],  # orthorhombic
            [3, 4, 5, 95, 90, 90],  # monoclinic a
            [3, 4, 5, 90, 95, 90],  # monoclinic b
            [3, 4, 5, 90, 90, 95],  # monoclinic c
            [3, 4, 5, 91, 92, 93],  # triclinic
        ):
            for sgno in (
                2,  # P-1,
                10,  # P2/m
                123,  # P4/mmm
                162,  # R-3m
                191,  # P6/mm
                221,  # Pm3m
            ):
                cell = unitcell(uc, sgno)
                UBs = [m.dot(cell.B) for m in Umats]
                # convert one at a time:
                o1 = [cell.get_orix_orien(m).to_matrix() for m in UBs]
                # convert as an array:
                o2 = [m.to_matrix() for m in cell.get_orix_orien(UBs)]
                self.assertTrue(np.allclose(o1, o2))


class TestIPFCol(unittest.TestCase):
    def test_cubic(self):
        # set up B matrix
        cubic_ucell_array = np.array([3.0, 3.0, 3.0, 90.0, 90.0, 90.0])
        cubic_ucell = unitcell(cubic_ucell_array, symmetry="F")

        B = cubic_ucell.B

        # set up reference unit cell
        ref_cubic_ucell_array = np.array([3, 3, 3, 90.0, 90.0, 90.0])
        ref_cubic_ucell = unitcell(ref_cubic_ucell_array, symmetry=225)

        # set up U matrix

        # we want:
        # a* to point in +Z lab
        # b* to point in +X lab
        # c* to point in +Y lab

        U = np.squeeze(R.from_euler("zyx", [-90, 0, -90], degrees=True).as_matrix())

        UB = np.dot(U, B)

        # each of the 100, 010, 001, etc directions should have red colours (cubic symmetry)

        axes = np.array(
            [[1.0, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1]]
        )

        red = np.array([[1.0, 0.0, 0.0]])

        for axis in axes:
            rgb = ref_cubic_ucell.get_ipf_colour(UBs=UB, axis=axis)
            self.assertTrue(np.allclose(rgb, red))

    def test_cubic_110(self):
        cubic_ucell_array = np.array([3.0, 3.0, 3.0, 90.0, 90.0, 90.0])
        cubic_ucell = unitcell(cubic_ucell_array, symmetry="F")

        B = cubic_ucell.B

        # set up reference unit cell
        ref_cubic_ucell_array = np.array([3, 3, 3, 90.0, 90.0, 90.0])
        ref_cubic_ucell = unitcell(ref_cubic_ucell_array, symmetry=225)

        # set up U matrix

        # we want:
        # some (110) normal to point along X
        # rotate 45 degrees in Z
        # should mean 0.5 e1 + 0.5 e2 + 0.0 e3

        U = np.squeeze(R.from_euler("z", [45.0], degrees=True).as_matrix())

        UB = np.dot(U, B)

        # 110 should be pointing along X
        # so in an IPF-X plot, the colour should be green

        rgb = ref_cubic_ucell.get_ipf_colour(UBs=UB, axis=np.array([1, 0, 0]))
        green = np.array([0.0, 1.0, 0.0])
        self.assertTrue(np.allclose(rgb, green, atol=0.01))

    def test_cubic_111(self):
        # set up B matrix
        cubic_ucell_array = np.array([3.0, 3.0, 3.0, 90.0, 90.0, 90.0])
        cubic_ucell = unitcell(cubic_ucell_array, symmetry="F")

        B = cubic_ucell.B

        # set up reference unit cell
        ref_cubic_ucell_array = np.array([3, 3, 3, 90.0, 90.0, 90.0])
        ref_cubic_ucell = unitcell(ref_cubic_ucell_array, symmetry=225)

        # set up U matrix

        # we want:
        # some (111) normal to point along X
        # rotate 45 degrees in Z

        # You need 111 along an axis, so the angle between 111 and the 110 axis is:
        cosangle = (
            np.dot((1, 1, 1), (1, 1, 0))
            / np.linalg.norm((1, 1, 1))
            / np.linalg.norm((1, 1, 0))
        )
        angle = np.degrees(np.arccos(cosangle))
        U = R.from_euler("ZYZ", (0, angle, -45), degrees=True).as_matrix()
        UB = np.dot(U, B)
        # 111 should be pointing along X : check this
        g111 = UB.dot((1, 1, 1))
        self.assertTrue(np.allclose((1, 0, 0), g111 / np.linalg.norm(g111)))
        # so in an IPF-X plot, the colour should be blue-ish
        rgb = ref_cubic_ucell.get_ipf_colour(UBs=UB, axis=np.array([1, 0, 0]))
        blue = np.array([0.0, 0.0, 1.0])
        self.assertTrue(np.allclose(rgb, blue, atol=0.01))

    def test_hexagonal_fast(self):
        # set up hexagonal B matrix
        hex_ucell_array = np.array([3, 3, 4, 90.0, 90.0, 120.0])
        hex_ucell = unitcell(hex_ucell_array, symmetry="P")

        B = hex_ucell.B

        # set up reference unit cell
        ref_hex_ucell_array = np.array([3, 3, 4.0, 90.0, 90.0, 120.0])
        ref_hex_ucell = unitcell(ref_hex_ucell_array, symmetry=194)

        # set up U matrix

        # we want:
        # a* to point in +Z lab
        # b* to point in +X lab
        # c* to point in +Y lab

        U = np.squeeze(R.from_euler("zyx", [-90, 0, -90], degrees=True).as_matrix())
        UB = np.dot(U, B)
        
        rgb_slow = ref_hex_ucell.get_ipf_colour(UBs=UB, axis=np.array([0, 0, 1]))
        rgb_fast = ref_hex_ucell.get_ipf_colour_fast(Us=U, axis=np.array([0, 0, 1]))
        self.assertTrue(np.allclose(rgb_slow, rgb_fast))
        

        # on an IPF-Z:
        # a* should be parallel to Z
        # (100) (10-10) should be along Z
        # on TSL, (10-10) is blue

        rgb_slow = ref_hex_ucell.get_ipf_colour(UBs=UB, axis=np.array([0, 0, 1]))
        rgb_fast = ref_hex_ucell.get_ipf_colour_fast(Us=U, axis=np.array([0, 0, 1]))
        blue = np.array([0.0, 0.0, 1.0])
        self.assertTrue(np.allclose(rgb_slow, blue, atol=0.01))
        self.assertTrue(np.allclose(rgb_fast, blue, atol=0.01))

        # on an IPF-Y, c* should be parallel to Y
        # c* is (0001) normal

        rgb_slow = ref_hex_ucell.get_ipf_colour(UBs=UB, axis=np.array([0, 1, 0]))
        rgb_fast = ref_hex_ucell.get_ipf_colour_fast(Us=U, axis=np.array([0, 1, 0]))
        red = np.array([1.0, 0.0, 0.0])

        self.assertTrue(np.allclose(rgb_slow, red, atol=0.01))
        self.assertTrue(np.allclose(rgb_fast, red, atol=0.01))

    def test_hexagonal(self):
        # set up hexagonal B matrix
        hex_ucell_array = np.array([3, 3, 4, 90.0, 90.0, 120.0])
        hex_ucell = unitcell(hex_ucell_array, symmetry="P")

        B = hex_ucell.B

        # set up reference unit cell
        ref_hex_ucell_array = np.array([3, 3, 4.0, 90.0, 90.0, 120.0])
        ref_hex_ucell = unitcell(ref_hex_ucell_array, symmetry=194)

        # set up U matrix

        # we want:
        # a* to point in +Z lab
        # b* to point in +X lab
        # c* to point in +Y lab

        U = np.squeeze(R.from_euler("zyx", [-90, 0, -90], degrees=True).as_matrix())

        UB = np.dot(U, B)

        # on an IPF-Z:
        # a* should be parallel to Z
        # (100) (10-10) should be along Z
        # on TSL, (10-10) is blue

        rgb = ref_hex_ucell.get_ipf_colour(UBs=UB, axis=np.array([0, 0, 1]))
        blue = np.array([0.0, 0.0, 1.0])
        self.assertTrue(np.allclose(rgb, blue, atol=0.01))

        # on an IPF-Y, c* should be parallel to Y
        # c* is (0001) normal

        rgb = ref_hex_ucell.get_ipf_colour(UBs=UB, axis=np.array([0, 1, 0]))
        red = np.array([1.0, 0.0, 0.0])

        self.assertTrue(np.allclose(rgb, red, atol=0.01))


class TestOrixAlignment(unittest.TestCase):
    """Symmetry-equivalent descriptions of a grain are the same orientation for orix (issue 674)"""

    MONO = [3, 4, 5, 90, 100, 90]
    HEX = [4.9, 4.9, 5.4, 90, 90, 120]
    RHOMB = [5, 5, 5, 70, 70, 70]
    CASES = [(MONO, 4), (MONO, 14), (MONO, 15), ([3, 4, 5, 90, 90, 90], 62),
             ([3, 3, 5, 90, 90, 90], 88), ([3, 3, 5, 90, 90, 90], 139),
             (HEX, 148), (HEX, 149), (HEX, 150), (HEX, 157), (HEX, 162), (HEX, 176), (HEX, 194),
             ([4.9, 4.9, 13, 90, 90, 120], 166), (RHOMB, 148), (RHOMB, 166),
             ([4, 4, 4, 90, 90, 90], 205), ([4, 4, 4, 90, 90, 90], 225), ([3, 4, 5, 80, 85, 95], 2)]

    def equivalents(self, cell):
        from ImageD11 import sym_u
        return sym_u.cartesian(sym_u.laue_rotations(cell.spacegroup, cell.cell_choice), cell.B)

    def test_same_colour(self):
        from orix.vector import Vector3d
        from orix.plot import IPFColorKeyTSL
        Us = R.random(20, random_state=42).as_matrix()
        for lp, sg in self.CASES:
            cell = unitcell(lp, sg)
            for axis in np.eye(3):
                key = IPFColorKeyTSL(cell.orix_phase.point_group, direction=Vector3d(axis))
                rgb = key.orientation2color(cell.get_orix_orien_fast(Us))
                for S in self.equivalents(cell):
                    rgb_S = key.orientation2color(cell.get_orix_orien_fast(np.matmul(Us, S)))
                    # orix's colour key has ~1e-5 noise (also cubic); the bug was ~0.9
                    self.assertTrue(np.allclose(rgb, rgb_S, atol=1e-3), (sg, cell.cell_choice))

    def test_same_orientation(self):
        # centrosymmetric point groups, where the Laue rotations are all in the point group
        Us = R.random(20, random_state=42).as_matrix()
        for lp, sg in [(self.MONO, 14), (self.MONO, 15), (self.HEX, 162), (self.RHOMB, 166)]:
            cell = unitcell(lp, sg)
            o = cell.get_orix_orien_fast(Us)
            for S in self.equivalents(cell):
                angle = o.angle_with(cell.get_orix_orien_fast(np.matmul(Us, S)))
                self.assertTrue(np.allclose(angle, 0, atol=1e-6), sg)

    def test_identity_where_orix_agrees(self):
        for lp, sg in self.CASES:
            cell = unitcell(lp, sg)
            changed = cell.spacegroup in (4, 14, 15, 149, 157, 162) or cell.cell_choice == "rhombohedral"
            self.assertEqual(not np.allclose(cell.orix_alignment, np.eye(3)), changed, sg)

    def test_convention(self):
        # the choice orix leaves open (unitcell.orix_alignment): a* on orix's x, then b to +z
        for lp, sg in [(self.MONO, 14), (self.HEX, 162), (self.HEX, 149)]:
            cell = unitcell(lp, sg)
            M = cell.B.dot(np.linalg.inv(cell.orix_phase.structure.lattice.recbase))
            astar = cell.orix_alignment.dot(M[0])
            self.assertTrue(np.allclose(astar / np.linalg.norm(astar), (1, 0, 0)), sg)
        cell = unitcell(self.MONO, 14)
        M = cell.B.dot(np.linalg.inv(cell.orix_phase.structure.lattice.recbase))
        b = cell.orix_alignment.dot(M.T.dot(np.linalg.inv(cell.B)[1]))
        self.assertTrue(np.allclose(b / np.linalg.norm(b), (0, 0, 1)))

    def test_rhombohedral_as_hexagonal(self):
        # an R group on rhombohedral axes gets the colours of its hexagonal (obverse) axes
        from orix.vector import Vector3d
        from orix.plot import IPFColorKeyTSL
        Us = R.random(20, random_state=42).as_matrix()
        a, al = 5., np.radians(70.)
        hexagonal = [2 * a * np.sin(al / 2), 2 * a * np.sin(al / 2), a * np.sqrt(3 * (1 + 2 * np.cos(al))), 90, 90, 120]
        P = np.array([[1., -1., 0.], [0., 1., -1.], [1., 1., 1.]])
        for sg in (148, 166):
            r, h = unitcell(self.RHOMB, sg), unitcell(hexagonal, sg)
            N = h.B.dot(P).dot(np.linalg.inv(r.B)).T     # U on hexagonal axes = U.N
            for axis in np.eye(3):
                key = IPFColorKeyTSL(r.orix_phase.point_group, direction=Vector3d(axis))
                self.assertTrue(np.allclose(key.orientation2color(r.get_orix_orien_fast(Us)),
                                            key.orientation2color(h.get_orix_orien_fast(np.matmul(Us, N))),
                                            atol=1e-3), sg)

    def test_slow_is_fast(self):
        # get_orix_orien(UB) takes U as grain.U does, also for strained UBs
        from ImageD11.grain import grain
        cell = unitcell(self.MONO, 14)
        rng = np.random.RandomState(42)
        Us = R.random(20, random_state=42).as_matrix()
        eps = 1e-3 * rng.normal(size=(20, 3, 3))
        UBs = np.matmul(np.matmul(np.eye(3) + eps + np.swapaxes(eps, 1, 2), Us), cell.B)
        Ug = np.array([grain(np.linalg.inv(UB)).U for UB in UBs])
        self.assertTrue(np.allclose(cell.get_orix_orien(UBs).to_matrix(),
                                    cell.get_orix_orien_fast(Ug).to_matrix()))


if __name__ == "__main__":
    unittest.main()
