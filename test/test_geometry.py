import sys

import numpy as np
import unittest

from ImageD11.sinograms import geometry


class TestSampleToLab(unittest.TestCase):
    def test_sample_to_lab(self):
        # set up a grain in the sample reference frame
        # such that assuming no alignment problems, rotating it by 30 degrees should bring it onto the y-axis
        # with (lx, ly) = (0, 1 + dty)

        sx = 0.5  # um
        sy = np.sqrt(3) / 2  # um
        omega = 30  # degrees

        # however we have an alignment problem with y0

        y0 = 101  # um - dty value where the rotation axis hits the beam
        dty = 100  # um

        # such that the grain actually ends up in the beam (ly = ybeam)

        desired_lx = 0.0  # um
        desired_ly = 0.0  # um

        lx, ly = geometry.sample_to_lab(sx, sy, y0, dty, omega)

        self.assertAlmostEqual(desired_lx, lx)
        self.assertAlmostEqual(desired_ly, ly)

    def test_lab_to_sample(self):
        # for now, assume no alignment problems
        y0 = 0  # um
        dty = 0  # um
        # set up a grain that's in the beam at 90 degrees
        lx = -1.0  # um
        ly = 0  # um
        omega = 90.0  # degrees
        # in the sample reference frame, should be (0, 1)
        desired_sx = 0.0
        desired_sy = 1.0

        sx, sy = geometry.lab_to_sample(lx, ly, y0, dty, omega)

        self.assertAlmostEqual(desired_sx, sx)
        self.assertAlmostEqual(desired_sy, sy)

        # now set up alignment problem
        # put dty at where we think the beam is
        # beam actually intersects at dty = y0
        # so rotation axis is slightly to the right (negative y)

        y0 = 36  # um
        dty = 35  # um

        # set up a grain that's in the beam at 90 degrees
        lx = -1.0  # um
        ly = 0.0  # um

        desired_sx = 1.0
        desired_sy = 1.0

        sx, sy = geometry.lab_to_sample(lx, ly, y0, dty, omega)

        self.assertAlmostEqual(desired_sx, sx)
        self.assertAlmostEqual(desired_sy, sy)

    def test_cycle(self):
        sx = 4.32
        sy = -6.49
        y0 = 42
        dty = -24
        omega = -456

        lx, ly = geometry.sample_to_lab(sx, sy, y0, dty, omega)
        sx_final, sy_final = geometry.lab_to_sample(lx, ly, y0, dty, omega)

        self.assertAlmostEqual(sx, sx_final)
        self.assertAlmostEqual(sy, sy_final)


class TestSampleToStep(unittest.TestCase):
    def test_sample_to_step(self):
        ystep = 2.0  # um/px
        x = 100.0  # um
        y = 200.0  # um
        desired_si = 50.0  # px
        desired_sj = -100.0  # px

        si, sj = geometry.sample_to_step(x, y, ystep)

        self.assertSequenceEqual((desired_si, desired_sj), (si, sj))

    def test_step_to_sample(self):
        ystep = 2.0  # um/px
        si = 50.0  # px
        sj = -100.0  # px
        desired_x = 100.0  # um
        desired_y = 200.0  # um

        sx, sy = geometry.step_to_sample(si, sj, ystep)

        self.assertSequenceEqual((desired_x, desired_y), (sx, sy))

    def test_cycle(self):
        ystep = 2.0  # um/px
        si = 50.0  # px
        sj = -100.0  # px

        sx, sy = geometry.step_to_sample(si, sj, ystep)
        si_out, sj_out = geometry.sample_to_step(sx, sy, ystep)

        self.assertSequenceEqual((si, sj), (si_out, sj_out))


class TestStepToRecon(unittest.TestCase):
    def test_step_to_recon(self):
        recon_shape = (200, 200)  # px
        si = 50  # px
        sj = 100  # px
        desired_ri = 150  # px
        desired_rj = 200  # px

        ri, rj = geometry.step_to_recon(si, sj, recon_shape)

        self.assertSequenceEqual((desired_ri, desired_rj), (ri, rj))

    def test_recon_to_step(self):
        recon_shape = (200, 200)  # px
        ri = 0.0  # px
        rj = 1.0  # px
        desired_si = -100.0  # px
        desired_sj = -99.0  # px

        si, sj = geometry.recon_to_step(ri, rj, recon_shape)

        self.assertSequenceEqual((desired_si, desired_sj), (si, sj))

    def test_cycle(self):
        recon_shape = (200, 200)  # px
        si = 50  # px
        sj = 100  # px

        ri, rj = geometry.step_to_recon(si, sj, recon_shape)
        si_out, sj_out = geometry.recon_to_step(ri, rj, recon_shape)

        self.assertSequenceEqual((si, sj), (si_out, sj_out))


class TestDtyValuesGrainInBeam(unittest.TestCase):
    def test_grain_position1(self):
        sx = 1.0  # um
        sy = 0.0  # um
        omega = 90.0  # degrees

        # a grain is at (1, 0) at 0 degrees
        # at 90 degrees, it rotates to (1, 0)
        # dty has to be -1 to bring it into the beam

        desired_dty = -1.0
        y0 = 0.0

        dty = geometry.dty_values_grain_in_beam(sx, sy, y0, omega)

        self.assertEqual(desired_dty, dty)

    def test_grain_position2(self):
        sx = 0.0  # um
        sy = 1.0  # um
        omega = 0.0  # degrees

        # a grain is at (0, 1) at 0 degrees
        # dty has to be -1 to bring it into the beam

        desired_dty = -1.0
        y0 = 0.0

        dty = geometry.dty_values_grain_in_beam(sx, sy, y0, omega)

        self.assertEqual(desired_dty, dty)

    def test_grain_position3(self):
        sx = -1.0  # um
        sy = 0.0  # um
        omega = 90.0  # degrees

        # a grain is at (-1, 0) at 0 degrees
        # at 90 degrees, it rotates to (0, -1)
        # dty has to be 1 to bring it into the beam

        desired_dty = 1.0
        y0 = 0.0

        dty = geometry.dty_values_grain_in_beam(sx, sy, y0, omega)

        self.assertEqual(desired_dty, dty)

    def test_grain_position4(self):
        sx = 0.0  # um
        sy = -1.0  # um
        omega = 0.0  # degrees

        # a grain is at (0, -1) at 0 degrees
        # dty has to be 1 to bring it into the beam

        desired_dty = 1.0
        y0 = 0.0

        dty = geometry.dty_values_grain_in_beam(sx, sy, y0, omega)

        self.assertEqual(desired_dty, dty)

    def test_grain_position_with_y0(self):
        sx = 0.0  # um
        sy = 1.0  # um
        omega = 0.0  # degrees

        # a grain is at (0, 1) at 0 degrees
        # if y0 = 0, dty has to be -1 to bring it into the beam
        # however y0 = 5
        # dty has to be at 5 to bring rotation axis into the beam
        # true desired value is 5 - 1 = 4

        desired_dty = 4.0
        y0 = 5.0

        dty = geometry.dty_values_grain_in_beam(sx, sy, y0, omega)

        self.assertEqual(desired_dty, dty)


class TestSxSyY0FromDtyOmega(unittest.TestCase):
    def test_sin(self):
        # define a sine wave
        omega = np.arange(0, 180, 1)  # degrees
        dty = np.sin(np.radians(omega))
        # starts at 0, 0
        # so grain is on rotation axis at that point
        # grain can only be in x somehow
        # omega increases as dty increases in first quadrant
        # therefore position should be (-1, 0)
        # as sample is translated in +y, need to rotate +ve to bring grain back into beam

        desired_sx, desired_sy = (-1, 0)
        desired_y0 = 0

        sx, sy, y0 = geometry.sx_sy_y0_from_dty_omega(dty, omega)

        self.assertAlmostEqual(desired_sx, sx)
        self.assertAlmostEqual(desired_sy, sy)
        self.assertAlmostEqual(desired_y0, y0)

    def test_cos(self):
        # define a cos wave
        omega = np.arange(0, 180, 1)  # degrees
        dty = np.cos(np.radians(omega))
        # dty starts at 1
        # so grain is at (0, -1) at omega = 0

        desired_sx, desired_sy = (0, -1)
        desired_y0 = 0

        sx, sy, y0 = geometry.sx_sy_y0_from_dty_omega(dty, omega)

        self.assertAlmostEqual(desired_sx, sx)
        self.assertAlmostEqual(desired_sy, sy)
        self.assertAlmostEqual(desired_y0, y0)

    def test_cos_with_y0(self):
        # define a cos wave
        omega = np.arange(0, 180, 1)  # degrees
        # we have a shift in the 'centre' of the cos wave
        desired_y0 = 5
        # shift cos graph in dty accordingly
        dty = desired_y0 + np.cos(np.radians(omega))
        # dty starts at 1
        # so grain is at (0, -1) at omega = 0
        # offset in dty should not be affected by y0

        desired_sx, desired_sy = (0, -1)

        sx, sy, y0 = geometry.sx_sy_y0_from_dty_omega(dty, omega)

        self.assertAlmostEqual(desired_sx, sx)
        self.assertAlmostEqual(desired_sy, sy)
        self.assertAlmostEqual(desired_y0, y0)


class TestDtyToDtyi(unittest.TestCase):
    def test_simple_sequence(self):
        ystep = 0.1  # um/px
        ymin = 13.5  # um
        ymax = 14.5  # um
        dty = np.arange(ymin, ymax + ystep, ystep)  # um

        desired_dtyi = np.arange(0, len(dty))
        dtyi = geometry.dty_to_dtyi(dty, ystep, ymin)

        self.assertTrue(np.allclose(desired_dtyi, dtyi))

        dty_back = geometry.dtyi_to_dty(dtyi, ystep, ymin)
        self.assertTrue(np.allclose(dty, dty_back))


class TestFitSamplePositionFromRecon(unittest.TestCase):
    def test_simple_recon(self):
        ni, nj = (100, 200)  # size of recon
        recon = np.zeros((ni, nj), dtype=float)  # blank reconstruction image
        recon_centre = (ni // 2, nj // 2)  # find centre of recon
        centre_offset = 10  # 10 px away from centre
        ri, rj = (recon_centre[0] + centre_offset, recon_centre[1] - centre_offset)  # recon position of hot pixel
        recon[ri, rj] = 1.0
        ystep = 2.0  # um/px
        desired_sx, desired_sy = 20.0, 20.0  # 20 px away from centre
        # remember y sign is flipped!
        # we are slightly up and to the left of centre in recon space
        # slightly +ve in both x and y
        # centre_offset px worth - that's 20 um
        # x, y should be (20, 20)
        sx, sy = geometry.fit_sample_position_from_recon(recon, ystep)

        self.assertSequenceEqual((desired_sx, desired_sy), (sx, sy))


class TestDtyiMaskFromSample(unittest.TestCase):
    def test_simple_value(self):
        y0 = 14.0
        sx = 1.0
        sy = 0.0
        # omega = 0 and omega = 180 should hit, no others
        # if we keep dty at 0
        omega = np.arange(0, 181, 1)
        dty = np.zeros_like(omega) + 14.0
        # make ystep small so only omega values very close to desired are accepted
        ystep = 0.01
        ymin = 14.0
        dtyi = geometry.dty_to_dtyi(dty, ystep, ymin)
        # true values should be where omega is a multiple of 180
        desired_mask = np.mod(omega, 180) == 0
        mask = geometry.dtyimask_from_sample(sx, sy, omega, dtyi, y0, ystep, ymin)

        self.assertTrue(np.allclose(desired_mask, mask))

    def test_weird_angle(self):
        y0 = 14.0
        sx = 1.0
        sy = 1.0
        # omega = 135 degrees should hit, no others
        # if we keep dty at 0
        omega = np.arange(0, 181, 1)
        dty = np.zeros_like(omega) + 14.0
        # make ystep small so only omega values very close to desired are accepted
        ystep = 0.01
        ymin = 14.0
        dtyi = geometry.dty_to_dtyi(dty, ystep, ymin)
        # true values should be where omega is a multiple of 180
        desired_mask = omega == 135
        mask = geometry.dtyimask_from_sample(sx, sy, omega, dtyi, y0, ystep, ymin)

        self.assertTrue(np.allclose(desired_mask, mask))


class TestDtyMask(unittest.TestCase):
    def setUp(self):
        om = np.arange(180)
        self.ystep = 0.1
        self.ymin = -5
        y = np.arange(-5, 5.1, self.ystep)
        self.shape = (len(y), len(om))
        self.omega = np.empty(self.shape, float)
        self.dty = np.empty(self.shape, float)
        self.omega[:] = om[np.newaxis, :]
        self.dty[:] = y[:, np.newaxis]
        self.dtyi = self.dty / self.ystep
        self.sinomega = np.sin(np.radians(self.omega))
        self.cosomega = np.cos(np.radians(self.omega))

    def test_cos_sin(self):
        x = 12.0
        y = 13.0
        for y0 in (-10.05, 0):
            for si in (-12, 0, 13):
                for sj in (-11, 0, 2):
                    m1 = geometry.dtyimask_from_step(si, sj, self.omega, self.dtyi, y0, self.ystep, self.ymin)
                    m2 = geometry.dtyimask_from_step_sincos(si, sj, self.sinomega, self.cosomega, self.dtyi, y0,
                                                            self.ystep, self.ymin)
                    self.assertTrue((m1 == m2).all())


@unittest.skipIf(int(sys.version_info.major) == 2, "can't test iradon on Python 2")
class TestFullLoop(unittest.TestCase):
    def test_full_loop(self):
        from ImageD11.sinograms.roi_iradon import run_iradon
        # tricky almost-half-acquisition scan
        # we scanned around 14 mm
        # but the true y0 is actually around 14.5 mm
        # also we didn't scan symmetrically
        y0 = 13.5 * 1000  # um

        ystep = 10.0
        ymin = 14 * 1000 - 750
        ymax = 14 * 1000 + 100

        yrange = ymax - ymin
        ny = int(yrange // ystep) + 1
        ybincens = np.linspace(ymin, ymax, ny)

        omega = np.arange(0, 361, 1)

        sx = 539.86
        sy = -510.25

        dty = geometry.dty_values_grain_in_beam(sx, sy, y0, omega)
        dtyi = geometry.dty_to_dtyi(dty, ystep, ybincens[0])

        # fill sinogram image
        sino = np.zeros((ny, len(omega)), dtype=float)

        for i in range(sino.shape[0]):
            for j in range(sino.shape[1]):
                this_dtyi = dtyi[j]
                sino[i, j] = 1 / (50 * np.cbrt((np.abs(i - this_dtyi))) + 0.01)

        shift, pad = geometry.sino_shift_and_pad(y0, ny, ymin, ystep)
        recon = run_iradon(sino, omega, pad=pad, shift=shift)
        ri_calc, rj_calc = np.array(np.where(recon == recon.max())).flatten()
        sx_calc, sy_calc = geometry.recon_to_sample(ri_calc, rj_calc, recon.shape, ystep)
        self.assertTrue(np.abs(sx - sx_calc) < ystep)
        self.assertTrue(np.abs(sy - sy_calc) < ystep)


if __name__ == "__main__":
    unittest.main()
