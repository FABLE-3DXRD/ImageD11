import numpy as np
import unittest

from ImageD11.sinograms import geometry


class TestSampleToLab(unittest.TestCase):
    def test_sample_to_lab(self):
        sx = 0.5  # um
        sy = np.sqrt(3) / 2  # um
        y0 = 35  # um
        dty = 34  # um
        omega = 30  # degrees
        desired_lx = 0.0  # um
        desired_ly = 0.0  # um

        # assuming no y0/dty problems, rotating by 30 degrees should bring the grain
        # onto the x-axis
        # with coords (1, 0)
        # however because y0 != dty
        # rotation axis is further negative y than expected
        # shift the sample to the right (negative y)
        # puts the grain on the origin

        lx, ly = geometry.sample_to_lab(sx, sy, y0, dty, omega)

        self.assertAlmostEqual(desired_lx, lx)
        self.assertAlmostEqual(desired_ly, ly)

    def test_lab_to_sample(self):
        # let's define a grain that's in the beam at 90 degrees
        lx = -1.0  # um
        ly = 0.0  # um
        y0 = 35  # um
        dty = 34  # um
        omega = 90.0  # degrees
        # dty is less than y0, so the rotation axis is further to the right (negative lab y) than expected
        desired_sx = 1.0
        desired_sy = 1.0

        # (-1, 0) in the lab
        # goes to (-1, 1) in the translated frame
        # then we rotate reference frame to match

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

        x, y = geometry.step_to_sample(si, sj, ystep)

        self.assertSequenceEqual((desired_x, desired_y), (x, y))

    def test_cycle(self):
        ystep = 2.0  # um/px
        si = 50.0  # px
        sj = -100.0  # px

        x, y = geometry.step_to_sample(si, sj, ystep)
        si_out, sj_out = geometry.sample_to_step(x, y, ystep)

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


class TestXYY0OmegaToDty(unittest.TestCase):
    def test_grain_position1(self):
        x = 1.0  # um
        y = 0.0  # um
        omega = 90.0  # degrees

        # a grain is at (1, 0) at 0 degrees
        # at 90 degrees, it rotates to (1, 0)
        # dty has to be -1 to bring it into the beam

        desired_dty = -1.0
        y0 = 0.0

        dty = geometry.x_y_y0_omega_to_dty(omega, x, y, y0)

        self.assertEqual(desired_dty, dty)

    def test_grain_position2(self):
        x = 0.0  # um
        y = 1.0  # um
        omega = 0.0  # degrees

        # a grain is at (0, 1) at 0 degrees
        # dty has to be -1 to bring it into the beam

        desired_dty = -1.0
        y0 = 0.0

        dty = geometry.x_y_y0_omega_to_dty(omega, x, y, y0)

        self.assertEqual(desired_dty, dty)

    def test_grain_position3(self):
        x = -1.0  # um
        y = 0.0  # um
        omega = 90.0  # degrees

        # a grain is at (-1, 0) at 0 degrees
        # at 90 degrees, it rotates to (0, -1)
        # dty has to be 1 to bring it into the beam

        desired_dty = 1.0
        y0 = 0.0

        dty = geometry.x_y_y0_omega_to_dty(omega, x, y, y0)

        self.assertEqual(desired_dty, dty)

    def test_grain_position4(self):
        x = 0.0  # um
        y = -1.0  # um
        omega = 0.0  # degrees

        # a grain is at (0, -1) at 0 degrees
        # dty has to be 1 to bring it into the beam

        desired_dty = 1.0
        y0 = 0.0

        dty = geometry.x_y_y0_omega_to_dty(omega, x, y, y0)

        self.assertEqual(desired_dty, dty)

    def test_grain_position_with_y0(self):
        x = 0.0  # um
        y = 1.0  # um
        omega = 0.0  # degrees

        # a grain is at (0, 1) at 0 degrees
        # if y0 = 0, dty has to be -1 to bring it into the beam
        # however y0 = 5
        # dty has to be at 5 to bring rotation axis into the beam
        # true desired value is 5 - 1 = 4

        desired_dty = 4.0
        y0 = 5.0

        dty = geometry.x_y_y0_omega_to_dty(omega, x, y, y0)

        self.assertEqual(desired_dty, dty)


class TestDtyOmegaToXYY0(unittest.TestCase):
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

        desired_x, desired_y = (-1, 0)
        desired_y0 = 0

        x, y, y0 = geometry.dty_omega_to_x_y_y0(dty, omega)

        self.assertAlmostEqual(desired_x, x)
        self.assertAlmostEqual(desired_y, y)
        self.assertAlmostEqual(desired_y0, y0)

    def test_cos(self):
        # define a cos wave
        omega = np.arange(0, 180, 1)  # degrees
        dty = np.cos(np.radians(omega))
        # dty starts at 1
        # so grain is at (0, -1) at omega = 0

        desired_x, desired_y = (0, -1)
        desired_y0 = 0

        x, y, y0 = geometry.dty_omega_to_x_y_y0(dty, omega)

        self.assertAlmostEqual(desired_x, x)
        self.assertAlmostEqual(desired_y, y)
        self.assertAlmostEqual(desired_y0, y0)

    def test_cos_with_y0(self):
        # define a cos wave
        omega = np.arange(0, 180, 1)  # degrees
        desired_y0 = 5
        dty = desired_y0 + np.cos(np.radians(omega))
        # dty starts at 1
        # so grain is at (0, -1) at omega = 0
        # offset in dty should not affect this

        desired_x, desired_y = (0, -1)

        x, y, y0 = geometry.dty_omega_to_x_y_y0(dty, omega)

        self.assertAlmostEqual(desired_x, x)
        self.assertAlmostEqual(desired_y, y)
        self.assertAlmostEqual(desired_y0, y0)


class TestDtyToDtyi(unittest.TestCase):
    def test_simple_sequence(self):
        dty = np.array([0., 2., 4., 6., 8])  # um
        ystep = 2.0  # um/px
        desired_dtyi = np.array([0, -1, -2, -3, -4])
        dtyi = geometry.dty_to_dtyi(dty, ystep)

        self.assertTrue(np.allclose(desired_dtyi, dtyi))


class TestDtyToDtyiForSinogram(unittest.TestCase):
    def test_simple_sequence(self):
        dty = np.array([0., 2., 4., 6., 8])  # um
        ystep = 2.0  # um/px
        ymin = -2.0  # um
        desired_dtyi = np.array([1, 2, 3, 4, 5])
        dtyi = geometry.dty_to_dtyi_for_sinogram(dty, ystep, ymin)

        self.assertTrue(np.allclose(desired_dtyi, dtyi))


class TestFitSamplePositionFromRecon(unittest.TestCase):
    def test_simple_recon(self):
        ni, nj = (100, 200)  # size of recon
        recon = np.zeros((ni, nj), dtype=float)  # blank reconstruction image
        recon_centre = (ni // 2, nj // 2)  # find centre of recon
        centre_offset = 10  # 10 px away from centre
        ri, rj = (recon_centre[0] + centre_offset, recon_centre[1] - centre_offset)  # recon position of hot pixel
        recon[ri, rj] = 1.0
        ystep = 2.0  # um/px
        desired_x, desired_y = 20.0, 20.0  # 20 px away from centre
        # remember y sign is flipped!
        # we are slightly up and to the left of centre in recon space
        # slightly +ve in both x and y
        # centre_offset px worth - that's 20 um
        # x, y should be (20, 20)
        x, y = geometry.fit_sample_position_from_recon(recon, ystep)

        self.assertSequenceEqual((desired_x, desired_y), (x, y))


class TestDtyMask(unittest.TestCase):
    def setUp(self):
        om = np.arange(180)
        self.ystep = 0.1
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
                    m1 = geometry.dtyimask_from_step(si, sj, self.omega, self.dtyi, y0, self.ystep)
                    m2 = geometry.dtyimask_from_sincos(si, sj, self.sinomega, self.cosomega, self.dtyi, y0, self.ystep)
                    self.assertTrue((m1 == m2).all())


if __name__ == "__main__":
    unittest.main()
