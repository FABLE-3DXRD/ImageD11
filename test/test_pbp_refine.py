import numpy as np
import unittest

import sys
if int(sys.version_info.major) == 2:
    raise unittest.SkipTest('Skipping PBP tests on Python 2')
else:
    from ImageD11.sinograms import pbp_refine
    from ImageD11 import transform


class TestDetectorRotationMatrix(unittest.TestCase):
    def test_random(self):
        # test 1000 times:
        for _ in range(1000):
            tilt_x = 0.01 * (np.random.random() * 2 - 1)
            tilt_y = 0.01 * (np.random.random() * 2 - 1)
            tilt_z = 0.01 * (np.random.random() * 2 - 1)

            mat_id11 = transform.detector_rotation_matrix(tilt_x, tilt_y, tilt_z)
            mat_numba = pbp_refine.detector_rotation_matrix(tilt_x, tilt_y, tilt_z)

            self.assertTrue(np.allclose(mat_id11, mat_numba))


class TestComputeXYZLab(unittest.TestCase):
    def test_random(self):
        # test 1000 times:
        for _ in range(1000):
            npeaks = 100
            sc = np.random.random(npeaks) * 2000
            fc = np.random.random(npeaks) * 2000
            y_center = np.random.random() * 1000
            z_center = np.random.random() * 1000
            y_size = 2048
            z_size = 2048
            tilt_x = 0.01 * (np.random.random() * 2 - 1)
            tilt_y = 0.01 * (np.random.random() * 2 - 1)
            tilt_z = 0.01 * (np.random.random() * 2 - 1)
            distance = np.random.random() * 1500
            
            peaks = np.vstack((sc, fc))
            
            rotvec_id11 = transform.compute_xyz_lab(peaks=peaks,
                                                    y_center=y_center,
                                                    z_center=z_center,
                                                    y_size=y_size,
                                                    z_size=z_size,
                                                    tilt_x=tilt_x,
                                                    tilt_y=tilt_y,
                                                    tilt_z=tilt_z,
                                                    distance=distance)
            
            rotvec_numba = pbp_refine.compute_xyz_lab(sc=sc,
                                                            fc=fc,
                                                            y_center=y_center,
                                                            z_center=z_center,
                                                            y_size=y_size,
                                                            z_size=z_size,
                                                            tilt_x=tilt_x,
                                                            tilt_y=tilt_y,
                                                            tilt_z=tilt_z,
                                                            distance=distance)

            self.assertTrue(np.allclose(rotvec_id11, rotvec_numba))

class TestComputeTthEtaFromXYZ(unittest.TestCase):
    def test_random(self):
        # test 1000 times:
        for _ in range(1000):
            npeaks = 100
            peaks_xyz = np.random.random((3, npeaks)) * 10000 - 5000
            wedge = np.random.random()
            chi = np.random.random()
            t_x = np.random.random()
            t_y = np.random.random()
            t_z = np.random.random()
            omega = np.random.random(npeaks) * 380 - 180
            
            tth_id11, eta_id11 = transform.compute_tth_eta_from_xyz(peaks_xyz=peaks_xyz, omega=omega, wedge=wedge, chi=chi, t_x=t_x, t_y=t_y, t_z=t_z)
            
            tth_numba, eta_numba = pbp_refine.compute_tth_eta_from_xyz(peaks_xyz=peaks_xyz, omega=omega, wedge=wedge, chi=chi, t_x=t_x, t_y=t_y, t_z=t_z)

            self.assertTrue(np.allclose(tth_id11, tth_numba))
            self.assertTrue(np.allclose(eta_id11, eta_numba))

class TestComputeTthEta(unittest.TestCase):
    def test_random(self):
        # test 1000 times:
        for _ in range(1000):
            npeaks = 100
            sc = np.random.random(npeaks) * 2000
            fc = np.random.random(npeaks) * 2000
            y_center = np.random.random() * 1000
            z_center = np.random.random() * 1000
            y_size = 2048
            z_size = 2048
            tilt_x = 0.01 * (np.random.random() * 2 - 1)
            tilt_y = 0.01 * (np.random.random() * 2 - 1)
            tilt_z = 0.01 * (np.random.random() * 2 - 1)
            distance = np.random.random() * 1500
            wedge = np.random.random()
            chi = np.random.random()
            t_x = np.random.random()
            t_y = np.random.random()
            t_z = np.random.random()
            omega = np.random.random(npeaks) * 380 - 180
            
            peaks = np.vstack((sc, fc))
            
            tth_id11, eta_id11 = transform.compute_tth_eta(peaks=peaks,
                                                    y_center=y_center,
                                                    z_center=z_center,
                                                    y_size=y_size,
                                                    z_size=z_size,
                                                    tilt_x=tilt_x,
                                                    tilt_y=tilt_y,
                                                    tilt_z=tilt_z,
                                                    distance=distance,
                                                    wedge=wedge,
                                                    chi=chi,
                                                    t_x=t_x,
                                                    t_y=t_y,
                                                    t_z=t_z,
                                                    omega=omega)
            
            tth_numba, eta_numba = pbp_refine.compute_tth_eta(sc=sc,
                                                            fc=fc,
                                                            y_center=y_center,
                                                            z_center=z_center,
                                                            y_size=y_size,
                                                            z_size=z_size,
                                                            tilt_x=tilt_x,
                                                            tilt_y=tilt_y,
                                                            tilt_z=tilt_z,
                                                            distance=distance,
                                                            wedge=wedge,
                                                            chi=chi,
                                                            t_x=t_x,
                                                            t_y=t_y,
                                                            t_z=t_z,
                                                            omega=omega)

            self.assertTrue(np.allclose(tth_id11, tth_numba))
            self.assertTrue(np.allclose(eta_id11, eta_numba))

            
class TestComputeKVectors(unittest.TestCase):
    def test_random(self):
        # test 1000 times:
        for _ in range(1000):
            npeaks = 100
            tth = np.random.random(npeaks) * 30
            eta = np.random.random(npeaks) * 360 - 180
            wvln = np.random.random()
            
            k_id11 = transform.compute_k_vectors(tth, eta, wvln)
            k_numba = pbp_refine.compute_k_vectors(tth, eta, wvln)
            
            self.assertTrue(np.allclose(k_id11, k_numba))

            
class TestComputeGFromK(unittest.TestCase):
    def test_random(self):
        # test 1000 times:
        for _ in range(1000):
            npeaks = 100
            k = np.random.random((3, npeaks))
            omega = np.random.random(npeaks) * 360 - 180
            wedge = np.random.random()
            chi = np.random.random()
            
            gve_id11 = transform.compute_g_from_k(k, omega, wedge, chi)
            gve_numba = pbp_refine.compute_g_from_k(k, omega, wedge, chi)
            
            self.assertTrue(np.allclose(gve_id11, gve_numba))

class TestComputeGVectors(unittest.TestCase):
    def test_random(self):
        # test 1000 times:
        for _ in range(1000):
            npeaks = 100
            tth = np.random.random(npeaks) * 30
            eta = np.random.random(npeaks) * 360 - 180
            omega = np.random.random(npeaks) * 360 - 180
            wvln = np.random.random()
            wedge = np.random.random()
            chi = np.random.random()
            
            gve_id11 = transform.compute_g_vectors(tth, eta, omega, wvln, wedge, chi)
            gve_numba = pbp_refine.compute_g_vectors(tth, eta, omega, wvln, wedge, chi)
            
            self.assertTrue(np.allclose(gve_id11, gve_numba))
            
if __name__ == "__main__":
    unittest.main()