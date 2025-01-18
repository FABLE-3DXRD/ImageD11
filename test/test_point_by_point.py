import numpy as np
import unittest

import sys
if int(sys.version_info.major) == 2:
    raise unittest.SkipTest('Skipping PBP tests on Python 2')
else:
    from ImageD11.sinograms import point_by_point
    from ImageD11 import transform


class TestDetectorRotationMatrix(unittest.TestCase):
    def test_random(self):
        # test 1000 times:
        for _ in range(1000):
            tilt_x = 0.01 * (np.random.random() * 2 - 1)
            tilt_y = 0.01 * (np.random.random() * 2 - 1)
            tilt_z = 0.01 * (np.random.random() * 2 - 1)

            mat_id11 = transform.detector_rotation_matrix(tilt_x, tilt_y, tilt_z)
            mat_numba = point_by_point.detector_rotation_matrix(tilt_x, tilt_y, tilt_z)

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
            
            rotvec_numba = point_by_point.compute_xyz_lab(sc=sc,
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
            
            tth_numba, eta_numba = point_by_point.compute_tth_eta_from_xyz(peaks_xyz=peaks_xyz, omega=omega, wedge=wedge, chi=chi, t_x=t_x, t_y=t_y, t_z=t_z)

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
            
            tth_numba, eta_numba = point_by_point.compute_tth_eta(sc=sc,
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
            k_numba = point_by_point.compute_k_vectors(tth, eta, wvln)
            
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
            gve_numba = point_by_point.compute_g_from_k(k, omega, wedge, chi)
            
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
            gve_numba = point_by_point.compute_g_vectors(tth, eta, omega, wvln, wedge, chi)
            
            self.assertTrue(np.allclose(gve_id11, gve_numba))
            
class TestCountUniquePeaks(unittest.TestCase):
    def test_random(self):
        # prepare some random input arrays
        npeaks = 1000
        # hkl integer with hmax 6
        hkli = np.random.randint(0, 6, size=(3, npeaks))
        # etasign (+1 or -1)
        etasign = np.random.randint(0, 2, size=(npeaks,))*2 - 1
        # dtyi integer -10 to +10 incl.
        dtyi = np.random.randint(-10, 11, size=(npeaks,))
        
        indices_numpy = np.unique(np.vstack((hkli, etasign, dtyi)), axis=1, return_inverse=True)[1]
        indices_numba = point_by_point.count_unique_peaks(hkli, etasign, dtyi)
        
        self.assertTrue(np.allclose(indices_numpy, indices_numba))

        
class TestCountUniquePeaksNoDtyi(unittest.TestCase):
    def test_random(self):
        # prepare some random input arrays
        npeaks = 1000
        # hkl integer with hmax 6
        hkli = np.random.randint(0, 6, size=(3, npeaks))
        # etasign (+1 or -1)
        etasign = np.random.randint(0, 2, size=(npeaks,))*2 - 1
        
        indices_numpy = np.unique(np.vstack((hkli, etasign)), axis=1, return_inverse=True)[1]
        indices_numba = point_by_point.count_unique_peaks_no_dtyi(hkli, etasign)
        
        self.assertTrue(np.allclose(indices_numpy, indices_numba))

        
class TestGveNorm(unittest.TestCase):
    def test_random(self):
        npeaks = 1000
        gves = np.random.random((3, npeaks))
        
        gve_norm_numpy = np.linalg.norm(gves, axis=0)
        gve_norm_numba = point_by_point.gve_norm(gves)
        
        self.assertTrue(np.allclose(gve_norm_numpy, gve_norm_numba))

        
class TestDivideWhere(unittest.TestCase):
    def test_random(self):
        arr1 = np.random.random((1000,1000))
        arr2 = np.random.random((1000,1000))
        out = np.ones_like(arr1)
        wherearr = np.random.randint(0, 2, (1000,1000))
        
        result_numpy = np.divide(arr1, arr2, out=out, where=wherearr!=0)
        result_numba = point_by_point.divide_where(arr1, arr2, out, wherearr)
        
        self.assertTrue(np.allclose(result_numpy, result_numba))
        
if __name__ == "__main__":
    unittest.main()
