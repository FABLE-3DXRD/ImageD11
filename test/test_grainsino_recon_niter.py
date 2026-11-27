import sys
import unittest

if int(sys.version_info.major) == 2:
    raise unittest.SkipTest("Skipping sinogram tests on Python 2")

import numpy as np
import ImageD11.grain
import ImageD11.sinograms.dataset
import ImageD11.sinograms.roi_iradon
import ImageD11.sinograms.sinogram as sinogram


class test_recon_niter(unittest.TestCase):
    """recon_niter reaches the iterative backends, astra included (#676)"""

    def setUp(self):
        self.calls = {}
        self.saved = sinogram.run_astra, ImageD11.sinograms.roi_iradon.run_mlem

        def fake(name):
            def f(sino, angles, **kwargs):
                self.calls[name] = kwargs
                return np.zeros((sino.shape[0], sino.shape[0]))
            return f
        sinogram.run_astra = fake("astra")
        ImageD11.sinograms.roi_iradon.run_mlem = fake("mlem")
        self.gs = sinogram.GrainSinogram(ImageD11.grain.grain(np.eye(3)),
                                         ImageD11.sinograms.dataset.DataSet())
        self.gs.ssino = np.ones((5, 4))
        self.gs.sinoangles = np.linspace(0, 180, 4, endpoint=False)

    def tearDown(self):
        sinogram.run_astra, ImageD11.sinograms.roi_iradon.run_mlem = self.saved

    def test_recon_niter(self):
        self.gs.update_recon_parameters(niter=7)
        for method in ("astra", "mlem"):
            self.gs.recon(method=method)
            self.assertEqual(self.calls[method]["niter"], 7, method)

    def test_explicit_niter_wins(self):
        self.gs.update_recon_parameters(niter=7)
        self.gs.recon(method="astra", niter=3)
        self.assertEqual(self.calls["astra"]["niter"], 3)

    def test_default(self):
        self.gs.recon(method="astra")
        self.assertNotIn("niter", self.calls["astra"])


if __name__ == "__main__":
    unittest.main()
