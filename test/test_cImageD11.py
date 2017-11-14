



from ImageD11 import cImageD11,  transform

import unittest
import numpy as np

TEST_FORTRAN=False
if TEST_FORTRAN:
    from ImageD11 import fImageD11

class test_compute_gv(unittest.TestCase):
    def setUp(self):
        N = 203

        self.XLYLZL = np.array([ np.linspace(0,2000,N),
                                 np.linspace(10,2010,N),
                                 np.ones(N) ])
        self.XLYLZLT = self.XLYLZL.T.copy()
        self.omega = np.linspace(-180,180,N)

    def runTest(self):
        tth, eta = transform.compute_tth_eta_from_xyz( self.XLYLZL,
                                                       self.omega*self.omegasign,
                                                       t_x=self.t_x,
                                                       t_y=self.t_y,
                                                       t_z=self.t_z,
                                                       wedge=self.wedge,
                                                       chi=self.chi)
        gve1 = transform.compute_g_vectors( tth, eta, self.omega*self.omegasign,
                                            self.wvln, wedge=self.wedge,
                                            chi=self.chi)
        t = np.array( (self.t_x, self.t_y, self.t_z))
        gve2 = np.empty( (len(tth), 3), np.float)
        cImageD11.compute_gv(self.XLYLZLT,
                             self.omega,self.omegasign,
                             self.wvln,
                             self.wedge,self.chi,
                             t,gve2)
        
        #        print t, self.wedge,self.chi
#        for i in range( len(tth)):
#            print i, gve1[:,i],gve1[:,i]-gve2[i],gve1[:,i]-gve3[:,i]
#        print "Test C"
        self.assertTrue( np.allclose(gve1, gve2.T ))
        #        print "Test Fortran"
        if TEST_FORTRAN:
            gve3 = np.empty( (3, len(tth)), np.float, order='F')
            fImageD11.compute_gv(self.XLYLZL,
                                 self.omega,self.omegasign,
                                 self.wvln,
                                 self.wedge,self.chi,
                                 t,
                                 gve3)

        
    def test_5_10(self):
        self.wedge = 5.
        self.chi = 10.
        self.omegasign=1.0
        self.wvln=0.125
        self.t_x=20.
        self.t_y=30.
        self.t_z=40.
        self.runTest()
        
    def test_5_0(self):
        self.wedge = 5.
        self.chi = 0.
        self.omegasign=1.0
        self.wvln=0.125
        self.t_x=20.
        self.t_y=30.
        self.t_z=40.
        self.runTest()

    def test_0_5(self):
        self.wedge = 0.
        self.chi = 5.
        self.omegasign=1.0
        self.wvln=0.125
        self.t_x=20.
        self.t_y=30.
        self.t_z=40.
        self.runTest()

        
    def test_5_10m(self):
        self.wedge = 5.
        self.chi = 10.
        self.omegasign=-1.0
        self.wvln=0.125
        self.t_x=200.
        self.t_y=300.
        self.t_z=40.
        self.runTest()
        
    def test_5_0_tt(self):
        self.wedge = 5.
        self.chi = 0.
        self.omegasign=1.0
        self.wvln=0.125
        self.t_x=200.
        self.t_y=300.
        self.t_z=-400.
        self.runTest()
        
    def test_0_0_t0(self):
        self.wedge = 0.
        self.chi = 0.
        self.omegasign=1.0
        self.wvln=0.125
        self.t_x=00.
        self.t_y=00.
        self.t_z=00.
        self.runTest()

    def test_0_0_t10(self):
        self.wedge = 0.
        self.chi = 0.
        self.omegasign=1.0
        self.wvln=0.125
        self.t_x=200.
        self.t_y=300.
        self.t_z=400.
        self.runTest()

        
if __name__ ==  "__main__":
    unittest.main()

