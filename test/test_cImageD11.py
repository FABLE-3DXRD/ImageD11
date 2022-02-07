



from ImageD11 import cImageD11,  transform, parameters

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
        gve2 = np.empty( (len(tth), 3), float)
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
            gve3 = np.empty( (3, len(tth)), float, order='F')
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



class test_compute_xlylzl(unittest.TestCase):
    def setUp(self):
        npk = 10240
        self.sc = np.random.random( npk )*2048
        self.fc = np.random.random( npk )*2048
        self.pars = parameters.parameters()
        self.pars.set('z_center', 980. )
        self.pars.set('y_center', 1010. )
        self.pars.set('z_size', 48. )
        self.pars.set('y_size', 49. )
        self.pars.set('distance', 100100. )

    def test1(self):
        pks = self.sc, self.fc
        p = self.pars
        testtilts = [-0.5, 0., 0.24]
        tilts = [(x,y,z) for x in testtilts for y in testtilts for z in testtilts]
        pars = np.array( ( p.get("z_center"),
                           p.get("y_center"),
                           p.get("z_size"),
                           p.get("y_size")))
        dist = np.array( (p.get("distance"),0.,0.))
        ok = 0
        for o11,o12,o21,o22 in [ [ 1,0,0,1],
                                 [-1,0,0,1],
                                 [-1,0,0,-1],
                                 [ 1,0,0,-1],
                                 [ 0,1,1,0],
                                 [ 0,-1,1,0],
                                 [ 0,-1,-1,0],
                                 [ 0,1,-1,0]]:
            for tx,ty,tz in tilts:
                p.parameters['tilt_x']=tx
                p.parameters['tilt_y']=ty
                p.parameters['tilt_z']=tz
                p.parameters['o11']=o11
                p.parameters['o12']=o12
                p.parameters['o21']=o21
                p.parameters['o22']=o22

                xlylzl = transform.compute_xyz_lab( pks,
                                  **p.parameters )

                dmat = transform.detector_rotation_matrix(
                    p.get('tilt_x'), p.get('tilt_y'), p.get('tilt_z'))
                fmat = np.array( [[ 1,   0,   0],
                                  [ 0, p.get('o22'), p.get('o21')],
                                  [ 0, p.get('o12'), p.get('o11')]])
                r = np.dot(dmat, fmat)

                outxyz = np.zeros( (len(self.sc), 3), float )
                cImageD11.compute_xlylzl( self.sc,
                                          self.fc,
                                          pars, r.ravel(), dist,
                                          outxyz)
                error = np.abs(outxyz - xlylzl.T)
                self.assertTrue( error.max() < 1e-6 )
                ok += 1
        print("tested xlylzl %d times"%(ok))




if __name__ ==  "__main__":
    unittest.main()

