

from ImageD11 import transform
import numpy as np
import unittest

class test_uncomputegv( unittest.TestCase ):

    def setUp(self):
        self.tth = np.array(  [60, 10, 15,  20,  25, 1.9555],float)
        self.wvln = 0.5
        self.eta = np.array(  [20, 10, 120, -20, 340, -73 ],float)
        self.omega = np.array([60, 90, 180,  60, 97,  131],float)
        self.np = len(self.tth)

    def test_5_10(self):
        self.w = 5.0
        self.c = 10.0
        self.runtest()

    def test_5_0(self):
        self.w = 5.0
        self.c = 0.0
        self.runtest()

    def test_0_10(self):
        self.w = 0.0
        self.c = 10.0
        self.runtest()

    def test_0_0(self):
        self.w = 0.0
        self.c = 0.0
        self.runtest()

    def runtest(self):
        g = transform.compute_g_vectors( self.tth,
                self.eta,
                self.omega,
                self.wvln,
                wedge = self.w,
                chi = self.c
                )
        tth, eta, omega =  transform.uncompute_g_vectors( g, 
                self.wvln, 
                wedge = self.w,
                chi = self.c )
#        print tth
#        print eta
#        print omega
        
        #print "#  w  c   i     tth      tth     etao     etaC    omegao   omegaC    etaC2    omegaC2"
        for i in range(len(tth)):
         #   print self.w, self.c, i, 
            best = np.argmin( np.abs(angmod(np.array(omega)[:,i] - self.omega[i] )))
            deta = angmod(np.array(eta)[best,i] - self.eta[i] )
            domega = angmod(np.array(omega)[best,i] - self.omega[i] )

         #   print ("%8.2f "*8)%(self.tth[i], tth[i],self.eta[i], eta[best][i],\
         #       self.omega[i], omega[best][i],eta[1-best][i],omega[1-best][i])
            

            self.assertAlmostEqual( self.tth[i], tth[i], 5 )
            self.assertAlmostEqual( deta, 0.0, 5 )
            self.assertAlmostEqual( domega, 0.0, 5 )

def angmod(x):
    return np.degrees(np.arctan2( np.sin(np.radians(x)), np.cos(np.radians(x)) ))

if __name__ ==  "__main__":
    unittest.main()


