

from ImageD11 import transform
import numpy as np
import unittest

class test_uncomputegv( unittest.TestCase ):

    def setUp(self):
        self.tth = np.array( [60, 10, 15, 20,25, 1.9555],np.float)
        self.wvln = 0.5
        self.eta = np.array( [90, 10,120,-20,340, -73 ],np.float)
        self.omega = np.array([60,90,180, 60,97,  131],np.float)
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
        #print tth
        #print eta
        #print omega
        for i in range(len(tth)):
            print
            print self.w, self.c,
            print i,self.tth[i], tth[i],
            best = np.argmin( np.array(eta)[:,i] - self.eta[i] )
            print i,self.tth[i], tth[i],self.eta[i], eta[best][i],\
            self.omega[i], omega[best][i]
            self.assertAlmostEqual( self.tth[i], tth[i], 6 )
            self.assertAlmostEqual( self.eta[i], eta[best][i], 6 )
            self.assertAlmostEqual( self.omega[i], omega[best][i], 6 )
            
if __name__ ==  "__main__":
    unittest.main()


