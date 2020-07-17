
from __future__ import print_function, division


import unittest, os
import numpy as np
from ImageD11.grain import read_grain_file



class test_eps_grain( unittest.TestCase):
    def setUp(self):
        pth = os.path.split(__file__)[0]
        self.testdata = np.loadtxt(os.path.join(pth,"test.out"), skiprows = 30)
        self.titles = list( 
            "cell__a cell__b cell__c cell_alpha cell_beta cell_gamma u11 u12 u13 u21 u22 u23 u31 u32 u33".split() + 
            "eps11_c eps22_c eps33_c eps12_c eps13_c eps23_c eps11_s eps22_s eps33_s eps12_s eps13_s eps23_s".split() + 
            "sig11_c sig22_c sig33_c sig12_c sig13_c sig23_c sig11_s sig22_s sig33_s sig12_s sig13_s sig23_s".split() 
            )
        self.grains = read_grain_file(os.path.join(pth,"CuAlBe_scan10.map"))
        a = 5.77
        self.dzero = ( a, a, a, 90., 90., 90. )

    def test_e6_sample(self):
        t11 = self.testdata[:,self.titles.index("eps11_s")]/ 100.0
        t22 = self.testdata[:,self.titles.index("eps22_s")]/ 100.0
        t33 = self.testdata[:,self.titles.index("eps33_s")]/ 100.0
        t12 = self.testdata[:,self.titles.index("eps12_s")]/ 100.0
        t23 = self.testdata[:,self.titles.index("eps23_s")]/ 100.0
        t13 = self.testdata[:,self.titles.index("eps13_s")]/ 100.0
        for i, g in enumerate(self.grains):
            e11,e12,e13,e22,e23,e33 = g.eps_sample( self.dzero )
            ok =  np.allclose( 
                ( e11,    e12,    e13,    e22,    e23,    e33 ),
                ( t11[i], t12[i], t13[i], t22[i], t23[i], t33[i] ), rtol=.1, atol=1e-4)
            if not ok:
                print(i)
                for x,y in zip( ( e11,    e12,    e13,    e22,    e23,    e33 ),
                                ( t11[i], t12[i], t13[i], t22[i], t23[i], t33[i] )):
                    print(x,y,(x-y),(x-y)/(x+y))
            self.assertTrue( ok )

    def test_e6_grain(self):
        t11 = self.testdata[:,self.titles.index("eps11_c")]/ 100.0
        t22 = self.testdata[:,self.titles.index("eps22_c")]/ 100.0
        t33 = self.testdata[:,self.titles.index("eps33_c")]/ 100.0
        t12 = self.testdata[:,self.titles.index("eps12_c")]/ 100.0
        t23 = self.testdata[:,self.titles.index("eps23_c")]/ 100.0
        t13 = self.testdata[:,self.titles.index("eps13_c")]/ 100.0
        for i, g in enumerate(self.grains):
            e11,e12,e13,e22,e23,e33 = g.eps_grain( self.dzero )
            ok =  np.allclose( 
                ( e11,    e12,    e13,    e22,    e23,    e33 ),
                ( t11[i], t12[i], t13[i], t22[i], t23[i], t33[i] ),rtol=0.1, atol=1e-4)
            if not ok:
                print(i)
                for x,y in zip( ( e11,    e12,    e13,    e22,    e23,    e33 ),
                                ( t11[i], t12[i], t13[i], t22[i], t23[i], t33[i] )):
                    print(x,y,(x-y)/(x+y))
            self.assertTrue( ok )



if __name__=="__main__":
    unittest.main()
