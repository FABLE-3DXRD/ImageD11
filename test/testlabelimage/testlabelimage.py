import unittest
import os
import numpy as np

from ImageD11 import labelimage
from ImageD11 import columnfile


class test_moments(unittest.TestCase):
    def setUp(self):
        self.data = np.array(
           [[20, 20,  0,  0,  2,  0,  0,  0],
            [20, 20,  0,  0,  0,  0,  6,  0],
            [ 0,  0,  0,  0,  0,  0,  0,  0],
            [ 0,  0,  0,  0,  4,  0,  0, 10],
            [ 1, 20,  0,  2,  4,  10,  0, 0],
            [ 1,  0,  0,  0,  5,  0,  0,  0],
            [ 1,  0,  0,  0,  10, 0,  0,  0],
            [ 1,  0,  0,  0,  0,  0,  0,  0]],
           np.float32
            )
        self.outfile = "momentstest.out"
        lio = labelimage.labelimage(self.data.shape , self.outfile)
        # print "go"
        lio.peaksearch(self.data, 0.1, 1.)
        lio.mergelast()
        lio.peaksearch(self.data, 0.1, 2.)
        lio.mergelast()
        lio.finalise()
        self.lio = lio


    def tearDown(self):
        try:
            os.remove("momentest.out")
        except:
            pass


    def testpyvsc(self):
        from ImageD11.cImageD11 import s_1, s_I, s_I2, \
            s_fI, s_ffI, s_sI, s_ssI, s_sfI, s_oI, s_ooI, s_foI, s_soI, \
            bb_mn_f, bb_mn_s, bb_mx_f, bb_mx_s, bb_mn_o, bb_mx_o, \
            mx_I, mx_I_f, mx_I_s, mx_I_o, \
            avg_i, f_raw, s_raw, o_raw, f_cen, s_cen, \
            m_ss, m_ff, m_oo, m_sf, m_so, m_fo
        from math import sqrt
        #print lio.res.shape, lio.lastres.shape
        #print "res",lio.res[:,m_sf]
        #print "lastres",lio.lastres[:,m_sf]
        for i in self.lio.res:
            #print i
            #import sys
            #sys.exit()
            if i[s_1] < 0.1:
                # Merged with another
                continue
            avg_i_p = i[s_I]/i[s_1]
            self.assertAlmostEqual( avg_i_p, i[avg_i] , 6 )
            # first moments
            fraw_p = i[s_fI]/i[s_I]
            self.assertAlmostEqual( fraw_p, i[f_raw] , 6 )

            sraw_p = i[s_sI]/i[s_I]
            self.assertAlmostEqual( sraw_p, i[s_raw] , 6 )

            oraw_p = i[s_oI]/i[s_I]
            self.assertAlmostEqual( oraw_p, i[o_raw] , 6 )

            # second moments - the plus 1 is the zero width = 1 pixel
            mss_p = sqrt( i[s_ssI]/i[s_I] - sraw_p*sraw_p + 1 )
            self.assertAlmostEqual( mss_p, i[m_ss] , 6 )

            mff_p = sqrt( i[s_ffI]/i[s_I] - fraw_p*fraw_p + 1 )
            self.assertAlmostEqual( mff_p, i[m_ff] , 6 )

            moo_p = sqrt( i[s_ooI]/i[s_I] - oraw_p*oraw_p + 1 )
            self.assertAlmostEqual( moo_p, i[m_oo] , 6 )

            msf_p = ( i[s_sfI]/i[s_I] - sraw_p*fraw_p )/mss_p/mff_p
            self.assertAlmostEqual( msf_p, i[m_sf] , 6 )

            mso_p = ( i[s_soI]/i[s_I] - sraw_p*oraw_p )/mss_p/moo_p
            self.assertAlmostEqual( mso_p, i[m_so] , 6 )

            mfo_p = ( i[s_foI]/i[s_I] - fraw_p*oraw_p )/mff_p/moo_p
            self.assertAlmostEqual( mfo_p, i[m_fo] , 6 )


class test_more_moments(test_moments):
    def setUp(self):
        self.data1 = np.array(
           [[20, 20,  0,  0,  2,  0,  0,  0],
            [20, 20,  0,  0,  0,  0,  6,  0],
            [ 0,  0,  0,  0,  0,  0,  0,  0],
            [ 0,  0,  0,  0,  4,  0,  0, 10],
            [ 1, 20,  0,  2,  4,  10,  0, 0],
            [ 1,  0,  0,  0,  5,  0,  0,  0],
            [ 1,  0,  0,  0,  10, 0,  0,  0],
            [ 1,  0,  0,  0,  0,  0,  0,  0]],
           np.float32
            )
        self.data2 = np.array(
           [[10, 00,  0,  0,  2,  0,  0,  0],
            [00, 20,  0,  0,  0,  0,  6,  0],
            [ 0,  0,  0,  0,  0,  0,  0,  0],
            [ 0,  0,  0,  0,  4,  0,  0, 10],
            [ 1, 20,  0,  2,  4,  10,  0, 0],
            [ 1,  0,  0,  0,  5,  0,  0,  0],
            [ 1,  0,  0,  0,  10, 0,  0,  0],
            [ 1,  0,  0,  0,  0,  0,  0,  0]],
           np.float32
            )
        self.data3 = np.array(
           [[00, 00,  0,  0,  2,  0,  0,  0],
            [00, 10,  0,  0,  0,  0,  6,  0],
            [ 0,  0,  0,  0,  0,  0,  0,  0],
            [ 0,  0,  0,  0,  4,  0,  0, 10],
            [ 1, 20,  0,  2,  4,  10,  0, 0],
            [ 1,  0,  0,  0,  5,  0,  0,  0],
            [ 1,  0,  0,  0,  10, 0,  0,  0],
            [ 1,  0,  0,  0,  0,  0,  0,  0]],
           np.float32
            )
        self.outfile = "momentstest.out"
        lio = labelimage.labelimage(self.data1.shape , self.outfile)
        # print "go"
        lio.peaksearch(self.data1, 0.1, 1.)
        lio.mergelast()
        lio.peaksearch(self.data2, 0.1, 2.)
        lio.mergelast()
        lio.peaksearch(self.data3, 0.1, 3.)
        lio.mergelast()
        lio.finalise()
        self.lio = lio


    def tearDown(self):
        try:
            os.remove("momentest.out")
        except:
            pass



class test_labelimage(unittest.TestCase):
    def setUp(self):
        self.dims = (200,300)
        self.data1 = np.zeros(self.dims, np.float32)
        self.data2 = np.zeros(self.dims, np.float32)
        self.data3 = np.zeros(self.dims, np.float32)
        self.data4 = np.zeros(self.dims, np.float32)
        self.data5 = np.zeros(self.dims, np.float32)
        self.outfile = "l2.out"


    def test_1(self):
        self.data1[20:41,120:141] = 1.
        self.data2[20:41,120:141] = 1.
        self.data3[20:41,120:141] = 1.
        self.data4[20:41,120:141] = 1.
        self.data5[20:41,120:141] = 1.
        self.data3[30,130] = 2. # in the middle
        self.data1[55,65] = 1.
        self.data2[55,65] = 2.
        self.data3[55,65] = 3.
        self.data4[55,65] = 2.
        self.data5[55,65] = 1.

        lio = labelimage.labelimage(self.dims , self.outfile)
        lio.peaksearch(self.data1, 0.1, 1.)
        lio.mergelast()
        lio.peaksearch(self.data2, 0.1, 2.)
        lio.mergelast()
        lio.peaksearch(self.data3, 0.1, 3.)
        lio.mergelast()
        lio.peaksearch(self.data4, 0.1, 4.)
        lio.mergelast()
        lio.peaksearch(self.data5, 0.1, 5.)
        lio.mergelast()
        lio.finalise()
        lio.outfile.close()

        co = columnfile.columnfile(self.outfile)



        self.assertEqual(co.nrows,2)
        # first peaks 20->41, 120->141
        self.assertAlmostEqual(co.fc[0],130.,6)
        self.assertAlmostEqual(co.sc[0],30.,6)
        self.assertAlmostEqual(co.omega[0],3.0,6)
        self.assertAlmostEqual(co.Number_of_pixels[0],21*21*5,6)
        self.assertAlmostEqual(co.avg_intensity[0],1,1)
        self.assertAlmostEqual(co.IMax_int[0], 2 , 6)
        self.assertAlmostEqual(co.IMax_f[0], 130., 6)
        self.assertAlmostEqual(co.IMax_s[0], 30.0, 6)
        self.assertAlmostEqual(co.IMax_o[0], 3.0, 6)

        self.assertAlmostEqual(co.fc[1],65.,6)
        self.assertAlmostEqual(co.sc[1],55.,6)
        self.assertAlmostEqual(co.omega[1],3.0,6)
        self.assertAlmostEqual(co.Number_of_pixels[1],5,6)
        self.assertAlmostEqual(co.avg_intensity[1],1.8,1)
        self.assertAlmostEqual(co.IMax_int[1], 3 , 6)
        self.assertAlmostEqual(co.IMax_f[1], 65., 6)
        self.assertAlmostEqual(co.IMax_s[1], 55., 6)
        self.assertAlmostEqual(co.IMax_o[1], 3.0, 6)


        mask = co.Number_of_pixels > 10
        co.filter(mask)
        self.assertEqual(co.nrows, 1)
        co.writefile("l1.out")


if __name__=="__main__":
    unittest.main()
