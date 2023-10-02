import unittest

import numpy as np

from ImageD11 import cImageD11


class test_connectedpixels(unittest.TestCase):
    def setUp(self):
        dims = (10,10)
        dark = np.ones(dims)
        flood = np.ones(dims)

    def test_1(self):
        for t in [np.uint8, np.int8, np.uint16, np.int16,
                  np.int32, np.uint32, np.float32, float]:
            data = np.array( [[ 0, 0, 0, 0, 0, 0, 0],
                             [ 0, 0, 0, 1, 0, 0, 0],
                             [ 0, 0, 0, 0, 0, 0, 0],
                             [ 0, 0, 0, 0, 0, 0, 0]],t)
            bl = np.zeros(data.shape,np.int32)
            cImageD11.connectedpixels(
                data.astype(np.float32),
                bl,
                0.1)  # threshold
            err = np.sum(np.ravel(data-bl))
            self.assertEqual(err, 0)

    def test_1_shape(self):
        for t in [np.uint8, np.int8, np.uint16, np.int16,
                  np.int32, np.uint32, np.float32, float]:
            data = np.array( [[ 1, 0, 1, 0, 1, 0, 1],
                             [ 1, 0, 1, 0, 1, 0, 1],
                             [ 1, 0, 1, 0, 0, 1, 0],
                             [ 1, 1, 1, 1, 1, 0, 1]],t)
            bl = np.zeros(data.shape,np.int32)
            cImageD11.connectedpixels(
                data.astype(np.float32),
                bl,
                0.1)  # threshold
            err = np.sum(np.ravel(data-bl))
            self.assertEqual(err, 0)

    def test_2_shapes(self):
        for t in [np.uint8, np.int8, np.uint16, np.int16,
                  np.int32, np.uint32, np.float32, float]:
            data = np.array( [[ 1, 0, 1, 0, 2, 0, 2],
                             [ 1, 0, 1, 0, 2, 0, 2],
                             [ 1, 0, 1, 0, 0, 2, 0],
                             [ 1, 1, 1, 0, 2, 0, 2]],t)
            bl = np.zeros(data.shape,np.int32)
            cImageD11.connectedpixels(
                data.astype(np.float32),
                bl,
                0.1)  # threshold
            err = np.sum(np.ravel(data-bl))
            self.assertEqual(err, 0)




class test_blobproperties(unittest.TestCase):

    def test_prop_names(self):
        names="""  s_1,       /* 1 Npix */
  s_I,       /* 2 Sum intensity */
  s_I2,      /* 3 Sum intensity^2 */
  s_fI,      /* 4 Sum f * intensity */
  s_ffI,     /* 5 Sum f * f* intensity */
  s_sI,      /* 6 Sum s * intensity */
  s_ssI,     /* 7 Sum s * s * intensity */
  s_sfI,     /* 8 Sum f * s * intensity */
  s_oI,         /* 9 sum omega * intensity */
  s_ooI,        /*  */
  s_soI,        /* 10 sum omega * s * intensity */
  s_foI,        /* 11 sum omega * f * intensity */
  mx_I,      /* 12  Max intensity */
  mx_I_f,    /* 13 fast at Max intensity */
  mx_I_s,    /* 14 slow at Max intensity */
  mx_I_o,    /* 15 omega at max I */
  bb_mx_f,      /* 16 max of f */
  bb_mx_s,      /* 17 max of s */
  bb_mx_o,      /* 18 max of omega */
  bb_mn_f,      /* 19 min of f */
  bb_mn_s,      /* 20 min of s */
  bb_mn_o,      /* 22 min of o */
  avg_i,   /* Average intensity */
  f_raw,       /* fast centre of mass */
  s_raw,       /* slow centre of mass */
  o_raw,       /* omega centre of mass */
  m_ss,   /* moments */
  m_ff,
  m_oo,
  m_sf,
  m_so,
  m_fo,
  f_cen,  /* Filled in elsewhere  - spatial correction */
  s_cen,  /* ditto */
  dety, /*Filled in elsewhere  - flip to HFPO book */
  detz, /*Filled in elsewhere  - flip to HFPO book */
  NPROPERTY ,   /* Number of properties if starting at 0 */ """
        nl = names.split("\n")
        namelist = [n.split(",")[0].lstrip().rstrip()
                    for n in nl]
        i = 0
        while i < len(namelist):
            self.assertEqual(i,getattr(cImageD11,namelist[i]))
            i += 1

    def test_find_max(self):
        for t in [np.uint8, np.int8, np.uint16, np.int16,
                  np.int32, np.uint32, np.float32, float]:
            data = np.array( [[ 1, 0, 1],
                             [ 1, 0, 1],
                             [ 1, 8, 1],
                             [ 1, 1, 1]],t)
            bl = np.zeros(data.shape,np.int32)
            npx = cImageD11.connectedpixels(
                data.astype(np.float32),bl,0.1)
            self.assertEqual(npx,1)
            err = np.sum(np.ravel(data-bl))
            self.assertEqual(err, 7) # 8-1
            res = cImageD11.blobproperties(data.astype(np.float32),
                                           bl,
                                           npx,
                                           omega=22.)
            from ImageD11.cImageD11 import s_1, s_I, s_I2, \
                s_fI, s_ffI, s_sI, s_ssI, s_sfI, \
                bb_mn_f, bb_mn_s, bb_mx_f, bb_mx_s,\
                bb_mn_o, bb_mx_o, \
                mx_I, mx_I_f, mx_I_s, mx_I_o
            #            print res,res.shape
            self.assertAlmostEqual(res[0][s_1],10)
            self.assertAlmostEqual(res[0][mx_I],8)
            self.assertAlmostEqual(res[0][mx_I_f],1)
            self.assertAlmostEqual(res[0][mx_I_s],2)
            self.assertAlmostEqual(res[0][mx_I_o],22)




class testbloboverlaps(unittest.TestCase):
    def test1(self):
        import sys
        data1 =  np.array([[ 1, 0, 1, 0, 0, 0, 0],
                          [ 1, 0, 1, 0, 0, 0, 0],
                          [ 1, 0, 1, 1, 0, 0, 0],
                           [ 1, 1, 1, 0, 0, 0, 0]],np.float32)
        bl1 = np.zeros(data1.shape, np.int32)
        np1 = cImageD11.connectedpixels(data1,bl1,0.1)
        data2 =  np.array([[ 0, 0, 0, 0, 2, 0, 2],
                          [ 0, 0, 0, 0, 2, 0, 2],
                          [ 0, 0, 0, 2, 0, 2, 0],
                           [ 0, 0, 0, 0, 2, 0, 2]],np.float32)
        bl2 = np.zeros(data2.shape, np.int32)
        np2 = cImageD11.connectedpixels(data2,bl2,0.1)
        r1 = cImageD11.blobproperties(data1, bl1, np1, omega=-10.0)
        r2 = cImageD11.blobproperties(data2, bl2, np2, omega=10.0)

        cImageD11.bloboverlaps(bl1,np1,r1,
                                     bl2,np2,r2, verbose=0)
        # check r1 is zeroed
        err = np.sum(np.ravel(r1))
        self.assertAlmostEqual(err,0.,6)
        from ImageD11.cImageD11 import s_1, s_I, s_I2, \
            s_fI, s_ffI, s_sI, s_ssI, s_sfI, \
            bb_mn_f, bb_mn_s, bb_mx_f, bb_mx_s,\
            bb_mn_o, bb_mx_o,\
            mx_I, mx_I_f, mx_I_s
        # check total pixels
        self.assertAlmostEqual(r2[0,s_1], 18.0, 6)
        self.assertAlmostEqual(r2[0,s_I], 26.0, 6)
        self.assertAlmostEqual(r2[0,bb_mn_f], 0.0, 6)
        self.assertAlmostEqual(r2[0,bb_mx_f], 6.0, 6)
        self.assertAlmostEqual(r2[0,bb_mn_s], 0.0, 6)
        self.assertAlmostEqual(r2[0,bb_mx_s], 3.0, 6)
        self.assertAlmostEqual(r2[0,bb_mn_o], -10.0, 6)
        self.assertAlmostEqual(r2[0,bb_mx_o],  10.0, 6)

    def test2(self):
        import sys
        data1 =  np.array([[ 1, 0, 0, 0, 0, 1, 1],
                          [ 1, 0, 2, 0, 0, 1, 1],
                          [ 1, 0, 2, 2, 0, 0, 0],
                           [ 0, 0, 2, 0, 0, 0, 0]],np.float32)
        bl1 = np.zeros(data1.shape, np.int32)
        np1 = cImageD11.connectedpixels(data1,bl1,0.1)
        data2 =  np.array([[ 0, 0, 0, 0, 2, 0, 0],
                          [ 0, 0, 0, 0, 2, 2, 0],
                          [ 0, 0, 0, 2, 0, 0, 0],
                           [ 0, 0, 0, 0, 0, 0, 0]],np.float32)
        bl2 = np.zeros(data2.shape,np.int32)
        np2 = cImageD11.connectedpixels(data2,bl2,0.1)
        r1 = cImageD11.blobproperties(data1, bl1, np1)
        r2 = cImageD11.blobproperties(data2, bl2, np2)

        cImageD11.bloboverlaps(bl1,np1,r1,
                                     bl2,np2,r2, verbose=0)
        # check which peaks are zeroed
        self.assertAlmostEqual(np.sum(r1[1]),0.,6)
        self.assertAlmostEqual(np.sum(r1[2]),0.,6)

        from ImageD11.cImageD11 import s_1, s_I, s_I2, \
            s_fI, s_ffI, s_sI, s_ssI, s_sfI, \
            bb_mn_f, bb_mn_s, bb_mx_f, bb_mx_s,\
            bb_mn_o, bb_mn_f,\
            mx_I, mx_I_f, mx_I_s
        self.assertAlmostEqual(r1[0,s_1], 3, 6)
        self.assertAlmostEqual(r1[0,s_I], 3, 6)

        # check total pixels
        self.assertAlmostEqual(r2[0,s_1], 12.0, 6)
        self.assertAlmostEqual(r2[0,s_I], 20.0, 6)

    def test3(self):
        import sys
        data1 =  np.array([[ 0, 0, 0, 0, 0, 0, 0],
                          [ 0, 0, 0, 0, 0, 0, 0],
                          [ 0, 1, 1, 1, 1, 1, 0],
                           [ 0, 0, 0, 0, 0, 0, 0]],np.float32)
        bl1 = np.zeros(data1.shape,np.int32)
        np1 = cImageD11.connectedpixels(data1,bl1,0.1)
        data2 =  np.array([[ 0, 0, 0, 0, 0, 0, 0],
                          [ 0, 2, 0, 0, 0, 2, 0],
                          [ 0, 2, 0, 0, 0, 2, 0],
                           [ 0, 2, 0, 0, 0, 2, 0]],np.float32)
        bl2 = np.zeros(data2.shape,np.int32)
        np2 = cImageD11.connectedpixels(data2,bl2,0.1)
        r1 = cImageD11.blobproperties(data1, bl1, np1)
        r2 = cImageD11.blobproperties(data2, bl2, np2)

        cImageD11.bloboverlaps(bl1,np1,r1,
                                     bl2,np2,r2, verbose=0)

        # check which peaks are zeroed
        self.assertAlmostEqual(np.sum(r1[0]),0.,6)
        # Results pile up all in r2[0]
        self.assertAlmostEqual(np.sum(r2[1]),0.,6)
        from ImageD11.cImageD11 import s_1, s_I, s_I2, \
            s_fI, s_ffI, s_sI, s_ssI, s_sfI, \
            bb_mn_f, bb_mn_s, bb_mx_f, bb_mx_s,\
            mx_I, mx_I_f, mx_I_s
        self.assertAlmostEqual(r2[0,s_1],11, 6)
        self.assertAlmostEqual(r2[0,s_I],17, 6)





if __name__=="__main__":
    unittest.main()
