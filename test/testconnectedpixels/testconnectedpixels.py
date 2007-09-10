
from ImageD11 import connectedpixels


import unittest

# keep use of Numeric in mind
import numpy.oldnumeric as n





class test_connectedpixels(unittest.TestCase):
    def setUp(self):
        dims = (10,10)
        dark = n.ones(dims)
        flood = n.ones(dims)

    def test_1(self):
        for t in [n.UInt8, n.Int8, n.UInt16, n.Int16,
                  n.Int32, n.UInt32, n.Float32, n.Float]:
            data = n.array( [[ 0, 0, 0, 0, 0, 0, 0],
                             [ 0, 0, 0, 1, 0, 0, 0],
                             [ 0, 0, 0, 0, 0, 0, 0],
                             [ 0, 0, 0, 0, 0, 0, 0]],t)
            bl = n.zeros(data.shape,n.Int32)
            connectedpixels.connectedpixels(
                data,
                bl,
                0.1)  # threshold 
            err = n.sum(n.ravel(data-bl))
            self.assertEqual(err, 0)

    def test_1_shape(self):
        for t in [n.UInt8, n.Int8, n.UInt16, n.Int16,
                  n.Int32, n.UInt32, n.Float32, n.Float]:
            data = n.array( [[ 1, 0, 1, 0, 1, 0, 1],
                             [ 1, 0, 1, 0, 1, 0, 1],
                             [ 1, 0, 1, 0, 0, 1, 0],
                             [ 1, 1, 1, 1, 1, 0, 1]],t)
            bl = n.zeros(data.shape,n.Int32)
            connectedpixels.connectedpixels(
                data,
                bl,
                0.1)  # threshold 
            err = n.sum(n.ravel(data-bl))
            self.assertEqual(err, 0)

    def test_2_shapes(self):
        for t in [n.UInt8, n.Int8, n.UInt16, n.Int16,
                  n.Int32, n.UInt32, n.Float32, n.Float]:
            data = n.array( [[ 1, 0, 1, 0, 2, 0, 2],
                             [ 1, 0, 1, 0, 2, 0, 2],
                             [ 1, 0, 1, 0, 0, 2, 0],
                             [ 1, 1, 1, 0, 2, 0, 2]],t)
            bl = n.zeros(data.shape,n.Int32)
            connectedpixels.connectedpixels(
                data,
                bl,
                0.1)  # threshold 
            err = n.sum(n.ravel(data-bl))
            self.assertEqual(err, 0)

    def test_2_transpose(self):
        for t in [n.UInt8, n.Int8, n.UInt16, n.Int16,
                  n.Int32, n.UInt32, n.Float32, n.Float]:
            data = n.array( [[ 1, 0, 1, 0, 2, 0, 2],
                             [ 1, 0, 1, 0, 2, 0, 2],
                             [ 1, 0, 1, 0, 0, 2, 0],
                             [ 1, 1, 1, 0, 2, 0, 2]],t)
            bl = n.zeros(data.shape,n.Int32)
            self.assertRaises(ValueError,
                              connectedpixels.connectedpixels,
                              *(n.transpose(data),bl,0.1))
            connectedpixels.connectedpixels(
                n.transpose(data),n.transpose(bl),0.1)
            err = n.sum(n.ravel(data-bl))
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
  bb_mn_o,      /* 21 min of o */  
  NPROPERTY,     /* Number of properties if starting at 0 */ """
        namelist = [n.split(",")[0].lstrip().rstrip() 
                    for n in names.split("\n")]
        i = 0
        while i < len(namelist):
            self.assertEqual(i,getattr(connectedpixels,namelist[i]))
            i += 1 

    def test_find_max(self):
        for t in [n.UInt8, n.Int8, n.UInt16, n.Int16,
                  n.Int32, n.UInt32, n.Float32, n.Float]:
            data = n.array( [[ 1, 0, 1],
                             [ 1, 0, 1],
                             [ 1, 8, 1],
                             [ 1, 1, 1]],t)
            bl = n.zeros(data.shape,n.Int32)
            self.assertRaises(ValueError,
                              connectedpixels.connectedpixels,
                              *(n.transpose(data),bl,0.1))
            np = connectedpixels.connectedpixels(
                n.transpose(data),n.transpose(bl),0.1)
            self.assertEqual(np,1)
            err = n.sum(n.ravel(data-bl))
            self.assertEqual(err, 7) # 8-1
            res = connectedpixels.blobproperties(data, bl, np)
            from ImageD11.connectedpixels import s_1, s_I, s_I2, \
                s_fI, s_ffI, s_sI, s_ssI, s_sfI, \
                bb_mn_f, bb_mn_s, bb_mx_f, bb_mx_s,\
                mx_I, mx_I_f, mx_I_s 
            #            print res,res.shape
            self.assertAlmostEqual(res[0][s_1],10)
            self.assertAlmostEqual(res[0][mx_I],8)
            self.assertAlmostEqual(res[0][mx_I_f],1)
            self.assertAlmostEqual(res[0][mx_I_s],2)



    def test_2_transpose(self):
        for t in [n.UInt8, n.Int8, n.UInt16, n.Int16,
                  n.Int32, n.UInt32, n.Float32, n.Float]:
            data = n.array( [[ 1, 0, 1, 0, 2, 0, 2],
                             [ 1, 0, 1, 0, 2, 0, 2],
                             [ 1, 0, 1, 0, 0, 2, 0],
                             [ 1, 1, 1, 0, 2, 0, 2]],t)
            bl = n.zeros(data.shape,n.Int32)
            self.assertRaises(ValueError,
                              connectedpixels.connectedpixels,
                              *(n.transpose(data),bl,0.1))
            np = connectedpixels.connectedpixels(
                n.transpose(data),n.transpose(bl),0.1)
            self.assertEqual(np,2)
            err = n.sum(n.ravel(data-bl))
            self.assertEqual(err, 0)
            res = connectedpixels.blobproperties(data, bl, np)
            from ImageD11.connectedpixels import s_1, s_I, s_I2, \
                s_fI, s_ffI, s_sI, s_ssI, s_sfI, \
                bb_mn_f, bb_mn_s, bb_mx_f, bb_mx_s,\
                mx_I, mx_I_f, mx_I_s 
            #            print res,res.shape
            self.assertAlmostEqual(res[0][s_1],9)
            self.assertAlmostEqual(res[1][s_1],7)
            self.assertAlmostEqual(res[0][s_I],9)
            self.assertAlmostEqual(res[1][s_I],14)
            self.assertAlmostEqual(res[0][s_I2],9)
            self.assertAlmostEqual(res[1][s_I2],28)
            # [[ 1, 0, 1, 0, 2, 0, 2],     --> Fast
            #  [ 1, 0, 1, 0, 2, 0, 2],     |
            #  [ 1, 0, 1, 0, 0, 2, 0],     |
            #  [ 1, 1, 1, 0, 2, 0, 2]],t)  V Slow
            # f*I:
            # f= 0, 1, 2, 3, 4, 5, 6
            # [[ 0, 0, 2, 0, 8, 0, 12],     --> Fast
            #  [ 0, 0, 2, 0, 8, 0, 12],     |
            #  [ 0, 0, 2, 0, 0,10, 0 ],     |
            #  [ 0, 1, 2, 0, 8, 0, 12]],t)  V Slow
            #         =9     =70
            self.assertAlmostEqual(res[0][s_fI],9)
            self.assertAlmostEqual(res[1][s_fI],70)
            # s*I:
            # s=
            # 0[[ 0, 0, 0, 0, 0, 0, 0],     --> Fast
            # 1 [ 1, 0, 1, 0, 2, 0, 2],     |
            # 2 [ 2, 0, 2, 0, 0, 4, 0],     |
            # 3 [ 3, 3, 3, 0, 6, 0, 6]],t)  V Slow
            #         =15     =20
            self.assertAlmostEqual(res[0][s_sI],15)
            self.assertAlmostEqual(res[1][s_sI],20)
            # Bounding box
            self.assertAlmostEqual(res[0][bb_mn_f],0)
            self.assertAlmostEqual(res[1][bb_mn_f],4)
            self.assertAlmostEqual(res[0][bb_mx_f],2)
            self.assertAlmostEqual(res[1][bb_mx_f],6)
            self.assertAlmostEqual(res[0][bb_mn_s],0)
            self.assertAlmostEqual(res[1][bb_mn_s],0)
            self.assertAlmostEqual(res[0][bb_mx_s],3)
            self.assertAlmostEqual(res[1][bb_mx_s],3)

            
class testbloboverlaps(unittest.TestCase):
    def test1(self):
        data1 =  n.array([[ 1, 0, 1, 0, 0, 0, 0],
                          [ 1, 0, 1, 0, 0, 0, 0],
                          [ 1, 0, 1, 1, 0, 0, 0],
                          [ 1, 1, 1, 0, 0, 0, 0]])
        bl1 = n.zeros(data1.shape)
        np1 = connectedpixels.connectedpixels(data1,bl1,0.1)
        data2 =  n.array([[ 0, 0, 0, 0, 2, 0, 2],
                          [ 0, 0, 0, 0, 2, 0, 2],
                          [ 0, 0, 0, 2, 0, 2, 0],
                          [ 0, 0, 0, 0, 2, 0, 2]])
        bl2 = n.zeros(data2.shape)
        np2 = connectedpixels.connectedpixels(data2,bl2,0.1)

        r1 = connectedpixels.blobproperties(data1, bl1, np1)

        r2 = connectedpixels.blobproperties(data2, bl2, np2)
        
        print r1
        print r2
        connectedpixels.bloboverlaps(bl1,np1,r1, 
                                     bl2,np2,r2, verbose=0)
        err = n.sum(n.ravel(r2))
        print r1
        print r2
        self.assertAlmostEqual(err,0.,6)
        
        

    
    
class test_roisum(unittest.TestCase):
    
    def test_avg(self):
        data = n.ones((10,10),n.Float)
        res = connectedpixels.roisum(data,0,9,0,9)
        self.assertAlmostEqual(1,res,0,6)

    def test_for_overflow(self):
        data = n.ones((10,10),n.UInt8)+255+128
        # Ensure the array type does indeed overflow
        self.assertEqual(data[0,0],128)
        # And that the sum does not
        res = connectedpixels.roisum(data,0,9,0,9)
        self.assertAlmostEqual(128.,res,0,6)

    def test_orient(self):
        # add a testcase to work out which way up is the roi...
        data = n.array( [[ 0, 0, 0, 1, 1, 1],
                         [ 0, 0, 0, 1, 1, 1],
                         [ 0, 0, 0, 1, 1, 1]], n.Float)
        self.assertAlmostEqual(connectedpixels.roisum(data, 0, 2, 0, 2),
                               0.)
        self.assertAlmostEqual(connectedpixels.roisum(data, 0, 2, 3, 5),
                               1.)
        self.assertRaises(ValueError,
                          connectedpixels.roisum,
                          *(data, 3, 5, 0, 2))
        self.assertRaises(ValueError,
                          connectedpixels.roisum,
                          *(data, 0, 7, 0, 2))
        self.assertRaises(ValueError,
                          connectedpixels.roisum,
                          *(data, 0, 2, 0, 7))
        self.assertRaises(ValueError,
                          connectedpixels.roisum,
                          *(data, 0, 2, -1, 2))
        self.assertRaises(ValueError,
                          connectedpixels.roisum,
                          *(data, -1, 2, 0, 2))


if __name__=="__main__":
    unittest.main()
