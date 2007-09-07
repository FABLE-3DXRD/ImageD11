
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
    def setUp(self):
        self.labels = ["s_1",
                       "s_I",
                       "s_I2",
                       "s_fI",
                       "s_ffI",
                       "s_sI",
                       "s_ssI",
                       "s_sfI",
                       "mx_I",
                       "mx_I_f",
                       "mx_I_s",
                       "bb_mx_f",
                       "bb_mx_s" ,
                       "bb_mn_f" ,
                       "bb_mn_s" ]

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
            # print bl,bl.dtype
            res = connectedpixels.blobproperties(data, bl, np)
            # print res
                
                                           
    

    
    
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
