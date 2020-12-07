

from ImageD11 import cImageD11

import unittest
import numpy as np

class test1(unittest.TestCase):
    def test_put(self):
        data = np.zeros(10,np.float32)
        ind  = np.arange(10).astype(np.intp)
        vals = np.ones(10,np.float32)
        cImageD11.put_incr( data, ind, vals )
        assert (data == vals).all()
        cImageD11.put_incr( data, ind, vals )
        assert (data == 2*vals).all()

    def test_put_twice(self):
        data = np.zeros(10,np.float32)
        ind  = np.ones(10,np.intp)
        vals = np.ones(10,np.float32)
        cImageD11.put_incr( data, ind, vals )
        assert (data == np.array( [0, 10] + [0]*8 , float)).all()
        cImageD11.put_incr( data, ind, vals )
        assert (data == np.array( [0, 20] + [0]*8 , float)).all()

    def test_as_flat(self):
        data = np.zeros( (10, 10), np.float32 )
        ind  = np.ones( 10 , np.intp)*50
        vals = np.ones( 10 , np.float32)
        cImageD11.put_incr( np.ravel(data) , ind, vals )
        assert ( np.ravel(data)[50] == 10 )
        assert ( np.ravel(data)[49] == 0 )
        assert ( np.ravel(data)[51] == 0 )

if __name__=="__main__":
    unittest.main()
        
