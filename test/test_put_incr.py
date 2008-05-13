

from ImageD11 import closest

import unittest
import numpy as np

class test1(unittest.TestCase):
    def test_put(self):
        data = np.zeros(10,np.float)
        ind  = np.arange(10).astype(np.intc)
        vals = np.ones(10,np.float)
        closest.put_incr( data, ind, vals )
        assert (data == vals).all()
        closest.put_incr( data, ind, vals )
        assert (data == 2*vals).all()

    def test_put_twice(self):
        data = np.zeros(10,np.float)
        ind  = np.ones(10,np.intc)
        vals = np.ones(10,np.float)
        closest.put_incr( data, ind, vals )
        assert (data == np.array( [0, 10] + [0]*8 , np.float)).all()
        closest.put_incr( data, ind, vals )
        assert (data == np.array( [0, 20] + [0]*8 , np.float)).all()

if __name__=="__main__":
    unittest.main()
        
