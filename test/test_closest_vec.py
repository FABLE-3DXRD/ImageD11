
from ImageD11 import cImageD11
import numpy as np, time
import unittest


def closest_pyvec( X ):
    nv = X.shape[0]
    nd = X.shape[1]
    res = np.zeros( nv, int )
    for i in range(nv):
        dif = X - X[i]
        err = (dif*dif).sum(axis=1)
        err[i] = 1e9
        res[i] = np.argmin( err )
    return res

class test_closest_vec( unittest.TestCase ):
    def setUp( self ):
        np.random.seed( 42 )
        self.data = np.random.random( (2048,12) )

    def test_same_as_python( self ):
        t1 = time.time()
        r1 = closest_pyvec( self.data )
        t2 = time.time()
        r2 = np.zeros( self.data.shape[0], 'i' )
        cImageD11.closest_vec( self.data, r2 )
        t3 = time.time()
        print( "python",t2-t1,"c",t3-t2)
        self.assertTrue( (r1 == r2).all() )

if __name__=="__main__":
    unittest.main()
            
