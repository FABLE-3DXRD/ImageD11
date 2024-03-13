
import unittest
import numpy as np
from ImageD11 import cImageD11

# print(cImageD11.cImageD11_avx.__file__)

def py_coo( msk ):
    i, j = np.mgrid[0:msk.shape[0], 0:msk.shape[1]]
    return i[msk].astype(np.uint16), j[msk].astype(np.uint16)


class test_clean_mask( unittest.TestCase ):
    def setUp( self ):
        shape = (9,11)
        data = np.zeros( shape, bool )
        data[3:5, 3:5] = 1 # keep
        self.target = data.copy()
        data[1,1] = 1 # isolated, to be removed
        data[3,0] = 1
        data[-1,-1] = 1
        data[-1,5] = 1
        data[5,-1] = 1
        data[6,0]  = 1 # adjacent to prev in mem but on edge
        self.src = data.copy()
    def test1(self):
        testmask = np.zeros( self.src.shape, np.int8 )
        npx = cImageD11.clean_mask( self.src.astype(np.int8), testmask )
        if 0:
            print('npx=',npx)
            print( self.src)
            print(self.src.view(dtype=np.int8))
            print(testmask)
            print(testmask.sum())
        self.assertTrue( self.target.sum() == npx )
        self.assertTrue( (testmask == self.target).all() )
    def test_coo(self):
        npx = self.target.sum()
        i = np.empty( npx, np.uint16)
        j = np.empty( npx, np.uint16)
        tmp = np.empty( self.target.shape[0], np.int32 )
        ok = cImageD11.mask_to_coo( self.target.astype(np.int8), i, j, tmp )
        self.assertTrue( ok == 0 )
        pi, pj = py_coo( self.target )
        self.assertTrue( ( pi == i ).all() )
        self.assertTrue( ( pj == j ).all() )
    def test_coo2(self):
        npx = self.src.sum()
        i = np.empty( npx, np.uint16)
        j = np.empty( npx, np.uint16)
        tmp = np.empty( self.target.shape[0], np.int32 )
        ok = cImageD11.mask_to_coo( self.src.astype(np.int8), i, j, tmp )
        self.assertTrue( ok == 0 )
        pi, pj = py_coo( self.src )
        self.assertTrue( ( pi == i ).all() )
        self.assertTrue( ( pj == j ).all() )
        

if __name__=="__main__":
    unittest.main()

        
