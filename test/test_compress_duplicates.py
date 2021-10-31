
import unittest
import numpy as np
import scipy.sparse
from ImageD11 import cImageD11

def _sum_duplicates(row, col, data):
    """ for testing - to be moved or removed """
    # From scipy.sparse.coo_matrix
    # Assumes (data, row, col) not in canonical format.
    if len(data) == 0:
        return row, col, data
    order = np.lexsort((col,row))
    row = row[order]
    col = col[order]
    data = data[order]
    unique_mask = ((row[1:] != row[:-1]) |
                   (col[1:] != col[:-1]))
    unique_mask = np.append(True, unique_mask)
    row = row[unique_mask]
    col = col[unique_mask]
    unique_inds, = np.nonzero(unique_mask)
    data = np.add.reduceat(data, unique_inds, dtype=data.dtype)
    return row, col, data            


class testpls(unittest.TestCase):
    N=1000
    def test1(self):
        np.random.seed(42)
        N = self.N
        a = np.random.random(N*2)
        i, j = (a*10).astype('i').reshape(2,N)
        oj = np.empty( N, 'i')
        oi = np.empty( N, 'i')
        tmp = np.empty( 20, 'i' )
        d = np.ones( N, int )
        r,c,d = _sum_duplicates( i, j, d )
        n = cImageD11.compress_duplicates( i, j, oi, oj, tmp )
        self.assertTrue( n == len(r) )
        self.assertTrue((r==i[:n]).all())
        self.assertTrue((c==j[:n]).all())
        self.assertTrue((d==oi[:n]).all())
    def test2(self):
        np.random.seed(42)
        N = self.N
        a = np.random.random(N*2)
        i, j = (a*10).astype('i').reshape(2,N)
        oj = np.empty( N, 'i')
        oi = np.empty( N, 'i')
        tmp = np.empty( 20, 'i' )
        d = np.ones( N, int )
        coo = scipy.sparse.coo_matrix((d,(i.copy(),j.copy())))
        csr = coo.tocsr()
        coo = csr.tocoo()
        #print("",i,"\n",j)
        n = cImageD11.compress_duplicates( i, j, oi, oj, tmp )
        self.assertEqual(cImageD11.sparse_is_sorted( i[:n].astype(np.uint16),
                                                     j[:n].astype(np.uint16)
                                                  ), 0)
        r = coo.row
        c = coo.col
        d = coo.data
        self.assertTrue( n == len(r) )
        self.assertTrue((r==i[:n]).all())
        self.assertTrue((c==j[:n]).all())
        self.assertTrue((d==oi[:n]).all())

    
if __name__=="__main__":
    unittest.main()
            
