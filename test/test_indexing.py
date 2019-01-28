
from __future__ import print_function

from ImageD11.indexing import ubi_fit_2pks
from ImageD11.unitcell import unitcell
import numpy as np
import time
import unittest, cProfile, pstats

def make_random_orientations( N ):
    """
    Generate random quaternions and convert to U matrices
    """
    q = np.random.standard_normal( (4, N) )
    s = 1/(q*q).sum( axis=0 )
    U = np.zeros( (N, 3, 3), float )
    r,i,j,k = 0,1,2,3
    U[:,0,0] = 1 - 2*s*(q[j]*q[j]+q[k]*q[k])
    U[:,0,1] = 2*s*(q[i]*q[j]-q[k]*q[r])
    U[:,0,2] = 2*s*(q[i]*q[k]+q[j]*q[r])
    U[:,1,0] = 2*s*(q[i]*q[j]+q[k]*q[r])    
    U[:,1,1] = 1 - 2*s*(q[i]*q[i]+q[k]*q[k])
    U[:,1,2] = 2*s*(q[j]*q[k]-q[i]*q[r])
    U[:,2,0] = 2*s*(q[i]*q[k]-q[j]*q[r])
    U[:,2,1] = 2*s*(q[j]*q[k]+q[i]*q[r])
    U[:,2,2] = 1 - 2*s*(q[i]*q[i]+q[j]*q[j])
    return U

class test_2pks( unittest.TestCase ):
    def setUp(self):
        self.a = 6
        self.c = 5
        self.U = make_random_orientations( 1 )[0]
        
    def test_tetragonal(self):
        a,c = self.a, self.c
        cell = unitcell( [ a*1.05, a*1.05, c*0.95, 90., 90.,90.], "P" )
        ub = np.dot( self.U, cell.B )
        h1 = (1,0,0)
        h2 = (0,0,1)
        g1 = np.dot( ub, h1 ) # a calc, wrong 
        g2 = np.dot( ub, h2 ) # c
        ideal = unitcell( [ a, a, c, 90., 90.,90.], "P" )
        ideal.makerings(0.3)
        ideal.orient( 0, g1, 1, g2 )
        gve = np.vstack( (g1 , g2) ).T
        for ubi in ideal.UBIlist:
            ufit = ubi_fit_2pks( ubi, g1, g2 )
            hold = np.dot( ubi, gve )
            drlvold = abs(hold - np.round(hold)).sum()
            hnew = np.dot( ufit, gve )
            drlvnew = abs(hnew - np.round(hnew)).sum()
            self.assertTrue( drlvold > 0 )
            self.assertAlmostEqual( drlvnew, 0 )

    def test_hexagonal(self):
        a,c = self.a, self.c
        cell = unitcell( [ a*1.05, a*1.05, c*0.95, 90., 90.,120.], "P" )
        ub = np.dot( self.U, cell.B )
        h1 = (1,0,0)
        h2 = (0,0,1)
        g1 = np.dot( ub, h1 ) # a calc, wrong 
        g2 = np.dot( ub, h2 ) # c
        ideal = unitcell( [ a, a, c, 90., 90.,120.], "P" )
        ideal.makerings(0.3)
        ideal.orient( 0, g1, 1, g2 )
        gve = np.vstack( (g1 , g2) ).T
        for ubi in ideal.UBIlist:
            ufit = ubi_fit_2pks( ubi, g1, g2 )
            hold = np.dot( ubi, gve )
            drlvold = abs(hold - np.round(hold)).sum()
            hnew = np.dot( ufit, gve )
            drlvnew = abs(hnew - np.round(hnew)).sum()
            self.assertTrue( drlvold > 0 )
            self.assertAlmostEqual( drlvnew, 0 )

        
        


if __name__=="__main__":
    unittest.main()
            
