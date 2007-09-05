


from ImageD11 import gv_general
from ImageD11 import transform

import unittest

import Numeric as n
import RandomArray
import math

def array_diff(a1,a2):
    if a1.shape != a2.shape:
        raise Exception("Shape mismatch")
    diff = n.ravel(a1) - n.ravel(a2)
    err = n.dot(diff, diff)
    #    print "array_diff",err
    return err

class test_g_to_k(unittest.TestCase):
    
    def setUp(self):
        self.np = 5
        self.tth = n.array([60,10,15,20,25],n.Float)
        self.wvln = 0.5
        self.eta = n.array([90,10,120,-20,340],n.Float)
        self.omega = n.array([60,90,180,60,97],n.Float)
        
    def test_0_0(self):
        """ wedge,chi = 0 """
        g = transform.compute_g_vectors(self.tth, 
                                            self.eta, 
                                            self.omega, 
                                            self.wvln, 
                                            wedge=0.0, 
                                            chi=0.0 )
        sol1, sol2, valid = gv_general.g_to_k( g,  #
                                               self.wvln,
                                               axis = n.array([0,0,-1], n.Float),
                                               pre = None,
                                               post = None )

        c0 = n.cos(transform.radians(self.omega))
        c1 = n.cos(transform.radians(sol1))
        c2 = n.cos(transform.radians(sol2))

        s0 = n.sin(transform.radians(self.omega))
        s1 = n.sin(transform.radians(sol1))
        s2 = n.sin(transform.radians(sol2))

        err1 = n.absolute(c1-c0) + n.absolute(s1-s0)
        err2 = n.absolute(c2-c0) + n.absolute(s2-s0)
        err = n.minimum(err1,err2)
        self.assertAlmostEqual( array_diff( err, n.zeros(self.np)), 0, 6)

# if False:
class dont_test_g_to_k_many(test_g_to_k):
    
    def setUp(self):
        self.np = 5
        self.tth = RandomArray.random(self.np)*2+1
        self.wvln = 0.154
        self.eta = RandomArray.random(self.np)*360.
        self.omega = RandomArray.random(self.np)*360.
        


class test_k_to_g(unittest.TestCase):
    
    def setUp(self):
        self.np = 6
        self.tth = RandomArray.random(self.np)*20
        self.wvln = 0.154
        self.eta = (RandomArray.random(self.np)-0.5)*180.
        self.omega = n.array([90,180,12,97,230,199],n.Float)
        # print "Called setup"
        
    def test_0_0_0(self):
        """ wedge, chi = 0 """
        SANITY = False
        om = n.zeros(self.np,n.Float)
        g_old = transform.compute_g_vectors(self.tth, 
                                            self.eta, 
                                            om, 
                                            self.wvln, 
                                            wedge=0.0, 
                                            chi=0.0 )
        if SANITY: print g_old.shape, g_old[:,0]
        k = transform.compute_k_vectors(self.tth, 
                                        self.eta, 
                                        self.wvln)
        if SANITY: print "k=",k[:,0]
        g_new = gv_general.k_to_g(k, om)
        if SANITY: print g_new.shape,g_new[:,0]
        if SANITY: print "end routine"
        self.assertAlmostEqual( array_diff( g_new, g_old ), 0, 6)
        
    def test_0_0(self):
        """ wedge, chi = 0 """
        SANITY = False
        if SANITY: print "*"*80
        om = n.zeros(self.np,n.Float)
        g_old = transform.compute_g_vectors(self.tth, 
                                            self.eta, 
                                            self.omega, 
                                            self.wvln, 
                                            wedge=0.0, 
                                            chi=0.0 )
        if SANITY: print g_old.shape, g_old[:,:3]
        k = transform.compute_k_vectors(self.tth, 
                                        self.eta, 
                                        self.wvln)
        if SANITY: print "k=",k[:,:3]
        g_new = gv_general.k_to_g(k, self.omega)
        if SANITY: print g_new.shape,g_new[:,:3]
        if SANITY: print "end routine"
        if SANITY: print "*"*80
        self.assertAlmostEqual( array_diff( g_new, g_old ), 0, 6)
        
class test_rotation_axis(unittest.TestCase):
    def test_Rz0(self):
        """ zero rotation around z give identity"""
        o = gv_general.rotation_axis( [ 0, 0, 1] , 0 )
        self.assertAlmostEqual( array_diff( o.to_matrix() , 
                                            n.identity(3, n.Float) ) ,
                                            0.0, 6)

    def test_Rz30(self):
        """ right handed, 30 degrees around z """
        o = gv_general.rotation_axis( [ 0, 0, 1], 30 )
        c30 = math.sqrt(3)/2. # cos(30)
        s30 = 0.5 # sin(30)
        r = n.transpose(n.array( [[ c30, s30, 0 ],
                                  [-s30, c30, 0 ],
                                  [ 0  , 0  , 1 ]]) )
        e = array_diff( o.to_matrix() , r ) 
        self.assertAlmostEqual(  e , 0.0 , 6)
        

    def test_from_matrix(self):
        """test_from_matrix: Check matrix construction"""
        c30 = math.sqrt(3)/2.
        s30 = math.sin(math.radians(30.))
        m = n.transpose( n.array( [[ c30, s30, 0 ],
                                   [-s30, c30, 0 ],
                                   [ 0  , 0  , 1 ]] )    )
        o = gv_general.axis_from_matrix(m)
        self.assertAlmostEqual( array_diff (o.to_matrix(), m) , 0, 6 )
        self.assertAlmostEqual( array_diff (o.direction , n.array([0,0,1])), 0, 6 )
        self.assertNotAlmostEqual( array_diff (o.direction , n.array([0,1,1])), 0, 6 )
        self.assertAlmostEqual( o.angle, 30.0 , 6 )
        self.assertNotAlmostEqual( o.angle, 31.0 , 6 )


    def test_rotate_vectors(self):
        """rotate 3x3 vectors, angs or not"""
        o = gv_general.rotation_axis( [0,0,1] , 90.0 )
        vecs = n.transpose(n.array([[ 0, 0, 1 ], 
                                    [ 0, 1, 0 ],
                                    [ 1, 0, 0 ],
                                    [ 1, 0, 0 ]], n.Float))
        res =  n.transpose(n.array([[ 0, 0, 1 ], 
                                    [-1, 0, 0 ],
                                    [ 0, 1, 0 ],
                                    [ 0, 1, 0 ]], n.Float))
        m1, m2 = o.rotate_vectors(vecs), res
        self.assertEqual(m1.shape, m2.shape)
        err = "\n\n"+str(m1.astype(n.Int))+" != \n\n"+str(m2)
        self.assertAlmostEqual(array_diff(m1, m2),0,6,err)
        angs = n.array([90.]*4,n.Float)
        m3, m4 = o.rotate_vectors(vecs,angs), res
        err = "\n\n"+str(m3.astype(n.Int))+" != \n\n"+str(m4)
        self.assertAlmostEqual(array_diff(m3, m4),0,6,err)
        
        
    def test_rotate_vectors_a(self):
        """rotate again 3x3, with different angles"""
        o = gv_general.rotation_axis( [0,0,1] , 90 )
        vecs = n.transpose(n.array([[ 0, 0, 1 ], 
                                    [ 0, 1, 0 ],
                                    [ 1, 0, 0 ]]))
        res =  n.transpose(n.array([[ 0, 0, 1 ], 
                                    [-1, 0, 0 ],
                                    [ 0, 1, 0 ]]))
        self.assertAlmostEqual(array_diff(
             o.rotate_vectors(vecs), res ),0,6)
        r2 = o.rotate_vectors(vecs, [90,90,90] )
        #print "\n",n.array_str(r2, suppress_small=1)
        #print res
        self.assertAlmostEqual(array_diff(
            r2 , res ),0,6)
        self.assertNotAlmostEqual(array_diff(
             o.rotate_vectors(vecs, [90,91,90] ), res ),0,6)
        
    def test_rotate_vectors1(self):
        """ rotate a single vector """
        o = gv_general.rotation_axis( [0,0,1] , 180 )
        vecs = n.transpose(n.array([[ 1, 1, 0 ]]))
        res =  n.transpose(n.array([[-1, -1, 0 ]]))
        self.assertAlmostEqual(array_diff(
          o.rotate_vectors(vecs), res ),0,6)
        self.assertAlmostEqual(array_diff(
           o.rotate_vectors(vecs, [180] ), res ),0,6)
        r1 = o.rotate_vectors(vecs, [181] )
        self.assertNotAlmostEqual(array_diff(r1 , res ),0,6)

        
        
    def test_rotate_vectors2(self):
        """ rotate more vectors """
        o = gv_general.rotation_axis( [0,0,1] , 180 )        
        vecs = n.transpose(n.array([[ 0, 0, 1 ], #1
                                    [ 0, 1, 0 ], #2
                                    [ 0, 1, 1 ], #3
                                    [ 1, 1, 0 ], #4
                                    [-1, 1, 0 ], #5
                                    [ 1, 0, 0 ], #6
                                    [ 0, 1, 1 ]],n.Float))
        res =  n.transpose(n.array([[ 0, 0, 1 ],  #1
                                    [ 0, -1, 0 ], #2
                                    [ 0, -1, 1 ], #3
                                    [-1, -1, 0 ], #4
                                    [ 1, -1, 0 ], #5
                                    [-1, 0 , 0 ], #6
                                    [0, -1 , 1]],n.Float))
        r1 = o.rotate_vectors(vecs)
        #print n.array_str(r1,precision=1,suppress_small=1)
        #print res
        #print n.array_str(r1-res,precision=1,suppress_small=0)
        self.assertAlmostEqual(array_diff( r1 , res) ,0 ,6 )
        self.assertAlmostEqual(array_diff( o.rotate_vectors(vecs, 
                               #1  2    3   4   5   6   7
                                 [180,180,180,180,180,180,180] ), res ), 0,5)


if __name__ ==  "__main__":
    unittest.main()

