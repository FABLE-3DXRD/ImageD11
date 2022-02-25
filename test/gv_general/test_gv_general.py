
from __future__ import print_function

from ImageD11 import gv_general
from ImageD11 import transform

import unittest

import numpy as np
import math

def array_diff(a1,a2):
    if a1.shape != a2.shape:
        raise Exception("Shape mismatch")
    diff = np.ravel(a1) - np.ravel(a2)
    err = np.dot(diff, diff)
    #    print "array_diff",err
    return err

class test_g_to_k(unittest.TestCase):

    def setUp(self):
        self.tth = np.array( [60, 10, 15, 20,25, 1.9555],float)
        self.wvln = 0.5
        self.eta = np.array( [90, 10,120,-20,340, -73 ],float)
        self.omega = np.array([60,90,180, 60,97,  131],float)
        self.np = len(self.tth)

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
                                               axis = np.array([0,0,-1],
                                                               float),
                                               pre = np.eye(3),
                                               post = np.eye(3) )

        c0 = np.cos(transform.radians(self.omega))
        c1 = np.cos(transform.radians(sol1))
        c2 = np.cos(transform.radians(sol2))

        s0 = np.sin(transform.radians(self.omega))
        s1 = np.sin(transform.radians(sol1))
        s2 = np.sin(transform.radians(sol2))

        err1 = np.absolute(c1-c0) + np.absolute(s1-s0)
        err2 = np.absolute(c2-c0) + np.absolute(s2-s0)
        err = np.minimum(err1,err2)
        self.assertAlmostEqual( array_diff( err, np.zeros(self.np)), 0, 6)

    def test_10_0(self):
        """ wedge,chi = 10,0 """
        w=10.0
        c=0
        g = transform.compute_g_vectors(self.tth,
                                            self.eta,
                                            self.omega,
                                            self.wvln,
                                            wedge=w,
                                            chi=c )
        post = gv_general.wedgechi( wedge=w, chi=c)

        sol1, sol2, valid = gv_general.g_to_k( g,  #
                                               self.wvln,
                                               axis = np.array([0,0,-1],
                                                               float),
                                               pre = None,
                                               post = post )

        c0 = np.cos(transform.radians(self.omega))
        c1 = np.cos(transform.radians(sol1))
        c2 = np.cos(transform.radians(sol2))

        s0 = np.sin(transform.radians(self.omega))
        s1 = np.sin(transform.radians(sol1))
        s2 = np.sin(transform.radians(sol2))

        err1 = np.absolute(c1-c0) + np.absolute(s1-s0)
        err2 = np.absolute(c2-c0) + np.absolute(s2-s0)
        err = np.minimum(err1,err2)
        self.assertAlmostEqual( array_diff( err, np.zeros(self.np)), 0, 6)

    def test_0_10(self):
        """ wedge,chi = 0,10 """
        w=0.
        c=10.
        g = transform.compute_g_vectors(self.tth,
                                            self.eta,
                                            self.omega,
                                            self.wvln,
                                            wedge=w,
                                            chi=c )

        post = gv_general.wedgechi( wedge=w, chi=c)

        sol1, sol2, valid = gv_general.g_to_k( g,  #
                                               self.wvln,
                                               axis = np.array([0,0,-1],
                                                               float),
                                               pre = None,
                                               post = post )

        c0 = np.cos(transform.radians(self.omega))
        c1 = np.cos(transform.radians(sol1))
        c2 = np.cos(transform.radians(sol2))

        s0 = np.sin(transform.radians(self.omega))
        s1 = np.sin(transform.radians(sol1))
        s2 = np.sin(transform.radians(sol2))

        err1 = np.absolute(c1-c0) + np.absolute(s1-s0)
        err2 = np.absolute(c2-c0) + np.absolute(s2-s0)
        err = np.minimum(err1,err2)
        self.assertAlmostEqual( array_diff( err, np.zeros(self.np)), 0, 6)
#        print "passed",w,c


    def test_5_10(self):
        """ wedge,chi = 5,10 """
        w,c = 5.0, 10.0
        g = transform.compute_g_vectors(self.tth,
                                        self.eta,
                                        self.omega,
                                        self.wvln,
                                        wedge=w,
                                        chi=c )

        post = gv_general.wedgechi( wedge=w, chi=c )

        sol1, sol2, valid = gv_general.g_to_k( g,  #
                                               self.wvln,
                                               axis = np.array([0,0,-1],
                                                               float),
                                               pre = None,
                                               post = post )

        c0 = np.cos(transform.radians(self.omega))
        c1 = np.cos(transform.radians(sol1))
        c2 = np.cos(transform.radians(sol2))

        s0 = np.sin(transform.radians(self.omega))
        s1 = np.sin(transform.radians(sol1))
        s2 = np.sin(transform.radians(sol2))

        err1 = np.absolute(c1-c0) + np.absolute(s1-s0)
        err2 = np.absolute(c2-c0) + np.absolute(s2-s0)
        err = np.minimum(err1,err2)
        #print "\ns1",sol1
        #print "s2",sol2
        #print "omega",self.omega,"\n"
        self.assertAlmostEqual( array_diff( err, np.zeros(self.np)), 0, 6)




if False:
 class dont_test_g_to_k_many(test_g_to_k):

    def setUp(self):
        self.np = 5
        self.tth = np.random.random(self.np)*2+1
        self.wvln = 0.154
        self.eta = np.random.random(self.np)*360.
        self.omega = np.random.random(self.np)*360.



class test_k_to_g(unittest.TestCase):

    def setUp(self):
        self.np = 6
        self.tth = np.random.random(self.np)*20
        self.wvln = 0.154
        self.eta = (np.random.random(self.np)-0.5)*180.
        self.omega = np.array([90,180,12,97,230,199],float)
        # print "Called setup"

    def test_0_0_0(self):
        """ wedge, chi = 0 """
        SANITY = False
        om = np.zeros(self.np,float)
        g_old = transform.compute_g_vectors(self.tth,
                                            self.eta,
                                            om,
                                            self.wvln,
                                            wedge=0.0,
                                            chi=0.0 )
        if SANITY: print (g_old.shape, g_old[:,0])
        k = transform.compute_k_vectors(self.tth,
                                        self.eta,
                                        self.wvln)
        if SANITY: print("k=",k[:,0])
        g_new = gv_general.k_to_g(k, om)
        if SANITY: print(g_new.shape,g_new[:,0])
        if SANITY: print( "end routine")
        self.assertAlmostEqual( array_diff( g_new, g_old ), 0, 6)

    def test_0_0(self):
        """ wedge, chi = 0 """
        SANITY = False
        if SANITY: print ("*"*80)
        om = np.zeros(self.np,float)
        g_old = transform.compute_g_vectors(self.tth,
                                            self.eta,
                                            self.omega,
                                            self.wvln,
                                            wedge=0.0,
                                            chi=0.0 )
        if SANITY: print (g_old.shape, g_old[:,:3])
        k = transform.compute_k_vectors(self.tth,
                                        self.eta,
                                        self.wvln)
        if SANITY: print( "k=",k[:,:3])
        g_new = gv_general.k_to_g(k, self.omega,
                                  axis=[0,0,-1] )
        if SANITY: print( g_new.shape,g_new[:,:3])
        if SANITY: print( "end routine")
        if SANITY: print( "*"*80)
        self.assertAlmostEqual( array_diff( g_new, g_old ), 0, 6)

    def test_25_30(self):
        """ wedge, chi = 25,30 """
        w,c=25,30
        SANITY = False
        if SANITY: print ("*"*80)
        om = np.zeros(self.np,float)
        g_old = transform.compute_g_vectors(self.tth,
                                            self.eta,
                                            self.omega,
                                            self.wvln,
                                            wedge=w,
                                            chi=c )
        if SANITY: print (g_old.shape,"\n", g_old[:,:3])
        k = transform.compute_k_vectors(self.tth,
                                        self.eta,
                                        self.wvln)
        if SANITY: print( "k=",k[:,:3])
        post = gv_general.chiwedge( chi=c, wedge=w )
        g_new = gv_general.k_to_g(k, self.omega,
                                  axis=[0,0,-1],
                                  post=post )
        if SANITY: print( g_new.shape,g_new[:,:3])
        if SANITY: print( "end routine")
        if SANITY: print( "*"*80)
        self.assertAlmostEqual( array_diff( g_new, g_old ), 0, 6)


class test_rotation_axis(unittest.TestCase):
    def test_Rz0(self):
        """ zero rotation around z give identity"""
        o = gv_general.rotation_axis( [ 0, 0, 1] , 0 )
        self.assertAlmostEqual( array_diff( o.to_matrix() ,
                                            np.identity(3, float) ) ,
                                            0.0, 6)

    def test_Rz30(self):
        """ right handed, 30 degrees around z """
        o = gv_general.rotation_axis( [ 0, 0, 1], 30 )
        c30 = math.sqrt(3)/2. # cos(30)
        s30 = 0.5 # sin(30)
        r = np.transpose(np.array( [[ c30, s30, 0 ],
                                   [-s30, c30, 0 ],
                                   [ 0  , 0  , 1 ]]) )
        e = array_diff( o.to_matrix() , r )
        self.assertAlmostEqual(  e , 0.0 , 6)


    def test_from_matrix(self):
        """test_from_matrix: Check matrix construction"""
        c30 = math.sqrt(3)/2.
        s30 = math.sin(math.radians(30.))
        m = np.transpose( np.array( [[ c30, s30, 0 ],
                                    [-s30, c30, 0 ],
                                    [ 0  , 0  , 1 ]] ) )
        o = gv_general.axis_from_matrix(m)
        self.assertAlmostEqual( array_diff (o.to_matrix(), m) , 0, 6 )
        self.assertAlmostEqual( array_diff (o.direction , np.array([0,0,1])),
                                0, 6 )
        self.assertNotAlmostEqual( array_diff (o.direction , np.array([0,1,1])),
                                   0, 6 )
        self.assertAlmostEqual( o.angle, 30.0 , 6 )
        self.assertNotAlmostEqual( o.angle, 31.0 , 6 )


    def test_rotate_vectors(self):
        """rotate 3x3 vectors, angs or not"""
        o = gv_general.rotation_axis( [0,0,1] , 90.0 )
        vecs = np.transpose(np.array([[ 0, 0, 1 ],
                                     [ 0, 1, 0 ],
                                     [ 1, 0, 0 ],
                                     [ 1, 0, 0 ]], float))
        res =  np.transpose(np.array([[ 0, 0, 1 ],
                                     [-1, 0, 0 ],
                                     [ 0, 1, 0 ],
                                     [ 0, 1, 0 ]], float))
        m1, m2 = o.rotate_vectors(vecs), res
        self.assertEqual(m1.shape, m2.shape)
        err = "\n\n"+str(m1.astype(int))+" != \n\n"+str(m2)
        self.assertAlmostEqual(array_diff(m1, m2),0,6,err)
        angs = np.array([90.]*4,float)
        m3, m4 = o.rotate_vectors(vecs,angs), res
        err = "\n\n"+str(m3.astype(int))+" != \n\n"+str(m4)
        self.assertAlmostEqual(array_diff(m3, m4),0,6,err)


    def test_rotate_vectors_a(self):
        """rotate again 3x3, with different angles"""
        o = gv_general.rotation_axis( [0,0,1] , 90 )
        vecs = np.transpose(np.array([[ 0, 0, 1 ],
                                     [ 0, 1, 0 ],
                                     [ 1, 0, 0 ]]))
        res =  np.transpose(np.array([[ 0, 0, 1 ],
                                     [-1, 0, 0 ],
                                     [ 0, 1, 0 ]]))
        self.assertAlmostEqual(array_diff(
             o.rotate_vectors(vecs), res ),0,6)
        r2 = o.rotate_vectors(vecs, [90,90,90] )
        #print "\n",np.array_str(r2, suppress_small=1)
        #print res
        self.assertAlmostEqual(array_diff(
            r2 , res ),0,6)
        self.assertNotAlmostEqual(array_diff(
             o.rotate_vectors(vecs, [90,91,90] ), res ),0,6)

    def test_rotate_vectors1(self):
        """ rotate a single vector """
        o = gv_general.rotation_axis( [0,0,1] , 180 )
        vecs = np.transpose(np.array([[ 1, 1, 0 ]]))
        res =  np.transpose(np.array([[-1, -1, 0 ]]))
        self.assertAlmostEqual(array_diff(
          o.rotate_vectors(vecs), res ),0,6)
        self.assertAlmostEqual(array_diff(
           o.rotate_vectors(vecs, [180] ), res ),0,6)
        r1 = o.rotate_vectors(vecs, [181] )
        self.assertNotAlmostEqual(array_diff(r1 , res ),0,6)



    def test_rotate_vectors2(self):
        """ rotate more vectors """
        o = gv_general.rotation_axis( [0,0,1] , 180 )
        vecs = np.transpose(np.array([[ 0, 0, 1 ], #1
                                     [ 0, 1, 0 ], #2
                                     [ 0, 1, 1 ], #3
                                     [ 1, 1, 0 ], #4
                                     [-1, 1, 0 ], #5
                                     [ 1, 0, 0 ], #6
                                     [ 0, 1, 1 ]],float))
        res =  np.transpose(np.array([[ 0, 0, 1 ],  #1
                                     [ 0, -1, 0 ], #2
                                     [ 0, -1, 1 ], #3
                                     [-1, -1, 0 ], #4
                                     [ 1, -1, 0 ], #5
                                     [-1, 0 , 0 ], #6
                                     [0, -1 , 1]],float))
        r1 = o.rotate_vectors(vecs)
        #print np.array_str(r1,precision=1,suppress_small=1)
        #print res
        #print np.array_str(r1-res,precision=1,suppress_small=0)
        self.assertAlmostEqual(array_diff( r1 , res) ,0 ,6 )
        self.assertAlmostEqual(array_diff( o.rotate_vectors(vecs,
                               #1  2    3   4   5   6   7
                                 [180,180,180,180,180,180,180] ), res ), 0,5)


if __name__ ==  "__main__":
    unittest.main()

