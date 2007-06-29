

from ImageD11 import gv_general
import unittest

import Numeric as n
import math

def array_diff(a1,a2):
    if a1.shape != a2.shape:
        raise Exception("Shape mismatch")
    diff = n.ravel(a1) - n.ravel(a2)
    return n.dot(diff, diff)

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
        e = array_diff( o.to_matrix() ,n.array( [[ c30, s30, 0 ],
                                                 [-s30, c30, 0 ],
                                                 [ 0  , 0  , 1 ]] ) ) 
        self.assertAlmostEqual(  e , 0.0 , 6)
        

    def test_from_matrix(self):
        """Check matrix construction"""
        c30 = math.sqrt(3)/2.
        s30 = math.sin(math.radians(30.))
        m =  n.array( [[ c30, s30, 0 ],
                       [-s30, c30, 0 ],
                       [ 0  , 0  , 1 ]] )    
        o = gv_general.axis_from_matrix(m)
        self.assertAlmostEqual( array_diff (o.to_matrix(), m) , 0, 6 )
        self.assertAlmostEqual( array_diff (o.direction , n.array([0,0,1])), 0, 6 )
        self.assertNotAlmostEqual( array_diff (o.direction , n.array([0,1,1])), 0, 6 )
        self.assertAlmostEqual( o.angle, 30.0 , 6 )
        self.assertNotAlmostEqual( o.angle, 31.0 , 6 )


    def test_rotate_vectors(self):
        """rotate 3x3 vectors"""
        o = gv_general.rotation_axis( [0,0,1] , 90.0 )
        vecs = [[ 0, 0, 1 ], 
                [ 0, 1, 0 ],
                [ 1, 0, 0 ],
                [ 1, 0, 0 ]]
        res =  n.array([[ 0, 0, 1 ], 
                        [-1, 0, 0 ],
                        [ 0, 1, 0 ],
                        [ 0, 1, 0 ]])
        m1, m2 = o.rotate_vectors(vecs), res 
        err = "\n\n"+str(m1.astype(n.Int))+" != \n\n"+str(m2)
        self.assertAlmostEqual(array_diff(m1,m2),0,6,err)
        
        
    def test_rotate_vectors_a(self):
        """rotate again 3x3, with different angles"""
        o = gv_general.rotation_axis( [0,0,1] , 90 )
        vecs = [[ 0, 0, 1 ], 
                [ 0, 1, 0 ],
                [ 1, 0, 0 ]]
        res =  n.array([[ 0, 0, 1 ], 
                        [-1, 0, 0 ],
                        [ 0, 1, 0 ]])
        self.assertAlmostEqual(array_diff(
             o.rotate_vectors(vecs), res ),0,6)
        self.assertAlmostEqual(array_diff(
             o.rotate_vectors(vecs, [90,90,90] ), res ),0,6)
        self.assertNotAlmostEqual(array_diff(
             o.rotate_vectors(vecs, [90,91,90] ), res ),0,6)
        
    def test_rotate_vectors1(self):
        """ rotate a single vector """
        o = gv_general.rotation_axis( [0,0,1] , 180 )
        vecs = [[ 1, 1, 0 ]]
        res =  n.array([[-1, -1, 0 ]])
        self.assertAlmostEqual(array_diff(
          o.rotate_vectors(vecs), res ),0,6)
        self.assertAlmostEqual(array_diff(
           o.rotate_vectors(vecs, [180] ), res ),0,6)
        r1 = o.rotate_vectors(vecs, [181] )
        self.assertNotAlmostEqual(array_diff(r1 , res ),0,6)

        
        
    def test_rotate_vectors2(self):
        """ rotate 5 vectors """
        o = gv_general.rotation_axis( [0,0,1] , 180 )
        vecs = [[ 0, 0, 1 ], 
                [ 0, 1, 0 ],
                [ 1, 1, 0 ],
                [-1, 1, 0 ],
                [ 1, 0, 0 ]]
        res =  n.array([[ 0, 0, 1 ], 
                        [ 0, -1, 0 ],
                        [-1, -1, 0 ],
                        [ 1, -1, 0 ],
                        [-1, 0 , 0 ]])
        self.assertEqual(o.rotate_vectors(vecs), res )
        self.assertEqual(o.rotate_vectors(vecs, [180,180,180,180,180] ), res )


if __name__ ==  "__main__":
    unittest.main()

