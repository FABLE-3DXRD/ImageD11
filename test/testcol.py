import unittest
from ImageD11 import columnfile

class testcol( unittest.TestCase):
    def setUp(self):
        open("test.col","w").write("""#  x  a  b  c
0 1 2 3
5 6 1 2
0 4 2 3
32 5 6 2
""")
    def test1(self):
        c  = columnfile.columnfile( "test.col" ) 
        self.assertEqual( list(c.b) , [ 2, 1, 2, 6] )
        c.removerows( "b" , [1] )
        self.assertEqual( c.nrows , 3)
        self.assertEqual( list(c.b)  ,  [ 2, 2, 6] )
    def test2(self):
        c  = columnfile.columnfile( "test.col" ) 
        self.assertEqual( list(c.b)  ,  [ 2, 1, 2, 6])
        c.removerows( "b" , [2] )
        self.assertEqual( c.nrows  ,  2)
        self.assertEqual( list(c.b)  ,  [ 1, 6])
    def test3(self):            
        c  = columnfile.columnfile( "test.col" ) 
        self.assertEqual( list(c.a)  ,  [ 1, 6, 4, 5])
        c.removerows( "a" , [1] )
        c.removerows( "a" , [4] )
        self.assertEqual( c.nrows  ,  2)
        self.assertEqual( list(c.b)  ,  [ 1, 6])
    def test4(self):
        c  = columnfile.columnfile( "test.col" ) 
        self.assertEqual( list(c.a)  ,  [ 1, 6, 4, 5])
        c.removerows( "a" , [1, 4] )
        self.assertEqual( c.nrows  ,  2)
        self.assertEqual( list(c.b)  ,  [ 1, 6])


if __name__=="__main__":
    unittest.main()
