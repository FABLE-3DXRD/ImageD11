
import unittest

""" Write some test cases for the columnfile stuff """

from ImageD11 import columnfile

class test1( unittest.TestCase ):
    def setUp( self ):
        f=open("test.col","w")
        f.write("""#  x  a  b  c
0 1 2 3
5 6 1 2
0 4 2 3
32 5 6 2
""")
        f.close()

    def test1(self):
        c  = columnfile.columnfile( "test.col" ) 
        self.assertEqual( list(c.b), [ 2, 1, 2, 6] )
        c.removerows( "b" , [1] )
        self.assertEqual( c.nrows , 3 )
        self.assertEqual( list(c.b) , [ 2, 2, 6] )
        
    def test2(self):
        c  = columnfile.columnfile( "test.col" ) 
        assert list(c.b) == [ 2, 1, 2, 6]
        c.removerows( "b" , [2] )
        assert c.nrows == 2
        assert list(c.b) == [ 1, 6]

    def test3(self):
        c  = columnfile.columnfile( "test.col" ) 
        assert list(c.a) == [ 1, 6, 4, 5]
        c.removerows( "a" , [1] )
        c.removerows( "a" , [4] )
        assert c.nrows == 2
        assert list(c.b) == [ 1, 6]

    def test4(self):
        c  = columnfile.columnfile( "test.col" ) 
        assert list(c.a) == [ 1, 6, 4, 5]
        c.removerows( "a" , [1, 4] )
        assert c.nrows == 2
        assert list(c.b) == [ 1, 6]

    def testsort(self):
        c  = columnfile.columnfile( "test.col" ) 
        c.sortby('a')
        assert( list(c.a.astype(int))  == [ 1, 4, 5, 6] )

    def testreorder(self):
        c  = columnfile.columnfile( "test.col" ) 
        c.reorder([3,2,1,0])
        assert( list(c.a)  == [ 5., 4., 6., 1.] )

    def testmissingfile(self):
        with self.assertRaises( (IOError, FileNotFoundError) ):
            columnfile.columnfile( "a_filename_that_does_not_exist.flt") 
        
if __name__ == '__main__':
    unittest.main()
