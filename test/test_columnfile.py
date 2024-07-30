
import unittest

""" Write some test cases for the columnfile stuff """

from ImageD11 import columnfile
import ImageD11.columnfile
import numpy as np

class testgeom( unittest.TestCase ):
    def setUp( self ):
        with open("testgeom.flt","w") as f:
            f.write("""
# chi = 0.0
# distance = 5881.14814327
# o11 = 1
# o12 = 0
# o21 = 0
# o22 = 1
# omegasign = 1.0
# t_x = 31.8613635188
# t_y = -28.5519205159
# t_z = 0.0
# tilt_x = -0.00319993715326
# tilt_y = -0.00227706383676
# tilt_z = 0.0172853900689
# wavelength = 0.2469981
# wedge = 0.0761143411257
# y_center = 1080.16471254
# y_size = 1.5
# z_center = 978.22872971
# z_size = 1.5
#  sc  fc  omega  Number_of_pixels
20  30  0  1
21  31  30  1
22  32  60  1
23  33  70  1
24  34  80  1
""")
    def test1( self ):
        c  = columnfile.columnfile("testgeom.flt")
        c.updateGeometry( )

    def testfilter(self):
        c = columnfile.columnfile("testgeom.flt")
        d = c.copy()
        d.filter( abs(d.sc-22) > .1 )
        self.assertTrue( d.nrows == 4 )
        d.filter( abs(d.sc-22) > .1 )
        self.assertTrue( d.nrows == 4 )
        self.assertTrue( c.nrows == 5 )
        d = c.copy()
        self.assertTrue( d.nrows == 5 )

    def testfilter2(self):
        c = columnfile.columnfile("testgeom.flt")
        d = c.copy()
        a = d.bigarray
        d.filter( abs(d.sc-22) > .1 )
        self.assertTrue( d.nrows == 4 )
        d.filter( abs(d.sc-22) > .1 )
        self.assertTrue( d.nrows == 4 )
        self.assertTrue( c.nrows == 5 )
        d = c.copy()
        self.assertTrue( d.nrows == 5 )

    def testgeom(self):
        cold = columnfile.columnfile("testgeom.flt")
        cnew = columnfile.columnfile("testgeom.flt")
        cols = ('xl','yl','zl',"tth", "eta", "ds", "gx", "gy", "gz")
        for t in ((0,0,0), (100,200,300), None):
            cold.updateGeometry( pars=None, translation=t, fast=True )
            cnew.updateGeometry( pars=None, translation=t, fast=False )
            for col in cols:
                e =  np.allclose( cold[col], cnew[col] ) 
                if not e:
                    print(col, cold[col], cnew[col], t)
                self.assertTrue(e)
                                


class test1( unittest.TestCase ):
    def setUp( self ):
        f=open("test.col","w")
        f.write("""#  x  a  b  c
0 1 2 3
5 6 1 2
0 4 2 3
32 5 6 2
33 7 9 0
""")
        f.close()
        
    def test_to_hdf( self ): 
        c = columnfile.columnfile( "test.col" )
        columnfile.colfile_to_hdf( c, "testcolfile.hdf" )
        h = columnfile.columnfile( "testcolfile.hdf" )
        for t in c.titles:
            assert( (c.getcolumn(t) == h.getcolumn(t)).all())
            assert t in h.titles

        
    def testaddcol1( self ):
        c = columnfile.columnfile( "test.col")
        c.addcolumn( [5,4,1,2,0], 'alice' )
        self.assertEqual( c.titles[-1] , 'alice')
        self.assertEqual( list(c.alice) , [5,4,1,2,0])
        self.assertEqual( c.ncols , 5)
        self.assertEqual( c.nrows , 5)

    def testaddcol2( self ):
        c = columnfile.columnfile( "test.col")
        c.addcolumn( [5,4,1,2,9], 'a' )
        self.assertEqual( c.titles.index("a") , 1 )
        self.assertEqual( list(c.a) , [5,4,1,2,9])
        self.assertEqual( c.ncols , 4)
        self.assertEqual( c.nrows , 5)

        
    def test1(self):
        c  = columnfile.columnfile( "test.col" ) 
        self.assertEqual( list(c.b), [ 2, 1, 2, 6, 9] )
        c.removerows( "b" , [1] )
        self.assertEqual( c.nrows , 4 )
        self.assertEqual( list(c.b) , [ 2, 2, 6, 9] )
        
    def test2(self):
        c  = columnfile.columnfile( "test.col" ) 
        assert list(c.b) == [ 2, 1, 2, 6, 9]
        c.removerows( "b" , [2] )
        assert c.nrows == 3
        assert list(c.b) == [ 1, 6, 9]

    def test3(self):
        c  = columnfile.columnfile( "test.col" ) 
        assert list(c.a) == [ 1, 6, 4, 5, 7]
        assert list(c.b) == [ 2, 1, 2, 6, 9]
        c.removerows( "a" , [1] )
        c.removerows( "a" , [4] )
        assert c.nrows == 3
        assert list(c.b) == [ 1, 6, 9]

    def test4(self):
        c  = columnfile.columnfile( "test.col" ) 
        assert list(c.a) == [ 1, 6, 4, 5, 7]
        c.removerows( "a" , [1, 4] )
        assert c.nrows == 3
        assert list(c.b) == [ 1, 6, 9]

    def testsort(self):
        c  = columnfile.columnfile( "test.col" ) 
        c.sortby('a')
        assert( list(c.a.astype(int))  == [ 1, 4, 5, 6, 7] )

    def testreorder(self):
        c  = columnfile.columnfile( "test.col" ) 
        c.reorder([4,3,2,1,0])
        assert( list(c.a)  == [ 7, 5, 4, 6, 1] )

    def testmissingfile(self):
        e = False
        try:
            columnfile.columnfile( "a_filename_that_does_not_exist.flt")        
        except:
            # got an exception
            e = True
        assert( e )


class test_issue289( unittest.TestCase ):  
    """
    https://github.com/FABLE-3DXRD/ImageD11/issues/289
    @AxelHenningsson
    """
    def setUp(self):
        self.colf = ImageD11.columnfile.colfile_from_dict( {'data': np.random.rand(4,)} )
        
    def testMutationInPlace(self):
        # mutation by in-place-operation - OK!
        self.colf.data *= 10
        self.assertTrue( np.allclose(self.colf['data'], self.colf.data) )
        
    def testMutationByRef(self):
        # mutation by reference to a scalar
        self.colf.data = 10
        self.assertTrue( np.allclose(self.colf['data'], self.colf.data) )
        self.assertEqual( self.colf['data'][0], 10) 
        
    def testAddAttr(self):
        self.colf.newthing = 12
        self.assertEqual( self.colf.newthing, 12 )
        
        
        
        
        
if __name__ == '__main__':
    unittest.main()
