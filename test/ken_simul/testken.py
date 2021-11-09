
# Create objects to manipulate - they hold your data
#
from ImageD11 import indexing
import unittest, sys, os

class fake():
    def write(self,*a,**k):
        pass
    def flush(self,*a,**k):
        pass

class testKen(unittest.TestCase):
    def test12(self):
        ss=sys.stdout
        sys.stdout = fake()
        try:
            myindexer = indexing.indexer()
            myindexer.readgvfile( os.path.join( os.path.split(__file__)[0] ,'simAl.gve' ) )

            myindexer.updateparameters( )
            myindexer.parameterobj.set_parameters(
                {'minpks': '80',
                 'ring_1': '0',
                 'ring_2': '1',
                } )
            myindexer.loadpars( )
            myindexer.assigntorings( )
            myindexer.find( )
            myindexer.scorethem( )
            # myindexer.saveubis(  'test.ubi' )
            myindexer.assigntorings( )
            # myindexer.saveindexing(  'test.idx' )
            self.assertEqual( len( myindexer.ubis ) , 20 )
        finally:
            sys.stdout=ss

if __name__=="__main__":
    unittest.main()
