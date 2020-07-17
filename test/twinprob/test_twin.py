#

import unittest
from ImageD11 import indexing
import sys, os

class fake():
    def write(self,*a,**k):
        pass
    def flush(self,*a,**k):
        pass

class testtwin( unittest.TestCase ):
    def testtwin(self):
        myindexer = indexing.indexer()
        myindexer.readgvfile(os.path.join(
            os.path.dirname(__file__), 'Al.gve' ))
        myindexer.updateparameters( )
        myindexer.parameterobj.set_parameters(
	{'minpks': 3, 'ring_1': 20, 'ring_2': 20} )
        myindexer.loadpars( )
        myindexer.assigntorings( )
        myindexer.find( )
        myindexer.scorethem( )
        self.assertEqual( len(myindexer.ubis), 1)

if __name__=="__main__":
    unittest.main()
