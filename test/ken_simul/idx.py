
import time
start = time.time()
# Create objects to manipulate - they hold your data
#
from ImageD11 import indexing
myindexer = indexing.indexer()
#
# Your work starts here:
#    
myindexer.readgvfile( 'simAl.gve' )
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
myindexer.saveubis(  'test.ubi' )
myindexer.assigntorings( )
print "Time:",time.time()-start
myindexer.saveindexing(  'test.idx' )
print "Time:", time.time()-start
