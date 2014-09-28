#
from ImageD11 import indexing
myindexer = indexing.indexer()
myindexer.readgvfile(  'Al.gve' )
myindexer.updateparameters( )
myindexer.parameterobj.set_parameters(
	{'minpks': 3, 'ring_1': 20, 'ring_2': 20} )
myindexer.loadpars( )
myindexer.assigntorings( )
myindexer.find( )
myindexer.scorethem( )
assert len(myindexer.ubis) == 1, "Should only have one grain"
print "Seems to be OK"

