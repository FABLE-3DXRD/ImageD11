

# Create objects to manipulate - they hold your data
#
from ImageD11 import indexing
myindexer = indexing.indexer()

#
# Your work starts here:
#
myindexer.readgvfile(  'Al1000/Al1000.gve' )
myindexer.updateparameters( )
myindexer.parameterobj.set_parameters(  {
    'minpks': 100,
    'max_grains': 1001,
    'ring_1': 5,
    'ring_2': 5,
    'cosine_tol': -0.001,
    } )
myindexer.loadpars( )
myindexer.assigntorings( )
myindexer.find( )
myindexer.scorethem( )
myindexer.saveubis(  'allgrid.map' )
