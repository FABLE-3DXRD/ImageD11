import sys
from ImageD11 import indexing
myindexer = indexing.indexer()
myindexer.parameterobj.set_parameters(  {
    'max_grains': 10000,
    'minpks' : 50,
} )
myindexer.loadpars( )
myindexer.readgvfile(  sys.argv[1] )
myindexer.score_all_pairs(10)
myindexer.saveubis(  sys.argv[2] )
assert len(myindexer.ubis)==1000
print("Got 1000 grains")
