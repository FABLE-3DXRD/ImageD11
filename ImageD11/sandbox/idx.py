# Create objects to manipulate - they hold your data
#
from ImageD11 import peakmerge, indexing, transformer
mypeakmerger = peakmerge.peakmerger()
mytransformer = transformer.transformer()
myindexer = indexing.indexer()
#
# Your work starts here:
#
import sys

pars = sys.argv[1]
flt  = sys.argv[2]
gve  = sys.argv[3]
ubi  = sys.argv[4]
mytransformer.loadfileparameters( pars )                               
mytransformer.loadfiltered( flt )
mytransformer.compute_tth_eta( )
mytransformer.addcellpeaks( )
mytransformer.computegv( )
mytransformer.savegv( gve )

myindexer.updateparameters( )
myindexer.parameterobj.set_parameters(  {
    'ds_tol': '0.02', 'minpks': '40', 'uniqueness': '0.5',
    'hkl_tol': '0.05', 'eta_range': '5',
    'ring_1': '2', 'wavelength': '0.26458', 'ring_2': '2', 'cosine_tol': '0.02'} )
myindexer.loadpars( )

myindexer.readgvfile( gve )

myindexer.assigntorings( )

myindexer.find( )
myindexer.scorethem( )
for r1, r2, npks in [
    ( 1,  1, 50 ),
    ( 1,  1, 40 ),
    ( 0,  1, 50 ),
    ( 2,  1, 50 ),
    ( 3,  1, 50 ),
    ( 4,  1, 50 ),
    ( 1 , 1, 50 ),
    ( 0,  1, 40 ),
    ( 2,  1, 40 ),
    ( 3,  1, 40 ),
    ( 4,  1, 40 ),
    ( 1 , 1, 30 ),
    ( 0,  1, 30 ),
    ( 2,  1, 30 ),
    ( 3,  1, 30 ),
    ( 4,  1, 30 ),
    ( 1 , 0, 30 ),
    ( 0,  0, 30 ),
    ( 2,  0, 30 ),
    ( 3,  0, 30 ),
    ( 4,  0, 30 ),
    ]:
    
    myindexer.updateparameters( )
    myindexer.parameterobj.set_parameters(   
    { 'ring_1' : r1,
      'ring_2' : r2,
      'minpks' : npks })
    myindexer.assigntorings( )
    myindexer.loadpars( )
    myindexer.find( )
    myindexer.scorethem( )


myindexer.saveubis( ubi )

