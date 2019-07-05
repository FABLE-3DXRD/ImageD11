
import logging
logging.basicConfig(level=logging.INFO)
# Create objects to manipulate - they hold your data
#
from ImageD11 import peakmerge, indexing, transformer
mypeakmerger = peakmerge.peakmerger()
mytransformer = transformer.transformer()
myindexer = indexing.indexer()
#
# Your work starts here:
#    
mypeakmerger.readpeaks(  'eu.pks' )
mypeakmerger.harvestpeaks(  numlim=(0.0, 1000.0)  , omlim=(-100.0, 100.0)  )
mypeakmerger.mergepeaks( )
mypeakmerger.filter( )
mypeakmerger.savepeaks(  'eu.flt' )
mytransformer.loadfiltered(  'eu.flt' )
mytransformer.loadfileparameters(  'eu2.pars' )
mytransformer.compute_tth_eta( )
mytransformer.addcellpeaks( )
mytransformer.fit( 0.0 , 14.0 )
mytransformer.saveparameters(  'eu3fitted.pars' )
mytransformer.computegv( )
mytransformer.savegv(  'eu3.gve' )
myindexer.readgvfile(  'eu3.gve' )
myindexer.updateparameters( )
p = {'ds_tol': '0.005', 'minpks': '100', 'uniqueness': '0.5', 'hkl_tol': '0.05', 'eta_range': '0.0', 'ring_1': '1', 'wavelength': '0.153684', 'ring_2': '1', 'cosine_tol': '0.002'} 
myindexer.parameterobj.set_parameters(  p )
myindexer.loadpars( )
myindexer.assigntorings( )
myindexer.find( )
myindexer.scorethem( )
myindexer.histogram_drlv_fit( )
myindexer.saveubis(  'eu3.ubi' )
myindexer.saveindexing(  'eu3.idx' )

