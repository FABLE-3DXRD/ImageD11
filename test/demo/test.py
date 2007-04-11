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
mytransformer.updateparameters( )
mytransformer.parameterobj.set_parameters(  {'cell_beta': '90', 'distance': '302287.142', 'omegasign': '1', 'cell__a': '4.8192', 'cell__c': '8.5156', 'cell__b': '7.6983', 'fit_tolerance': '0.02', 'chi': '0.0', 'z_center': '502.19212', 'cell_lattice_[P,A,B,C,I,F,R]': 'I', 'cell_alpha': '90', 'cell_gamma': '90', 'y_size': '96.14045', 'tilt_z': '0.001618', 'tilt_y': '-0.002322', 'wedge': '0.0', 'wavelength': '0.153684', 'z_size': '93.54977', 'y_center': '502.618176'} )
mytransformer.compute_tth_eta( )
mytransformer.addcellpeaks( )
mytransformer.fit(  14.0 )
mytransformer.saveparameters(  'eu3.pars' )
mytransformer.computegv( )
mytransformer.savegv(  'eu3.gve' )
myindexer.readgvfile(  'eu3.gve' )
myindexer.updateparameters( )
myindexer.parameterobj.set_parameters(  {'ds_tol': '0.005', 'minpks': '100', 'uniqueness': '0.5', 'hkl_tol': '0.05', 'eta_range': '0.0', 'ring_1': '0', 'wavelength': '0.153684', 'ring_2': '0', 'cosine_tol': '0.002'} )
myindexer.loadpars( )
myindexer.assigntorings( )
myindexer.find( )
myindexer.scorethem( )
myindexer.histogram_drlv_fit( )
myindexer.saveubis(  'eu3.ubi' )
myindexer.saveindexing(  'eu3.idx' )

