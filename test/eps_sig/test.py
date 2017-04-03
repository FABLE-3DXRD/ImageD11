# Create objects to manipulate - they hold your data
#
from ImageD11 import peakmerge, indexing, transformer, eps_sig_solver
mypeakmerger = peakmerge.peakmerger()
mytransformer = transformer.transformer()
myindexer = indexing.indexer()
mysolver  = eps_sig_solver.solver()
#
# Your work starts here:
#
mysolver.loadmap(  u'CuAlBe_scan10.map' )
mysolver.loadpars(  u'mypar.par' )
mysolver.updateparameters( )
mysolver.compute_write_eps_sig(  u'test.out' )

