# Create objects to manipulate - they hold your data
#
from ImageD11 import transformer

mytransformer = transformer.transformer()

#
# Your work starts here:
#
import sys
mytransformer.loadfiltered( sys.argv[1] )
mytransformer.loadfileparameters(  sys.argv[2] )
mytransformer.compute_tth_eta( )
mytransformer.addcellpeaks( )
mytransformer.computegv( )
mytransformer.savegv( sys.argv[3] )
