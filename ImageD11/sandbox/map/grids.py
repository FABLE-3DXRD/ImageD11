
from __future__ import print_function
# Create objects to manipulate - they hold your data
#
import numpy as np
from ImageD11 import peakmerge, indexing, transformer
mypeakmerger = peakmerge.peakmerger()
mytransformer = transformer.transformer()
myindexer = indexing.indexer()
#
# Your work starts here:
#



# Peaks file, tolerance, number of peaks


#flts = [ ("90-269_t50.flt",     0.05, 50 ),
#         ("90-269_t50_gt2.flt", 0.05, 50 ),
#         ("com_t50_gt2.flt",    0.10, 50),
#         ] 


def grids( flt, tol, npk ):
    import os
    mapfile = open(flt.replace("flt","maptxt"),"a")
    mapfile.write("# %s %f %d\n"%( flt, tol, npk ))
    dir = flt.replace(".flt","_mapdir")
    os.mkdir(dir)
    mytransformer.loadfiltered(  flt )
    mytransformer.loadfileparameters(  'krill.prm' )
    mytransformer.updateparameters( )

    first = True

    for tx in range(-750, 751, 25):
        for ty in range(-750, 751, 25):  
            mytransformer.parameterobj.set_parameters(  {
              't_x': tx,
              't_y': ty })
            
            mytransformer.compute_tth_eta( )
            mytransformer.computegv( )
            if first:
                mytransformer.savegv(  'map.gve' )
                myindexer = indexing.indexer()
                myindexer.readgvfile(  'map.gve' )
                first = False
            else:
                # Do what readgvfile does...
                # Set grain assignments to blanks
                # Copy in the g-vectors
                gv = np.array([ 
                    mytransformer.getcolumn("gx"),
                    mytransformer.getcolumn("gy"),
                    mytransformer.getcolumn("gz") ] ).T
                myindexer.gv = gv
                myindexer.ubis = []
                myindexer.ga = -1*np.ones( (gv.shape[0],)  ,int)
                myindexer.gvflat = np.reshape(np.fromstring(
                    myindexer.gv.tostring(),float), myindexer.gv.shape)
                # ds for rings
                print(gv.shape)
                myindexer.ds = np.sqrt( gv[:,0]*gv[:,0] +
                                        gv[:,1]*gv[:,1] +
                                        gv[:,2]*gv[:,2] )
                myindexer.ds = myindexer.ds.clip( 0, 1.25 )
                myindexer.scores = []

                
            myindexer.updateparameters( )
            myindexer.parameterobj.set_parameters(  {
              'hkl_tol': tol,
              'minpks': npk,
              
              'ds_tol': '0.01',
              'max_grains': '100',
              'uniqueness': '0.5',
              'eta_range': '0.0',
              'ring_1': '3',
              'wavelength': '0.26464',
              'ring_2': '5',
              'cosine_tol': '0.01'} )
            myindexer.loadpars( )
            myindexer.assigntorings( )
            myindexer.find( )
            myindexer.scorethem( )
            mapfile.write("%d %d "%(tx,ty))
            if len(myindexer.scores)>0:
                myindexer.saveubis(  os.path.join(dir, 'x%d_y%d.ubi'%(tx,ty) ))
              
                mapfile.write("%d %d %s\n"%(len(myindexer.scores),
                                          max(myindexer.scores),
                                          str(myindexer.scores)))
            else:
                mapfile.write("0  0  []\n")
            mapfile.flush()
              

import sys

flt = sys.argv[1]
tol = float(sys.argv[2])
npk = int(sys.argv[3])

                                          
grids( flt, tol, npk )
