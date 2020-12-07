
from __future__ import print_function, division
import numpy as np
import pylab as pl

"""
Look for a new figure of merit to replace drlv

[hkl] = [[ubi]] . [gve]

"""

def norm3(vec):
    return np.sqrt( (vec*vec).sum(axis=0) )

def pkfom( gr, gvectors, axis=(0,0,1) ):
    """
    gr is a grain carrying ub and ubi
    gvector is an array of g-gvectors (shape (n,3))

    projections onto a gvector-axis co-ordinate system
    Normalised to length of vector (strain + 2 angles in radians)
    
    Critique : not pixels and image number
    """
    assert np.dot(axis, axis) == 1
    assert gvectors.shape[0] == 3
    #    assert (np.dot(gr.ub, gr.ubi) == np.eye(3)).all()
    hrkrlr = np.dot( gr.ubi, gvectors )
    hikili = np.round( hrkrlr )
    gcalc  = np.dot( gr.ub, hikili )    #  units A-1
    gerror = gvectors - gcalc           #  A-1
    modg   = norm3( gvectors )          #  A-1
    gen    = gerror / modg              # relative error in g (pure number)
    normg  = gvectors / modg            #    pure number
    omdir   = np.cross( normg.T, axis )   # ~rotation (pure number)
    etadir  = np.cross( omdir, normg.T )  # ~azimuth  (pure number)
    e0 = ( gen * normg    ).sum(axis=0)  # ~parallel (pure number, strain)
    e1 = ( gen * etadir.T ).sum(axis=0)  #           (pure number, angle rad)
    e2 = ( gen * omdir.T  ).sum(axis=0)  # Rotation  (pure number, angle rad)
    return e0,e1,e2
    

def test_eu():
    from ImageD11.grain import read_grain_file
    from ImageD11.columnfile import columnfile
    from ImageD11.parameters import read_par_file
    gr  = read_grain_file("./eu3.map")[0]
    par = read_par_file("./0.par")
    flt = columnfile( "../demo/eu.flt" )
    flt.updateGeometry(par)
    gve = np.array( [flt.gx, flt.gy, flt.gz], float)
    e0,e1,e2 = pkfom( gr, gve )
    flt.addcolumn( e0, "e0" )
    flt.addcolumn( e1, "e1" )
    flt.addcolumn( e2, "e2" )
    pl.subplot(131)
    pl.scatter( flt.tth, flt.e0, c=flt.eta, edgecolor='none')
    pl.subplot(132)
    pl.scatter( flt.tth, flt.e1, c=flt.eta, edgecolor='none')
    pl.subplot(133)
    pl.scatter( flt.tth, flt.e2, c=flt.eta, edgecolor='none')
    pl.show()
    return gr, flt, par
    
    
    
if __name__=="__main__":
    test_eu()
    


