
from __future__ import print_function


#
from ImageD11 import transform, transformer, indexing, parameters
from ImageD11.cImageD11 import score
import sys

    


def gridgrains( ul, flt, pars, 
                minx=-750, maxx=750, stepx=25,
                miny=-750, maxy=750, stepy=25,
                tol = 0.05,
                ) :
    trn = transformer.transformer()
    trn.loadfiltered( flt )
    trn.parameterobj = parameters.parameters( **pars )

    peaks = [trn.getcolumn(trn.xname),
             trn.getcolumn(trn.yname) ]
    peaks_xyz =  transform.compute_xyz_lab( peaks,
                       **trn.parameterobj.get_parameters() )
    
    omega = trn.getcolumn( trn.omeganame )
    trn.updateparameters()
    tx = minx - stepx
    n = 0
    sys.stderr.write("Using tol = %f\n"%(tol))
    while tx <= maxx:
         tx = tx + stepx
         ty = miny - stepy
         while ty <= maxy:
             ty = ty + stepy      
             trn.parameterobj.set_parameters(  {
                 't_x': tx,
                 't_y': ty })
             pars = trn.parameterobj.get_parameters()
             tth, eta = transform.compute_tth_eta_from_xyz(
                 peaks_xyz,
                 omega,
                 **pars)
             if 'omegasign' in pars:
                 om_sgn = float(pars["omegasign"])
             else:
                om_sgn =  1.0
                
             gv = transform.compute_g_vectors(
                 tth, eta, omega*om_sgn, **pars)
             print(tx, ty, end=' ')
             for ubi in ul:
                 ns = score( ubi, gv.T, tol )
                 print(ns, end=' ')
                 n += ns
             print()
    return n


if __name__ == "__main__":
    flt = sys.argv[1]
    par = sys.argv[2]
    grains = sys.argv[3]
    tol = float(sys.argv[4])

    ul = indexing.readubis( grains )

    po = parameters.parameters()
    po.loadparameters( par )
    pars = po.get_parameters()
    
    gridgrains( ul, flt, pars , tol = tol)
                
