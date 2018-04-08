import pylab as pl
import sys
import numpy as np

from ImageD11 import transform, grain, columnfile, parameters



def mod360( angs, targets ):
    diff = np.asarray( angs ) - targets
    turns = np.round( diff / 360.0 )
    return angs - turns * 360.

def calccolumns( p, c, g, pks = None):
    """
    p  = parameters
    c  = columnfile
    g  = grain
    pks = selection of which peaks to use for the grain. If none we do all.
    """
    t = g.translation
    p.set('t_x', t[0] )
    p.set('t_y', t[1] )
    p.set('t_z', t[2] )
    # tth, eta of the spots given this origin position
    if pks is None:
        pks = ~np.isnan( c.sc )
    tth, eta = transform.compute_tth_eta( [c.sc[pks], c.fc[pks]],
                                          omega = c.omega[pks],
                                          **p.parameters )
    # g-vectors given this origin
    gve = transform.compute_g_vectors( tth, eta, c.omega[pks],
                                       p.get("wavelength"),
                                       p.get("wedge"),
                                       p.get("chi") )
    # Find hkl indices
    hklo = np.dot( g.ubi, gve )
    hkli = np.round( hklo )
    diff = hklo - hkli
    # Now uncompute to get ideal spot positions in each case:
    gcalc = np.dot( g.ub, hkli )
    tthcalc, etacalc, omegacalc = transform.uncompute_g_vectors(
        gcalc, p.get("wavelength"), p.get("wedge"), p.get("chi"))
    # which eta to choose - there are two possibilities
    #   ... eta should always be in -180, 180 as it is arctan2
    smatch = np.sign(eta) == np.sign( etacalc[0] )
    etachosen = np.where( smatch , etacalc[0], etacalc[1] )
    omegachosen = np.where( smatch , omegacalc[0], omegacalc[1] )
    # deal with the mod360 for omega
    omegachosen = mod360( omegachosen, c.omega[pks] )
    fcc , scc = transform.compute_xyz_from_tth_eta(tthcalc,
                                                   etachosen,
                                                   omegachosen,
                                                   **p.parameters)
    # how close to reciprocal lattice point (could throw out some peaks here...)
    drlv = np.sqrt((diff*diff).sum(axis=0))
    # save arrays to colfile:
    c.drlv[pks] = drlv
    c.tth_per_grain[pks] = tth
    c.eta_per_grain[pks] = eta
    c.tthcalc[pks] = tthcalc
    c.etacalc[pks] = etachosen
    c.omegacalc[pks] = omegachosen
    c.sccalc[pks] = scc
    c.fccalc[pks] = fcc
    return c
    

def addcols( p, c, gl ):
    c.addcolumn(np.zeros(c.nrows), "omegacalc")
    c.addcolumn(np.zeros(c.nrows), "etacalc")
    c.addcolumn(np.zeros(c.nrows), "drlv")
    c.addcolumn(np.zeros(c.nrows), "tthcalc")
    c.addcolumn(np.zeros(c.nrows), "sccalc")
    c.addcolumn(np.zeros(c.nrows), "fccalc")
    c.addcolumn(np.zeros(c.nrows)-1, "tth_per_grain")
    c.addcolumn(np.zeros(c.nrows)-1, "eta_per_grain")    
    for i,g in enumerate(gl):
        pks = c.labels == i
        c = calccolumns( p, c, g, pks)
        print(i,pks.sum(),(c.tth_per_grain>0).sum())
    return c

if __name__=="__main__":
    p = parameters.read_par_file( sys.argv[1] )
    gl = grain.read_grain_file( sys.argv[2] )
    c = columnfile.columnfile( sys.argv[3] )
    c = addcols( p, c, gl )
    c.writefile(sys.argv[4])

    
