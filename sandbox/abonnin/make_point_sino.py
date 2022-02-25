
from __future__ import print_function

import numpy as np, pylab as pl


def read_log(fname):
    grains = []
    spots = []
    for line in open(fname).readlines():
        if line.find("Grain") == 0 and line.find(",")>0:
            if len(spots)>0:
                grains.append(spots)
            spots = []
        if line[0]!=" ":
            continue
        vals = line.split()
        if len(vals) == 22:
            spots.append(int(vals[2]))

    if len(spots)>0:
        grains.append(spots)
    return grains

def make_point_sino( cf, points ):
    indices = []
    for p in points:
        try:
            i = cf.indices.index(p)
        except:
            print(p)
            raise
    indices = [cf.indices.index( p ) for p in points]
    
    om = cf.omega.take( indices )
    xpos =cf.xpos.take( indices )
    sum_intensity = cf.sum_intensity.take( indices )
    return om, xpos, sum_intensity

def sino_func( pars, omega, obs ):
    r = np.pi * omega/ 180.0
    calc = pars[0] + np.sin( r )*pars[1] + np.cos( r )*pars[2]
    # print pars, obs.shape, r.shape
    return obs - calc 

def fit_sino( omega, obs, debug=False ):
    from scipy.optimize import leastsq
    pars = np.array([ 15.0, 1.0, 1.0 ])
    mask = np.ones( len(omega) )
    indice = list(range( len(omega)))
    for i in range(10):
        omeganow = np.compress( mask, omega )
        obsnow   = np.compress( mask, obs)
        indice   = np.compress( mask, list(range(len(omega))))
        print("# npks",len(obsnow), len(omeganow))
        newpars, status = leastsq( sino_func, pars, (omeganow, obsnow) )
        diff = sino_func( newpars, omeganow, obsnow )
        error = (diff-np.median(diff))/diff.std()
        bad = np.compress( abs(error)>2.5, indice )
        if len(bad) == 0:
            break
        np.put( mask , bad, 0 )
        if debug:
            pl.ion()
            pl.clf()
            pl.plot(omega, obs, "o")
            pl.plot(omeganow, obsnow, "o" )
            pl.plot(omeganow, obsnow-diff, "+")
            pl.show()
            input("Hit a key to continue")
    print("#  Constant  A*sin(om)   B*cos(om)")
    print("%f %f %f"%tuple(newpars))
    return omeganow, obsnow, obsnow-diff, indice


if __name__=="__main__":
    import sys
    grains = read_log( sys.argv[1] )
    print("# Number of grains", len(grains))
    from ImageD11.columnfile import columnfile
#    import matplotlib.pylab as pl
    c = columnfile( sys.argv[2] )
    c.indices = list(c.spot3d_id.astype(int))
    ng = 0
    for g in grains:
        print("\n\n#",len(g))
        om, xpos, intens = make_point_sino( c, g)
        omused, xposused, calc, indice = fit_sino( om, xpos )
        usedpks = np.take(g, indice)
        intens = np.take( intens, indice)
        print("#  spot3d_id  omega  xpos  intensity")
        for i in range(len(usedpks)):
            print("%-7d  %8.1f %9.5f %15g"%(usedpks[i], omused[i], xposused[i], intens[i]))
        # if 0:
        #     pl.clf()
        #     pl.title("Grain number %d"%(ng))
        #     pl.plot(om, xpos,"o")
        #     pl.plot(omused, xposused,"o")
        #     pl.plot(omused, calc,"+")
        #     pl.ion()
        #     pl.show()
        #     pl.savefig("grainfit_%d.png"%(ng))
        ng += 1
