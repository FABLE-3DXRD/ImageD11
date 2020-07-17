
import numpy as np
from ImageD11 import indexing, fft_index_refac

o = indexing.indexer()
o.readgvfile("test.vecs")
g = fft_index_refac.grid( 128, 1., 5.)
g.gv_to_grid_new( o.gv )
g.fft()
g.props()
g.peaksearch( "pairFFT.pks" )
g.read_peaks( "pairFFT.pks" )
tol=0.1
tol2=tol*tol

def refine_vector( v, gv, tol=0.25, ncycles=25, precision=1e-6 ):
    """
    Input
    v  :  real space lattice vector
    gv : scattering vectors (no errors)
    tol : tolerance for initial assignment to peak
    ncycles : how much to refine
    precision : when to stop refining in v += s when |s|/|v| < precision
    Returns
    vref : refined real space lattice vector

    Uses:
    hr = v[0]*gv[i,0] + v[0]*gv[i,1] + v[1]*gv[i,2]  1 per peak
    hi = round(hr) = calc                            1 per peak
    diff = calc - hr                                 1 per peak
    gradient = d(diff)/dv[i] = gv[j,i]               3 per peak
    matrix = gradient[i].gradient[j]                 3x3 peaks^2
    rhs = gradient . diff                            3

    Solve for shifts and iterates reducing tolerance as 1/(ncycles+1)
    Peaks that change index during refinement are removed
    Stops on : 
       no peaks reassigned
       ncycles runs out
       |shift| / |vec| < precision
    """
    assert gv.shape[1] == 3
    vref = v.copy()
    mg = (gv*gv).sum(axis=1)
    mold = None
    wt = 1./mg # + mg.max())
    for i in range(ncycles):
        hr = np.dot( vref, gv.T )   # npks
        hi = np.round( hr )         # npks
        diff = hi - hr              # npks
        m = abs(diff) < tol/(i+1)   #  keep or not ?
        diff = (diff*wt)[m]
        grad = (gv.T*wt).T[m]
        rhs = np.dot( grad.T, diff )
        mat = np.dot( grad.T, grad )
        U,s,V = np.linalg.svd( mat )
        one_over_s = np.where( s/s.max() < 1e-6, 1, 1./s )
        mati = np.dot( U, np.dot(np.diag( one_over_s ), V ) )
        vshft = np.dot( mati, rhs )
        vref = vref + vshft
        if mold is not None:
            changed = mold ^ m   # XOR
            nch = changed.sum()
            # Down weight peaks that are different to before for next round
            if nch == 0:
                break
            wt[ changed ] = 0.
            #print nch,i
        mold = m
        slen = np.sqrt((vshft * vshft).sum())
        vlen = np.sqrt((vref*vref).sum())
#        print "shift = ",slen,vlen,slen/vlen, tol/(i+1)
        if abs(slen/vlen) < precision:
            break
    return vref
    
def scr(v): 
    h=np.dot(v, o.gv.T)
    return h-np.round(h)

    

idx = []
drl = []
vs = []
for v in g.UBIALL:
    v = refine_vector( v, o.gv )
    vs.append(v)
    hr = np.dot( v, o.gv.T )
    hi = np.round( hr )
    dr = (hr - hi)
    drlv2 = dr*dr
    drl.append( drlv2 )
    m = drlv2 < tol2
    idx.append( m )
    print( v, m.sum() )

idx = np.array( idx ) # dims vecs x npk
pair = np.zeros( (len(idx),len(idx)))
for i in range(len(g.UBIALL)):
    for j in range(i, len(g.UBIALL)):
        pair[i,j] = pair[j,i] = (idx[i]*idx[j]).sum()
        

import pylab as pl
pl.imshow( pair, interpolation='nearest' )
pl.show()

