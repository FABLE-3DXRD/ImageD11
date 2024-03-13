
from __future__ import print_function, division
import sys, os
import numpy as np, pylab as pl
from ImageD11 import columnfile, unitcell

def find_pair( c,
               omegatol = 0.25,
               etatol = 0.25,
               tthtol = 0.1):
    omi = np.round( c.omega / omegatol ).astype( int )
    f180 = 180./omegatol
    indsc = np.arange(c.nrows, dtype=int)
    pi = []
    pj = []
    for oi in range( omi.min(), omi.max()+1):
        i = omi == oi
        j = abs(omi - (oi+f180)) < 1.5
        if i.sum() > 0 and j.sum() > 0:
            # potential pairs
            ei = c.eta[i]
            ej = c.eta[j]
            eji = np.where( ej < 0, -180-ej, 180-ej )
            deta = abs(ei - eji[:,np.newaxis]) # ni x nj
            dtth = abs( c.tth[i] - c.tth[j][:,np.newaxis] )
            #print(deta, deta.shape, i.sum(),j.sum())
            #print(dtth)
            m = (deta < etatol)&(dtth < tthtol)
            inds = np.arange( len(ei)*len(ej), dtype=np.int ) [ m.ravel() ]
            ip = inds % len(ei)
            jp = inds // len(ei)
            #print("Found",m.ravel().sum())
            
            #for k in range(len(inds)):
                #print(k, ip[k], jp[k], m[jp[k],ip[k]], i.sum(), j.sum())
            ic = np.compress( i, indsc)[ ip ]
            jc = np.compress( j, indsc)[ jp ]
            pi.append( ic )
            pj.append( jc )
                # print( k,ip[k], jp[k], c.tth[ic],c.tth[jc],c.eta[ic],c.eta[jc],c.omega[ic],c.omega[jc])
        #print(len(pi))
    pi = np.concatenate( pi )
    pj = np.concatenate( pj )
    print(pi.shape,pj.shape,c.nrows)
    return pi, pj

def calc_tth_eta( c, pi, pj ):
    dX = c.xl[pi] + c.xl[pj]
    dY = c.yl[pi] + c.yl[pj]
    dZ = c.zl[pi] - c.zl[pj]
    r = np.sqrt(dY*dY + dZ*dZ)
    tth = np.degrees( np.arctan2( r, dX )  )
    eta = np.degrees(np.arctan2( dY, dZ ))
    return tth, eta
    
def tth_ds( ds, w ):
    return np.degrees( 2 * np.arcsin( w * ds / 2. ) )

def strain( tthvals, etavals ):
    # cot(theta) . d(theta)
    th = np.radians(tthvals / 2)
    th0 = th.mean()
    strain = (th0 - th) / np.tan( th0 )
    # Fit this as being linear
    s2e = np.sin( np.radians(etavals) )**2
    p = np.polyfit( s2e, strain, 1 )
    scalc = np.polyval(p , s2e )
    return strain, th0, p, scalc, s2e

def d_tth( tth, w ):
    return w  / 2 / np.sin(np.radians(tth/2))

def fit_hkl( c, pi, pj, tth, eta, stem ):
    u = unitcell.unitcell_from_parameters( c.parameters )
    u.makerings( c.ds.max()*1.1, tol=0.001 )
    w = c.parameters.get( "wavelength" )
    tthv = [ tth_ds( ds, w ) for ds in u.ringds ]
    print("(h, k, l)  tthcalc  <tthobs>  Gradient_sin2p Intercept_sin2p")
    for tth0,d in zip(tthv[:46], u.ringds):
        print( u.ringhkls[d][0], end=" " )
        m = abs(tth - tth0) < 0.025 # tolerance in tth
        try:
            s,th0,p, sc, s2e = strain( tth[m], eta[m] )
            print( tth0, np.degrees(th0*2), p[0], p[1] )
            sys.stdout.flush()
        except:
            print()
            continue
        try:
            f=pl.figure(2, dpi=200, figsize=(20,15))
            f.clf()
            f.add_subplot(2,3,1)
            #pl.subplot(231)
            pl.plot( eta[m], tth[m],"+" )
            pl.ylabel("tth")
            f.add_subplot(2,3,2)
                        #pl.subplot(23)
            pl.plot( eta[m], d_tth(tth[m],w),"+" )
            pl.ylabel("d-spacing")
            #pl.subplot(233)
            f.add_subplot(2,3,3)
            pl.plot( eta[m], s, "+")
            pl.plot( eta[m], sc, "x")
            #pl.subplot(234)
            f.add_subplot(2,3,4)
            pl.plot( s2e,tth[m],"+" )
            pl.ylabel("tth")
            #pl.subplot(235)
            f.add_subplot(2,3,5)
            pl.plot( s2e,d_tth(tth[m],w),"+" )
            pl.ylabel("d-spacing")
            f.add_subplot(2,3,6)
            pl.plot( s2e, s, "+")
            pl.plot( s2e, sc, "x")
            pl.ylabel("strain")
            hkl = "_%d_%d_%d.png"%(tuple(u.ringhkls[d][0]))
            f.savefig( stem + hkl)
#            pl.show()
            

        except:
            print("Fail")
            raise

    
if __name__=="__main__":
    c = columnfile.columnfile(sys.argv[1])
    stem = os.path.join("png", sys.argv[1].split(".")[0])
    c.parameters.loadparameters( sys.argv[2] )
    c.updateGeometry()
    d=c.copy()
    pi, pj = find_pair(c)
    tth, eta = calc_tth_eta( c, pi, pj )
    pl.subplot(221)
    pl.plot(d.tth, d.eta,"r.")
    pl.plot( tth, eta, "g.")
    pl.subplot(222)
    pl.plot(d.tth, d.eta,"r.")
    pl.plot( tth, eta, "g.")
    pl.xlim( 3,15.1)
    pl.subplot(223)
    pl.plot(d.tth, d.eta,"r.")
    pl.plot( tth, eta, "g.")
    pl.xlim( 11,11.75)
    pl.subplot(224)
    pl.plot(d.tth, d.eta,"r.")
    pl.plot( tth, eta, "g.")
    pl.xlim( 14.3,14.85)
#    pl.show()
    pl.savefig( stem+"_all.png")
    f = open( sys.argv[1].replace(".hdf",".tth_eta"), "w")
    f.write("# From %s\n"%(sys.argv[1]))
    f.write("#  tth  eta \n")
    for i in range(len(tth)):
        f.write("%f %f\n"%(tth[i], eta[i]))
    f.close()
    fit_hkl( c, pi, pj, tth, eta, stem )
    
