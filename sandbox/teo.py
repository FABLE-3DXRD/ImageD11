

# Try fitting using theta/eta/omega instead of x/y/omega or gx/gy/gz

# Each peak comes with sc, fc, omega
#  convert this to tth, eta, omega and dtth/dt, deta/dt

# omega depends on UB (not translations)

# sc,fc depends on UB and translations

# tth depends on B, translations but not U

# eta depends on UB and translations

# variables to fit: U(3), B(6), t(3)
#  B0 B
# write U as Unow + Ushift...


# TODO :
# ... replace numeric derivatives by analytical
# ... reformulate UB as rx,ry,rz,a*,b*,c*,alpha* etc
# ... refactor LSQ to allow symmetry constraints



from ImageD11 import columnfile, grain, parameters, transform, gv_general, indexing
import sys
import numpy as np

TESTING=0

if TESTING:
    import pylab as pl
    pl.ion()

def Dcalc( UB, hkls, wedge, chi, wavelength, etao):
    """
    ysign determines which solution (there are usually two)
    """
    g = np.dot( UB  , hkls )
    if 0:
        print g.shape
        print hkls.shape
        for i in range(10):
            print g[:,i],
            print hkls[:,i],
            print np.dot(np.linalg.inv( UB ), g[:,i] )

    wc = gv_general.wedgechi( wedge=wedge, chi=chi )
    omega1, omega2, valid = gv_general.g_to_k(
        g, wavelength,axis=[0,0,-1], pre=None, post=wc )

    cwT = gv_general.chiwedge( wedge=wedge, chi=chi ).T
    
    k1 = gv_general.k_to_g( g, omega1, axis=[0,0,1],
                            pre = cwT, post=None)
    k2 = gv_general.k_to_g( g, omega2, axis=[0,0,1],
                            pre = cwT, post=None)
    #
    # k[1,:] = -ds*c*sin(eta)
    # ------    -------------   .... tan(eta) = -k1/k2
    # k[2,:] =  ds*c*cos(eta)
    #
    eta_one = np.degrees( np.arctan2(-k1[1, :], k1[2, :]) )
    eta_two = np.degrees( np.arctan2(-k2[1, :], k2[2, :]) )

    # pick the solution which is closest to matching in eta/omega
    d1 = eta_one - etao
    d2 = eta_two - etao
    etac   = np.where( abs(d1) < abs(d2), eta_one, eta_two)
    omegac = np.where( abs(d1) < abs(d2), omega1 , omega2 )

    mod_g = np.sqrt((g*g).sum(axis=0))
    tthc = np.degrees( 2*np.arcsin(wavelength*mod_g/2) )
    # Dtth / dB
    # Dtth dT == 0, this is the computed value
    return  tthc, etac, omegac

def Dobs( t, sc, fc, omega, pars ):
    mypars = pars.parameters.copy()
    mypars['t_x'] = t[0]
    mypars['t_y'] = t[1]
    mypars['t_z'] = t[2]
    tth0, eta0 = transform.compute_tth_eta( (sc, fc), omega=omega, **mypars )
    dtth=np.zeros( (3,len(sc)) )
    deta=np.zeros( (3,len(sc)) )
    s = 1.0
    for i,p in enumerate(("t_x","t_y","t_z")):
        # reset
        mypars['t_x'] = t[0]
        mypars['t_y'] = t[1]
        mypars['t_z'] = t[2]
        # shift + 
        mypars[p] = t[i] + s/2.0
        tth1, eta1 = transform.compute_tth_eta( (sc, fc),
                                                omega=omega, **mypars )
        mypars[p] = t[i] - s/2.0
        tth2, eta2 = transform.compute_tth_eta( (sc, fc),
                                                omega=omega, **mypars )
        dtth[ i ] = (tth1 - tth2)/s
        deta[ i ] = (eta1 - eta2)/s # check wraps at 180 !!!
    return tth0, eta0, dtth, deta

def angmod( a ):
    ar = np.radians(a)
    return np.degrees( np.arctan2( np.sin( ar ), np.cos( ar ) ) )

def fitone( UB, t, sc, fc, omega, hkls, pars):
    npks = len(omega)
    ttho, etao, dtth, deta =  Dobs( t, sc, fc, omega, pars )
    tthc, etac, omegac = Dcalc( UB, hkls,
                                pars.get('wedge'), pars.get('chi'),
                                pars.get('wavelength'),
                                etao)

    s = 1e-5
    # Gradient arrays (move to Dcalc)
    dtthUB = np.zeros((3,3,npks), np.float)
    detaUB = np.zeros((3,3,npks), np.float)
    domUB = np.zeros((3,3,npks), np.float)
    for i in range(3):
        for j in range(3):
            UBc = UB.copy()
            UBc[i,j] -= s/2
            tths1, etas1, omegas1 = Dcalc( UBc, hkls,
                                        pars.get('wedge'), pars.get('chi'),
                                        pars.get('wavelength'),
                                        etao)
            UBc = UB.copy()
            UBc[i,j] += s/2
            tths2, etas2, omegas2 = Dcalc( UBc, hkls,
                                        pars.get('wedge'), pars.get('chi'),
                                        pars.get('wavelength'),
                                        etao)
            dtthUB[i,j] = (tths1 - tths2)
            detaUB[i,j] = (etas1 - etas2)
            domUB[i,j] = (omegas1 - omegas2)
            if TESTING:
                pl.subplot(3,3,i*3+j+1)
                pl.plot(detaUB[i,j])

    # Differences
    errtth = tthc - ttho
    erreta = etac - etao
    erromega = angmod( omegac - omega )

    if TESTING:
        pl.figure()
        pl.subplot(311)
        pl.plot(etac,errtth, "+")
        pl.subplot(312)
        pl.plot(etac, erromega, "+")
        pl.subplot(313)
        pl.plot(etac, erreta, "+")
        pl.show()
        raw_input()

    
    # Least squares - move to subroutine
    m = np.zeros( (12,12), np.float )
    r = np.zeros( (12,), np.float )
    for i in range(3):
        r[i]   += np.dot( dtth[i], errtth )
        r[i]   += np.dot( deta[i], erreta )
        for j in range(3):
            m[i,j] += np.dot( dtth[i], dtth[j] )
            m[i,j] += np.dot( deta[i], deta[j] )
    # ub
    for i in range(3):
        for j in range(3): # loop over UB
            k = i*3+j+3
            r[k] += np.dot( dtthUB[i,j], errtth )
            r[k] += np.dot( detaUB[i,j], erreta )
            r[k] += np.dot( domUB[i,j], erromega )
            for l in range(3):
                for o in range(3):
                    n = l*3+o+3
                    m[k,n] += np.dot( dtthUB[i,j], dtthUB[l,o] )
                    m[k,n] += np.dot( detaUB[i,j], detaUB[l,o] )
                    m[k,n] += np.dot( domUB[i,j], domUB[l,o] )
    try:
        im = np.linalg.inv(m)
        sh = np.dot( im, r )
    except:
        print "Failed"
        sh = np.ones(12)
        pass
    # End Least squares subroutine ?
    #
    # Apply shifts
    UBnew = np.reshape( UB.ravel() + s*sh[-9:], (3,3) )
    return t+sh[:3], UBnew
    

def refit_makemap( colf, pars, grains ):
    for i,g in enumerate( grains ):
        d = colf.copy()
        d.filter( d.labels == i )
        sc = d.sc
        fc = d.fc
        om = d.omega
        hkls = np.array( ( d.h, d.k, d.l ) )
        #print "Translation",g.translation
        #print "Cell",indexing.ubitocellpars(np.linalg.inv(g.ub))        
        t, UB = fitone( g.ub, g.translation, sc, fc, om, hkls, pars)
        #print "Translation",t
        #print "Cell",indexing.ubitocellpars(np.linalg.inv(UB))
        for i in range(3):
            t, UB = fitone( UB, t, sc, fc, om, hkls, pars)
            if TESTING:
                print "translation",t
                print "Cell",indexing.ubitocellpars(np.linalg.inv(UB))
        print "Translation",t
        print "Cell",indexing.ubitocellpars(np.linalg.inv(UB))
        g.translation = t
        g.set_ubi = np.linalg.inv( UB )
    return grains

if __name__=="__main__":
    colfile, parfile, grainsfile, newgrainsfile = sys.argv[1:5]
    c = columnfile.columnfile( colfile )
    p = parameters.read_par_file( parfile )
    g = grain.read_grain_file( grainsfile )
    
    grain.write_grain_file( newgrainsfile, refit_makemap( c, p, g) )


    # python teo.py ../test/simul_1000_grains/g187.flt ../test/simul_1000_grains/Al1000/Al1000.par ../test/simul_1000_grains/g187.map teo187fitted.map
