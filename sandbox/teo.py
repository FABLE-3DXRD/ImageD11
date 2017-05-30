

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

def Dcalc( UB, hkls, wedge, chi, wavelength, etao, omega):
    """
    ysign determines which solution (there are usually two)
    """
    g = np.dot( UB  , hkls )
    assert g.dtype == np.float
    if 0:
        print g.shape
        print hkls.shape
        for i in range(10):
            print g[:,i],
            print hkls[:,i],
            print np.dot(np.linalg.inv( UB ), g[:,i] )

    wc = gv_general.wedgechi( wedge=wedge, chi=chi )
    assert wc.dtype == np.float
    omega1, omega2, valid = gv_general.g_to_k(
        g, wavelength,axis=[0,0,-1], pre=None, post=wc )
    assert omega1.dtype == np.float
    cwT = gv_general.chiwedge( wedge=wedge, chi=chi ).T
    
    k1 = gv_general.k_to_g( g, omega1, axis=[0,0,1],
                            pre = cwT, post=None)
    k2 = gv_general.k_to_g( g, omega2, axis=[0,0,1],
                            pre = cwT, post=None)
    assert k1.dtype == np.float
    #
    # k[1,:] = -ds*c*sin(eta)
    # ------    -------------   .... tan(eta) = -k1/k2
    # k[2,:] =  ds*c*cos(eta)
    #
    eta_one = np.degrees( np.arctan2(-k1[1, :], k1[2, :]) )
    eta_two = np.degrees( np.arctan2(-k2[1, :], k2[2, :]) )

    # pick the solution which is closest to matching in eta/omega
    d1 = angmod(eta_one - etao)
    d2 = angmod(eta_two - etao)
    etac   = np.where( abs(d1) < abs(d2), eta_one, eta_two)
    omegac = np.where( abs(d1) < abs(d2), omega1 , omega2 )

    if TESTING > 2:
        err = abs(etac - etao)
        for i in range(len(err)):
            if err[i] < 1:
                continue
            print i, eta_one[i], eta_two[i], omega1[i], omega2[i]
            print "\t", etao[i], omega[i], etac[i], omegac[i]
        

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
        deta[ i ] = angmod(eta1 - eta2)/s 
    return tth0, eta0, dtth, deta

def angmod( a ):
    ar = np.radians(a)
    return np.degrees( np.arctan2( np.sin( ar ), np.cos( ar ) ) )

class SVDLSQ(object):
    def __init__(self, parnames):
        """
        List of parameter names for fitting.
        """
        self.parnames = parnames
        self.gradients = {}
        self.differences = []
        self.weights = []
        for name in self.parnames:
            self.gradients[name] = []

    def addobs(self, difference, gradients, weights):
        """
        Add in a series of observations with weights
        Records these ready to solve later
        """
        for name in self.parnames:
            self.gradients[name].append( gradients[name] )
        self.differences.append( difference )
        self.weights.append( weights )

    def solve( self ):
        """
        Solve the least squares problem using an SVD
        TODO ... fill in some residual and error information
        """
        N = len(self.parnames)
        M = 0
        for diff in self.differences:
            M += len(diff)
        A = np.zeros( (N, M), np.float )
        b = np.zeros( M, np.float )
        iM = 0
        #print N, M, [len(x) for x in self.differences]
        for i, (diff, wt) in enumerate(zip(self.differences, self.weights)):
            b[iM:iM+len(diff)] = diff*wt
            for j, name in enumerate(self.parnames):
                A[j,iM:iM+len(diff)] = self.gradients[name][i]*wt
            iM += len(diff)
        #assert (b == np.concatenate( self.differences )).all()
        try:
            U, s, V = np.linalg.svd( A, full_matrices = False, compute_uv = True )
        except np.linalg.LinAlgError:
            print "SVD failure!!!"
            raise
        # 
        invS = np.where( s > s.max()*1e-9, 1/s, 0 )
        self.imat = np.dot(  U*invS, V )
        solution = np.dot( self.imat , b )
        self.U = U
        self.s = s
        self.V = V
        self.errmat = np.dot( U*invS*invS, U.T )
        if TESTING > 0:
            lsqmat = np.dot( A, A.T )
            errmat = np.linalg.inv( lsqmat )
            print solution
            print errmat
            print svderrmat
            assert np.allclose( self.errmat, errmat )
        self.esds = np.sqrt( np.diag( self.errmat ) )
        self.condition = s.max() / s.min() # sqrt ?
        self.sh_esd = solution / self.esds
        return solution



def fitone( UB, t, sc, fc, omega, hkls, pars):
    npks = len(omega)

    pnames = "t0 t1 t2".split() + ["UB%d%d"%(i,j) for i in range(3) for j in range(3)]
    
    S = SVDLSQ( pnames )
    # To do : Allocate LSQ problem gradient array here and keep S in grain (?)
    gradients = {}


    
    ttho, etao, dtth, deta =  Dobs( t, sc, fc, omega, pars )
    tthc, etac, omegac = Dcalc( UB, hkls,
                                pars.get('wedge'), pars.get('chi'),
                                pars.get('wavelength'),
                                etao, omega)

    # Differences
    errtth = tthc - ttho
    erreta = etac - etao
    erromega = angmod( omegac - omega )


    if TESTING>1:
        pl.figure()
        pl.subplot(311)
        pl.plot(etac,errtth, "+")
        pl.subplot(312)
        pl.plot(etac, erromega, "+")
        pl.subplot(313)
        pl.plot(etac, erreta, "+")
        pl.show()
        raw_input()

    
    s = 1e-7
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
                                           etao, omega)
            UBc = UB.copy()
            UBc[i,j] += s/2
            tths2, etas2, omegas2 = Dcalc( UBc, hkls,
                                        pars.get('wedge'), pars.get('chi'),
                                        pars.get('wavelength'),
                                           etao, omega)
            dtthUB[i,j] = (tths1 - tths2)
            detaUB[i,j] = (etas1 - etas2)
            domUB[i,j] = (omegas1 - omegas2)
            if TESTING>1:
                pl.subplot(3,3,i*3+j+1)
                pl.plot(detaUB[i,j])
                

    
    for i in range(3):
        gradients["t%d"%(i)] = dtth[i]
        for j in range(3):
            gradients["UB%d%d"%(i,j)] = dtthUB[i,j]
    S.addobs( errtth, gradients, 10. )
    for i in range(3):
        gradients["t%d"%(i)] = deta[i]
        for j in range(3):
            gradients["UB%d%d"%(i,j)] = detaUB[i,j]
    S.addobs( erreta, gradients, 1. )
    for i in range(3):
        gradients["t%d"%(i)] = 0.
        for j in range(3):
            gradients["UB%d%d"%(i,j)] = domUB[i,j]
    S.addobs( erromega, gradients, 1. )
    
    sh = S.solve()
    
    UBnew = np.reshape( UB.ravel() + s*sh[-9:], (3,3) )
    
    return t+sh[:3], UBnew, S
    

def refit_makemap( colf, pars, grains ):
    global TESTING
    for i,g in enumerate( grains ):
        d = colf.copy()
        d.filter( d.labels == i )
        sc = d.sc
        fc = d.fc
        om = d.omega
        hkls = np.array( ( d.h, d.k, d.l ) )
#        print "Translation",g.translation
#        print "Cell",indexing.ubitocellpars(np.linalg.inv(g.ub))
        if i == 1004:
            TESTING=3
        t, UB, S = fitone( g.ub, g.translation, sc, fc, om, hkls, pars)
        for ncycle in range(5):
            t, UB, S = fitone( UB, t, sc, fc, om, hkls, pars)
            if S.sh_esd.max() < 1e-10:
                break
            if TESTING:
                print "translation",t
                print "Cell",indexing.ubitocellpars(np.linalg.inv(UB))
        print "Translation %.5f %.5f %.5f"%tuple(t)
        print "Cell %.7f %.7f %.7f %.8f %.8f %.8f"%(indexing.ubitocellpars(np.linalg.inv(UB)))
        g.translation = t
        g.set_ubi(np.linalg.inv( UB ))
    return grains

# TODO:
#
# - Least squares as SVD problem
# - Allow a constraint matrix to be used
# - Compute weights from experimental data (sig/cov etc)
# - Compute error estimates
# - Allow fitting of other geometrical parameters
# - Allow refinement in terms of other parameterisation of UB (e.g. U, B)
# - Documentation


if __name__=="__main__":
    colfile, parfile, grainsfile, newgrainsfile = sys.argv[1:5]
    c = columnfile.columnfile( colfile )
    p = parameters.read_par_file( parfile )
    g = grain.read_grain_file( grainsfile )    
    grain.write_grain_file( newgrainsfile, refit_makemap( c, p, g) )


    # python teo.py ../test/simul_1000_grains/g187.flt ../test/simul_1000_grains/Al1000/Al1000.par ../test/simul_1000_grains/g187.map teo187fitted.map
