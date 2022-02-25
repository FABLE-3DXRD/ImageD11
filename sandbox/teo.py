
from __future__ import print_function

from six.moves import input

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
    Finds "C"alculated 2thetac, etac, omegac for a grain
    (Note that these do not depend on the grain position)

    UB = 3x3 UB matrix
    hkls = 3xn list of peaks (normally integers)
    wedge = tilt around y under rotation axis
    chi = tilt around beam under rotation axis
    wavelength = of the radiation, scales UB
    etao = observed azimuth values (angle on powder ring)
    omega = observed(?) sample rotation angles

    ???    ysign determines which solution (there are usually two)
    """
    # computed g-vectors
    g = np.dot( UB  , hkls )
    assert g.dtype == float
    if 0:
        print(g.shape)
        print(hkls.shape)
        for i in range(10):
            print(g[:,i], end=' ')
            print(hkls[:,i], end=' ')
            print(np.dot(np.linalg.inv( UB ), g[:,i] ))
    
    # wedge/chi matrix for axis orientation
    wc = gv_general.wedgechi( wedge=wedge, chi=chi )
    assert wc.dtype == float
    # computed omega angles for reflections
    omega1, omega2, valid = gv_general.g_to_k(
        g, wavelength,axis=[0,0,-1], pre=None, post=wc )
    assert omega1.dtype == float
    # hum - isn't this wc-1 ? 
    cwT = gv_general.chiwedge( wedge=wedge, chi=chi ).T
    # k vectors, scattering in laboratory at angle omega1/2
    k1 = gv_general.k_to_g( g, omega1, axis=[0,0,1],
                            pre = cwT, post=None)
    k2 = gv_general.k_to_g( g, omega2, axis=[0,0,1],
                            pre = cwT, post=None)
    assert k1.dtype == float
    # Computed azimuth angles
    #
    # k[1,:] = -ds*c*sin(eta)
    # ------    -------------   .... tan(eta) = -k1/k2
    # k[2,:] =  ds*c*cos(eta)
    #
    eta_one = np.degrees( np.arctan2(-k1[1, :], k1[2, :]) )
    eta_two = np.degrees( np.arctan2(-k2[1, :], k2[2, :]) )
    # 
    # We select the one matching the observed data
    # FIXME - doesn't make sense to do this here
    # pick the solution which is closest to matching in eta/omega
    d1 = angmod(eta_one - etao)
    d2 = angmod(eta_two - etao)
    # Selecting which is the computed value for peaks
    etac   = np.where( abs(d1) < abs(d2), eta_one, eta_two)
    omegac = np.where( abs(d1) < abs(d2), omega1 , omega2 )

    if TESTING > 2:
        err = abs(etac - etao)
        for i in range(len(err)):
            if err[i] < 1:
                continue
            print(i, eta_one[i], eta_two[i], omega1[i], omega2[i])
            print("\t", etao[i], omega[i], etac[i], omegac[i])
        
    # Two theta angles
    mod_g = np.sqrt((g*g).sum(axis=0))
    tthc = np.degrees( 2*np.arcsin(wavelength*mod_g/2) )
    # Return computed values (no derivatives here)
    return  tthc, etac, omegac

def Dobs( t, sc, fc, omega, pars ):
    """
    Compute tth0, eta0 and derivatives given:
    t = (tx,ty,tz) grain position
    sc,fc = pixel positions of spots
    omega = sample rotation angle
    pars = ImageD11 parameters (dict) with wavelength, distance, tilts, etc

    returns tth, eta, dtth, deta
    """
    mypars = pars.parameters.copy()
    mypars['t_x'] = t[0]
    mypars['t_y'] = t[1]
    mypars['t_z'] = t[2]
    # tth, eta for the spots
    tth0, eta0 = transform.compute_tth_eta( (sc, fc), omega=omega, **mypars )
    # arrays for dtth/dt_{xyz} and deta/dt_{xyz}
    dtth=np.zeros( (3,len(sc)) )
    deta=np.zeros( (3,len(sc)) )
    # Step size for numerical derivative
    s = 1.0
    for i,p in enumerate(("t_x","t_y","t_z")):
        # reset
        mypars['t_x'] = t[0]
        mypars['t_y'] = t[1]
        mypars['t_z'] = t[2]
        # shift + and -
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
    """
    Brings angle into range of arctan2 (normally -180, 180)
    """
    ar = np.radians(a)
    return np.degrees( np.arctan2( np.sin( ar ), np.cos( ar ) ) )

class SVDLSQ(object):
    def __init__(self, parnames):
        """
        parnames = List of parameter names for fitting.
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

        difference, gradients, weights = block addins
        difference = [nobs]
        gradients = dict[ parnames ][nobs]
        weights = [nobs]
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
        N = len(self.parnames)        # Number of parameters
        M = 0                         # Differences
        for diff in self.differences: # block addins
            M += len(diff)
        # Least squares problem matrix A.x=b
        A = np.zeros( (N, M), float )
        b = np.zeros( M, float )
        iM = 0
        #print N, M, [len(x) for x in self.differences]
        # FIXME: This loop is doing a numpy.concatenate
        for i, (diff, wt) in enumerate(zip(self.differences, self.weights)):
            b[iM:iM+len(diff)] = diff*wt
            for j, name in enumerate(self.parnames):
                A[j,iM:iM+len(diff)] = self.gradients[name][i]*wt
            iM += len(diff)
        #assert (b == np.concatenate( self.differences )).all()
        # Compute the SVD of A:
        try:
            U, s, V = np.linalg.svd( A, full_matrices = False, compute_uv = True )
        except np.linalg.LinAlgError:
            print("SVD failure!!!")
            raise
        # Clip out small singular values:
        invS = np.where( s > s.max()*1e-9, 1/s, 0 )
        # Pseudo inverse matrix
        self.imat = np.dot(  U*invS, V )
        solution = np.dot( self.imat , b )
        # Save some stuff
        self.U = U
        self.s = s
        self.V = V
        self.errmat = np.dot( U*invS*invS, U.T )
        if TESTING > 0:
            lsqmat = np.dot( A, A.T )
            errmat = np.linalg.inv( lsqmat )
            print(solution)
            print(errmat)
#            print(errmat)
            assert np.allclose( self.errmat, errmat )
        self.esds = np.sqrt( np.diag( self.errmat ) )
        self.condition = s.max() / s.min() # sqrt ?
        self.sh_esd = solution / self.esds
        # These are shifts
        return solution



def fitone( UB, t, sc, fc, omega, hkls, pars):
    """
    Fit a grain
    """
    npks = len(omega)
    pnames = "t0 t1 t2".split() + ["UB%d%d"%(i,j) for i in range(3) 
                                for j in range(3)]
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
    erreta = etac - etao # Why not angmod here ? 
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
        input()

    
    s = 1e-7
    # FIXME: Gradient arrays (move to Dcalc)
    dtthUB = np.zeros((3,3,npks), float)
    detaUB = np.zeros((3,3,npks), float)
    domUB = np.zeros((3,3,npks), float)
    
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
                

    # Three block addins: twotheta, eta, omega:
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
    
    # Find shifts
    sh = S.solve()
    
    # new UB
    UBnew = np.reshape( UB.ravel() + s*sh[-9:], (3,3) )
    
    # returns new translation and new UB
    return t+sh[:3], UBnew, S
    

def refit_makemap( colf, pars, grains ):
    """
    Takes a columnfile, parameters and list of grains
    and repeats the refinement of each grain, keeping
    the same grain assignments
    Returns : new grains
    """ 
    global TESTING
    for i,g in enumerate( grains ):
        d = colf.copy()
        d.filter( d.labels == i )
        # Take peaks belonging to this grain only
        sc = d.sc
        fc = d.fc
        om = d.omega
        hkls = np.array( ( d.h, d.k, d.l ) )
        if i == 1004:
            TESTING=3
        # Get new t, UB, S
        t, UB, S = fitone( g.UB, g.translation, sc, fc, om, hkls, pars)
        for ncycle in range(5):
            t, UB, S = fitone( UB, t, sc, fc, om, hkls, pars)
            if S.sh_esd.max() < 1e-10:
                break
            if TESTING:
                print("translation",t)
                print("Cell",indexing.ubitocellpars(np.linalg.inv(UB)))
        print("Translation %.5f %.5f %.5f"%tuple(t))
        print("Cell %.7f %.7f %.7f %.8f %.8f %.8f"%(indexing.ubitocellpars(np.linalg.inv(UB))))
        g.translation = t
        g.set_ubi(np.linalg.inv( UB ))
    return grains

# TODO:
#
# DONE : Least squares as SVD problem
# - Allow a constraint matrix to be used
# - Compute weights from experimental data (sig/cov etc)
# DONE (but not weighted or printed) : Compute error estimates
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
