
from __future__ import print_function

from ImageD11.columnfile import columnfile
from ImageD11.parameters import read_par_file
from ImageD11.unitcell   import unitcell_from_parameters
from ImageD11.grain      import read_grain_file, write_grain_file
from ImageD11 import transform
import numpy as np

# idea:
#    Rotate XL,YL,ZL into laboratory frame
#    Rotate beam into laboratory frame
#    Make a new compute-g-vectors where rotations were already applied
#    Fit grain position and orientation in this frame


def getCxyz(colf, pars):
    """
    Find spot positions in crystal frame
    colf = columnfile (or other object) holding .sc, .fc, .omega
    pars = parameters for ImageD11 transform module
    """
    wedge = pars.get("wedge")
    chi   = pars.get("chi")
    peaks_xyz = transform.compute_xyz_lab([ colf.sc,
                                            colf.fc ],
                                          **pars.parameters)
    fun = transform.compute_g_from_k
    peaks_Cxyz = fun( peaks_xyz, colf.omega, 
                      wedge=wedge, 
                      chi=chi)
    b = np.zeros(peaks_Cxyz.shape)
    b[0,:] = 1.0/pars.get("wavelength")
    beam_Cxyz =  fun( b, colf.omega, 
                      wedge=wedge, 
                      chi=chi)
    return peaks_Cxyz, beam_Cxyz

def compute_Cgve(txyz, peaks_Cxyz, beam_Cxyz, wavelength):
    """
    Computes g-vectors from spots and beam in crystal frame
    txyz = (3,) origin position
    peaks_Cxyz = (3,n) spot position on rotated detector
    beam_Cxyz  = (3,n) vectors of rotated incident beam
                       beam is already normalised to wavelength
    wavelength = float radation, normalises length
    returns scattering vectors (g-vectors)
    """
    # Output ray - full length
    r_f = (peaks_Cxyz.T - txyz).T
    # Normalise to unit vector and scale for wavelength
    fac =  wavelength * np.sqrt( (r_f*r_f).sum(axis=0) )
    np.divide( r_f, fac, r_f )
    # Compute scattering vector
    np.subtract( r_f, beam_Cxyz, r_f )
    return r_f

def compute_dgdt( txyz, peaks_Cxyz, beam_Cxyz, wavelength):
    """
    Compute the gvectors and the derivatives with respect to tx,ty,tz
    txyz = (3,) origin position
    peaks_Cxyz = (3,n) spot position on rotated detector
    beam_Cxyz  = (3,n) vectors of rotated incident beam
                       beam is already normalised to wavelength
    wavelength = float radation, normalises length
    returns scattering vectors (g-vectors) and d(g)/d(txyz)
    """
    r = (peaks_Cxyz.T - txyz).T # 3xn
    # Normalise to unit vector and scale for wavelength
    f = wavelength * np.sqrt( (r*r).sum(axis=0) )
    # d|r| -> r/|r| and f/wvln = |r|
    dfdT = - wavelength * wavelength * r / f # 3xn
    # over-writes in place
    np.divide( r, f, r )
    # quotient rule d(r/f) : (r'.f - r.f')/f.f
    drdT = np.zeros((3,3,len(f)))
    # d( r/f )/dT -> r is a vector, r[0], r[1], r[2]
    #                f is a scalar
    # o = np.ones(len(f))
    # please derive/explain this... note r is already r/f
    drdT[0] = ( (dfdT[0] * r).T + (1,0,0) ).T / f
    drdT[1] = ( (dfdT[1] * r).T + (0,1,0) ).T / f
    drdT[2] = ( (dfdT[2] * r).T + (0,0,1) ).T / f
    np.subtract( r, beam_Cxyz, r )
    return r, drdT


def fit_ub_t( ub, translation, hkl, peaks_Cxyz, beam_Cxyz, wavelength):
    """
    Fits the ub and grain origin to a list of assigned peaks
    All unit cell and orientations parameters are free
    Runs 2 cycles (empirically this converges)
    
    ub = (3,3) input UB matrix (so that h ~= UB.g)
    translation = (3,) grain origin
    hkl = (3,n) integer peak assignments
    peaks_Cxyz = (3,n) spot positions in crystal frame (getCxyz)
    beam_Cxyz = (3,n) beam directions in crystal frame (getCxyz)
    wavelength = float, radiation, normalises length
    returns fitted UB and translation
    """
    npk = len(hkl[0])
    dg = np.zeros( (12,3,npk) )
    for i in range(3):
        for j in range(3):
            dg[i*3+j,i] = hkl[j] # ub in 0->8
            # Why? :
            # d(gcalc3)/d(ub9) = h
            # dgdub = np.zeros((3,3,3,c.nrows))
            # g0 = ub00.h0 + ub01.h1 + ub02.h2
            # g1 = ub10.h0 + ub11.h1 + ub12.h2
            # g2 = ub20.h0 + ub21.h1 + ub22
            # only 9 out of 27 non-zero. Should roll this another way. (outer h,h.T)
            # ubi ubj gk
            # print dgdub[0,0,0].shape, hi[0].shape
            # dgdub[0,0,0]=hi[0]
            # dgdub[0,1,0]=hi[1]
            # dgdub[0,2,0]=hi[2]
            # dgdub[1,0,1]=hi[0]
            # dgdub[1,1,1]=hi[1]
            # dgdub[1,2,1]=hi[2]
            # dgdub[2,0,2]=hi[0]
            # dgdub[2,1,2]=hi[1]
            # dgdub[2,2,2]=hi[2]
            # gobs = g0 + dgdT.t
            # gobs - gcalc = dgdT.t
            # t = (gobs - gcalc).
    tnew = translation.copy()
    ubnew = ub.copy()
    # empirically it converges to 3 decimal places in 1 cycle
    # ...since: dgobsdt seems to depend on t we run a couple of cycles
    for _ in range(2):
        gobs, dgobsdt = compute_dgdt( tnew, peaks_Cxyz, beam_Cxyz, wavelength )
        # Note dgdub=h does not change here
        for i in range(3):
            dg[i+9] = dgobsdt[i]       # translation at 9,10,11
        gcalc =  np.dot( ubnew , hkl ) # 3xn
        gdiff = gcalc - gobs
        #    print((gdiff*gdiff).ravel().sum(),tnew)
        dg.shape = 12,3*npk
        mat = np.dot( dg, dg.T )
        rhs = np.dot( dg, gdiff.ravel() )
        imat = np.linalg.inv( mat )
        shifts = np.dot (imat, rhs )
        dg.shape =  (12,3,npk) 
        ubnew = ubnew  - np.reshape(shifts[:9],(3,3))
        tnew  = tnew - shifts[9:]
    return ubnew, tnew

def fitagrain( gr, pars ):
    """
    """
    t = gr.translation.copy()
    ub = gr.ub.copy()
    pks, beam = getCxyz( gr, pars )
    # This is a little ugly. Can the wavelength live in the B matrix?
    hi = np.round( np.dot( gr.ubi, compute_Cgve( t, pks, beam, pars.get("wavelength") ) ) )
    ubnew, tnew = fit_ub_t( ub, t, hi, pks, beam, pars.get("wavelength"))
    return ubnew, tnew


def main():
    import sys, time
    c = columnfile(sys.argv[1])
    p = read_par_file(sys.argv[2])
    u = unitcell_from_parameters(p)
    gl = read_grain_file(sys.argv[3])
    if gl[0].translation is None:
        gl[0].translation = np.array((0.,0.,0.))
    start = time.time()
    # Setup and assign hkls
    w = p.get("wavelength")
    peaks_Cxyz, beam_Cxyz = getCxyz( c, p )
    t = gl[0].translation.copy()
    ub = gl[0].ub.copy()
    ubi = gl[0].ubi.copy()
    gve = compute_Cgve(t, peaks_Cxyz, beam_Cxyz, w)
    hi = np.round( np.dot( ubi, gve ) )
    lastgof = 1e9
    ubn, tn = fit_ub_t( ub, t, hi, peaks_Cxyz, beam_Cxyz, w)
    print("Before\nt=",gl[0].translation)
    print("UB=",gl[0].ub)
    gl[0].set_ubi( np.linalg.inv( ubn ) )
    gl[0].translation = tn    
    dt = time.time() - start
    print("time calculating",dt,"gps",1/dt)
    print("After\nt=",gl[0].translation)
    print("UB=",gl[0].ub)
    write_grain_file(sys.argv[4],gl)


def main2():
    import sys
    c = columnfile(sys.argv[1])
    p = read_par_file(sys.argv[2])
    gl = read_grain_file(sys.argv[3])
    for i,g in enumerate(gl):
        mask = c.labels == i
        g.sc = np.compress( mask, c.sc )
        g.fc = np.compress( mask, c.fc )
        g.omega = np.compress( mask, c.omega )
        ubnew, tnew = fitagrain( g, p )
        g.set_ubi( np.linalg.inv( ubnew ) )
        g.translation[:] = tnew
        print(i, len(g.sc), tnew)
    write_grain_file( sys.argv[4], gl )
        




if __name__=="__main__":
    main2()






# python rotdex.py g10.flt cubic.par g10.ubi g10.t1
