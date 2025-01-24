
import sys
import numpy as np, pylab as pl
from ImageD11 import refinegrains, unitcell, transform


# Need:



parfile = sys.argv[1]
peaksfile = sys.argv[2]
grainsfile = sys.argv[3]
tolerance = 0.05
OmFloat = True
OmSlop = 0.5
nrings = 20

# one way ...
rgr = refinegrains.refinegrains( tolerance=tolerance,
                                 OmFloat=OmFloat,
                                 OmSlop=OmSlop)
rgr.loadparameters( parfile )
rgr.loadfiltered( peaksfile )
rgr.readubis( grainsfile )
rgr.generate_grains( )
print("Assigning peak labels")
rgr.assignlabels()

# Compute a histogram of omega

def anorm(x):
    o = np.radians(x)    
    return np.arctan2( np.sin( o ), np.cos( o ) )

omr = anorm( rgr.scandata[peaksfile].omega )
print(omr.min(), omr.max())
osl = np.radians( OmSlop )
obins = np.arange( -np.pi - osl/2, np.pi+osl, osl )
omega_histogram, _ = np.histogram( omr, bins = obins )

pl.figure()
pl.plot( (obins[1:]+obins[:-1])/2, omega_histogram, "-")
pl.title("Omega histogram")

print('Grain  #peaks')

for k in list(rgr.grains.keys()):
    g = rgr.grains[k]
    rgr.set_translation( *k )
    rgr.compute_gv( g, update_columns=True )
    print(k, g.ind.shape)
    hkl_real = np.dot( g.ubi, rgr.gv.T )
    mod_g = np.sqrt( rgr.gv * rgr.gv ).sum( axis = 1 )
    hkl_int = np.round( hkl_real ).astype(int)
    g_calc = np.dot( g.UB, hkl_int )
    g_error = rgr.gv.T - g_calc
    mod_g_error = np.sqrt( ( g_error * g_error ).sum(axis=0) )
    sign_eta = np.sign( rgr.scandata[peaksfile].eta_per_grain[ g.ind ] ).astype(int)
    # Got h,k,l,sign_eta
    #    pl.figure(2)
    #    pl.cla()
    #    pl.hist( mod_g_error, bins=np.linspace(0,0.01,100) )
    h,k,l = hkl_int

    observed_peaks = set( zip( h,k,l,sign_eta ) )
    print( len(observed_peaks) )

    # Lattice / symmetry come in here !
    uc = unitcell.unitcell( g.unitcell, rgr.parameterobj.get('cell_lattice_[P,A,B,C,I,F,R]') )
    uc.makerings( mod_g.max()*1.05 )

    wvln, wedg, chi = (rgr.parameterobj.get(p) for p in ('wavelength', 'wedge', 'chi'))
    
    for ds in uc.ringds[:nrings][::-1]:
        hkls = uc.ringhkls[ds]
        gcalc = np.dot( g.UB, np.transpose(hkls))
        tth, eta12, omega12 = transform.uncompute_g_vectors( gcalc, wvln, wedg, chi )
        m = 0
        hit = 0
        expected = 0
        for j,(h,k,l) in enumerate(hkls):
            for i in range(2):
                s = np.sign( eta12[i][j] ).astype(int)
                if (h,k,l,s) in observed_peaks:
                    hit += 1
                m += 1
                obin = np.round(( anorm( omega12[i][j] ) - obins[0])/osl).astype(int)
                if omega_histogram[ obin ] == 0 or tth[j] < 1e-3:
                    continue
                # tth/eta map?
                expected += 1

        if expected > 0:
            p = 100 * hit / expected
        else:
            p = 0

        print("( %4d %4d % 4d )  %4d    [ %4d / %4d ] = %6.2f %%"%(hkls[-1][0],
                                                                 hkls[-1][1],
                                                                 hkls[-1][2],
                                                                 m, hit, expected, p) )


        


# TODO:
#    Needs a space group symbol per grain
#    Needs the list of rings (or hkls) to be tested
#    Needs to check detector position for module gaps
#    Ideally taking account of error in position and/or error in intensity
