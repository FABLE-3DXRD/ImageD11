
"""WARNING: work in progress"""


# given a dataset + parameters + unit cell (+ rings) + sparsepixels : create polefigues
import numpy as np
import pylab as pl
import os
import ImageD11.sinograms.dataset
import ImageD11.sinograms.geometry
import ImageD11.sparseframe
import ImageD11.transform
import ImageD11.unitcell
import ImageD11.parameters


def ring_hkls(ucell, dsmax, n_rings=None):
    """
    List the powder rings of ucell out to dsmax (in 1/d), with their
    multiplicities.

    hkls come back ring-by-ring, so a per-ring block of hkls can be
    recovered from mults (e.g. with np.cumsum(mults)). Pass n_rings to
    keep only the first n_rings low-angle rings (dsmax must still be
    large enough to reach them).
    """
    ucell.makerings(dsmax)
    dstars = ucell.ringds if n_rings is None else ucell.ringds[:n_rings]
    hkls, mults = [], []
    for dstar in dstars:
        hklring = ucell.ringhkls[dstar]
        mults.append(len(hklring))
        hkls.extend(hklring)
    return np.array(hkls), mults


def calc_gve_angles(pmap, hkls, pars):
    """
    For every indexed voxel of a point-by-point map (pmap), use its best
    UBI matrix to predict where each of hkls should diffract to.

    Returns voxel_mask, selecting the indexed voxels out of pmap's 2D
    voxel grid, together with the calculated (tth, eta, omega). Each hkl
    has two solutions (a Friedel pair), so eta and omega each come back
    as a pair of arrays (eta1, eta2) / (omega1, omega2), shaped
    (n_indexed_voxels, len(hkls)).
    """
    voxel_mask = ~np.isnan(pmap.best_ubi[:, :, 0, 0])
    ub = np.linalg.inv(pmap.best_ubi[voxel_mask])
    # gve[xyz, voxel, hkl]
    gve = np.einsum('vij,hj->ivh', ub, hkls, order='C')
    gve = gve.reshape(3, -1, copy=False)
    tth, (eta1, eta2), (omega1, omega2) = ImageD11.transform.uncompute_g_vectors(
        gve, pars.get('wavelength')
    )
    shape = ub.shape[0], len(hkls)
    tth = tth.reshape(shape, copy=False)
    eta1 = eta1.reshape(shape, copy=False)
    eta2 = eta2.reshape(shape, copy=False)
    omega1 = omega1.reshape(shape, copy=False)
    omega2 = omega2.reshape(shape, copy=False)
    return voxel_mask, tth, eta1, eta2, omega1, omega2


def calc_dty(ds, voxel_mask, omega1, omega2, y0):
    """
    Convert the calculated omega values into the dty motor position that
    would have put each voxel on the beam at that omega, so they can be
    overlaid on the measured (omega, dty) sinogram.
    """
    N = voxel_mask.shape[0]
    si, sj = np.mgrid[0:N, 0:N]
    si = si[voxel_mask] - N // 2
    sj = sj[voxel_mask] - N // 2
    sx, sy = ImageD11.sinograms.geometry.step_to_sample(si, sj, ystep=ds.ystep)
    dty1 = ImageD11.sinograms.geometry.dtycalc(omega1, sx[:, np.newaxis], sy[:, np.newaxis], y0)
    dty2 = ImageD11.sinograms.geometry.dtycalc(omega2, sx[:, np.newaxis], sy[:, np.newaxis], y0)
    return dty1, dty2


def plot_ring_polefigures(ref_ucell, mults, eta1, eta2, omega1, omega2, cf_2d, ds_tol=0.005, bin_deg=1.0):
    """
    For each ring described by mults, compare the (eta, omega) pole figure
    forward-simulated from the point-by-point map (every indexed voxel's
    best UBI, projected through hkl) against the measured 2D peaks lying
    close to that ring's d-spacing.

    bin_deg sets the (eta, omega) histogram bin width in degrees - a finer
    bin_deg gives sharper pole figure spots but needs more indexed voxels
    per spot to show up. dataset.ostep (the scan's angular step) is a
    natural alternative choice.
    """
    start = 0
    ebins = np.arange(-180, 180, bin_deg)
    obins = np.arange(-90, 90, bin_deg)
    for ring, m in enumerate(mults):
        end = start + m
        # both Friedel-pair solutions go in the same histogram - two separate
        # hist2d calls on one axes would paint over each other
        eta_calc = np.concatenate((eta1[:, start:end].ravel(), eta2[:, start:end].ravel()))
        omega_calc = np.concatenate((omega1[:, start:end].ravel(), omega2[:, start:end].ravel()))
        f, a = pl.subplots(2, 1, constrained_layout=True, figsize=(9, 7))
        a[0].hist2d(eta_calc, omega_calc, bins=(ebins, obins), norm='log')
        title = 'Ring %d %s' % (ring, str(ref_ucell.ringhkls[ref_ucell.ringds[ring]][-1]))
        a[0].set(xlabel='eta', ylabel='omega', title='Calc: ' + title)
        cfr = cf_2d.copyrows(abs(cf_2d.ds - ref_ucell.ringds[ring]) < ds_tol)
        a[1].hist2d(cfr.eta, cfr.omega, bins=(ebins, obins), weights=cfr.sum_intensity, norm='log')
        a[1].set(xlabel='eta', ylabel='omega', title='Obs: ' + title)
        start = end


def plot_ring_sinograms(ds, cf_2d, ref_ucell, mults, omega1, omega2, dty1, dty2,
                         ds_tol=0.005, omega_step=5):
    """
    For each ring described by mults, compare the (omega, dty) sinogram of
    the calculated peaks (from the point-by-point map) against the
    measured 2D peaks lying close to that ring's d-spacing.
    """
    start = 0
    b = ds.obinedges[::omega_step]
    yb = ds.ybinedges
    for ring, m in enumerate(mults):
        end = start + m
        # both Friedel-pair solutions go in the same histogram - two separate
        # hist2d calls on one axes would paint over each other
        omega_calc = np.concatenate((omega1[:, start:end].ravel(), omega2[:, start:end].ravel()))
        dty_calc = np.concatenate((dty1[:, start:end].ravel(), dty2[:, start:end].ravel()))
        f, a = pl.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
        a[0].hist2d(omega_calc, dty_calc, bins=(b, yb), norm='log')
        title = 'Ring %d %s' % (ring, str(ref_ucell.ringhkls[ref_ucell.ringds[ring]][-1]))
        a[0].set(xlabel='Omega', ylabel='dty', title='Calc: ' + title)
        cfr = cf_2d.copyrows(abs(cf_2d.ds - ref_ucell.ringds[ring]) < ds_tol)
        a[1].hist2d(cfr.omega, cfr.dty, bins=(b, yb), norm='log')
        a[1].set(xlabel='Omega', ylabel='dty', title='Obs: ' + title)
        start = end


def main( pars,
          dataset,
          
        ):
    pass


class splatter:
    """ splats pixels onto polefigures """
    def __init__(self, 
                 parfile,
                 dxfile = '/data/id11/nanoscope/Eiger/spatial_20210415_JW/e2dx.edf',
                 dyfile = '/data/id11/nanoscope/Eiger/spatial_20210415_JW/e2dy.edf'):
        self.pardict = ImageD11.parameters.parameters(parfile).parameters
        self.pardict['dxfile'] = dxfile
        self.pardict['dyfile'] = dyfile
        self.pLUT = ImageD11.transform.pixelLUT( self.pardict )
        
    def process(self, dsfile):
        ds = ImageD11.sinograms.dataset( dsfile )
        


    

if 0:
    import numba
    ds = ImageD11.sinograms.dataset.load('/data/visitor/blc14570/id11/20230425/PROCESSED_DATA/ds_Martensite500C_DTz50.h5')
    ds.omega.shape, len(ds.sparsefiles)
    chosen = ds.omega.shape[0]//2 + 1
    sps = ImageD11.sparseframe.SparseScan( os.path.join(ds.analysispath, ds.sparsefiles[chosen]), ds.limapath )
    pars = ImageD11.parameters.read_par_file( '/data/visitor/blc14570/id11/20230425/PROCESSED_DATA/FeBCC.par' ).parameters
    pars['dxfile'] = '/data/id11/nanoscope/Eiger/spatial_20210415_JW/e2dx.edf'
    pars['dyfile'] = '/data/id11/nanoscope/Eiger/spatial_20210415_JW/e2dy.edf'
    plut = ImageD11.transform.PixelLUT( pars )
    ds = 2 * np.sin( np.radians( plut.tth / 2 ) ) / pars['wavelength']
    a0 = 2.86
    uc = ImageD11.unitcell.unitcell( [a0,a0,a0,90,90,90,],'I')
    uc.makerings( ds.max() )
    labels = np.zeros( ds.shape, dtype=np.uint32 )
    for i, dsr in enumerate( uc.ringds ):
        tol = 4e-3 * ds + 1e-2
        m = abs( ds - dsr ) < tol
        labels[m] = i + 1
    # eta values
    etabin = np.round( ( plut.eta % 360 ) / 0.1 ).astype(int)
    etabin.min(), etabin.max()
    output_shape = (labels.max()+1, sps.shape[0], 3601)
    pfs = np.zeros( output_shape, int)
    ipf = labels[ sps.row, sps.col ]
    ifrm = np.zeros( len(sps.row), int )
    for i in range( len(sps.ipt) - 1):
        ifrm[ sps.ipt[i]: sps.ipt[i+1] ] = i
    ieta = etabin[ sps.row, sps.col ]        
    tth_step = plut.tth.max() / 3000
    @numba.njit
    def accumulate( ipf, ifrm, ieta, intensity, output ):
        for j in range(ipf.size):
            output[ ipf.flat[j], ifrm.flat[j], ieta.flat[j] ] += intensity[j]
    accumulate( ipf, ifrm, ieta, sps.intensity, pfs )
    
    tthpf = np.zeros( output_shape, int)
    tthbin = np.round( plut.tth / tth_step ).astype(int)
    assert tthbin.max() <= 3001
    itth = tthbin[ sps.row, sps.col ]
    accumulate( ipf, ifrm, itth, sps.intensity, tthpf )
    
    f,a = pl.subplots(3,7, figsize=(21,7))
    for ax, p, j in zip(a.ravel(), pfs[1:], range(len(pfs[1:]))):
        vmx = max( p.max(), 1 )
        ax.imshow(p,  norm=pl.matplotlib.colors.LogNorm(vmin=0.1, vmax=vmx))
        hkls = uc.ringhkls[ uc.ringds[j] ]
        title=str(hkls[-1]) + "  M=%d"%(len(hkls))
        ax.set( xlabel= 'eta', ylabel='omega', xticks = [], yticks = [], title=title )
        
    