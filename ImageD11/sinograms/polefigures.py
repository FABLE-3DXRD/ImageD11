
"""WARNING: work in progress"""


# given a dataset + parameters + unit cell (+ rings) + sparsepixels : create polefigues
import numpy as np
import pylab as pl
import os
import ImageD11.sinograms.dataset
import ImageD11.sparseframe
import ImageD11.transform
import ImageD11.unitcell
import ImageD11.parameters


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
        
    