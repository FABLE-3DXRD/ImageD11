
from __future__ import print_function


from ImageD11.columnfile import columnfile
from ImageD11 import transform, parameters, unitcell, grain, gv_general
import numpy as np, time
from matplotlib import pylab
import numba

def k_to_g(k, cosomega, sinomega, chiwedge, out = None ):
    """ Rotate k vector to g vector, e.g. from lab to crystal """
    if out is None:
        g = np.empty_like( k )
    else:
        g = out 
    assert (g.shape == k.shape) and (g.dtype == float)
    cwk = np.dot( chiwedge, k )
    g[0,:] =  cosomega*cwk[0,:] + sinomega*cwk[1,:]
    g[1,:] = -sinomega*cwk[0,:] + cosomega*cwk[1,:]
    g[2,:] = cwk[2,:]
    return g

def XSXB_to_gv( XS, XB, tx, ty, tz, wavelength):
    """ 
    XS - position of spot on detector (rotated)
    XB - direction of X-ray beam (rotated)
    tx, ty, tz - grain position
    wavelength
    Computes g-vector
    """
    xst = (XS.T - [tx,ty,tz]).T
    modxst = np.sqrt( (xst*xst).sum(axis=0) )
    g = (xst/modxst - XB)/wavelength
    return g

# derivatives:

def d_XSXB_to_gv( XS, XB, tx, ty, tz, wavelength):
    """
    As XSXB_to_gv but with derivative on tx, ty, tz
    """
    dgdt = np.zeros( (3, 3, XS.shape[1]), dtype=float )
    rscale = 1/wavelength
    # Xspot - X0 = scattered output rays
    xst = XS.copy()
    xst[0] -= tx
    xst[1] -= ty
    xst[2] -= tz
    #    xst = XS - np.array((tx, ty, tz))[:,np.newaxis]
    # dxst_dt = -1  # derivatice w.r.t translation
    # reciprocal length of vector squared
    rmodxst2 = 1/(xst[0]*xst[0] + xst[1]*xst[1] + xst[2]*xst[2])
    # dmodxst2 = 2
    rmodxst = np.sqrt( rmodxst2 )
    # dmodxst_dt = dmodxst_dxst . dxst_dt
    # dmodxst2_dt = -2*xst[0] , -2*xst[1], -2*xst[2]
    # d(sqrt(f))/dx = 1/2 1/sqrt(f) . df/dx
    dmodxst_dt = -xst*rmodxst
    # From XSXB_to_gv
    #    g = np.array((    (xst[0]/modxst - XB[0])/wavelength,\
    #                      (xst[1]/modxst - XB[1])/wavelength,\
    #                      (xst[2]/modxst - XB[2])/wavelength ) )
    g = rscale * (xst * rmodxst - XB )
    # g[0] = (xst[0]/modxst - XB[0])/wavelength
    # g[1] = (xst[1]/modxst - XB[1])/wavelength
    # g[2] = (xst[2]/modxst - XB[2])/wavelength
    # Want dg/dt
    # dg/dt = (dxst/dt*1/modxst + xst.dmodxst_dt/modxst2)/w
    for i in range(3):
        for j in range(3):
            if i==j:
                dgdt[i,j] = -rmodxst-xst[i]*dmodxst_dt[j]*rmodxst2
            else:
                dgdt[i,j] = -xst[i]*dmodxst_dt[j]*rmodxst2
    np.multiply( dgdt, rscale, dgdt )
    return g, dgdt

# TODO:
#    test this derivative to check whether it is correct (seems OK on simul_1000_grains)
#    make the 'units' of position match the 'units' of UB (GSAS normalising diagonal)
#    error bars, sh/esd for convergence
#    weights on the data... set up 3x3 weight matrices
#    impact of using omega obs?
#    



class fittrans(object):
    """ Fits grain translations """
    def __init__(self, fname, parfile):
        print("Reading", fname)
        self.colfile = columnfile( fname )
        self.colfile.parameters.loadparameters( parfile )
        print("Setting up")
        self.compute_XLYLZL()

    def refine(self, gr, tol = 0.05):
        tx,ty,tz = gr.translation
        wvln = float(self.colfile.parameters.get("wavelength"))
        if hasattr(gr, "pks") and gr.pks is not None:
            #print "Got pks"
            pks = gr.pks
            XS = self.XS[:,pks]  # 3xn
            XB = self.XB[:,pks]  # 3xn
        else:
            #print "New pks"
            XS = self.XS
            XB = self.XB
            pks = np.arange( len(XS[0]) )
        npks = len(pks)
        # Compute the g-vectors and derivatives with respect to translations
        gv, dgdt = d_XSXB_to_gv( XS, XB, tx, ty, tz, wvln)
        hklr = np.dot( gr.ubi, gv )
        hkli = np.round(hklr)
        hkle = hklr - hkli
        scor = (hkle*hkle).sum(axis=0)
        #print "before",len(pks),pks        
        if tol is not None:
            msk = scor < tol*tol
            wt = np.where( msk, 1., 0.)
            use = np.compress( scor < tol*tol,  np.arange(npks) )
            #print "after",len(pks),pks
            gr.pks = pks[use]
        else:
            wt = 1
            use = np.arange( 0, len(gv[0]), dtype=int)
        #print "score = ", scor[pks].sum()/len(pks), len(pks), tol
        # peaks to use are those with scor OK
        #
        #  UB.h = gvcalc
        #  dg/dUB = h
        #  dg/dT  = found above
        gcalc = np.dot( gr.UB,  hkli ) # 3x3, 3xn
        diff = np.take(gv - gcalc, use, axis=1) # 3xn
        # print diff.shape, pks.shape
        # gv[0],[1],[2] = 3
        # tx, ty, ty    = 3
        # UB00 ... UB22 = 9
        # want derivative of diff w.r.t each variable
        grads = np.zeros( (12, 3, npks ) )   # 12 vars, 3xn
        #for j in range(3): # tx, ty, tz
        #    for i in range(3):                #print 1+j*3+i,
        #        #print ret[1+j*3+i].shape
        #        grads[j, i] = dgdt[j*3+i][use]
        grads[:3, :] = dgdt[:,:]
        
        for i in range(3):
            for j in range(3):
                # gx = 0h + 1k + 2l
                # gy = 3h + 4k + 5l
                # gz = 6h + 7k + 8l
                # i is gx, gy, gz
                # j is ub elements
                grads[3+j+i*3,i] = hkli[j]
                #     print grads[3+j+i*3,i]
        #       ... there are a lot of zero here. 
        #       ... dgx / dub[yz] == 0 ... of course.
        #
        # 
        # future: put some weights in. Not necessarily diagonal
        #
        # M.A.M.s = M.A.d
        # M = gradients
        # A = weight matrix
        # s = shifts
        # d = differences
        # 
        grd = grads.reshape( (12, -1) )
        # (12,n),(n,12) -> 12,12
        # weights will be 3x3 blocked per peak
        # wmat = np.dot( grd, wt, grd.T )
        # wrhs = np.dot( grd, wt, diff  )
        mat  = np.dot( grd, grd.T )
        rhs  = np.dot( grd, diff.ravel() )
        shifts, residuals, rank, s = np.linalg.lstsq( mat, rhs )
        tx = tx - shifts[0]
        ty = ty - shifts[1]
        tz = tz - shifts[2]
        gr.translation = [tx,ty,tz]
        ub = gr.UB.ravel()
        np.add(ub, shifts[3:], ub)
        gr.set_ubi( np.linalg.inv( np.reshape(ub,(3,3))))
        gr.npks = len(use)
        #1/0
        return gr

    def compute_XLYLZL(self):
        self.peaks_xyz = transform.compute_xyz_lab( [self.colfile.sc,
                                                    self.colfile.fc],
                                                  **self.colfile.parameters.parameters)
        # These are the peak positions in the laboratory frame
        #
        # Now get the beam and peaks in the crystal frame
        wedge = float(self.colfile.parameters.get('wedge'))
        chi   = float(self.colfile.parameters.get('chi'))
        self.chiwedge = gv_general.chiwedge( chi, wedge )
        omega = np.radians( self.colfile.omega )
        cosomega = np.cos( omega )
        sinomega = np.sin( omega )
        beam = np.zeros( (3,len(omega)), float)
        beam[0,:] = 1.
        self.XB = k_to_g( beam, cosomega, sinomega, self.chiwedge)
        self.XS = k_to_g( self.peaks_xyz, cosomega, sinomega, self.chiwedge )



if __name__=="__main__":
    import sys
    from ImageD11.indexing import ubitocellpars
    o = fittrans( sys.argv[1], sys.argv[2] )
    gl = grain.read_grain_file( sys.argv[3] )
    gref = gl[0]
    import time
    start = time.time()
    gfl = []
    ng = 0
    # Take old peak assignments:
    if 1:
        print("Using existing peak assignments")
        inds = np.arange(o.colfile.nrows,dtype=int)
        for i,gref in enumerate(gl):
            gref.pks = np.compress( o.colfile.labels == i, inds )
        tols = [None,]*3
    else:
        tols = [0.05,0.01,0.0075]
    for gref in gl:
        print("BEFORE:            ",
              3*"% 7.3f "%tuple(gref.translation), 6*"%.6f "%(ubitocellpars(gref.ubi)))
        for ii,tol in enumerate(tols):
            #print gref.translation
            gref = o.refine(gref, tol=tol)
            #print i,gref.npks
#        gref.pks = None
        # re-assign after convergence
#        gref = o.refine( gref, tol=0.0075) 
        gfl.append( gref )
        ng += 1
        print(ng,ng*100.0/len(gl), "%.2f   "%(time.time()-start), gref.npks, end=' ')
        print((3*"%.4f ")%tuple(gref.translation), end=' ')
        print((6*"%.6f ")%ubitocellpars(gref.ubi))
    print(time.time()-start)
    grain.write_grain_file(sys.argv[4], gfl)
#    1/0
#    for i in range(2):
#        o.find_triplets( i*17 )
#    pylab.show()

# ~/ImageD11/test/simul_1000_grains % python ../../sandbox/fittrans.py Al1000/Al1000.flt Al1000/Al1000.par allgrid.map allgridfittrans.map
