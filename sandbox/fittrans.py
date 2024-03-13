
from __future__ import print_function


from ImageD11.columnfile import columnfile
from ImageD11 import transform, parameters, unitcell, grain, gv_general
import numpy as np, time
from matplotlib import pylab


def k_to_g(k, omega, wedge, chi ):
    """ Rotate k vector to g vector """
    #print "k_to_g",k,omega,wedge,chi
    om = np.radians( omega )
    g = np.zeros( k.shape )
    cw = gv_general.chiwedge( chi, wedge )
    cwk = np.dot( cw, k )
    co = np.cos( np.radians( omega )) # 1D
    so = np.sin( np.radians( omega )) # 1D
    g[0,:] =  co*cwk[0,:] + so*cwk[1,:]
    g[1,:] = -so*cwk[0,:] + co*cwk[1,:]
    g[2,:] = cwk[2,:]
    return g

def g_to_k(g, omega, wedge, chi):
    """ Rotate g vector to k vector """
    co = np.cos(np.radians(omega))
    so = np.sin(np.radians(omega))
    k = np.zeros( (3,len(omega)))
    k[0,:] = co * g[0,:] - so * g[1,:]
    k[1,:] = so * g[0,:] + co * g[1,:]
    k[2,:] = g[2,:]
    cw = gv_general.chiwedge( chi, wedge )
    k = np.dot( cw, k )
    return k


def grain_x0( omega, wedge, chi, tx, ty, tz ):
    """ Grain origin positions at Omega rotation """
    co = np.cos(np.radians(omega))
    so = np.sin(np.radians(omega))
    x0 = np.zeros(( 3, omega.shape[0]) )
    x0[0] = co*tx - so*ty
    x0[1] = so*tx + co*ty
    x0[2] = tz
    cw = gv_general.chiwedge( chi, wedge )
    x0 = np.dot( cw.T, x0 )
    return x0


def XL_to_gv_old( omega, wedge, chi, XL, wavelength, tx, ty, tz):
    """ Historic way to compute gv from peak position """
    beam = np.array( [ 1/wavelength, 0, 0] )
    x0 = grain_x0( omega, wedge, chi, tx, ty, tz)
    s1 = XL - x0
    mods = np.sqrt( (s1*s1).sum(axis=0) ) * wavelength
    k = ((s1/mods).T - beam).T
    g = k_to_g( k, omega, wedge, chi )
    return g

def XL_to_gv( omega, wedge, chi, XL, wavelength, tx, ty, tz):
    """ Compute gv by rotation of experiment around crystal """
    beam = np.zeros( (3,len(omega)), float)
    beam[0,:] = 1.0
    XB = k_to_g( beam, omega, wedge, chi)
    #print "XB",XB
    XS = k_to_g( XL, omega, wedge, chi )
    #print "XS",XS
    np.subtract( XS[0], tx, XS[0] )
    np.subtract( XS[1], ty, XS[1] )
    np.subtract( XS[2], tz, XS[2] )
    modxs = np.sqrt( (XS*XS).sum(axis=0))
    g = (XS/modxs - XB)/wavelength
    return g

def test_XL_to_gv():
    """ Check the old and new way give the samy gv """
    omega = np.array([0,20,60,127,243])
    wedge = 0
    chi = 0
    XL = np.ones((15))
    XL.shape=3,5
    wavelength = 1.2345
    T = np.array( [0,0,0] )
    gnew = XL_to_gv( omega, wedge, chi, XL, wavelength, T[0],T[1],T[2])
    gold = XL_to_gv_old( omega, wedge, chi, XL, wavelength, T[0],T[1],T[2])
    print(gnew-gold)

#test_XL_to_gv()


def XSXB_to_gv( XS, XB, tx, ty, tz, wavelength):
    """ 
    XS - position of spot on detector
    XB - direction of X-ray beam
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
    xst = XS[0] - tx, XS[1] - ty, XS[2] - tz
    # dxst_dt = -1
    modxst2 = xst[0]*xst[0] + xst[1]*xst[1] + xst[2]*xst[2]
    # dmodxst2 = 2
    modxst = np.sqrt( modxst2 )
    # dmodxst_dt = dmodxst_dxst . dxst_dt
    # dmodxst2_dt = -2*xst[0] , -2*xst[1], -2*xst[2]
    # d(sqrt(f))/dx = 1/2 1/sqrt(f) . df/dx
    dmodxst_dt = -xst[0]/modxst, -xst[1]/modxst , -xst[2]/modxst
    # From XSXB_to_gv
    #    g = np.array((    (xst[0]/modxst - XB[0])/wavelength,\
    #                      (xst[1]/modxst - XB[1])/wavelength,\
    #                      (xst[2]/modxst - XB[2])/wavelength ) )
    g = (xst / modxst - XB ) /wavelength
    #
    # g[0] = (xst[0]/modxst - XB[0])/wavelength
    # g[1] = (xst[1]/modxst - XB[1])/wavelength
    # g[2] = (xst[2]/modxst - XB[2])/wavelength
    #
    # Want dg/dt
    #
    # dg/dt = (dxst/dt*1/modxst + xst.dmodxst_dt/modxst2)/w
    dg0_dtx = -(1./modxst + xst[0]*dmodxst_dt[0]/modxst2)/wavelength
    dg1_dtx = -(xst[0]*dmodxst_dt[1]/modxst2)/wavelength
    dg2_dtx = -(xst[0]*dmodxst_dt[2]/modxst2)/wavelength
    dg0_dty = -(xst[1]*dmodxst_dt[0]/modxst2)/wavelength
    dg1_dty = -(1./modxst + xst[1]*dmodxst_dt[1]/modxst2)/wavelength
    dg2_dty = -(xst[1]*dmodxst_dt[2]/modxst2)/wavelength
    dg0_dtz = -(xst[2]*dmodxst_dt[0]/modxst2)/wavelength
    dg1_dtz = -(xst[2]*dmodxst_dt[1] /modxst2)/wavelength
    dg2_dtz = -(1./modxst + xst[2]*dmodxst_dt[2]/modxst2)/wavelength
    #
           #   1       2       3       4       5       6       7       8       9
    return g,dg0_dtx,dg1_dtx,dg2_dtx,dg0_dty,dg1_dty,dg2_dty,dg0_dtz,dg1_dtz,dg2_dtz


class fittrans(object):

    def __init__(self, fname, parfile):
        self.colfile = columnfile( fname )
        self.pars = parameters.read_par_file( parfile )
        
        self.ds_tol = 0.005
        self.bins=np.arange(0,0.5,0.005)
        print("Setting up")
        self.compute_XLYLZL()
        self.computegv()

    def computegv(self):
        sign = float(self.pars.get('omegasign'))
        tx  = float(self.pars.get('t_x'))
        ty  = float(self.pars.get('t_y'))
        tz  = float(self.pars.get('t_z'))
        wvln= float(self.pars.get('wavelength'))
        self.gv = XSXB_to_gv( self.XS, self.XB, tx, ty, tz, wvln)

    def refine(self, gr, tol=0.05):
        tx,ty,tz = gr.translation
        wvln = float(self.pars.get("wavelength"))
        if hasattr(gr, "pks") and gr.pks is not None:
            #print "Got pks"
            pks = gr.pks
            XS = self.XS[:,pks]
            XB = self.XB[:,pks]
        else:
            #print "New pks"
            XS = self.XS
            XB = self.XB
            pks = np.arange( len(XS[0]) ) 
        ret = d_XSXB_to_gv( XS, XB, tx, ty, tz, wvln)
        gv = ret[0]
        dg0_dt = ret[1],ret[4],ret[7]
        dg1_dt = ret[2],ret[5],ret[8]
        dg2_dt = ret[3],ret[6],ret[9]
        hklr = np.dot( gr.ubi, gv )
        hkli = np.round(hklr)
        hkle = hklr - hkli
        scor = np.sqrt((hkle*hkle).sum(axis=0))
        #print "before",len(pks),pks
        if tol is not None:
            use = np.compress( scor < tol,  np.arange(len(gv[0]) ))
            #print "after",len(pks),pks
            gr.pks = pks[use]
        else:
            use = np.arange( len(gr.pks), dtype=int )
        #print "score = ", scor[pks].sum()/len(pks), len(pks), tol
        # peaks to use are those with scor OK
        #
        #  UB.h = gvcalc
        #  dg/dUB = h
        #  dg/dT  = found above
        gcalc = np.dot( gr.UB,  hkli )
        diff = np.take(gv - gcalc, use, axis=1)
        # print diff.shape, pks.shape
        # gv[0],[1],[2] = 3
        # tx, ty, ty    = 3
        # UB00 ... UB22 = 9
        # want derivative of diff w.r.t each variable
        grads = np.zeros( (12, 3, len(use)) )
        for i in range(3):
            for j in range(3): # tx, ty, tz
                #print 1+j*3+i,
                #print ret[1+j*3+i].shape
                grads[j, i] = ret[1+j*3+i][use]
           #     print grads[j,i]
            for j in range(3):
                # gx = 0h + 1k + 2l
                # gy = 3h + 4k + 5l
                # gz = 6h + 7k + 8l
                # i is gx, gy, gz
                # j is ub elements
                grads[3+j+i*3 , i] = hkli[j][use]
           #     print grads[3+j+i*3,i]
        # grains = 12, xyz, pks
        mtx = np.zeros( (12,12) )
        for i in range(12):
            for j in range(i,12):
                for k in range(3):
                    mtx[i,j] += (grads[i,k,:] * grads[j,k,:]).sum()
                if j!=i:
                    mtx[j,i] = mtx[i,j]

        #    mtx = np.dot( grads, grads.T) # vector, outer, ?
        rhs = np.zeros( 12 )
        for i in range(12):
            for k in range(3):
                rhs[i] += (grads[i,k]*diff[k]).sum()
        #print mtx
        # print rhs
        imt = np.linalg.inv( mtx )
        shifts = np.dot( imt, rhs )
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
                                                  **self.pars.parameters)
        # These are the peak positions in the laboratory frame
        #
        # Now get the beam and peaks in the crystal frame
        wedge = float(self.pars.get('wedge'))
        chi   = float(self.pars.get('chi'))
        omega = self.colfile.omega
        beam = np.zeros( (3,len(omega)))
        beam[0,:] = 1
        self.XB = k_to_g( beam, omega, wedge, chi)
        self.XS = k_to_g( self.peaks_xyz, omega, wedge, chi )



    def makecell(self):
        p = self.pars.parameters
        c =[ float(v) for v in [
                p['cell__a'],
                p['cell__b'],
                p['cell__c'],
                p['cell_alpha'],
                p['cell_beta'],
                p['cell_gamma'] ]]
        self.cell = unitcell.unitcell(c, p['cell_lattice_[P,A,B,C,I,F,R]'])
        self.cell.makerings( limit = self.ds.max(), tol = self.ds_tol )

    def ringassign(self, ds=None):
        if ds is None:
            ds = self.colfile.ds
        diffs = np.ones( self.colfile.nrows, np.float32 )
        ra    = np.ones( self.colfile.nrows, np.int ) - 2
        for i in range(len(self.cell.ringds)):
            testdiff = abs( ds - self.cell.ringds[i] )
            msk = testdiff < diffs
            diffs = np.where( msk, testdiff, diffs )
            ra = np.where( msk, i, ra )
        self.ra = ra
        self.diffs = diffs
        pylab.figure()
        vol = np.sqrt(np.linalg.det(self.cell.g))
        rscp = pow(vol, 1.0/3.0)
        pylab.hist( diffs*rscp, bins = np.arange(0,0.5,0.005) )
        print((diffs*rscp < 0.05).sum())
        return ra, diffs

    def find_triplets(self, peak):
        gx = self.colfile.gx - self.colfile.gx[peak]
        gy = self.colfile.gy - self.colfile.gy[peak]
        gz = self.colfile.gz - self.colfile.gz[peak]
        ds = np.sqrt( gx*gx + gy*gy + gz*gz )
        ra, diffs = self.ringassign( ds )


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
        for ii,tol in enumerate(tols):
            #print gref.translation
            gref = o.refine(gref, tol=tol)
            #print ii, gref.translation, gref.npks,
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
