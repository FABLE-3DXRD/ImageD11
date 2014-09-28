

from ImageD11.columnfile import columnfile
from ImageD11 import transform, parameters, unitcell, grain
import numpy as np, time
from matplotlib import pylab

def wedgematrix( wedge ):
    cw = np.cos(np.radians(wedge))
    sw = np.sin(np.radians(wedge))
    w = np.array( [ [ cw, 0, sw], [0, 1, 0], [-sw, 0, cw]] )
    return w

def chimatrix( chi ):
    cc = np.cos(np.radians(chi))
    sc = np.sin(np.radians(chi))
    c = np.array( [ [ 1, 0, 0], [ 0, cc, sc],[ 0, -sc, cc]] )
    return c


def k_to_g(k, omega, wedge, chi ):
    #print "k_to_g",k,omega,wedge,chi
    om = np.radians( omega)
    g = np.zeros( k.shape )
    cw = np.dot( chimatrix( chi ), wedgematrix( wedge ) )
    cwk = np.dot( cw, k ) 
    co = np.cos( np.radians( omega )) # 1D
    so = np.sin( np.radians( omega )) # 1D
    g[0,:] =  co*cwk[0,:] + so*cwk[1,:]
    g[1,:] = -so*cwk[0,:] + co*cwk[1,:]
    g[2,:] = cwk[2,:]
    return g

def g_to_k(g, omega, wedge, chi):
    co = np.cos(np.radians(omega))
    so = np.sin(np.radians(omega))
    k = np.zeros( (3,len(omega)))
    k[0,:] = co * g[0,:] - so * g[1,:]
    k[1,:] = so * g[0,:] + co * g[1,:]
    k[2,:] = g[2,:]
    cw = np.dot( chimatrix( chi ), wedgematrix( wedge ) )
    k = np.dot( cw, k )
    return k


def grain_x0( omega, wedge, chi, tx, ty, tz ):
    co = np.cos(np.radians(omega))
    so = np.sin(np.radians(omega))
    x0 = np.zeros(( 3, omega.shape[0]) )
    x0[0] = co*tx - so*ty
    x0[1] = so*tx + co*ty
    x0[2] = tz
    cw = np.dot( chimatrix( chi ), wedgematrix( wedge ) )
    x0 = np.dot( cw.T, x0 )
    return x0


def XL_to_gv_old( omega, wedge, chi, XL, wavelength, tx, ty, tz):
    beam = np.array( [ 1/wavelength, 0, 0] )
    x0 = grain_x0( omega, wedge, chi, tx, ty, tz)
    s1 = XL - x0
    mods = np.sqrt( (s1*s1).sum(axis=0) ) * wavelength
    k = ((s1/mods).T - beam).T
    g = k_to_g( k, omega, wedge, chi )
    return g
    
def XL_to_gv( omega, wedge, chi, XL, wavelength, tx, ty, tz):
    beam = np.zeros( (3,len(omega)))
    beam[0,:] = 1
    XB = k_to_g( beam, omega, wedge, chi)
    #print "XB",XB
    XS = k_to_g( XL, omega, wedge, chi )
    #print "XS",XS
    np.subtract(XS[0], tx, XS[0] )
    np.subtract(XS[1], ty, XS[1] )
    np.subtract(XS[2], tz, XS[2] )
    modxs = np.sqrt( (XS*XS).sum(axis=0))
    g = (XS/modxs - XB)/wavelength
    return g

def test_XL_to_gv():
    omega = np.array([0,20,60,127,243])
    wedge = 0
    chi = 0
    XL = np.ones((15))
    XL.shape=3,5
    wavelength = 1.2345
    T = np.array( [0,0,0] )
    gnew = XL_to_gv( omega, wedge, chi, XL, wavelength, T[0],T[1],T[2])
    gold = XL_to_gv_old( omega, wedge, chi, XL, wavelength, T[0],T[1],T[2])
    print gnew-gold

#test_XL_to_gv()


def XSXB_to_gv( XS, XB, tx, ty, tz, wavelength):
    xst = (XS.T - [tx,ty,tz]).T
    modxst = np.sqrt( (xst*xst).sum(axis=0) )
    g = (xst/modxst - XB)/wavelength
    return g

# derivatives:
def d_XSXB_to_gv( XS, XB, tx, ty, tz, wavelength):
    debug = False
    if debug:
        print "in deriv, txtytz",tx,ty,tz
        print XS.shape, XB.shape
    xst = XS[0] - tx, XS[1] - ty, XS[2] - tz
    # dxst_dt = -1
    modxst2 = xst[0]*xst[0] + xst[1]*xst[1] + xst[2]*xst[2]
    # dmodxst2 = 2
    modxst = np.sqrt( modxst2 )
    # dmodxst_dt = dmodxst_dxst . dxst_dt
    # dmodxst2_dt = -2*xst[0] , -2*xst[1], -2*xst[2]
    # d(sqrt(f))/dx = 1/2 1/sqrt(f) . df/dx
    dmodxst_dt = -xst[0]/modxst, -xst[1]/modxst , -xst[2]/modxst  
    if 0: # OK
        # [f(x+h) - f(x)] / h
        t = (XS.T - [tx+1,ty,tz]).T
        dtx = np.sqrt( ( t * t ).sum(axis=0) ) -  modxst
        t = (XS.T - [tx,ty+1,tz]).T
        dty = np.sqrt( ( t * t ).sum(axis=0) ) - modxst
        t = (XS.T - [tx,ty,tz+1]).T
        dtz = np.sqrt( ( t * t ).sum(axis=0) ) - modxst
        print dtx[:3]
        print dty[:3]
        print dtz[:3]
        pylab.plot(dtx, dmodxst_dt[0],"o")
        pylab.plot(dty, dmodxst_dt[1],"o")
        pylab.plot(dtz, dmodxst_dt[2],"o")
        pylab.show()
        1/0
    #
    g =     (xst[0]/modxst - XB[0])/wavelength,\
            (xst[1]/modxst - XB[1])/wavelength,\
            (xst[2]/modxst - XB[2])/wavelength
    #
    # g[0] = (xst[0]/modxst - XB[0])/wavelength
    # g[1] = (xst[1]/modxst - XB[1])/wavelength
    # g[2] = (xst[2]/modxst - XB[2])/wavelength
    #
    # Want dg/dt
    # 
    # dg/dt = (dxst/dt*1/modxst + xst.dmodxst_dt/modxst2)/w
    dg0_dtx = -(1/modxst + xst[0]*dmodxst_dt[0]/modxst2)/wavelength 
    dg1_dtx = -(xst[0]*dmodxst_dt[1]/modxst2)/wavelength 
    dg2_dtx = -(xst[0]*dmodxst_dt[2]/modxst2)/wavelength 
    dg0_dty = -(xst[1]*dmodxst_dt[0]/modxst2)/wavelength 
    dg1_dty = -(1/modxst + xst[1]*dmodxst_dt[1]/modxst2)/wavelength 
    dg2_dty = -(xst[1]*dmodxst_dt[2]/modxst2)/wavelength 
    dg0_dtz = -(xst[2]*dmodxst_dt[0]/modxst2)/wavelength 
    dg1_dtz = -(xst[2]*dmodxst_dt[1] /modxst2)/wavelength 
    dg2_dtz = -(1/modxst + xst[2]*dmodxst_dt[2]/modxst2)/wavelength 
    # 
    if 0:
        delta=1
        t = (XS.T - [tx+delta,ty,tz]).T
        gt = ( t/np.sqrt( (t*t).sum(axis=0)) - XB) / wavelength
        numeric = (gt - g)/delta
        pylab.subplot(331)
        pylab.plot(numeric[0], dg0_dtx,"o" )
        pylab.subplot(332)
        pylab.plot(numeric[1], dg1_dtx,"o" )
        pylab.subplot(333)
        pylab.plot(numeric[2], dg2_dtx,"o" )
        delta=1
        t = (XS.T - [tx,ty+delta,tz]).T
        gt = ( t/np.sqrt( (t*t).sum(axis=0)) - XB) / wavelength
        numeric = (gt - g)/delta
        pylab.subplot(334)
        pylab.plot(numeric[0], dg0_dty,"o" )
        pylab.subplot(335)
        pylab.plot(numeric[1], dg1_dty,"o" )
        pylab.subplot(336)
        pylab.plot(numeric[2], dg2_dty,"o" )
        delta=1
        t = (XS.T - [tx,ty,tz+delta]).T
        gt = ( t/np.sqrt( (t*t).sum(axis=0)) - XB) / wavelength
        numeric = (gt - g)/delta
        pylab.subplot(337)
        pylab.plot(numeric[0], dg0_dtz,"o" )
        pylab.subplot(338)
        pylab.plot(numeric[1], dg1_dtz,"o" )
        pylab.subplot(339)
        pylab.plot(numeric[2], dg2_dtz,"o" )
        pylab.show()
        1/0
           #   1       2       3       4       5       6       7       8       9
    return g,dg0_dtx,dg1_dtx,dg2_dtx,dg0_dty,dg1_dty,dg2_dty,dg0_dtz,dg1_dtz,dg2_dtz


class newindexer(object):

    def __init__(self, fname, parfile):
        self.colfile = columnfile( fname ) 
        self.pars = parameters.read_par_file( parfile )

        self.ds_tol = 0.005
        self.bins=np.arange(0,0.5,0.005)

        self.compute_XLYLZL()
        self.computegv()
        self.makecell()
        self.ringassign()

    def computegv(self):
        sign = float(self.pars.get('omegasign'))
        tx  = float(self.pars.get('t_x'))
        ty  = float(self.pars.get('t_y'))
        tz  = float(self.pars.get('t_z'))
        wvln= float(self.pars.get('wavelength')) 

        start = time.clock()
        self.gv = XSXB_to_gv( self.XS, self.XB, tx, ty, tz, wvln)
        #print "time for gv only",time.clock()-start 
        #start = time.clock()
        ret = d_XSXB_to_gv( self.XS, self.XB, tx, ty, tz, wvln)
        #print "time for derivs too",time.clock()-start 
        self.deriv_g_t = ret[1:]
        if 0:
            tth, eta = transform.compute_tth_eta_from_xyz(
                             self.peaks_xyz, self.colfile.omega*sign,
                             **self.pars.parameters)
            gvtest  = transform.compute_g_vectors( tth, eta, self.colfile.omega*sign,
                        float(self.pars.get('wavelength')),
                        float(self.pars.get('wedge')),
                        float(self.pars.get('chi')))
            assert (abs(self.gv - gvtest)).max() < 1e-12, "gvector error"

            ret = d_XSXB_to_gv( self.XS, self.XB, tx, ty, tz, wvln)
            g = ret[0]
            dg0_dt = ret[1],ret[4],ret[7]
            dg1_dt = ret[2],ret[5],ret[8]
            dg2_dt = ret[3],ret[6],ret[9]
            assert (abs(g - gvtest)).max() < 1e-12, "gvector deriv error"
            delta = 1.0
            dgtx = XSXB_to_gv( self.XS, self.XB, tx+delta, ty, tz, wvln)-self.gv
            dgty = XSXB_to_gv( self.XS, self.XB, tx, ty+delta, tz, wvln)-self.gv
            dgtz = XSXB_to_gv( self.XS, self.XB, tx, ty, tz+delta, wvln)-self.gv
            def pe(a,b): 
                return abs(a-b)/abs(a+b)
            print abs((dg0_dt[0] - dgtx[0])).max()
            print abs((dg0_dt[1] - dgty[0])).max()
            print abs((dg0_dt[2] - dgtz[0])).max()
            print abs((dg1_dt[0] - dgtx[1])).max()
            print abs((dg1_dt[1] - dgty[1])).max()
            print abs((dg1_dt[2] - dgtz[1])).max()
            print abs((dg2_dt[0] - dgtx[2])).max()
            print abs((dg2_dt[1] - dgty[2])).max()
            print abs((dg2_dt[2] - dgtz[2])).max()

        self.ds = np.sqrt( (self.gv*self.gv).sum(axis=0) )
        self.colfile.addcolumn( self.gv[0], 'gx')
        self.colfile.addcolumn( self.gv[1], 'gy')
        self.colfile.addcolumn( self.gv[2], 'gz')
        self.colfile.addcolumn( self.ds, 'ds')

    def refine(self, gr, tol=0.05):
        tx,ty,tz = gr.translation
        wvln = float(self.pars.get("wavelength"))
        ret = d_XSXB_to_gv( self.XS, self.XB, tx, ty, tz, wvln)
        gv = ret[0]
        dg0_dt = ret[1],ret[4],ret[7]
        dg1_dt = ret[2],ret[5],ret[8]
        dg2_dt = ret[3],ret[6],ret[9]
        gv = XSXB_to_gv( self.XS, self.XB, tx, ty, tz, wvln )
        hklr = np.dot( gr.ubi, gv )
        hkli = np.round(hklr)
        hkle = hklr - hkli
        scor = np.sqrt((hkle*hkle).sum(axis=0))
        pks = np.compress( scor < tol,  np.arange( len(gv[0]) ) )
        print "score = ", np.take(scor, pks).sum()/len(pks), len(pks)
        if 0:
            pylab.hist(scor, bins=128)
            pylab.show()
        # peaks to use are those with scor OK
        #
        #    From Paciorek et al Acta A55 543 (1999)
        #  UB = R H-1
        #    R = sum_n r_n h_n^t
        #    H = sum_n h_n h_n^t
        #
        #  UB.h = gvcalc
        #  dg/dUB = h
        #  dg/dT  = found above
        gcalc = np.dot( gr.ub,  hkli )
        diff = np.take(gv - gcalc, pks, axis=1)
        # print diff.shape, pks.shape
        # gv[0],[1],[2] = 3
        # tx, ty, ty    = 3
        # UB00 ... UB22 = 9
        # want derivative of diff w.r.t each variable
        grads = np.empty( (12, 3, len(pks )) )
        for i in range(3):
            for j in range(3): # tx, ty, tz
                #print 1+j*3+i,
                #print ret[1+j*3+i].shape
                grads[j, i] = np.take( ret[1+j*3+i], pks )
            for j in range(3):
                # gx = 0h + 1k + 2l 
                # gy = 3h + 4k + 5l
                # gz = 6h + 7k + 8l
                # i is gx, gy, gz
                # j is ub elements
                #print j%3
                grads[3+j+i*3 , i] = np.take(hkli[j], pks)
        #       pylab.imshow(grads, aspect='auto' , interpolation='nearest')
        #       pylab.show()
        # grains = 12, xyz, pks
        mtx = np.zeros( (12,12) )
        for i in range(12):
            for j in range(i,12):
                for k in range(3):
                    mtx[i,j] += (grads[i,k,:] * grads[j,k,:]).sum()
                if j!=i:
                    mtx[j,i]=mtx[i,j]
    #    mtx = np.dot( grads, grads.T) # vector, outer, ?
        rhs = np.zeros( 12 )
        for i in range(12):
            for k in range(3):
                rhs[i] += (grads[i,k]*diff[k]).sum()
        #print rhs.shape
        #print mtx.shape
        #print mtx
        #pylab.clf()
        #pylab.imshow(np.log(np.abs(mtx)+1e-5),aspect='auto', interpolation='nearest' )
        #pylab.show()
        imt = np.linalg.inv( mtx )
        shifts = np.dot( imt, rhs )
        #print shifts
        tx = tx - shifts[0]
        ty = ty - shifts[1]
        tz = tz - shifts[2]
        gr.translation = [tx,ty,tz]
        ub = gr.ub.ravel()
        np.add(ub, shifts[3:], ub)
        gr.set_ubi( np.linalg.inv( np.reshape(ub,(3,3))))
        return gr
        

    def compute_XLYLZL(self):
        self.peaks_xyz = transform.compute_xyz_lab( [self.colfile.sc,
                                                    self.colfile.fc],
                                                  **self.pars.parameters)
        self.colfile.addcolumn( self.peaks_xyz[0], 'XL')
        self.colfile.addcolumn( self.peaks_xyz[1], 'YL')
        self.colfile.addcolumn( self.peaks_xyz[2], 'ZL')
        # These are the peak positions in the laboratory frame
        #
        # Now get the beam and peaks in the crystal frame
        wedge = float(self.pars.get('wedge'))
        chi   = float(self.pars.get('chi'))
        omega = self.colfile.omega
        beam = np.zeros( (3,len(omega)))
        beam[0,:] = 1
        self.XB = k_to_g( beam, omega, wedge, chi)
        #print "XB",XB
        self.colfile.addcolumn( self.XB[0], 'XB')
        self.colfile.addcolumn( self.XB[1], 'YB')
        self.colfile.addcolumn( self.XB[2], 'ZB')
        self.XS = k_to_g( self.peaks_xyz, omega, wedge, chi )
        self.colfile.addcolumn( self.XS[0], 'XS')
        self.colfile.addcolumn( self.XS[1], 'YS')
        self.colfile.addcolumn( self.XS[2], 'ZS')
        


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
        print (diffs*rscp < 0.05).sum()
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
    o = newindexer( sys.argv[1], sys.argv[2] )
    gl = grain.read_grain_file( sys.argv[3] )
    gref = gl[0]
    for i in range(5):
        print (3*"%.4f ")%tuple(gref.translation),
        print (6*"%.6f ")%ubitocellpars(gref.ubi)
        gref = o.refine(gref)
    print (3*"%.4f ")%tuple(gref.translation),
    print (6*"%.6f ")%ubitocellpars(gref.ubi)
    1/0
    for i in range(2):
        o.find_triplets( i*17 )
    pylab.show()

