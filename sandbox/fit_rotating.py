
from __future__ import print_function

from scipy.spatial.transform import Rotation, RotationSpline
from ImageD11 import grain, columnfile, cImageD11, transform
import sys, numpy as np, pylab as pl, os

folder = os.path.split(os.getcwd())[-1]

c = columnfile.columnfile( folder + "_t500.flt")# sys.argv[1] )
c.sortby('omega')
c.parameters.loadparameters( "../e2.par" ) # sys.argv[2] )
c.parameters.parameters['t_x']=0.
c.parameters.parameters['t_y']=0.
c.parameters.parameters['t_z']=0.
c.updateGeometry()
c.filter( abs( np.sin( np.radians( c.eta ))) > 0.05 )

g = grain.read_grain_file( folder+".map" )[0]# sys.argv[3] )[0] 

try:
    NROT = int(sys.argv[1])
except:
    NROT = 1
try:
    JUMP = float(sys.argv[2]) #-116.5
except:
    JUMP = None



gv = np.array( (c.gx, c.gy, c.gz) ).T.copy()

for name in 'h,k,l'.split(','):
    c.addcolumn(np.zeros(c.nrows),name)

step = c.nrows//10 + 1
start = 0

results = []

while start < c.nrows:
    end = min(start + step, c.nrows)
    gsel = gv[start:end].copy()
    order = np.argsort( c.ds[start:end] )
    gfit = gsel[order]
    ubi = g.ubi.copy()
    print(start,end,c.omega[start:end].mean(),end=" " )
    cImageD11.score_and_refine( ubi, gfit[:len(gfit)//3], 0.5)
    cImageD11.score_and_refine( ubi, gfit[:len(gfit)//3], 0.5)
    cImageD11.score_and_refine( ubi, gfit, 0.5)
    cImageD11.score_and_refine( ubi, gfit, 0.5)
    n, e = cImageD11.score_and_refine( ubi, gfit, 0.25)
    print(e,n,len(gfit))
    results.append(( c.omega[start:end].mean(), n, e, end-start, ubi ))
    hkl = np.round( np.dot( ubi, gsel.T ) )
    c.h[start:end]=hkl[0]
    c.k[start:end]=hkl[1]
    c.l[start:end]=hkl[2]
    start = end

# c.writefile("assigned.flt")

c.filter( (abs(c.h)+abs(c.k)+abs(c.l)) > 0)


class calcgrain:

    def __init__(self, hkl, sign_eta):
        assert hkl.shape[0] == 3
        self.hkl = hkl # peaks to be used
        self.sign_eta = sign_eta
        
    def gcalc( self, ub_or_spline ):
        """ send and average orientation of average plus rotation spline """
        if isinstance( ub_or_spline[0], RotationSpline ):
            rs, ub = ub_or_spline
            self.gve = np.dot( ub, self.hkl )
            bins = np.round( np.linspace( 0, len(self.gve[0]), NROT*10) ).astype(int)
            for i in range(len(bins)-1):
                bi = (bins[i]+bins[i+1])/2
                self.gve[:, bins[i]:bins[i+1]] = np.dot(rs(bi).as_dcm(),
                                                        self.gve[:, bins[i]:bins[i+1]])
        else:
            ub = ub_or_spline
            self.gve = np.dot( ub, self.hkl )
        return self.gve
    
    def angcalc( self, ub, wvln, wedge=0., chi=0.):
        gc = self.gcalc( ub )
        self.tth, (e0,e1), (o0,o1) = transform.uncompute_g_vectors(
            gc, wvln, wedge=wedge, chi=chi )
        assert e0.max() > 0 and e0.min() >= 0.
        assert e1.max() <=0 and e1.min() < 0   # negatives
        self.eta    = np.where( self.sign_eta > 0, e0, e1 )
        self.omega  = np.where( self.sign_eta > 0, o0, o1 )
        return self.tth, self.eta, self.omega

detparskeys = ( 'y_center', 'y_size', 'z_center', 'z_size',
                'tilt_x', 'tilt_y', 'tilt_z', 'wavelength',
                'distance', 'o11', 'o12', 'o21', 'o22' )

class calcdata:

    def __init__(self, scfc):
        self.sc, self.fc = scfc
        self.xlylzl=None
        
    def getxlylzl( self, **detpars ):
        self.xlylzl = transform.compute_xyz_lab(
            (self.sc, self.fc),
            **detpars )
        
    def getttheta( self, trans, omega, wedge=0., chi=0. ):
        #print("got trans",trans)
        self.tth, self.eta = transform.compute_tth_eta_from_xyz(
            self.xlylzl, omega,
            t_x=trans[0], t_y=trans[1], t_z=trans[2],
            wedge=wedge, chi=chi )
        return self.tth, self.eta

og = calcgrain( np.array( (c.h, c.k, c.l) ), np.sign(c.eta) )
od = calcdata ( np.array( (c.sc, c.fc)) )
od.getxlylzl( **c.parameters.parameters )

ubiavg = np.array( [r[-1] for r in results] ).mean(axis=0)
ub = np.linalg.inv( ubiavg )


tthc, etac, omegac = og.angcalc( ub, c.parameters.get('wavelength') )
tx,ty,tz=0,0,0
ttho, etao = od.getttheta( (tx, ty, tz), omegac ) # wedge, chi = 0 here

def mod360( x ):
    return np.mod( x + 180, 360 ) - 180


    

import lmfit 

def guess_weight( diff, sigma_cut=5 ):
    s = diff.std()
    c = sigma_cut * s
    m = abs(diff) < c
    s = diff[m].std()
    c = sigma_cut * s
    m = abs(diff) < c
    return m, diff/s

def guess_flat( diff, slop=0.475 ):
    d = mod360( diff )
    m0 = abs(d) < (slop*10)
    m =  abs(d) < slop
    e = np.where( d > slop, d-slop,
                  np.where( d < -slop, d+slop, 0 ) )  
    return m0, e
    

def B( a, b, c, al, be, ga):
    ca, cb, cg = [ np.cos( np.radians( x )) for x in (al, be, ga)]
    g = np.array( (( a*a, a*b*cg, a*c*cb),
                   ( a*b*cg, b*b, b*c*ca),
                   ( a*c*cb, b*c*ca, c*c)) )
    gi = np.linalg.inv( g )
    return np.linalg.cholesky( gi ).T
    
    
def U( rx, ry, rz ):
    u = Rotation.from_euler( 'XYZ', [np.radians(x) for x in (rx,ry,rz)]).as_dcm()
    return u

def getUB(p):
    x= np.dot( U0,
               B( p['a'].value,p['b'].value,p['c'].value,
                  p['al'].value,p['be'].value,p['ga'].value ) )
    return x

def residual( pars, x, cgr, cda, omega_obs, plot=False ):
    """
    pars = lmfit.Parameters
    x = ?
    cgr = calc grain object
    cda = calc data object
    """
#    aub = np.array( [[pars['ub%d%d'%(i,j)].value for j in range(3)] for i in range(3)] )
    aub = getUB( pars )
    #Uxyz = U( pars['rx'].value, pars['ry'].value, pars['rz'].value )
    #aub = np.dot( Uxyz, aub )
    if NROT == 1:
        aub = np.dot( Rotation.from_euler('XYZ',
                                          (pars['rx0'].value,
                                           pars['ry0'].value,
                                           pars['rz0'].value),
                                          degrees=True).as_dcm(), aub )
    else:
        times = np.linspace( 0, c.nrows, NROT)
        if JUMP is not None:
            ijump = np.searchsorted( c.omega, JUMP ) 
            ji = np.searchsorted( times, ijump )
            times[ji] = ijump+1
            times[ji-1] = ijump-1
        angles = [ ( pars['rx%d'%(i)].value, pars['ry%d'%(i)].value, pars['rz%d'%(i)].value)
               for i in range(NROT) ]
        rs = RotationSpline( times , Rotation.from_euler('XYZ', angles, degrees=True) )
        aub = rs,aub
    # FIXMEHERE SEND IN THE ROTATIONSPLINE ...
#    print(aub)
 #   1/0
    tthc, etac, omegac = cgr.angcalc( aub, pars['wavelength'],
                                      wedge=pars['wedge'].value )
    arg = {k:pars[k].value for k in detparskeys} 
    cda.getxlylzl( **arg )
    ttho, etao = cda.getttheta( (pars['tx'].value, pars['ty'].value, pars['tz'].value),
                                omegac, wedge=pars['wedge'].value ) # wedge, chi = 0 here
    mtth, etth = guess_weight( ttho - tthc )
    meta, eeta = guess_weight( etao - etac )
    mome, eome = guess_flat( omega_obs- omegac)
#    m = np.ones( ttho.shape, dtype=bool)
    m = mtth & meta & mome
    e= np.concatenate( (np.where( m, etth,0),
                        np.where( m, eeta,0),
                        np.where( m, eome,0) ) )

    if plot:
        print('in residual', m.sum(), (e*e).sum())
        pl.figure()
        pl.subplot(231)
        pl.plot( tthc, tthc-ttho, "." )
        pl.plot( tthc[~m], (tthc-ttho)[~m], "+" )
        pl.title(' tth ')
        pl.subplot(232)
        pl.plot( etac, etac-etao, ".")
        pl.plot( etac[~m], (etac-etao)[~m], "+")
        pl.title(' eta ');pl.ylabel('deg')
        pl.subplot(233)
        pl.plot( c.omega, mod360(omegac-c.omega), ".")
        pl.plot( c.omega[~m], mod360(omegac-c.omega)[~m], "+")
        pl.title(' omega ')
        pl.subplot(234)
        pl.plot( tthc, etth*m, "." )
        pl.title(' tth ')
        pl.subplot(235)
        pl.plot( omegac, eeta*m, ".")
        pl.title(' eta ')
        pl.subplot(236)
        pl.plot( omegac, eome*m, ".")
        pl.title(' omega ')
        if NROT > 1:        
            pl.figure()
            j = np.round(np.linspace( 0, len(c.omega), NROT*10)).astype(int)
            oj = [c.omega[j[i]:j[i+1]].mean() for i in range(len(j)-1) ]
            r = [rs((j[i]+j[i+1])/2) for i in range(len(j)-1)]
            pl.plot( oj, [x.as_euler('XYZ',degrees=True)[0] for x in r], 'o-', label='rx')
            pl.plot( oj, [x.as_euler('XYZ',degrees=True)[1] for x in r], 'o-', label='ry')
            pl.plot( oj, [x.as_euler('XYZ',degrees=True)[2] for x in r], 'o-', label='rz')
        pl.figure()
        if NROT == 1:
            ubi = np.linalg.inv( aub )
        else:
            ub = np.dot( aub[0](0).as_dcm(), aub[1] )
            ubi = np.linalg.inv( ub )
        gve = transform.compute_g_vectors( ttho, etao, c.omega,
                                           pars['wavelength'].value,
                                           wedge = pars['wedge'].value )
        h,k,l = np.dot( ubi, gve )
        

        pl.scatter( abs(h-np.round(h)), abs(k-np.round(k)), c=abs(l-np.round(l)), s = pow( c.sum_intensity, 1/3))
        pl.axes().set_aspect('equal')
        pl.xlim(-0.6,0.6)
        pl.ylim(-0.6,0.6)
        pl.xticks(np.linspace(-0.5,0.5,3))
        pl.yticks(np.linspace(-0.5,0.5,3))
        pl.grid(True)
        pl.colorbar()
        pl.show()
    return e
    
    
    


print(ub)
print("check")


p = lmfit.Parameters()
for k in detparskeys:
    p.add( k, c.parameters.get( k ), vary = False )

p.add('tx',0, vary=True)
p.add('ty',0, vary=True)
p.add('tz',0, vary=True)
p.add('wedge',0,vary=False)

import ImageD11.indexing
U0 = ImageD11.indexing.ubi_to_u( ubiavg )
for name, val in zip( 'a,b,c,al,be,ga'.split(','),
                      ImageD11.indexing.ubitocellpars( ubiavg) ):
    p.add( name, val, vary= True)

    
for i in range(NROT):
    p.add('rx%d'%(i),0,vary=True)
    p.add('ry%d'%(i),0,vary=True)
    p.add('rz%d'%(i),0,vary=True)

#for i in range(3):
#    for j in range(3):
#        p.add( "ub%d%d"%(i,j), ub[i,j], vary=True)

x = np.arange(c.nrows*3, dtype=np.int )


#e0 = residual(  p, x, og, od, c.omega , plot=False)
#p['ub00'].value+=1e-8
#e1 = residual(  p, x, og, od, c.omega , plot=False)
#pl.plot(  (e0-e1)*1e8,".")
#pl.show()
#1/0

resids={}
def iter_cb( p, i, r, *a, **k):
    resids[i] = r.copy()
    print(i, (r*r).sum())#, end=" " )
#    for k in resids.keys():
#        if i == k : continue
#        if (resids[i] == r).all():
#            print("==",k,end=" ")
#    print()

#import pdb;pdb.set_trace()
out = lmfit.minimize( residual, p, args=(x, og, od, c.omega ),
                      nan_policy='raise',
                      iter_cb=iter_cb,
                     )
print(lmfit.fit_report( out ))



p = out.params
residual(  p, x, og, od, c.omega , plot=True)
sys.exit()
p['tilt_x'].vary = True
p['tilt_y'].vary = True
p['tilt_z'].vary = True
p['wedge'].vary = True
p['y_center'].vary=True
p['distance'].vary=True
#p['z_center'].vary=True

out = lmfit.minimize( residual, p, args=(x, og, od, c.omega ) )
print(lmfit.fit_report( out ))
p=out.params
residual(  p, x, og, od, c.omega , plot=True)
p['distance'].vary = False
out = lmfit.minimize( residual, p, args=(x, og, od, c.omega ) )
print(lmfit.fit_report( out ))
p=out.params
residual(  p, x, og, od, c.omega , plot=True)

for k in detparskeys:
    print(k,p[k].value)
print('wedge',p['wedge'].value)
