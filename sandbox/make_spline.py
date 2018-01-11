
from __future__ import print_function

from ImageD11.transformer import transformer
import numpy as np, sys, os
from scipy.interpolate import bisplrep, bisplev
import pyFAI.spline
import pylab 

def dtthdpix( tr ):
    p0=tr.parameterobj.get_parameters().copy()
    t0,_ = tr.compute_tth_eta()
    p1 = p0.copy()
    p1['z_center'] += 1
    tr.parameterobj.set_parameters(p1)
    t1,_ = tr.compute_tth_eta()
    p1 = p0.copy()
    p1['y_center'] += 1
    tr.parameterobj.set_parameters(p1)
    t2,_ = tr.compute_tth_eta()
    d1 = t1-t0
    d2 = t2-t0
    d12 = np.sqrt(d1*d1+d2*d2)
    return d12


def write_spline(rets, retf, fname):
    s = pyFAI.spline.Spline()
    s.xmin=0
    s.xmax=4096
    s.ymin=0
    s.ymax=4096
    s.grid=1000
    s.pixelSize=(100.,100.)
    s.xSplineKnotsX = rets[0]
    s.xSplineKnotsY = rets[1]
    s.xSplineCoeff = rets[2]
    s.ySplineKnotsX = retf[0]
    s.ySplineKnotsY = retf[1]
    s.ySplineCoeff = retf[2]
    s.write(fname)
    
    

def fitspline(fltfile, parfile, splinefile):
    
    tr = transformer()
    tr.loadfiltered( fltfile )
    tr.loadfileparameters( parfile )
    print("fitting to assign peaks")
    tr.fit()
    tr.loadfileparameters( parfile )
    tthobs, eta = tr.compute_tth_eta()
    tthcalc = tthobs * 0.0 - 1
    for tthc, inds in zip( tr.tthc, tr.indices): # assignments
        tthcalc[inds] = tthc
    # project the error in the two directions
    sineta = np.sin(np.radians(eta))
    coseta = np.cos(np.radians(eta))
    # convert differences to pixels
    r  = np.sqrt( tr.colfile.yl**2 + tr.colfile.zl**2 )
    ps = tr.parameterobj.get("y_size")
    pix_obs = r / ps
    pix_calc = r * tthcalc/tthobs / ps
    diffs = pix_calc - pix_obs
    print(diffs)
    mask = tthcalc > 0
    xvals  = np.compress( mask, tr.colfile.s_raw )
    yvals  = np.compress( mask, tr.colfile.f_raw )
    dsvals = np.compress( mask, sineta * diffs )
    dfvals = np.compress( mask, coseta * diffs )
    w = np.ones(len(yvals))
    ss = 0.25
    for i in range(2):
        m = len(yvals)
        s = (m-np.sqrt(2*m))*ss
        print("s=",s)
        rets = bisplrep( yvals, xvals, dsvals, w=w, kx=3, ky=3, xb=0, xe=4096, 
                         yb=0, ye=4096, full_output=0, s = s, task = 0  )
        retf = bisplrep( yvals, xvals, dfvals, w=w, kx=3, ky=3, xb=0, xe=4096, 
                         yb=0, ye=4096, full_output=0, s = s,  task = 0  )
        print(rets, retf)
        dscalc = [ bisplev( y, x, rets ) for y,x in zip(yvals, xvals)]
        dfcalc = [ bisplev( y, x, retf ) for y,x in zip(yvals, xvals)]
        es = dscalc - dsvals
        ef = dfcalc - dfvals
        pylab.ion()
        pylab.figure(1)
        pylab.clf()
        pylab.subplot(121)
        pylab.xlim(0,4096)
        pylab.ylim(0,4096)
        pylab.scatter( xvals,yvals,c=es ,edgecolors='none')
        pylab.colorbar()
        pylab.subplot(122)
        pylab.scatter( xvals,yvals,c=ef ,edgecolors='none' )
        pylab.xlim(0,4096)
        pylab.ylim(0,4096)
        pylab.colorbar()
        pylab.figure(2)
        pylab.clf()
        pylab.hist( es, bins=128)
        pylab.hist( ef, bins=128)
        pylab.show()
        print("stddev=",((np.std(es) + np.std(ef))/2))
#        co = ((np.std(es) + np.std(ef))/2)*3
        co = float(input("cutoff"))
        print("Using cutoff",co)
        m = (np.abs( es ) < co ) & (np.abs( ef ) < co ) 
        yvals = np.compress(m, yvals)
        xvals = np.compress(m, xvals)
        dsvals = np.compress(m, dsvals)
        dfvals = np.compress(m, dfvals)
        w = np.ones(len(yvals))/((np.std(es) + np.std(ef))/2)
        print("w avg = ",w.mean())

    m = len(yvals)
    s = (m-np.sqrt(2*m))*ss
    print("s=",s)
    rets = bisplrep( yvals, xvals, dsvals, w=w, kx=3, ky=3, xb=0, xe=4096, 
                     yb=0, ye=4096, full_output=0, s = s, task = 0  )
    retf = bisplrep( yvals, xvals, dfvals, w=w, kx=3, ky=3, xb=0, xe=4096, 
                     yb=0, ye=4096, full_output=0, s = s,  task = 0  )

    print(rets, retf)
    write_spline( rets, retf, splinefile )
    pylab.show()
    input("End?")
    
    
def testspline( fltfile, parfile, splinefile):
#    print "Yop"
    sys.stdout.flush()
    os.system("rm test.flt")
    os.system("fix_spline.py %s test.flt %s"%(fltfile, splinefile))
    tr = transformer()
    tr.loadfiltered( "test.flt" )
    tr.loadfileparameters( parfile )
    tr.fit()


if __name__=="__main__":
    fitspline( sys.argv[1], sys.argv[2], sys.argv[3]) 
    testspline( sys.argv[1], sys.argv[2], sys.argv[3]) 
