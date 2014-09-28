"""
Interface to the bispev bit of fitpack.


fitpack (dierckx in netlib) --- A Python-C wrapper to FITPACK (by P. Dierckx).
        FITPACK is a collection of FORTRAN programs for CURVE and SURFACE
        FITTING with SPLINES and TENSOR PRODUCT SPLINES.

See
 http://www.cs.kuleuven.ac.be/cwis/research/nalag/research/topics/fitpack.html
or
 http://www.netlib.org/dierckx/index.html


(I think) This code was copied from fitpack.py and modified by J. P. Wright
so as to avoid having to install all of scipy. Hence the following notice:



Copyright 2002 Pearu Peterson all rights reserved,
Pearu Peterson <pearu@cens.ioc.ee>
Permission to use, modify, and distribute this software is given under the
terms of the SciPy (BSD style) license.  See LICENSE.txt that came with
scipy for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.

Pearu Peterson

"""



from ImageD11 import _splines
import numpy as np

def myasarray(a):
    if type(a) in [type(1.0),type(1L),type(1),type(1j)]:
        return np.asarray([a])
    elif hasattr(a, "shape") and len(a.shape)==0:
        # Takes care of mapping array(number) to array([number])
        return np.asarray([a])
    else:
        return np.asarray(a)

def bisplev(x,y,tck,dx=0,dy=0):
    """Evaluate a bivariate B-spline and its derivatives.
    Description:
      Return a rank-2 array of spline function values (or spline derivative
      values) at points given by the cross-product of the rank-1 arrays x and y.
      In special cases, return an array or just a float if either x or y or
      both are floats.
    Inputs:
      x, y -- Rank-1 arrays specifying the domain over which to evaluate the
              spline or its derivative.
      tck -- A sequence of length 5 returned by bisplrep containing the knot
             locations, the coefficients, and the degree of the spline:
             [tx, ty, c, kx, ky].
      dx, dy -- The orders of the partial derivatives in x and y respectively.
    Outputs: (vals, )
      vals -- The B-pline or its derivative evaluated over the set formed by
              the cross-product of x and y.
    """
    tx,ty,c,kx,ky=tck
    if not (0<=dx<kx): raise ValueError,"0<=dx=%d<kx=%d must hold"%(dx,kx)
    if not (0<=dy<ky): raise ValueError,"0<=dy=%d<ky=%d must hold"%(dy,ky)
    x,y=map(myasarray,[x,y])
    if (len(x.shape) != 1) or (len(y.shape) != 1):
        raise ValueError, "First two entries should be rank-1 arrays."
    z,ier=_splines._bispev(tx,ty,c,kx,ky,x,y,dx,dy)
    if ier==10: raise ValueError,"Invalid input data"
    if ier: raise TypeError,"An error occurred"
    z.shape=len(x),len(y)
    if len(z)>1: return z
    if len(z[0])>1: return z[0]
    return z[0][0]



_surfit_cache = {'tx': np.array([],'d'),'ty': np.array([],'d'),
                 'wrk': np.array([],'d'), 'iwrk':np.array([],'i')}

def bisplrep(x,y,z,w=None,xb=None,xe=None,yb=None,ye=None,kx=3,ky=3,task=0,s=None,
             eps=1e-16,tx=None,ty=None,full_output=0,nxest=None,nyest=None,quiet=1):
    """Find a bivariate B-spline representation of a surface.

    Description:

      Given a set of data points (x[i], y[i], z[i]) representing a surface
      z=f(x,y), compute a B-spline representation of the surface.

    Inputs:

      x, y, z -- Rank-1 arrays of data points.
      w -- Rank-1 array of weights. By default w=ones(len(x)).
      xb, xe -- End points of approximation interval in x.
      yb, ye -- End points of approximation interval in y.
                By default xb, xe, yb, ye = x[0], x[-1], y[0], y[-1]
      kx, ky -- The degrees of the spline (1 <= kx, ky <= 5).  Third order
                (kx=ky=3) is recommended.
      task -- If task=0, find knots in x and y and coefficients for a given
                smoothing factor, s.
              If task=1, find knots and coefficients for another value of the
                smoothing factor, s.  bisplrep must have been previously called
                with task=0 or task=1.
              If task=-1, find coefficients for a given set of knots tx, ty.
      s -- A non-negative smoothing factor.  If weights correspond to the inverse
           of the standard-deviation of the errors in z, then a good s-value
           should be found in the range (m-sqrt(2*m),m+sqrt(2*m)) where m=len(x)
      eps -- A threshold for determining the effective rank of an over-determined
             linear system of equations (0 < eps < 1) --- not likely to need
             changing.
      tx, ty -- Rank-1 arrays of the knots of the spline for task=-1
      full_output -- Non-zero to return optional outputs.
      nxest, nyest -- Over-estimates of the total number of knots.  If None then
                      nxest = max(kx+sqrt(m/2),2*kx+3),
                      nyest = max(ky+sqrt(m/2),2*ky+3)
      quiet -- Non-zero to suppress printing of messages.

    Outputs: (tck, {fp, ier, msg})

      tck -- A list [tx, ty, c, kx, ky] containing the knots (tx, ty) and
             coefficients (c) of the bivariate B-spline representation of the
             surface along with the degree of the spline.

      fp -- The weighted sum of squared residuals of the spline approximation.
      ier -- An integer flag about splrep success.  Success is indicated if
             ier<=0. If ier in [1,2,3] an error occurred but was not raised.
             Otherwise an error is raised.
      msg -- A message corresponding to the integer flag, ier.

    Remarks:

      SEE bisplev to evaluate the value of the B-spline given its tck
      representation.
    """
    x,y,z=map(myasarray,[x,y,z])
    x,y,z=map(np.ravel,[x,y,z])  # ensure 1-d arrays.
    m=len(x)
    if not (m==len(y)==len(z)): raise TypeError, 'len(x)==len(y)==len(z) must hold.'
    if w is None: w=np.ones(m,'d')
    else: w=myasarray(w)
    if not len(w) == m: raise TypeError,' len(w)=%d is not equal to m=%d'%(len(w),m)
    if xb is None: xb=x[0]
    if xe is None: xe=x[-1]
    if yb is None: yb=y[0]
    if ye is None: ye=y[-1]
    if not (-1<=task<=1): raise TypeError, 'task must be either -1,0, or 1'
    if s is None: s=m-np.sqrt(2*m)
    if tx is None and task==-1: raise TypeError, 'Knots_x must be given for task=-1'
    if tx is not None: _surfit_cache['tx']=myasarray(tx)
    nx=len(_surfit_cache['tx'])
    if ty is None and task==-1: raise TypeError, 'Knots_y must be given for task=-1'
    if ty is not None: _surfit_cache['ty']=myasarray(ty)
    ny=len(_surfit_cache['ty'])
    if task==-1 and nx<2*kx+2:
        raise TypeError, 'There must be at least 2*kx+2 knots_x for task=-1'
    if task==-1 and ny<2*ky+2:
        raise TypeError, 'There must be at least 2*ky+2 knots_x for task=-1'
    if not ((1<=kx<=5) and (1<=ky<=5)):
        raise TypeError, \
       'Given degree of the spline (kx,ky=%d,%d) is not supported. (1<=k<=5)'%(kx,ky)
    if m<(kx+1)*(ky+1): raise TypeError, 'm>=(kx+1)(ky+1) must hold'
    if nxest is None: nxest=kx+np.sqrt(m/2)
    if nyest is None: nyest=ky+np.sqrt(m/2)
    nxest,nyest=max(nxest,2*kx+3),max(nyest,2*ky+3)
    if task>=0 and s==0:
        nxest=int(kx+np.sqrt(3*m))
        nyest=int(ky+np.sqrt(3*m))
    if task==-1:
        _surfit_cache['tx']=myasarray(tx)
        _surfit_cache['ty']=myasarray(ty)
    tx,ty=_surfit_cache['tx'],_surfit_cache['ty']
    wrk=_surfit_cache['wrk']
    iwrk=_surfit_cache['iwrk']
    u,v,km,ne=nxest-kx-1,nyest-ky-1,max(kx,ky)+1,max(nxest,nyest)
    bx,by=kx*v+ky+1,ky*u+kx+1
    b1,b2=bx,bx+v-ky
    if bx>by: b1,b2=by,by+u-kx
    lwrk1=u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
    lwrk2=u*v*(b2+1)+b2
    tx,ty,c,o = _splines._surfit(x,y,z,w,xb,xe,yb,ye,kx,ky,task,s,eps,
                                   tx,ty,nxest,nyest,wrk,lwrk1,lwrk2)
    _surfit_cache['tx']=tx
    _surfit_cache['ty']=ty
    _surfit_cache['wrk']=o['wrk']
    ier,fp=o['ier'],o['fp']
    tck=[tx,ty,c,kx,ky]
    if ier<=0 and not quiet:
        import logging
        logging.debug("_surfit "+ str(_iermess2[ier][0]))
        logging.debug("\tkx,ky=%d,%d nx,ny=%d,%d m=%d fp=%f s=%f"%(kx,ky,len(tx),
                                                           len(ty),m,fp,s))
    ierm=min(11,max(-3,ier))
    if ierm>0 and not full_output:
        if ier in [1,2,3,4,5]:
            import logging
            logging.debug("_surfit Warning: "+str(_iermess2[ierm][0]))
            logging.debug("\tkx,ky=%d,%d nx,ny=%d,%d m=%d fp=%f s=%f"%(kx,ky,len(tx),
                                                           len(ty),m,fp,s))
        else:
            try:
                logging.error(str((_iermess2[ierm][1],_iermess2[ierm][0])))
                raise _iermess2[ierm][1],_iermess2[ierm][0]
            except KeyError:
                logging.error(str((_iermess2['unknown'][1],_iermess2['unknown'][0])))
                raise _iermess2['unknown'][1],_iermess2['unknown'][0]
    if full_output:
        try:
            return tck,fp,ier,_iermess2[ierm][0]
        except KeyError:
            return tck,fp,ier,_iermess2['unknown'][0]
    else:
        return tck
