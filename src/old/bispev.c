/* bispev.f -- translated by f2c (version 20041007).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"


/* Routines from netlib */
/* Copyright by Paul Dierckx */


/* Subroutine */ int bispev_(doublereal *tx, integer *nx, doublereal *ty, 
	integer *ny, doublereal *c__, integer *kx, integer *ky, doublereal *x,
	 integer *mx, doublereal *y, integer *my, doublereal *z__, doublereal 
	*wrk, integer *lwrk, integer *iwrk, integer *kwrk, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, iw, lwest;
    extern /* Subroutine */ int fpbisp_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

/*  subroutine bispev evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,... */
/*  ,my a bivariate spline s(x,y) of degrees kx and ky, given in the */
/*  b-spline representation. */

/*  calling sequence: */
/*     call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk, */
/*    * iwrk,kwrk,ier) */

/*  input parameters: */
/*   tx    : real array, length nx, which contains the position of the */
/*           knots in the x-direction. */
/*   nx    : integer, giving the total number of knots in the x-direction */
/*   ty    : real array, length ny, which contains the position of the */
/*           knots in the y-direction. */
/*   ny    : integer, giving the total number of knots in the y-direction */
/*   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the */
/*           b-spline coefficients. */
/*   kx,ky : integer values, giving the degrees of the spline. */
/*   x     : real array of dimension (mx). */
/*           before entry x(i) must be set to the x co-ordinate of the */
/*           i-th grid point along the x-axis. */
/*           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx. */
/*   mx    : on entry mx must specify the number of grid points along */
/*           the x-axis. mx >=1. */
/*   y     : real array of dimension (my). */
/*           before entry y(j) must be set to the y co-ordinate of the */
/*           j-th grid point along the y-axis. */
/*           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my. */
/*   my    : on entry my must specify the number of grid points along */
/*           the y-axis. my >=1. */
/*   wrk   : real array of dimension lwrk. used as workspace. */
/*   lwrk  : integer, specifying the dimension of wrk. */
/*           lwrk >= mx*(kx+1)+my*(ky+1) */
/*   iwrk  : integer array of dimension kwrk. used as workspace. */
/*   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my. */

/*  output parameters: */
/*   z     : real array of dimension (mx*my). */
/*           on succesful exit z(my*(i-1)+j) contains the value of s(x,y) */
/*           at the point (x(i),y(j)),i=1,...,mx;j=1,...,my. */
/*   ier   : integer error flag */
/*    ier=0 : normal return */
/*    ier=10: invalid input data (see restrictions) */

/*  restrictions: */
/*   mx >=1, my >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my */
/*   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx */
/*   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my */

/*  other subroutines required: */
/*    fpbisp,fpbspl */

/*  references : */
/*    de boor c : on calculating with b-splines, j. approximation theory */
/*                6 (1972) 50-62. */
/*    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths */
/*                applics 10 (1972) 134-149. */
/*    dierckx p. : curve and surface fitting with splines, monographs on */
/*                 numerical analysis, oxford university press, 1993. */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1987 */

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
    /* Parameter adjustments */
    --tx;
    --ty;
    --c__;
    --x;
    --z__;
    --y;
    --wrk;
    --iwrk;

    /* Function Body */
    *ier = 10;
    lwest = (*kx + 1) * *mx + (*ky + 1) * *my;
    if (*lwrk < lwest) {
	goto L100;
    }
    if (*kwrk < *mx + *my) {
	goto L100;
    }
    if ((i__1 = *mx - 1) < 0) {
	goto L100;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L10;
    }
L10:
    i__1 = *mx;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (x[i__] < x[i__ - 1]) {
	    goto L100;
	}
/* L20: */
    }
L30:
    if ((i__1 = *my - 1) < 0) {
	goto L100;
    } else if (i__1 == 0) {
	goto L60;
    } else {
	goto L40;
    }
L40:
    i__1 = *my;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (y[i__] < y[i__ - 1]) {
	    goto L100;
	}
/* L50: */
    }
L60:
    *ier = 0;
    iw = *mx * (*kx + 1) + 1;
    fpbisp_(&tx[1], nx, &ty[1], ny, &c__[1], kx, ky, &x[1], mx, &y[1], my, &
	    z__[1], &wrk[1], &wrk[iw], &iwrk[1], &iwrk[*mx + 1]);
L100:
    return 0;
} /* bispev_ */

/* Subroutine */ int fpbisp_(doublereal *tx, integer *nx, doublereal *ty, 
	integer *ny, doublereal *c__, integer *kx, integer *ky, doublereal *x,
	 integer *mx, doublereal *y, integer *my, doublereal *z__, doublereal 
	*wx, doublereal *wy, integer *lx, integer *ly)
{
    /* System generated locals */
    integer wx_dim1, wx_offset, wy_dim1, wy_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static doublereal h__[6];
    static integer i__, j, l, m, i1, j1, l1, l2;
    static doublereal tb, te, sp;
    static integer kx1, ky1;
    static doublereal arg;
    static integer nkx1, nky1;
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  ..subroutine references.. */
/*    fpbspl */
/*  .. */
    /* Parameter adjustments */
    --tx;
    --ty;
    --c__;
    --lx;
    wx_dim1 = *mx;
    wx_offset = 1 + wx_dim1;
    wx -= wx_offset;
    --x;
    --ly;
    wy_dim1 = *my;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    --z__;
    --y;

    /* Function Body */
    kx1 = *kx + 1;
    nkx1 = *nx - kx1;
    tb = tx[kx1];
    te = tx[nkx1 + 1];
    l = kx1;
    l1 = l + 1;
    i__1 = *mx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	arg = x[i__];
	if (arg < tb) {
	    arg = tb;
	}
	if (arg > te) {
	    arg = te;
	}
L10:
	if (arg < tx[l1] || l == nkx1) {
	    goto L20;
	}
	l = l1;
	l1 = l + 1;
	goto L10;
L20:
	fpbspl_(&tx[1], nx, kx, &arg, &l, h__);
	lx[i__] = l - kx1;
	i__2 = kx1;
	for (j = 1; j <= i__2; ++j) {
	    wx[i__ + j * wx_dim1] = h__[j - 1];
/* L30: */
	}
/* L40: */
    }
    ky1 = *ky + 1;
    nky1 = *ny - ky1;
    tb = ty[ky1];
    te = ty[nky1 + 1];
    l = ky1;
    l1 = l + 1;
    i__1 = *my;
    for (i__ = 1; i__ <= i__1; ++i__) {
	arg = y[i__];
	if (arg < tb) {
	    arg = tb;
	}
	if (arg > te) {
	    arg = te;
	}
L50:
	if (arg < ty[l1] || l == nky1) {
	    goto L60;
	}
	l = l1;
	l1 = l + 1;
	goto L50;
L60:
	fpbspl_(&ty[1], ny, ky, &arg, &l, h__);
	ly[i__] = l - ky1;
	i__2 = ky1;
	for (j = 1; j <= i__2; ++j) {
	    wy[i__ + j * wy_dim1] = h__[j - 1];
/* L70: */
	}
/* L80: */
    }
    m = 0;
    i__1 = *mx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = lx[i__] * nky1;
	i__2 = kx1;
	for (i1 = 1; i1 <= i__2; ++i1) {
	    h__[i1 - 1] = wx[i__ + i1 * wx_dim1];
/* L90: */
	}
	i__2 = *my;
	for (j = 1; j <= i__2; ++j) {
	    l1 = l + ly[j];
	    sp = 0.f;
	    i__3 = kx1;
	    for (i1 = 1; i1 <= i__3; ++i1) {
		l2 = l1;
		i__4 = ky1;
		for (j1 = 1; j1 <= i__4; ++j1) {
		    ++l2;
		    sp += c__[l2] * h__[i1 - 1] * wy[j + j1 * wy_dim1];
/* L100: */
		}
		l1 += nky1;
/* L110: */
	    }
	    ++m;
	    z__[m] = sp;
/* L120: */
	}
/* L130: */
    }
    return 0;
} /* fpbisp_ */

/* Subroutine */ int fpbspl_(doublereal *t, integer *n, integer *k, 
	doublereal *x, integer *l, doublereal *h__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal f;
    static integer i__, j;
    static doublereal hh[5];
    static integer li, lj;
    static doublereal one;

/*  subroutine fpbspl evaluates the (k+1) non-zero b-splines of */
/*  degree k at t(l) <= x < t(l+1) using the stable recurrence */
/*  relation of de boor and cox. */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  .. */
    /* Parameter adjustments */
    --t;
    --h__;

    /* Function Body */
    one = 1.;
    h__[1] = one;
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    hh[i__ - 1] = h__[i__];
/* L10: */
	}
	h__[1] = 0.;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    li = *l + i__;
	    lj = li - j;
	    f = hh[i__ - 1] / (t[li] - t[lj]);
	    h__[i__] += f * (t[li] - *x);
	    h__[i__ + 1] = f * (*x - t[lj]);
/* L20: */
	}
    }
    return 0;
} /* fpbspl_ */

/* Subroutine */ int parder_(doublereal *tx, integer *nx, doublereal *ty, 
	integer *ny, doublereal *c__, integer *kx, integer *ky, integer *nux, 
	integer *nuy, doublereal *x, integer *mx, doublereal *y, integer *my, 
	doublereal *z__, doublereal *wrk, integer *lwrk, integer *iwrk, 
	integer *kwrk, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, m, l1, l2, m0, m1;
    static doublereal ak;
    static integer nc, lx, ly, kx1, ky1;
    static doublereal fac;
    static integer kkx, kky, iwx, iwy, nxx, nyy, nkx1, nky1, lwest;
    extern /* Subroutine */ int fpbisp_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

/*  subroutine parder evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,... */
/*  ,my the partial derivative ( order nux,nuy) of a bivariate spline */
/*  s(x,y) of degrees kx and ky, given in the b-spline representation. */

/*  calling sequence: */
/*     call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk, */
/*    * iwrk,kwrk,ier) */

/*  input parameters: */
/*   tx    : real array, length nx, which contains the position of the */
/*           knots in the x-direction. */
/*   nx    : integer, giving the total number of knots in the x-direction */
/*   ty    : real array, length ny, which contains the position of the */
/*           knots in the y-direction. */
/*   ny    : integer, giving the total number of knots in the y-direction */
/*   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the */
/*           b-spline coefficients. */
/*   kx,ky : integer values, giving the degrees of the spline. */
/*   nux   : integer values, specifying the order of the partial */
/*   nuy     derivative. 0<=nux<kx, 0<=nuy<ky. */
/*   x     : real array of dimension (mx). */
/*           before entry x(i) must be set to the x co-ordinate of the */
/*           i-th grid point along the x-axis. */
/*           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx. */
/*   mx    : on entry mx must specify the number of grid points along */
/*           the x-axis. mx >=1. */
/*   y     : real array of dimension (my). */
/*           before entry y(j) must be set to the y co-ordinate of the */
/*           j-th grid point along the y-axis. */
/*           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my. */
/*   my    : on entry my must specify the number of grid points along */
/*           the y-axis. my >=1. */
/*   wrk   : real array of dimension lwrk. used as workspace. */
/*   lwrk  : integer, specifying the dimension of wrk. */
/*           lwrk >= mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1) */
/*   iwrk  : integer array of dimension kwrk. used as workspace. */
/*   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my. */

/*  output parameters: */
/*   z     : real array of dimension (mx*my). */
/*           on succesful exit z(my*(i-1)+j) contains the value of the */
/*           specified partial derivative of s(x,y) at the point */
/*           (x(i),y(j)),i=1,...,mx;j=1,...,my. */
/*   ier   : integer error flag */
/*    ier=0 : normal return */
/*    ier=10: invalid input data (see restrictions) */

/*  restrictions: */
/*   mx >=1, my >=1, 0 <= nux < kx, 0 <= nuy < ky, kwrk>=mx+my */
/*   lwrk>=mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1), */
/*   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx */
/*   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my */

/*  other subroutines required: */
/*    fpbisp,fpbspl */

/*  references : */
/*    de boor c : on calculating with b-splines, j. approximation theory */
/*                6 (1972) 50-62. */
/*   dierckx p. : curve and surface fitting with splines, monographs on */
/*                numerical analysis, oxford university press, 1993. */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1989 */

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
    /* Parameter adjustments */
    --tx;
    --ty;
    --c__;
    --x;
    --z__;
    --y;
    --wrk;
    --iwrk;

    /* Function Body */
    *ier = 10;
    kx1 = *kx + 1;
    ky1 = *ky + 1;
    nkx1 = *nx - kx1;
    nky1 = *ny - ky1;
    nc = nkx1 * nky1;
    if (*nux < 0 || *nux >= *kx) {
	goto L400;
    }
    if (*nuy < 0 || *nuy >= *ky) {
	goto L400;
    }
    lwest = nc + (kx1 - *nux) * *mx + (ky1 - *nuy) * *my;
    if (*lwrk < lwest) {
	goto L400;
    }
    if (*kwrk < *mx + *my) {
	goto L400;
    }
    if ((i__1 = *mx - 1) < 0) {
	goto L400;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L10;
    }
L10:
    i__1 = *mx;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (x[i__] < x[i__ - 1]) {
	    goto L400;
	}
/* L20: */
    }
L30:
    if ((i__1 = *my - 1) < 0) {
	goto L400;
    } else if (i__1 == 0) {
	goto L60;
    } else {
	goto L40;
    }
L40:
    i__1 = *my;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (y[i__] < y[i__ - 1]) {
	    goto L400;
	}
/* L50: */
    }
L60:
    *ier = 0;
    nxx = nkx1;
    nyy = nky1;
    kkx = *kx;
    kky = *ky;
/*  the partial derivative of order (nux,nuy) of a bivariate spline of */
/*  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy. */
/*  we calculate the b-spline coefficients of this spline */
    i__1 = nc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wrk[i__] = c__[i__];
/* L70: */
    }
    if (*nux == 0) {
	goto L200;
    }
    lx = 1;
    i__1 = *nux;
    for (j = 1; j <= i__1; ++j) {
	ak = (doublereal) kkx;
	--nxx;
	l1 = lx;
	m0 = 1;
	i__2 = nxx;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ++l1;
	    l2 = l1 + kkx;
	    fac = tx[l2] - tx[l1];
	    if (fac <= 0.f) {
		goto L90;
	    }
	    i__3 = nyy;
	    for (m = 1; m <= i__3; ++m) {
		m1 = m0 + nyy;
		wrk[m0] = (wrk[m1] - wrk[m0]) * ak / fac;
		++m0;
/* L80: */
	    }
L90:
	    ;
	}
	++lx;
	--kkx;
/* L100: */
    }
L200:
    if (*nuy == 0) {
	goto L300;
    }
    ly = 1;
    i__1 = *nuy;
    for (j = 1; j <= i__1; ++j) {
	ak = (doublereal) kky;
	--nyy;
	l1 = ly;
	i__2 = nyy;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ++l1;
	    l2 = l1 + kky;
	    fac = ty[l2] - ty[l1];
	    if (fac <= 0.f) {
		goto L220;
	    }
	    m0 = i__;
	    i__3 = nxx;
	    for (m = 1; m <= i__3; ++m) {
		m1 = m0 + 1;
		wrk[m0] = (wrk[m1] - wrk[m0]) * ak / fac;
		m0 += nky1;
/* L210: */
	    }
L220:
	    ;
	}
	++ly;
	--kky;
/* L230: */
    }
    m0 = nyy;
    m1 = nky1;
    i__1 = nxx;
    for (m = 2; m <= i__1; ++m) {
	i__2 = nyy;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ++m0;
	    ++m1;
	    wrk[m0] = wrk[m1];
/* L240: */
	}
	m1 += *nuy;
/* L250: */
    }
/*  we partition the working space and evaluate the partial derivative */
L300:
    iwx = nxx * nyy + 1;
    iwy = iwx + *mx * (kx1 - *nux);
    i__1 = *nx - (*nux << 1);
    i__2 = *ny - (*nuy << 1);
    fpbisp_(&tx[*nux + 1], &i__1, &ty[*nuy + 1], &i__2, &wrk[1], &kkx, &kky, &
	    x[1], mx, &y[1], my, &z__[1], &wrk[iwx], &wrk[iwy], &iwrk[1], &
	    iwrk[*mx + 1]);
L400:
    return 0;
} /* parder_ */


/* Subroutine */ int surfit_(integer *iopt, integer *m, doublereal *x, 
	doublereal *y, doublereal *z__, doublereal *w, doublereal *xb, 
	doublereal *xe, doublereal *yb, doublereal *ye, integer *kx, integer *
	ky, doublereal *s, integer *nxest, integer *nyest, integer *nmax, 
	doublereal *eps, integer *nx, doublereal *tx, integer *ny, doublereal 
	*ty, doublereal *c__, doublereal *fp, doublereal *wrk1, integer *
	lwrk1, doublereal *wrk2, integer *lwrk2, integer *iwrk, integer *kwrk,
	 integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, la, lf, ki, lh, kn, lq, ib1, jb1, ib3, km1, km2, kx1, 
	    ky1, lff, lco, nek, lfp, lbx, lby;
    static doublereal tol;
    static integer nxk, nyk, nmx, nmy, lsx, lsy, nreg, kmax, nest, ncest, 
	    maxit, nminx, nminy, nrint, kwest, lwest;
    extern /* Subroutine */ int fpsurf_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *);

/* given the set of data points (x(i),y(i),z(i)) and the set of positive */
/* numbers w(i),i=1,...,m, subroutine surfit determines a smooth bivar- */
/* iate spline approximation s(x,y) of degrees kx and ky on the rect- */
/* angle xb <= x <= xe, yb <= y <= ye. */
/* if iopt = -1 surfit calculates the weighted least-squares spline */
/* according to a given set of knots. */
/* if iopt >= 0 the total numbers nx and ny of these knots and their */
/* position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic- */
/* ally by the routine. the smoothness of s(x,y) is then achieved by */
/* minimalizing the discontinuity jumps in the derivatives of s(x,y) */
/* across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1). */
/* the amounth of smoothness is determined by the condition that f(p) = */
/* sum ((w(i)*(z(i)-s(x(i),y(i))))**2) be <= s, with s a given non-neg- */
/* ative constant, called the smoothing factor. */
/* the fit is given in the b-spline representation (b-spline coefficients */
/* c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval- */
/* uated by means of subroutine bispev. */

/* calling sequence: */
/*     call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest, */
/*    *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier) */

/* parameters: */
/*  iopt  : integer flag. on entry iopt must specify whether a weighted */
/*          least-squares spline (iopt=-1) or a smoothing spline (iopt=0 */
/*          or 1) must be determined. */
/*          if iopt=0 the routine will start with an initial set of knots */
/*          tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i= */
/*          1,...,ky+1. if iopt=1 the routine will continue with the set */
/*          of knots found at the last call of the routine. */
/*          attention: a call with iopt=1 must always be immediately pre- */
/*                     ceded by another call with iopt=1 or iopt=0. */
/*          unchanged on exit. */
/*  m     : integer. on entry m must specify the number of data points. */
/*          m >= (kx+1)*(ky+1). unchanged on exit. */
/*  x     : real array of dimension at least (m). */
/*  y     : real array of dimension at least (m). */
/*  z     : real array of dimension at least (m). */
/*          before entry, x(i),y(i),z(i) must be set to the co-ordinates */
/*          of the i-th data point, for i=1,...,m. the order of the data */
/*          points is immaterial. unchanged on exit. */
/*  w     : real array of dimension at least (m). before entry, w(i) must */
/*          be set to the i-th value in the set of weights. the w(i) must */
/*          be strictly positive. unchanged on exit. */
/*  xb,xe : real values. on entry xb,xe,yb and ye must specify the bound- */
/*  yb,ye   aries of the rectangular approximation domain. */
/*          xb<=x(i)<=xe,yb<=y(i)<=ye,i=1,...,m. unchanged on exit. */
/*  kx,ky : integer values. on entry kx and ky must specify the degrees */
/*          of the spline. 1<=kx,ky<=5. it is recommended to use bicubic */
/*          (kx=ky=3) splines. unchanged on exit. */
/*  s     : real. on entry (in case iopt>=0) s must specify the smoothing */
/*          factor. s >=0. unchanged on exit. */
/*          for advice on the choice of s see further comments */
/*  nxest : integer. unchanged on exit. */
/*  nyest : integer. unchanged on exit. */
/*          on entry, nxest and nyest must specify an upper bound for the */
/*          number of knots required in the x- and y-directions respect. */
/*          these numbers will also determine the storage space needed by */
/*          the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1). */
/*          in most practical situation nxest = kx+1+sqrt(m/2), nyest = */
/*          ky+1+sqrt(m/2) will be sufficient. see also further comments. */
/*  nmax  : integer. on entry nmax must specify the actual dimension of */
/*          the arrays tx and ty. nmax >= nxest, nmax >=nyest. */
/*          unchanged on exit. */
/*  eps   : real. */
/*          on entry, eps must specify a threshold for determining the */
/*          effective rank of an over-determined linear system of equat- */
/*          ions. 0 < eps < 1.  if the number of decimal digits in the */
/*          computer representation of a real number is q, then 10**(-q) */
/*          is a suitable value for eps in most practical applications. */
/*          unchanged on exit. */
/*  nx    : integer. */
/*          unless ier=10 (in case iopt >=0), nx will contain the total */
/*          number of knots with respect to the x-variable, of the spline */
/*          approximation returned. if the computation mode iopt=1 is */
/*          used, the value of nx should be left unchanged between sub- */
/*          sequent calls. */
/*          in case iopt=-1, the value of nx should be specified on entry */
/*  tx    : real array of dimension nmax. */
/*          on succesful exit, this array will contain the knots of the */
/*          spline with respect to the x-variable, i.e. the position of */
/*          the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the */
/*          position of the additional knots tx(1)=...=tx(kx+1)=xb and */
/*          tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat. */
/*          if the computation mode iopt=1 is used, the values of tx(1), */
/*          ...,tx(nx) should be left unchanged between subsequent calls. */
/*          if the computation mode iopt=-1 is used, the values tx(kx+2), */
/*          ...tx(nx-kx-1) must be supplied by the user, before entry. */
/*          see also the restrictions (ier=10). */
/*  ny    : integer. */
/*          unless ier=10 (in case iopt >=0), ny will contain the total */
/*          number of knots with respect to the y-variable, of the spline */
/*          approximation returned. if the computation mode iopt=1 is */
/*          used, the value of ny should be left unchanged between sub- */
/*          sequent calls. */
/*          in case iopt=-1, the value of ny should be specified on entry */
/*  ty    : real array of dimension nmax. */
/*          on succesful exit, this array will contain the knots of the */
/*          spline with respect to the y-variable, i.e. the position of */
/*          the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the */
/*          position of the additional knots ty(1)=...=ty(ky+1)=yb and */
/*          ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat. */
/*          if the computation mode iopt=1 is used, the values of ty(1), */
/*          ...,ty(ny) should be left unchanged between subsequent calls. */
/*          if the computation mode iopt=-1 is used, the values ty(ky+2), */
/*          ...ty(ny-ky-1) must be supplied by the user, before entry. */
/*          see also the restrictions (ier=10). */
/*  c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1). */
/*          on succesful exit, c contains the coefficients of the spline */
/*          approximation s(x,y) */
/*  fp    : real. unless ier=10, fp contains the weighted sum of */
/*          squared residuals of the spline approximation returned. */
/*  wrk1  : real array of dimension (lwrk1). used as workspace. */
/*          if the computation mode iopt=1 is used the value of wrk1(1) */
/*          should be left unchanged between subsequent calls. */
/*          on exit wrk1(2),wrk1(3),...,wrk1(1+(nx-kx-1)*(ny-ky-1)) will */
/*          contain the values d(i)/max(d(i)),i=1,...,(nx-kx-1)*(ny-ky-1) */
/*          with d(i) the i-th diagonal element of the reduced triangular */
/*          matrix for calculating the b-spline coefficients. it includes */
/*          those elements whose square is less than eps,which are treat- */
/*          ed as 0 in the case of presumed rank deficiency (ier<-2). */
/*  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of */
/*          the array wrk1 as declared in the calling (sub)program. */
/*          lwrk1 must not be too small. let */
/*            u = nxest-kx-1, v = nyest-ky-1, km = max(kx,ky)+1, */
/*            ne = max(nxest,nyest), bx = kx*v+ky+1, by = ky*u+kx+1, */
/*            if(bx.le.by) b1 = bx, b2 = b1+v-ky */
/*            if(bx.gt.by) b1 = by, b2 = b1+u-kx  then */
/*          lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1 */
/*  wrk2  : real array of dimension (lwrk2). used as workspace, but */
/*          only in the case a rank deficient system is encountered. */
/*  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of */
/*          the array wrk2 as declared in the calling (sub)program. */
/*          lwrk2 > 0 . a save upper boundfor lwrk2 = u*v*(b2+1)+b2 */
/*          where u,v and b2 are as above. if there are enough data */
/*          points, scattered uniformly over the approximation domain */
/*          and if the smoothing factor s is not too small, there is a */
/*          good chance that this extra workspace is not needed. a lot */
/*          of memory might therefore be saved by setting lwrk2=1. */
/*          (see also ier > 10) */
/*  iwrk  : integer array of dimension (kwrk). used as workspace. */
/*  kwrk  : integer. on entry kwrk must specify the actual dimension of */
/*          the array iwrk as declared in the calling (sub)program. */
/*          kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1). */
/*  ier   : integer. unless the routine detects an error, ier contains a */
/*          non-positive value on exit, i.e. */
/*   ier=0  : normal return. the spline returned has a residual sum of */
/*            squares fp such that abs(fp-s)/s <= tol with tol a relat- */
/*            ive tolerance set to 0.001 by the program. */
/*   ier=-1 : normal return. the spline returned is an interpolating */
/*            spline (fp=0). */
/*   ier=-2 : normal return. the spline returned is the weighted least- */
/*            squares polynomial of degrees kx and ky. in this extreme */
/*            case fp gives the upper bound for the smoothing factor s. */
/*   ier<-2 : warning. the coefficients of the spline returned have been */
/*            computed as the minimal norm least-squares solution of a */
/*            (numerically) rank deficient system. (-ier) gives the rank. */
/*            especially if the rank deficiency which can be computed as */
/*            (nx-kx-1)*(ny-ky-1)+ier, is large the results may be inac- */
/*            curate. they could also seriously depend on the value of */
/*            eps. */
/*   ier=1  : error. the required storage space exceeds the available */
/*            storage space, as specified by the parameters nxest and */
/*            nyest. */
/*            probably causes : nxest or nyest too small. if these param- */
/*            eters are already large, it may also indicate that s is */
/*            too small */
/*            the approximation returned is the weighted least-squares */
/*            spline according to the current set of knots. */
/*            the parameter fp gives the corresponding weighted sum of */
/*            squared residuals (fp>s). */
/*   ier=2  : error. a theoretically impossible result was found during */
/*            the iteration proces for finding a smoothing spline with */
/*            fp = s. probably causes : s too small or badly chosen eps. */
/*            there is an approximation returned but the corresponding */
/*            weighted sum of squared residuals does not satisfy the */
/*            condition abs(fp-s)/s < tol. */
/*   ier=3  : error. the maximal number of iterations maxit (set to 20 */
/*            by the program) allowed for finding a smoothing spline */
/*            with fp=s has been reached. probably causes : s too small */
/*            there is an approximation returned but the corresponding */
/*            weighted sum of squared residuals does not satisfy the */
/*            condition abs(fp-s)/s < tol. */
/*   ier=4  : error. no more knots can be added because the number of */
/*            b-spline coefficients (nx-kx-1)*(ny-ky-1) already exceeds */
/*            the number of data points m. */
/*            probably causes : either s or m too small. */
/*            the approximation returned is the weighted least-squares */
/*            spline according to the current set of knots. */
/*            the parameter fp gives the corresponding weighted sum of */
/*            squared residuals (fp>s). */
/*   ier=5  : error. no more knots can be added because the additional */
/*            knot would (quasi) coincide with an old one. */
/*            probably causes : s too small or too large a weight to an */
/*            inaccurate data point. */
/*            the approximation returned is the weighted least-squares */
/*            spline according to the current set of knots. */
/*            the parameter fp gives the corresponding weighted sum of */
/*            squared residuals (fp>s). */
/*   ier=10 : error. on entry, the input data are controlled on validity */
/*            the following restrictions must be satisfied. */
/*            -1<=iopt<=1, 1<=kx,ky<=5, m>=(kx+1)*(ky+1), nxest>=2*kx+2, */
/*            nyest>=2*ky+2, 0<eps<1, nmax>=nxest, nmax>=nyest, */
/*            xb<=x(i)<=xe, yb<=y(i)<=ye, w(i)>0, i=1,...,m */
/*            lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1 */
/*            kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1) */
/*            if iopt=-1: 2*kx+2<=nx<=nxest */
/*                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe */
/*                        2*ky+2<=ny<=nyest */
/*                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye */
/*            if iopt>=0: s>=0 */
/*            if one of these conditions is found to be violated,control */
/*            is immediately repassed to the calling program. in that */
/*            case there is no approximation returned. */
/*   ier>10 : error. lwrk2 is too small, i.e. there is not enough work- */
/*            space for computing the minimal least-squares solution of */
/*            a rank deficient system of linear equations. ier gives the */
/*            requested value for lwrk2. there is no approximation re- */
/*            turned but, having saved the information contained in nx, */
/*            ny,tx,ty,wrk1, and having adjusted the value of lwrk2 and */
/*            the dimension of the array wrk2 accordingly, the user can */
/*            continue at the point the program was left, by calling */
/*            surfit with iopt=1. */

/* further comments: */
/*  by means of the parameter s, the user can control the tradeoff */
/*   between closeness of fit and smoothness of fit of the approximation. */
/*   if s is too large, the spline will be too smooth and signal will be */
/*   lost ; if s is too small the spline will pick up too much noise. in */
/*   the extreme cases the program will return an interpolating spline if */
/*   s=0 and the weighted least-squares polynomial (degrees kx,ky)if s is */
/*   very large. between these extremes, a properly chosen s will result */
/*   in a good compromise between closeness of fit and smoothness of fit. */
/*   to decide whether an approximation, corresponding to a certain s is */
/*   satisfactory the user is highly recommended to inspect the fits */
/*   graphically. */
/*   recommended values for s depend on the weights w(i). if these are */
/*   taken as 1/d(i) with d(i) an estimate of the standard deviation of */
/*   z(i), a good s-value should be found in the range (m-sqrt(2*m),m+ */
/*   sqrt(2*m)). if nothing is known about the statistical error in z(i) */
/*   each w(i) can be set equal to one and s determined by trial and */
/*   error, taking account of the comments above. the best is then to */
/*   start with a very large value of s ( to determine the least-squares */
/*   polynomial and the corresponding upper bound fp0 for s) and then to */
/*   progressively decrease the value of s ( say by a factor 10 in the */
/*   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the */
/*   approximation shows more detail) to obtain closer fits. */
/*   to choose s very small is strongly discouraged. this considerably */
/*   increases computation time and memory requirements. it may also */
/*   cause rank-deficiency (ier<-2) and endager numerical stability. */
/*   to economize the search for a good s-value the program provides with */
/*   different modes of computation. at the first call of the routine, or */
/*   whenever he wants to restart with the initial set of knots the user */
/*   must set iopt=0. */
/*   if iopt=1 the program will continue with the set of knots found at */
/*   the last call of the routine. this will save a lot of computation */
/*   time if surfit is called repeatedly for different values of s. */
/*   the number of knots of the spline returned and their location will */
/*   depend on the value of s and on the complexity of the shape of the */
/*   function underlying the data. if the computation mode iopt=1 */
/*   is used, the knots returned may also depend on the s-values at */
/*   previous calls (if these were smaller). therefore, if after a number */
/*   of trials with different s-values and iopt=1, the user can finally */
/*   accept a fit as satisfactory, it may be worthwhile for him to call */
/*   surfit once more with the selected value for s but now with iopt=0. */
/*   indeed, surfit may then return an approximation of the same quality */
/*   of fit but with fewer knots and therefore better if data reduction */
/*   is also an important objective for the user. */
/*   the number of knots may also depend on the upper bounds nxest and */
/*   nyest. indeed, if at a certain stage in surfit the number of knots */
/*   in one direction (say nx) has reached the value of its upper bound */
/*   (nxest), then from that moment on all subsequent knots are added */
/*   in the other (y) direction. this may indicate that the value of */
/*   nxest is too small. on the other hand, it gives the user the option */
/*   of limiting the number of knots the routine locates in any direction */
/*   for example, by setting nxest=2*kx+2 (the lowest allowable value for */
/*   nxest), the user can indicate that he wants an approximation which */
/*   is a simple polynomial of degree kx in the variable x. */

/*  other subroutines required: */
/*    fpback,fpbspl,fpsurf,fpdisc,fpgivs,fprank,fprati,fprota,fporde */

/*  references: */
/*   dierckx p. : an algorithm for surface fitting with spline functions */
/*                ima j. numer. anal. 1 (1981) 267-283. */
/*   dierckx p. : an algorithm for surface fitting with spline functions */
/*                report tw50, dept. computer science,k.u.leuven, 1980. */
/*   dierckx p. : curve and surface fitting with splines, monographs on */
/*                numerical analysis, oxford university press, 1993. */

/*  author: */
/*    p.dierckx */
/*    dept. computer science, k.u. leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  creation date : may 1979 */
/*  latest update : march 1987 */

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fpsurf */
/*  .. */
/*  we set up the parameters tol and maxit. */
    /* Parameter adjustments */
    --w;
    --z__;
    --y;
    --x;
    --c__;
    --ty;
    --tx;
    --wrk1;
    --wrk2;
    --iwrk;

    /* Function Body */
    maxit = 20;
    tol = .001f;
/*  before starting computations a data check is made. if the input data */
/*  are invalid,control is immediately repassed to the calling program. */
    *ier = 10;
    if (*eps <= 0.f || *eps >= 1.f) {
	goto L70;
    }
    if (*kx <= 0 || *kx > 5) {
	goto L70;
    }
    kx1 = *kx + 1;
    if (*ky <= 0 || *ky > 5) {
	goto L70;
    }
    ky1 = *ky + 1;
    kmax = max(*kx,*ky);
    km1 = kmax + 1;
    km2 = km1 + 1;
    if (*iopt < -1 || *iopt > 1) {
	goto L70;
    }
    if (*m < kx1 * ky1) {
	goto L70;
    }
    nminx = kx1 << 1;
    if (*nxest < nminx || *nxest > *nmax) {
	goto L70;
    }
    nminy = ky1 << 1;
    if (*nyest < nminy || *nyest > *nmax) {
	goto L70;
    }
    nest = max(*nxest,*nyest);
    nxk = *nxest - kx1;
    nyk = *nyest - ky1;
    ncest = nxk * nyk;
    nmx = *nxest - nminx + 1;
    nmy = *nyest - nminy + 1;
    nrint = nmx + nmy;
    nreg = nmx * nmy;
    ib1 = *kx * nyk + ky1;
    jb1 = *ky * nxk + kx1;
    ib3 = kx1 * nyk + 1;
    if (ib1 <= jb1) {
	goto L10;
    }
    ib1 = jb1;
    ib3 = ky1 * nxk + 1;
L10:
    lwest = ncest * (ib1 + 2 + ib3) + (nrint + nest * km2 + *m * km1 << 1) + 
	    ib3;
    kwest = *m + nreg;
    if (*lwrk1 < lwest || *kwrk < kwest) {
	goto L70;
    }
    if (*xb >= *xe || *yb >= *ye) {
	goto L70;
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (w[i__] <= 0.f) {
	    goto L70;
	}
	if (x[i__] < *xb || x[i__] > *xe) {
	    goto L70;
	}
	if (y[i__] < *yb || y[i__] > *ye) {
	    goto L70;
	}
/* L20: */
    }
    if (*iopt >= 0) {
	goto L50;
    }
    if (*nx < nminx || *nx > *nxest) {
	goto L70;
    }
    nxk = *nx - kx1;
    tx[kx1] = *xb;
    tx[nxk + 1] = *xe;
    i__1 = nxk;
    for (i__ = kx1; i__ <= i__1; ++i__) {
	if (tx[i__ + 1] <= tx[i__]) {
	    goto L70;
	}
/* L30: */
    }
    if (*ny < nminy || *ny > *nyest) {
	goto L70;
    }
    nyk = *ny - ky1;
    ty[ky1] = *yb;
    ty[nyk + 1] = *ye;
    i__1 = nyk;
    for (i__ = ky1; i__ <= i__1; ++i__) {
	if (ty[i__ + 1] <= ty[i__]) {
	    goto L70;
	}
/* L40: */
    }
    goto L60;
L50:
    if (*s < 0.f) {
	goto L70;
    }
L60:
    *ier = 0;
/*  we partition the working space and determine the spline approximation */
    kn = 1;
    ki = kn + *m;
    lq = 2;
    la = lq + ncest * ib3;
    lf = la + ncest * ib1;
    lff = lf + ncest;
    lfp = lff + ncest;
    lco = lfp + nrint;
    lh = lco + nrint;
    lbx = lh + ib3;
    nek = nest * km2;
    lby = lbx + nek;
    lsx = lby + nek;
    lsy = lsx + *m * km1;
    fpsurf_(iopt, m, &x[1], &y[1], &z__[1], &w[1], xb, xe, yb, ye, kx, ky, s, 
	    nxest, nyest, eps, &tol, &maxit, &nest, &km1, &km2, &ib1, &ib3, &
	    ncest, &nrint, &nreg, nx, &tx[1], ny, &ty[1], &c__[1], fp, &wrk1[
	    1], &wrk1[lfp], &wrk1[lco], &wrk1[lf], &wrk1[lff], &wrk1[la], &
	    wrk1[lq], &wrk1[lbx], &wrk1[lby], &wrk1[lsx], &wrk1[lsy], &wrk1[
	    lh], &iwrk[ki], &iwrk[kn], &wrk2[1], lwrk2, ier);
L70:
    return 0;
} /* surfit_ */

/* Subroutine */ int fpback_(doublereal *a, doublereal *z__, integer *n, 
	integer *k, doublereal *c__, integer *nest)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l, m, i1, k1;
    static doublereal store;

/*  subroutine fpback calculates the solution of the system of */
/*  equations a*c = z with a a n x n upper triangular matrix */
/*  of bandwidth k. */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
    /* Parameter adjustments */
    --c__;
    --z__;
    a_dim1 = *nest;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    k1 = *k - 1;
    c__[*n] = z__[*n] / a[*n + a_dim1];
    i__ = *n - 1;
    if (i__ == 0) {
	goto L30;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	store = z__[i__];
	i1 = k1;
	if (j <= k1) {
	    i1 = j - 1;
	}
	m = i__;
	i__2 = i1;
	for (l = 1; l <= i__2; ++l) {
	    ++m;
	    store -= c__[m] * a[i__ + (l + 1) * a_dim1];
/* L10: */
	}
	c__[i__] = store / a[i__ + a_dim1];
	--i__;
/* L20: */
    }
L30:
    return 0;
} /* fpback_ */

/* Subroutine */ int fpsurf_(integer *iopt, integer *m, doublereal *x, 
	doublereal *y, doublereal *z__, doublereal *w, doublereal *xb, 
	doublereal *xe, doublereal *yb, doublereal *ye, integer *kxx, integer 
	*kyy, doublereal *s, integer *nxest, integer *nyest, doublereal *eta, 
	doublereal *tol, integer *maxit, integer *nmax, integer *km1, integer 
	*km2, integer *ib1, integer *ib3, integer *nc, integer *intest, 
	integer *nrest, integer *nx0, doublereal *tx, integer *ny0, 
	doublereal *ty, doublereal *c__, doublereal *fp, doublereal *fp0, 
	doublereal *fpint, doublereal *coord, doublereal *f, doublereal *ff, 
	doublereal *a, doublereal *q, doublereal *bx, doublereal *by, 
	doublereal *spx, doublereal *spy, doublereal *h__, integer *index, 
	integer *nummer, doublereal *wrk, integer *lwrk, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, bx_dim1, bx_offset, by_dim1, 
	    by_offset, spx_dim1, spx_offset, spy_dim1, spy_offset, i__1, i__2,
	     i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, l, n;
    static doublereal p, f1, f2, f3;
    static integer i1, i2, i3, j1, l1, l2, n1;
    static doublereal p1, p2, p3, x0, x1, y0, y1;
    static integer la, ii, lf, lh, in;
    static doublereal wi, rn, hx[6], zi, sq;
    static integer kx, ky, lx, ly, nx, ny;
    static doublereal hy[6];
    static integer kx1, kx2, ky1, ky2;
    static doublereal acc;
    static integer ibb;
    static doublereal arg, one, cos__, ten, eps, hxi, sin__;
    static integer nxe, nye;
    static doublereal piv;
    static integer num;
    static doublereal fac1, fac2;
    static integer jxy, nxx, nyy, ich1, ich3;
    static doublereal con1, con4, con9;
    static integer num1, nk1x, nk1y;
    static doublereal half;
    static integer ncof;
    static doublereal dmax__;
    static integer nreg, rank, iter;
    static doublereal fpms, pinv;
    static integer irot, jrot, iband;
    static doublereal sigma, fpmax;
    static integer nminx, nminy;
    static doublereal store;
    static integer nrint, iband1, lwest, iband3, iband4;
    extern /* Subroutine */ int fpback_(doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *);
    static integer ichang;
    extern /* Subroutine */ int fpdisc_(doublereal *, integer *, integer *, 
	    doublereal *, integer *), fporde_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *), 
	    fprank_(doublereal *, doublereal *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal fprati_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *), fprota_(doublereal *, 
	    doublereal *, doublereal *, doublereal *), fpgivs_(doublereal *, 
	    doublereal *, doublereal *, doublereal *);

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota */
/*  .. */
/*  set constants */
    /* Parameter adjustments */
    --nummer;
    --w;
    --z__;
    --y;
    --x;
    --ty;
    --tx;
    spy_dim1 = *m;
    spy_offset = 1 + spy_dim1;
    spy -= spy_offset;
    spx_dim1 = *m;
    spx_offset = 1 + spx_dim1;
    spx -= spx_offset;
    by_dim1 = *nmax;
    by_offset = 1 + by_dim1;
    by -= by_offset;
    bx_dim1 = *nmax;
    bx_offset = 1 + bx_dim1;
    bx -= bx_offset;
    --h__;
    q_dim1 = *nc;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    a_dim1 = *nc;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ff;
    --f;
    --c__;
    --coord;
    --fpint;
    --index;
    --wrk;

    /* Function Body */
    one = 1.f;
    con1 = .1f;
    con9 = .9f;
    con4 = .04f;
    half = .5f;
    ten = 10.f;
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 1: determination of the number of knots and their position.     c */
/* ****************************************************************     c */
/* given a set of knots we compute the least-squares spline sinf(x,y),  c */
/* and the corresponding weighted sum of squared residuals fp=f(p=inf). c */
/* if iopt=-1  sinf(x,y) is the requested approximation.                c */
/* if iopt=0 or iopt=1 we check whether we can accept the knots:        c */
/*   if fp <=s we will continue with the current set of knots.          c */
/*   if fp > s we will increase the number of knots and compute the     c */
/*      corresponding least-squares spline until finally  fp<=s.        c */
/* the initial choice of knots depends on the value of s and iopt.      c */
/*   if iopt=0 we first compute the least-squares polynomial of degree  c */
/*     kx in x and ky in y; nx=nminx=2*kx+2 and ny=nminy=2*ky+2.        c */
/*     fp0=f(0) denotes the corresponding weighted sum of squared       c */
/*     residuals                                                        c */
/*   if iopt=1 we start with the knots found at the last call of the    c */
/*     routine, except for the case that s>=fp0; then we can compute    c */
/*     the least-squares polynomial directly.                           c */
/* eventually the independent variables x and y (and the corresponding  c */
/* parameters) will be switched if this can reduce the bandwidth of the c */
/* system to be solved.                                                 c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  ichang denotes whether(1) or not(-1) the directions have been inter- */
/*  changed. */
    ichang = -1;
    x0 = *xb;
    x1 = *xe;
    y0 = *yb;
    y1 = *ye;
    kx = *kxx;
    ky = *kyy;
    kx1 = kx + 1;
    ky1 = ky + 1;
    nxe = *nxest;
    nye = *nyest;
    eps = sqrt(*eta);
    if (*iopt < 0) {
	goto L20;
    }
/*  calculation of acc, the absolute tolerance for the root of f(p)=s. */
    acc = *tol * *s;
    if (*iopt == 0) {
	goto L10;
    }
    if (*fp0 > *s) {
	goto L20;
    }
/*  initialization for the least-squares polynomial. */
L10:
    nminx = kx1 << 1;
    nminy = ky1 << 1;
    nx = nminx;
    ny = nminy;
    *ier = -2;
    goto L30;
L20:
    nx = *nx0;
    ny = *ny0;
/*  main loop for the different sets of knots. m is a save upper bound */
/*  for the number of trials. */
L30:
    i__1 = *m;
    for (iter = 1; iter <= i__1; ++iter) {
/*  find the position of the additional knots which are needed for the */
/*  b-spline representation of s(x,y). */
	l = nx;
	i__2 = kx1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    tx[i__] = x0;
	    tx[l] = x1;
	    --l;
/* L40: */
	}
	l = ny;
	i__2 = ky1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ty[i__] = y0;
	    ty[l] = y1;
	    --l;
/* L50: */
	}
/*  find nrint, the total number of knot intervals and nreg, the number */
/*  of panels in which the approximation domain is subdivided by the */
/*  intersection of knots. */
	nxx = nx - (kx1 << 1) + 1;
	nyy = ny - (ky1 << 1) + 1;
	nrint = nxx + nyy;
	nreg = nxx * nyy;
/*  find the bandwidth of the observation matrix a. */
/*  if necessary, interchange the variables x and y, in order to obtain */
/*  a minimal bandwidth. */
	iband1 = kx * (ny - ky1) + ky;
	l = ky * (nx - kx1) + kx;
	if (iband1 <= l) {
	    goto L130;
	}
	iband1 = l;
	ichang = -ichang;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    store = x[i__];
	    x[i__] = y[i__];
	    y[i__] = store;
/* L60: */
	}
	store = x0;
	x0 = y0;
	y0 = store;
	store = x1;
	x1 = y1;
	y1 = store;
	n = min(nx,ny);
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    store = tx[i__];
	    tx[i__] = ty[i__];
	    ty[i__] = store;
/* L70: */
	}
	n1 = n + 1;
	if ((i__2 = nx - ny) < 0) {
	    goto L80;
	} else if (i__2 == 0) {
	    goto L120;
	} else {
	    goto L100;
	}
L80:
	i__2 = ny;
	for (i__ = n1; i__ <= i__2; ++i__) {
	    tx[i__] = ty[i__];
/* L90: */
	}
	goto L120;
L100:
	i__2 = nx;
	for (i__ = n1; i__ <= i__2; ++i__) {
	    ty[i__] = tx[i__];
/* L110: */
	}
L120:
	l = nx;
	nx = ny;
	ny = l;
	l = nxe;
	nxe = nye;
	nye = l;
	l = nxx;
	nxx = nyy;
	nyy = l;
	l = kx;
	kx = ky;
	ky = l;
	kx1 = kx + 1;
	ky1 = ky + 1;
L130:
	iband = iband1 + 1;
/*  arrange the data points according to the panel they belong to. */
	fporde_(&x[1], &y[1], m, &kx, &ky, &tx[1], &nx, &ty[1], &ny, &nummer[
		1], &index[1], &nreg);
/*  find ncof, the number of b-spline coefficients. */
	nk1x = nx - kx1;
	nk1y = ny - ky1;
	ncof = nk1x * nk1y;
/*  initialize the observation matrix a. */
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    f[i__] = 0.f;
	    i__3 = iband;
	    for (j = 1; j <= i__3; ++j) {
		a[i__ + j * a_dim1] = 0.f;
/* L140: */
	    }
	}
/*  initialize the sum of squared residuals. */
	*fp = 0.f;
/*  fetch the data points in the new order. main loop for the */
/*  different panels. */
	i__3 = nreg;
	for (num = 1; num <= i__3; ++num) {
/*  fix certain constants for the current panel; jrot records the column */
/*  number of the first non-zero element in a row of the observation */
/*  matrix according to a data point of the panel. */
	    num1 = num - 1;
	    lx = num1 / nyy;
	    l1 = lx + kx1;
	    ly = num1 - lx * nyy;
	    l2 = ly + ky1;
	    jrot = lx * nk1y + ly;
/*  test whether there are still data points in the panel. */
	    in = index[num];
L150:
	    if (in == 0) {
		goto L250;
	    }
/*  fetch a new data point. */
	    wi = w[in];
	    zi = z__[in] * wi;
/*  evaluate for the x-direction, the (kx+1) non-zero b-splines at x(in). */
	    fpbspl_(&tx[1], &nx, &kx, &x[in], &l1, hx);
/*  evaluate for the y-direction, the (ky+1) non-zero b-splines at y(in). */
	    fpbspl_(&ty[1], &ny, &ky, &y[in], &l2, hy);
/*  store the value of these b-splines in spx and spy respectively. */
	    i__2 = kx1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		spx[in + i__ * spx_dim1] = hx[i__ - 1];
/* L160: */
	    }
	    i__2 = ky1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		spy[in + i__ * spy_dim1] = hy[i__ - 1];
/* L170: */
	    }
/*  initialize the new row of observation matrix. */
	    i__2 = iband;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		h__[i__] = 0.f;
/* L180: */
	    }
/*  calculate the non-zero elements of the new row by making the cross */
/*  products of the non-zero b-splines in x- and y-direction. */
	    i1 = 0;
	    i__2 = kx1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		hxi = hx[i__ - 1];
		j1 = i1;
		i__4 = ky1;
		for (j = 1; j <= i__4; ++j) {
		    ++j1;
		    h__[j1] = hxi * hy[j - 1] * wi;
/* L190: */
		}
		i1 += nk1y;
/* L200: */
	    }
/*  rotate the row into triangle by givens transformations . */
	    irot = jrot;
	    i__2 = iband;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		++irot;
		piv = h__[i__];
		if (piv == 0.f) {
		    goto L220;
		}
/*  calculate the parameters of the givens transformation. */
		fpgivs_(&piv, &a[irot + a_dim1], &cos__, &sin__);
/*  apply that transformation to the right hand side. */
		fprota_(&cos__, &sin__, &zi, &f[irot]);
		if (i__ == iband) {
		    goto L230;
		}
/*  apply that transformation to the left hand side. */
		i2 = 1;
		i3 = i__ + 1;
		i__4 = iband;
		for (j = i3; j <= i__4; ++j) {
		    ++i2;
		    fprota_(&cos__, &sin__, &h__[j], &a[irot + i2 * a_dim1]);
/* L210: */
		}
L220:
		;
	    }
/*  add the contribution of the row to the sum of squares of residual */
/*  right hand sides. */
L230:
/* Computing 2nd power */
	    d__1 = zi;
	    *fp += d__1 * d__1;
/*  find the number of the next data point in the panel. */
/* L240: */
	    in = nummer[in];
	    goto L150;
L250:
	    ;
	}
/*  find dmax, the maximum value for the diagonal elements in the reduced */
/*  triangle. */
	dmax__ = 0.f;
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (a[i__ + a_dim1] <= dmax__) {
		goto L260;
	    }
	    dmax__ = a[i__ + a_dim1];
L260:
	    ;
	}
/*  check whether the observation matrix is rank deficient. */
	sigma = eps * dmax__;
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (a[i__ + a_dim1] <= sigma) {
		goto L280;
	    }
/* L270: */
	}
/*  backward substitution in case of full rank. */
	fpback_(&a[a_offset], &f[1], &ncof, &iband, &c__[1], nc);
	rank = ncof;
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + q_dim1] = a[i__ + a_dim1] / dmax__;
/* L275: */
	}
	goto L300;
/*  in case of rank deficiency, find the minimum norm solution. */
/*  check whether there is sufficient working space */
L280:
	lwest = ncof * iband + ncof + iband;
	if (*lwrk < lwest) {
	    goto L780;
	}
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    ff[i__] = f[i__];
	    i__2 = iband;
	    for (j = 1; j <= i__2; ++j) {
		q[i__ + j * q_dim1] = a[i__ + j * a_dim1];
/* L290: */
	    }
	}
	lf = 1;
	lh = lf + ncof;
	la = lh + iband;
	fprank_(&q[q_offset], &ff[1], &ncof, &iband, nc, &sigma, &c__[1], &sq,
		 &rank, &wrk[la], &wrk[lf], &wrk[lh]);
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    q[i__ + q_dim1] /= dmax__;
/* L295: */
	}
/*  add to the sum of squared residuals, the contribution of reducing */
/*  the rank. */
	*fp += sq;
L300:
	if (*ier == -2) {
	    *fp0 = *fp;
	}
/*  test whether the least-squares spline is an acceptable solution. */
	if (*iopt < 0) {
	    goto L820;
	}
	fpms = *fp - *s;
	if (abs(fpms) <= acc) {
	    if (*fp <= 0.) {
		goto L815;
	    } else {
		goto L820;
	    }
	}
/*  test whether we can accept the choice of knots. */
	if (fpms < 0.f) {
	    goto L430;
	}
/*  test whether we cannot further increase the number of knots. */
	if (ncof > *m) {
	    goto L790;
	}
	*ier = 0;
/*  search where to add a new knot. */
/*  find for each interval the sum of squared residuals fpint for the */
/*  data points having the coordinate belonging to that knot interval. */
/*  calculate also coord which is the same sum, weighted by the position */
/*  of the data points considered. */
/* L310: */
	i__2 = nrint;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    fpint[i__] = 0.f;
	    coord[i__] = 0.f;
/* L320: */
	}
	i__2 = nreg;
	for (num = 1; num <= i__2; ++num) {
	    num1 = num - 1;
	    lx = num1 / nyy;
	    l1 = lx + 1;
	    ly = num1 - lx * nyy;
	    l2 = ly + 1 + nxx;
	    jrot = lx * nk1y + ly;
	    in = index[num];
L330:
	    if (in == 0) {
		goto L360;
	    }
	    store = 0.f;
	    i1 = jrot;
	    i__3 = kx1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		hxi = spx[in + i__ * spx_dim1];
		j1 = i1;
		i__4 = ky1;
		for (j = 1; j <= i__4; ++j) {
		    ++j1;
		    store += hxi * spy[in + j * spy_dim1] * c__[j1];
/* L340: */
		}
		i1 += nk1y;
/* L350: */
	    }
/* Computing 2nd power */
	    d__1 = w[in] * (z__[in] - store);
	    store = d__1 * d__1;
	    fpint[l1] += store;
	    coord[l1] += store * x[in];
	    fpint[l2] += store;
	    coord[l2] += store * y[in];
	    in = nummer[in];
	    goto L330;
L360:
	    ;
	}
/*  find the interval for which fpint is maximal on the condition that */
/*  there still can be added a knot. */
L370:
	l = 0;
	fpmax = 0.f;
	l1 = 1;
	l2 = nrint;
	if (nx == nxe) {
	    l1 = nxx + 1;
	}
	if (ny == nye) {
	    l2 = nxx;
	}
	if (l1 > l2) {
	    goto L810;
	}
	i__2 = l2;
	for (i__ = l1; i__ <= i__2; ++i__) {
	    if (fpmax >= fpint[i__]) {
		goto L380;
	    }
	    l = i__;
	    fpmax = fpint[i__];
L380:
	    ;
	}
/*  test whether we cannot further increase the number of knots. */
	if (l == 0) {
	    goto L785;
	}
/*  calculate the position of the new knot. */
	arg = coord[l] / fpint[l];
/*  test in what direction the new knot is going to be added. */
	if (l > nxx) {
	    goto L400;
	}
/*  addition in the x-direction. */
	jxy = l + kx1;
	fpint[l] = 0.f;
	fac1 = tx[jxy] - arg;
	fac2 = arg - tx[jxy - 1];
	if (fac1 > ten * fac2 || fac2 > ten * fac1) {
	    goto L370;
	}
	j = nx;
	i__2 = nx;
	for (i__ = jxy; i__ <= i__2; ++i__) {
	    tx[j + 1] = tx[j];
	    --j;
/* L390: */
	}
	tx[jxy] = arg;
	++nx;
	goto L420;
/*  addition in the y-direction. */
L400:
	jxy = l + ky1 - nxx;
	fpint[l] = 0.f;
	fac1 = ty[jxy] - arg;
	fac2 = arg - ty[jxy - 1];
	if (fac1 > ten * fac2 || fac2 > ten * fac1) {
	    goto L370;
	}
	j = ny;
	i__2 = ny;
	for (i__ = jxy; i__ <= i__2; ++i__) {
	    ty[j + 1] = ty[j];
	    --j;
/* L410: */
	}
	ty[jxy] = arg;
	++ny;
/*  restart the computations with the new set of knots. */
L420:
	;
    }
/*  test whether the least-squares polynomial is a solution of our */
/*  approximation problem. */
L430:
    if (*ier == -2) {
	goto L830;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 2: determination of the smoothing spline sp(x,y)                c */
/* *****************************************************                c */
/* we have determined the number of knots and their position. we now    c */
/* compute the b-spline coefficients of the smoothing spline sp(x,y).   c */
/* the observation matrix a is extended by the rows of a matrix,        c */
/* expressing that sp(x,y) must be a polynomial of degree kx in x and   c */
/* ky in y. the corresponding weights of these additional rows are set  c */
/* to 1./p.  iteratively we than have to determine the value of p       c */
/* such that f(p)=sum((w(i)*(z(i)-sp(x(i),y(i))))**2) be = s.           c */
/* we already know that the least-squares polynomial corresponds to     c */
/* p=0  and that the least-squares spline corresponds to p=infinity.    c */
/* the iteration process which is proposed here makes use of rational   c */
/* interpolation. since f(p) is a convex and strictly decreasing        c */
/* function of p, it can be approximated by a rational function r(p)=   c */
/* (u*p+v)/(p+w). three values of p(p1,p2,p3) with corresponding values c */
/* of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the c */
/* new value of p such that r(p)=s. convergence is guaranteed by taking c */
/* f1 > 0 and f3 < 0.                                                   c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    kx2 = kx1 + 1;
/*  test whether there are interior knots in the x-direction. */
    if (nk1x == kx1) {
	goto L440;
    }
/*  evaluate the discotinuity jumps of the kx-th order derivative of */
/*  the b-splines at the knots tx(l),l=kx+2,...,nx-kx-1. */
    fpdisc_(&tx[1], &nx, &kx2, &bx[bx_offset], nmax);
L440:
    ky2 = ky1 + 1;
/*  test whether there are interior knots in the y-direction. */
    if (nk1y == ky1) {
	goto L450;
    }
/*  evaluate the discontinuity jumps of the ky-th order derivative of */
/*  the b-splines at the knots ty(l),l=ky+2,...,ny-ky-1. */
    fpdisc_(&ty[1], &ny, &ky2, &by[by_offset], nmax);
/*  initial value for p. */
L450:
    p1 = 0.f;
    f1 = *fp0 - *s;
    p3 = -one;
    f3 = fpms;
    p = 0.f;
    i__1 = ncof;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p += a[i__ + a_dim1];
/* L460: */
    }
    rn = (doublereal) ncof;
    p = rn / p;
/*  find the bandwidth of the extended observation matrix. */
    iband3 = kx1 * nk1y;
    iband4 = iband3 + 1;
    ich1 = 0;
    ich3 = 0;
/*  iteration process to find the root of f(p)=s. */
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
	pinv = one / p;
/*  store the triangularized observation matrix into q. */
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ff[i__] = f[i__];
	    i__3 = iband;
	    for (j = 1; j <= i__3; ++j) {
		q[i__ + j * q_dim1] = a[i__ + j * a_dim1];
/* L470: */
	    }
	    ibb = iband + 1;
	    i__3 = iband4;
	    for (j = ibb; j <= i__3; ++j) {
		q[i__ + j * q_dim1] = 0.f;
/* L480: */
	    }
	}
	if (nk1y == ky1) {
	    goto L560;
	}
/*  extend the observation matrix with the rows of a matrix, expressing */
/*  that for x=cst. sp(x,y) must be a polynomial in y of degree ky. */
	i__3 = nk1y;
	for (i__ = ky2; i__ <= i__3; ++i__) {
	    ii = i__ - ky1;
	    i__2 = nk1x;
	    for (j = 1; j <= i__2; ++j) {
/*  initialize the new row. */
		i__4 = iband;
		for (l = 1; l <= i__4; ++l) {
		    h__[l] = 0.f;
/* L490: */
		}
/*  fill in the non-zero elements of the row. jrot records the column */
/*  number of the first non-zero element in the row. */
		i__4 = ky2;
		for (l = 1; l <= i__4; ++l) {
		    h__[l] = by[ii + l * by_dim1] * pinv;
/* L500: */
		}
		zi = 0.f;
		jrot = (j - 1) * nk1y + ii;
/*  rotate the new row into triangle by givens transformations without */
/*  square roots. */
		i__4 = ncof;
		for (irot = jrot; irot <= i__4; ++irot) {
		    piv = h__[1];
/* Computing MIN */
		    i__5 = iband1, i__6 = ncof - irot;
		    i2 = min(i__5,i__6);
		    if (piv == 0.f) {
			if (i2 <= 0) {
			    goto L550;
			} else {
			    goto L520;
			}
		    }
/*  calculate the parameters of the givens transformation. */
		    fpgivs_(&piv, &q[irot + q_dim1], &cos__, &sin__);
/*  apply that givens transformation to the right hand side. */
		    fprota_(&cos__, &sin__, &zi, &ff[irot]);
		    if (i2 == 0) {
			goto L550;
		    }
/*  apply that givens transformation to the left hand side. */
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
			l1 = l + 1;
			fprota_(&cos__, &sin__, &h__[l1], &q[irot + l1 * 
				q_dim1]);
/* L510: */
		    }
L520:
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
			h__[l] = h__[l + 1];
/* L530: */
		    }
		    h__[i2 + 1] = 0.f;
/* L540: */
		}
L550:
		;
	    }
	}
L560:
	if (nk1x == kx1) {
	    goto L640;
	}
/*  extend the observation matrix with the rows of a matrix expressing */
/*  that for y=cst. sp(x,y) must be a polynomial in x of degree kx. */
	i__2 = nk1x;
	for (i__ = kx2; i__ <= i__2; ++i__) {
	    ii = i__ - kx1;
	    i__3 = nk1y;
	    for (j = 1; j <= i__3; ++j) {
/*  initialize the new row */
		i__4 = iband4;
		for (l = 1; l <= i__4; ++l) {
		    h__[l] = 0.f;
/* L570: */
		}
/*  fill in the non-zero elements of the row. jrot records the column */
/*  number of the first non-zero element in the row. */
		j1 = 1;
		i__4 = kx2;
		for (l = 1; l <= i__4; ++l) {
		    h__[j1] = bx[ii + l * bx_dim1] * pinv;
		    j1 += nk1y;
/* L580: */
		}
		zi = 0.f;
		jrot = (i__ - kx2) * nk1y + j;
/*  rotate the new row into triangle by givens transformations . */
		i__4 = ncof;
		for (irot = jrot; irot <= i__4; ++irot) {
		    piv = h__[1];
/* Computing MIN */
		    i__5 = iband3, i__6 = ncof - irot;
		    i2 = min(i__5,i__6);
		    if (piv == 0.f) {
			if (i2 <= 0) {
			    goto L630;
			} else {
			    goto L600;
			}
		    }
/*  calculate the parameters of the givens transformation. */
		    fpgivs_(&piv, &q[irot + q_dim1], &cos__, &sin__);
/*  apply that givens transformation to the right hand side. */
		    fprota_(&cos__, &sin__, &zi, &ff[irot]);
		    if (i2 == 0) {
			goto L630;
		    }
/*  apply that givens transformation to the left hand side. */
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
			l1 = l + 1;
			fprota_(&cos__, &sin__, &h__[l1], &q[irot + l1 * 
				q_dim1]);
/* L590: */
		    }
L600:
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
			h__[l] = h__[l + 1];
/* L610: */
		    }
		    h__[i2 + 1] = 0.f;
/* L620: */
		}
L630:
		;
	    }
	}
/*  find dmax, the maximum value for the diagonal elements in the */
/*  reduced triangle. */
L640:
	dmax__ = 0.f;
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (q[i__ + q_dim1] <= dmax__) {
		goto L650;
	    }
	    dmax__ = q[i__ + q_dim1];
L650:
	    ;
	}
/*  check whether the matrix is rank deficient. */
	sigma = eps * dmax__;
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (q[i__ + q_dim1] <= sigma) {
		goto L670;
	    }
/* L660: */
	}
/*  backward substitution in case of full rank. */
	fpback_(&q[q_offset], &ff[1], &ncof, &iband4, &c__[1], nc);
	rank = ncof;
	goto L675;
/*  in case of rank deficiency, find the minimum norm solution. */
L670:
	lwest = ncof * iband4 + ncof + iband4;
	if (*lwrk < lwest) {
	    goto L780;
	}
	lf = 1;
	lh = lf + ncof;
	la = lh + iband4;
	fprank_(&q[q_offset], &ff[1], &ncof, &iband4, nc, &sigma, &c__[1], &
		sq, &rank, &wrk[la], &wrk[lf], &wrk[lh]);
L675:
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + q_dim1] /= dmax__;
/* L680: */
	}
/*  compute f(p). */
	*fp = 0.f;
	i__3 = nreg;
	for (num = 1; num <= i__3; ++num) {
	    num1 = num - 1;
	    lx = num1 / nyy;
	    ly = num1 - lx * nyy;
	    jrot = lx * nk1y + ly;
	    in = index[num];
L690:
	    if (in == 0) {
		goto L720;
	    }
	    store = 0.f;
	    i1 = jrot;
	    i__2 = kx1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		hxi = spx[in + i__ * spx_dim1];
		j1 = i1;
		i__4 = ky1;
		for (j = 1; j <= i__4; ++j) {
		    ++j1;
		    store += hxi * spy[in + j * spy_dim1] * c__[j1];
/* L700: */
		}
		i1 += nk1y;
/* L710: */
	    }
/* Computing 2nd power */
	    d__1 = w[in] * (z__[in] - store);
	    *fp += d__1 * d__1;
	    in = nummer[in];
	    goto L690;
L720:
	    ;
	}
/*  test whether the approximation sp(x,y) is an acceptable solution. */
	fpms = *fp - *s;
	if (abs(fpms) <= acc) {
	    goto L820;
	}
/*  test whether the maximum allowable number of iterations has been */
/*  reached. */
	if (iter == *maxit) {
	    goto L795;
	}
/*  carry out one more step of the iteration process. */
	p2 = p;
	f2 = fpms;
	if (ich3 != 0) {
	    goto L740;
	}
	if (f2 - f3 > acc) {
	    goto L730;
	}
/*  our initial choice of p is too large. */
	p3 = p2;
	f3 = f2;
	p *= con4;
	if (p <= p1) {
	    p = p1 * con9 + p2 * con1;
	}
	goto L770;
L730:
	if (f2 < 0.f) {
	    ich3 = 1;
	}
L740:
	if (ich1 != 0) {
	    goto L760;
	}
	if (f1 - f2 > acc) {
	    goto L750;
	}
/*  our initial choice of p is too small */
	p1 = p2;
	f1 = f2;
	p /= con4;
	if (p3 < 0.f) {
	    goto L770;
	}
	if (p >= p3) {
	    p = p2 * con1 + p3 * con9;
	}
	goto L770;
L750:
	if (f2 > 0.f) {
	    ich1 = 1;
	}
/*  test whether the iteration process proceeds as theoretically */
/*  expected. */
L760:
	if (f2 >= f1 || f2 <= f3) {
	    goto L800;
	}
/*  find the new value of p. */
	p = fprati_(&p1, &f1, &p2, &f2, &p3, &f3);
L770:
	;
    }
/*  error codes and messages. */
L780:
    *ier = lwest;
    goto L830;
L785:
    *ier = 5;
    goto L830;
L790:
    *ier = 4;
    goto L830;
L795:
    *ier = 3;
    goto L830;
L800:
    *ier = 2;
    goto L830;
L810:
    *ier = 1;
    goto L830;
L815:
    *ier = -1;
    *fp = 0.f;
L820:
    if (ncof != rank) {
	*ier = -rank;
    }
/*  test whether x and y are in the original order. */
L830:
    if (ichang < 0) {
	goto L930;
    }
/*  if not, interchange x and y once more. */
    l1 = 1;
    i__1 = nk1x;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l2 = i__;
	i__3 = nk1y;
	for (j = 1; j <= i__3; ++j) {
	    f[l2] = c__[l1];
	    ++l1;
	    l2 += nk1x;
/* L840: */
	}
    }
    i__3 = ncof;
    for (i__ = 1; i__ <= i__3; ++i__) {
	c__[i__] = f[i__];
/* L850: */
    }
    i__3 = *m;
    for (i__ = 1; i__ <= i__3; ++i__) {
	store = x[i__];
	x[i__] = y[i__];
	y[i__] = store;
/* L860: */
    }
    n = min(nx,ny);
    i__3 = n;
    for (i__ = 1; i__ <= i__3; ++i__) {
	store = tx[i__];
	tx[i__] = ty[i__];
	ty[i__] = store;
/* L870: */
    }
    n1 = n + 1;
    if ((i__3 = nx - ny) < 0) {
	goto L880;
    } else if (i__3 == 0) {
	goto L920;
    } else {
	goto L900;
    }
L880:
    i__3 = ny;
    for (i__ = n1; i__ <= i__3; ++i__) {
	tx[i__] = ty[i__];
/* L890: */
    }
    goto L920;
L900:
    i__3 = nx;
    for (i__ = n1; i__ <= i__3; ++i__) {
	ty[i__] = tx[i__];
/* L910: */
    }
L920:
    l = nx;
    nx = ny;
    ny = l;
L930:
    if (*iopt < 0) {
	goto L940;
    }
    *nx0 = nx;
    *ny0 = ny;
L940:
    return 0;
} /* fpsurf_ */

/* Subroutine */ int fpdisc_(doublereal *t, integer *n, integer *k2, 
	doublereal *b, integer *nest)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal h__[12];
    static integer i__, j, k, l, k1;
    static doublereal an;
    static integer ik, jk, lj, lk, lp, nk1;
    static doublereal fac;
    static integer lmk;
    static doublereal prod;
    static integer nrint;

/*  subroutine fpdisc calculates the discontinuity jumps of the kth */
/*  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1) */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local array.. */
/*  .. */
    /* Parameter adjustments */
    --t;
    b_dim1 = *nest;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    k1 = *k2 - 1;
    k = k1 - 1;
    nk1 = *n - k1;
    nrint = nk1 - k;
    an = (doublereal) nrint;
    fac = an / (t[nk1 + 1] - t[k1]);
    i__1 = nk1;
    for (l = *k2; l <= i__1; ++l) {
	lmk = l - k1;
	i__2 = k1;
	for (j = 1; j <= i__2; ++j) {
	    ik = j + k1;
	    lj = l + j;
	    lk = lj - *k2;
	    h__[j - 1] = t[l] - t[lk];
	    h__[ik - 1] = t[l] - t[lj];
/* L10: */
	}
	lp = lmk;
	i__2 = *k2;
	for (j = 1; j <= i__2; ++j) {
	    jk = j;
	    prod = h__[j - 1];
	    i__3 = k;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		++jk;
		prod = prod * h__[jk - 1] * fac;
/* L20: */
	    }
	    lk = lp + k1;
	    b[lmk + j * b_dim1] = (t[lk] - t[lp]) / prod;
	    ++lp;
/* L30: */
	}
/* L40: */
    }
    return 0;
} /* fpdisc_ */

/* Subroutine */ int fpgivs_(doublereal *piv, doublereal *ww, doublereal *
	cos__, doublereal *sin__)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal dd, one, store;

/*  subroutine fpgivs calculates the parameters of a givens */
/*  transformation . */
/*  .. */
/*  ..scalar arguments.. */
/*  ..local scalars.. */
/*  ..function references.. */
/*  .. */
    one = 1.f;
    store = abs(*piv);
    if (store >= *ww) {
/* Computing 2nd power */
	d__1 = *ww / *piv;
	dd = store * sqrt(one + d__1 * d__1);
    }
    if (store < *ww) {
/* Computing 2nd power */
	d__1 = *piv / *ww;
	dd = *ww * sqrt(one + d__1 * d__1);
    }
    *cos__ = *ww / dd;
    *sin__ = *piv / dd;
    *ww = dd;
    return 0;
} /* fpgivs_ */

/* Subroutine */ int fprank_(doublereal *a, doublereal *f, integer *n, 
	integer *m, integer *na, doublereal *tol, doublereal *c__, doublereal 
	*sq, integer *rank, doublereal *aa, doublereal *ff, doublereal *h__)
{
    /* System generated locals */
    integer a_dim1, a_offset, aa_dim1, aa_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, i1, i2, j1, j2, j3, m1, ii, ij, jj, kk, nl;
    static doublereal yi, fac, cos__, sin__, piv, stor1, stor2, stor3, store;
    extern /* Subroutine */ int fprota_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), fpgivs_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/*  subroutine fprank finds the minimum norm solution of a least- */
/*  squares problem in case of rank deficiency. */

/*  input parameters: */
/*    a : array, which contains the non-zero elements of the observation */
/*        matrix after triangularization by givens transformations. */
/*    f : array, which contains the transformed right hand side. */
/*    n : integer,wich contains the dimension of a. */
/*    m : integer, which denotes the bandwidth of a. */
/*  tol : real value, giving a threshold to determine the rank of a. */

/*  output parameters: */
/*    c : array, which contains the minimum norm solution. */
/*   sq : real value, giving the contribution of reducing the rank */
/*        to the sum of squared residuals. */
/* rank : integer, which contains the rank of matrix a. */

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fpgivs,fprota */
/*  .. */
    /* Parameter adjustments */
    --ff;
    --c__;
    --f;
    --h__;
    aa_dim1 = *n;
    aa_offset = 1 + aa_dim1;
    aa -= aa_offset;
    a_dim1 = *na;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    m1 = *m - 1;
/*  the rank deficiency nl is considered to be the number of sufficient */
/*  small diagonal elements of a. */
    nl = 0;
    *sq = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (a[i__ + a_dim1] > *tol) {
	    goto L90;
	}
/*  if a sufficient small diagonal element is found, we put it to */
/*  zero. the remainder of the row corresponding to that zero diagonal */
/*  element is then rotated into triangle by givens rotations . */
/*  the rank deficiency is increased by one. */
	++nl;
	if (i__ == *n) {
	    goto L90;
	}
	yi = f[i__];
	i__2 = m1;
	for (j = 1; j <= i__2; ++j) {
	    h__[j] = a[i__ + (j + 1) * a_dim1];
/* L10: */
	}
	h__[*m] = 0.f;
	i1 = i__ + 1;
	i__2 = *n;
	for (ii = i1; ii <= i__2; ++ii) {
/* Computing MIN */
	    i__3 = *n - ii;
	    i2 = min(i__3,m1);
	    piv = h__[1];
	    if (piv == 0.f) {
		goto L30;
	    }
	    fpgivs_(&piv, &a[ii + a_dim1], &cos__, &sin__);
	    fprota_(&cos__, &sin__, &yi, &f[ii]);
	    if (i2 == 0) {
		goto L70;
	    }
	    i__3 = i2;
	    for (j = 1; j <= i__3; ++j) {
		j1 = j + 1;
		fprota_(&cos__, &sin__, &h__[j1], &a[ii + j1 * a_dim1]);
		h__[j] = h__[j1];
/* L20: */
	    }
	    goto L50;
L30:
	    if (i2 == 0) {
		goto L70;
	    }
	    i__3 = i2;
	    for (j = 1; j <= i__3; ++j) {
		h__[j] = h__[j + 1];
/* L40: */
	    }
L50:
	    h__[i2 + 1] = 0.f;
/* L60: */
	}
/*  add to the sum of squared residuals the contribution of deleting */
/*  the row with small diagonal element. */
L70:
/* Computing 2nd power */
	d__1 = yi;
	*sq += d__1 * d__1;
L90:
	;
    }
/*  rank denotes the rank of a. */
    *rank = *n - nl;
/*  let b denote the (rank*n) upper trapezoidal matrix which can be */
/*  obtained from the (n*n) upper triangular matrix a by deleting */
/*  the rows and interchanging the columns corresponding to a zero */
/*  diagonal element. if this matrix is factorized using givens */
/*  transformations as  b = (r) (u)  where */
/*    r is a (rank*rank) upper triangular matrix, */
/*    u is a (rank*n) orthonormal matrix */
/*  then the minimal least-squares solution c is given by c = b' v, */
/*  where v is the solution of the system  (r) (r)' v = g  and */
/*  g denotes the vector obtained from the old right hand side f, by */
/*  removing the elements corresponding to a zero diagonal element of a. */
/*  initialization. */
    i__1 = *rank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    aa[i__ + j * aa_dim1] = 0.f;
/* L100: */
	}
    }
/*  form in aa the upper triangular matrix obtained from a by */
/*  removing rows and columns with zero diagonal elements. form in ff */
/*  the new right hand side by removing the elements of the old right */
/*  hand side corresponding to a deleted row. */
    ii = 0;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (a[i__ + a_dim1] <= *tol) {
	    goto L120;
	}
	++ii;
	ff[ii] = f[i__];
	aa[ii + aa_dim1] = a[i__ + a_dim1];
	jj = ii;
	kk = 1;
	j = i__;
/* Computing MIN */
	i__1 = j - 1;
	j1 = min(i__1,m1);
	if (j1 == 0) {
	    goto L120;
	}
	i__1 = j1;
	for (k = 1; k <= i__1; ++k) {
	    --j;
	    if (a[j + a_dim1] <= *tol) {
		goto L110;
	    }
	    ++kk;
	    --jj;
	    aa[jj + kk * aa_dim1] = a[j + (k + 1) * a_dim1];
L110:
	    ;
	}
L120:
	;
    }
/*  form successively in h the columns of a with a zero diagonal element. */
    ii = 0;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	++ii;
	if (a[i__ + a_dim1] > *tol) {
	    goto L200;
	}
	--ii;
	if (ii == 0) {
	    goto L200;
	}
	jj = 1;
	j = i__;
/* Computing MIN */
	i__1 = j - 1;
	j1 = min(i__1,m1);
	i__1 = j1;
	for (k = 1; k <= i__1; ++k) {
	    --j;
	    if (a[j + a_dim1] <= *tol) {
		goto L130;
	    }
	    h__[jj] = a[j + (k + 1) * a_dim1];
	    ++jj;
L130:
	    ;
	}
	i__1 = *m;
	for (kk = jj; kk <= i__1; ++kk) {
	    h__[kk] = 0.f;
/* L140: */
	}
/*  rotate this column into aa by givens transformations. */
	jj = ii;
	i__1 = ii;
	for (i1 = 1; i1 <= i__1; ++i1) {
/* Computing MIN */
	    i__3 = jj - 1;
	    j1 = min(i__3,m1);
	    piv = h__[1];
	    if (piv != 0.f) {
		goto L160;
	    }
	    if (j1 == 0) {
		goto L200;
	    }
	    i__3 = j1;
	    for (j2 = 1; j2 <= i__3; ++j2) {
		j3 = j2 + 1;
		h__[j2] = h__[j3];
/* L150: */
	    }
	    goto L180;
L160:
	    fpgivs_(&piv, &aa[jj + aa_dim1], &cos__, &sin__);
	    if (j1 == 0) {
		goto L200;
	    }
	    kk = jj;
	    i__3 = j1;
	    for (j2 = 1; j2 <= i__3; ++j2) {
		j3 = j2 + 1;
		--kk;
		fprota_(&cos__, &sin__, &h__[j3], &aa[kk + j3 * aa_dim1]);
		h__[j2] = h__[j3];
/* L170: */
	    }
L180:
	    --jj;
	    h__[j3] = 0.f;
/* L190: */
	}
L200:
	;
    }
/*  solve the system (aa) (f1) = ff */
    ff[*rank] /= aa[*rank + aa_dim1];
    i__ = *rank - 1;
    if (i__ == 0) {
	goto L230;
    }
    i__2 = *rank;
    for (j = 2; j <= i__2; ++j) {
	store = ff[i__];
/* Computing MIN */
	i__1 = j - 1;
	i1 = min(i__1,m1);
	k = i__;
	i__1 = i1;
	for (ii = 1; ii <= i__1; ++ii) {
	    ++k;
	    stor1 = ff[k];
	    stor2 = aa[i__ + (ii + 1) * aa_dim1];
	    store -= stor1 * stor2;
/* L210: */
	}
	stor1 = aa[i__ + aa_dim1];
	ff[i__] = store / stor1;
	--i__;
/* L220: */
    }
/*  solve the system  (aa)' (f2) = f1 */
L230:
    ff[1] /= aa[aa_dim1 + 1];
    if (*rank == 1) {
	goto L260;
    }
    i__2 = *rank;
    for (j = 2; j <= i__2; ++j) {
	store = ff[j];
/* Computing MIN */
	i__1 = j - 1;
	i1 = min(i__1,m1);
	k = j;
	i__1 = i1;
	for (ii = 1; ii <= i__1; ++ii) {
	    --k;
	    stor1 = ff[k];
	    stor2 = aa[k + (ii + 1) * aa_dim1];
	    store -= stor1 * stor2;
/* L240: */
	}
	stor1 = aa[j + aa_dim1];
	ff[j] = store / stor1;
/* L250: */
    }
/*  premultiply f2 by the transpoze of a. */
L260:
    k = 0;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	store = 0.f;
	if (a[i__ + a_dim1] > *tol) {
	    ++k;
	}
	j1 = min(i__,*m);
	kk = k;
	ij = i__ + 1;
	i__1 = j1;
	for (j = 1; j <= i__1; ++j) {
	    --ij;
	    if (a[ij + a_dim1] <= *tol) {
		goto L270;
	    }
	    stor1 = a[ij + j * a_dim1];
	    stor2 = ff[kk];
	    store += stor1 * stor2;
	    --kk;
L270:
	    ;
	}
	c__[i__] = store;
/* L280: */
    }
/*  add to the sum of squared residuals the contribution of putting */
/*  to zero the small diagonal elements of matrix (a). */
    stor3 = 0.f;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (a[i__ + a_dim1] > *tol) {
	    goto L310;
	}
	store = f[i__];
/* Computing MIN */
	i__1 = *n - i__;
	i1 = min(i__1,m1);
	if (i1 == 0) {
	    goto L300;
	}
	i__1 = i1;
	for (j = 1; j <= i__1; ++j) {
	    ij = i__ + j;
	    stor1 = c__[ij];
	    stor2 = a[i__ + (j + 1) * a_dim1];
	    store -= stor1 * stor2;
/* L290: */
	}
L300:
	fac = a[i__ + a_dim1] * c__[i__];
	stor1 = a[i__ + a_dim1];
	stor2 = c__[i__];
	stor1 *= stor2;
	stor3 += stor1 * (stor1 - store - store);
L310:
	;
    }
    fac = stor3;
    *sq += fac;
    return 0;
} /* fprank_ */

doublereal fprati_(doublereal *p1, doublereal *f1, doublereal *p2, doublereal 
	*f2, doublereal *p3, doublereal *f3)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal p, h1, h2, h3;

/*  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati */
/*  gives the value of p such that the rational interpolating function */
/*  of the form r(p) = (u*p+v)/(p+w) equals zero at p. */
/*  .. */
/*  ..scalar arguments.. */
/*  ..local scalars.. */
/*  .. */
    if (*p3 > 0.f) {
	goto L10;
    }
/*  value of p in case p3 = infinity. */
    p = (*p1 * (*f1 - *f3) * *f2 - *p2 * (*f2 - *f3) * *f1) / ((*f1 - *f2) * *
	    f3);
    goto L20;
/*  value of p in case p3 ^= infinity. */
L10:
    h1 = *f1 * (*f2 - *f3);
    h2 = *f2 * (*f3 - *f1);
    h3 = *f3 * (*f1 - *f2);
    p = -(*p1 * *p2 * h3 + *p2 * *p3 * h1 + *p3 * *p1 * h2) / (*p1 * h1 + *p2 
	    * h2 + *p3 * h3);
/*  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0. */
L20:
    if (*f2 < 0.f) {
	goto L30;
    }
    *p1 = *p2;
    *f1 = *f2;
    goto L40;
L30:
    *p3 = *p2;
    *f3 = *f2;
L40:
    ret_val = p;
    return ret_val;
} /* fprati_ */

/* Subroutine */ int fprota_(doublereal *cos__, doublereal *sin__, doublereal 
	*a, doublereal *b)
{
    static doublereal stor1, stor2;

/*  subroutine fprota applies a givens rotation to a and b. */
/*  .. */
/*  ..scalar arguments.. */
/* ..local scalars.. */
/*  .. */
    stor1 = *a;
    stor2 = *b;
    *b = *cos__ * stor2 + *sin__ * stor1;
    *a = *cos__ * stor1 - *sin__ * stor2;
    return 0;
} /* fprota_ */

/* Subroutine */ int fporde_(doublereal *x, doublereal *y, integer *m, 
	integer *kx, integer *ky, doublereal *tx, integer *nx, doublereal *ty,
	 integer *ny, integer *nummer, integer *index, integer *nreg)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, l, k1, l1, im;
    static doublereal xi, yi;
    static integer kx1, ky1, num, nyy, nk1x, nk1y;

/*  subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m */
/*  according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong */
/*  to. for each panel a stack is constructed  containing the numbers */
/*  of data points lying inside; index(j),j=1,2,...,nreg points to the */
/*  first data point in the jth panel while nummer(i),i=1,2,...,m gives */
/*  the number of the next data point in the panel. */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
    /* Parameter adjustments */
    --nummer;
    --y;
    --x;
    --tx;
    --ty;
    --index;

    /* Function Body */
    kx1 = *kx + 1;
    ky1 = *ky + 1;
    nk1x = *nx - kx1;
    nk1y = *ny - ky1;
    nyy = nk1y - *ky;
    i__1 = *nreg;
    for (i__ = 1; i__ <= i__1; ++i__) {
	index[i__] = 0;
/* L10: */
    }
    i__1 = *m;
    for (im = 1; im <= i__1; ++im) {
	xi = x[im];
	yi = y[im];
	l = kx1;
	l1 = l + 1;
L20:
	if (xi < tx[l1] || l == nk1x) {
	    goto L30;
	}
	l = l1;
	l1 = l + 1;
	goto L20;
L30:
	k = ky1;
	k1 = k + 1;
L40:
	if (yi < ty[k1] || k == nk1y) {
	    goto L50;
	}
	k = k1;
	k1 = k + 1;
	goto L40;
L50:
	num = (l - kx1) * nyy + k - *ky;
	nummer[im] = index[num];
	index[num] = im;
/* L60: */
    }
    return 0;
} /* fporde_ */

