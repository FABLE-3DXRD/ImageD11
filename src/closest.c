
#include <stdio.h>
#include <math.h>
typedef double vec[3];


#ifdef _MSC_VER
#define restrict __restrict
#define inline __inline
#endif

/* 
# ImageD11_v1.x Software for beamline ID11
# Copyright (C) 2005-2017  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/



inline
int conv_double_to_int_fast(double);

int conv_double_to_int_safe(double);

int inverse3x3(double A[3][3]);

/* Utils */

int inverse3x3(double H[3][3])
{
    double det, inverse[3][3];
    int i, j;
    /*
       # | a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
       # | a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
       # | a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |

       # DET=a11   (a33     a22    -a32     a23)-  
       a21   (a33      a12   -a32     a13)+      
       a31   (a23     a12    -a22     a13)
     */

    det = H[0][0] * (H[2][2] * H[1][1] - H[2][1] * H[1][2]) -
	H[1][0] * (H[2][2] * H[0][1] - H[2][1] * H[0][2]) +
	H[2][0] * (H[1][2] * H[0][1] - H[1][1] * H[0][2]);

    if (det != 0.) {
	inverse[0][0] = (H[2][2] * H[1][1] - H[2][1] * H[1][2]) / det;
	inverse[0][1] = -(H[2][2] * H[0][1] - H[2][1] * H[0][2]) / det;
	inverse[0][2] = (H[1][2] * H[0][1] - H[1][1] * H[0][2]) / det;
	inverse[1][0] = -(H[2][2] * H[1][0] - H[2][0] * H[1][2]) / det;
	inverse[1][1] = (H[2][2] * H[0][0] - H[2][0] * H[0][2]) / det;
	inverse[1][2] = -(H[1][2] * H[0][0] - H[1][0] * H[0][2]) / det;
	inverse[2][0] = (H[2][1] * H[1][0] - H[2][0] * H[1][1]) / det;
	inverse[2][1] = -(H[2][1] * H[0][0] - H[2][0] * H[0][1]) / det;
	inverse[2][2] = (H[1][1] * H[0][0] - H[1][0] * H[0][1]) / det;

	for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++)
		H[i][j] = inverse[i][j];

	return 0;

    } else {
	return -1;
    }
}

int conv_double_to_int_safe(double x)
{
    int a;
    a = (int)floor(x + 0.5);
    return a;
}

typedef union {
    int i;
    double d;
} a_union;


inline
int conv_double_to_int_fast(double x)
{
    /*return conv_double_to_int_safe(x); */
    /* This was benched as about eight times faster than the safe mode!! */
    /* Put in the reference for where this was found on the web TODO */
    const int p = 52;
    const double c_p1 = (1L << (p / 2));
    const double c_p2 = (1L << (p - p / 2));
    const double c_mul = c_p1 * c_p2;
    const double cs = 1.5 * c_mul;
    /* Hopefully this notes the aliasing
     * perhaps better to just return the safe version
     * not clear it is really faster any more?
     */
    a_union t;
    t.d = x + cs;
    return t.i;
    /* x += cs
       // const int a = *(int *)(&x);
       // return (a); */
}

void closest_vec( double x[], int dim, int nv, int closest[] ){
  /*
   * For each x it finds the closest neighbor
   *   this will grow as n^2, which means it rapidly becomes slow
   */
  int i, j, k, ib;
  double scor, best, t;

#pragma omp parallel for private(i,j,k,ib,scor,best,t)
  for(i=0; i<nv; i++){ /* source vector */
    /* init with something */
    j = (i+1)%nv;
    best = 0.;
    for( k=0; k<dim; k++){
      t = x[i*dim+k]-x[j*dim+k];
      best += t*t;
    }
    ib = j;
    /* now check all the others */
    for(j=0; j<nv; j++){
      if(i==j) continue;
      scor = 0.;
      for( k=0; k<dim; k++){
	t = x[i*dim+k]-x[j*dim+k];
	scor += t*t;
      }
      if( scor < best){
	ib = j;
	best = scor;
      }
    }
    closest[i] = ib;
  }
}

void closest(double x[], double v[], int *ribest, double *rbest, int nx, int nv)
{
    /* 
     * Finds value and index in x closest to a value in v
     */
    int i, j, ibest;
    double best;
    best = 99.;
    ibest = 0;
    for (i = 0; i < nx; i++) {
	for (j = 0; j < nv; j++) {
	    if (fabs(x[i] - v[j]) < best) {
		best = fabs(x[i] - v[j]);
		ibest = i;
	    }
	}
    }
    *ribest = ibest;
    *rbest = best;
}

int score(vec ubi[3], vec gv[], double tol, int ng)
{
    /*
     * Counts g-vectors indexed by ubi within tol
     */
    double sumsq, h, atol;
    int n, k, j;
    n = 0;
    atol = tol * tol;
    for (k = 0; k < ng; k++) {
	sumsq = 0.;
	for (j = 0; j < 3; j++) {
	    h = ubi[j][0] * gv[k][0] + ubi[j][1] * gv[k][1] +
		ubi[j][2] * gv[k][2];
	    h -= conv_double_to_int_fast(h);
	    sumsq += h * h;
	}
	if (sumsq < atol) {
	    n = n + 1;
	}
    }
    return n;
}

void score_and_refine(vec ubi[3], vec gv[], double tol,
		      int *n_arg, double *sumdrlv2_arg, int ng)
{
  /* ng = number of g vectors */
    double h0, h1, h2, t0, t1, t2;
    double sumsq, tolsq, sumdrlv2;
    double R[3][3], H[3][3], UB[3][3];
    int n, k, i, j, l, ih[3];
    /* Zero some stuff for refinement */
    for (i = 0; i < 3; i++) {
	ih[i] = 0;
	for (j = 0; j < 3; j++) {
	    R[i][j] = 0.;
	    H[i][j] = 0.;
	    UB[i][j] = 0.;
	}
    }
    tolsq = tol * tol;
    n = 0;
    sumdrlv2 = 0.;
    /* Test peaks */
    for (k = 0; k < ng; k++) {
	h0 = ubi[0][0] * gv[k][0] + ubi[0][1] * gv[k][1] + ubi[0][2] * gv[k][2];
	h1 = ubi[1][0] * gv[k][0] + ubi[1][1] * gv[k][1] + ubi[1][2] * gv[k][2];
	h2 = ubi[2][0] * gv[k][0] + ubi[2][1] * gv[k][1] + ubi[2][2] * gv[k][2];
	t0 = h0 - conv_double_to_int_fast(h0);
	t1 = h1 - conv_double_to_int_fast(h1);
	t2 = h2 - conv_double_to_int_fast(h2);
	sumsq = t0 * t0 + t1 * t1 + t2 * t2;
	if (sumsq < tolsq) {	/* Add into lsq problem */
	    n = n + 1;
	    sumdrlv2 += sumsq;
	    /*   From Paciorek et al Acta A55 543 (1999)
	     *   UB = R H-1
	     *   where:
	     *   R = sum_n r_n h_n^t
	     *   H = sum_n h_n h_n^t
	     *   r = g-vectors
	     *   h = hkl indices
	     *   The hkl integer indices are: */
	    ih[0] = conv_double_to_int_fast(h0);
	    ih[1] = conv_double_to_int_fast(h1);
	    ih[2] = conv_double_to_int_fast(h2);
	    /* The g-vector is: gv[k][012] */
	    for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
		    /* Robust weight factor, fn(tol), would go here */
		    R[i][j] = R[i][j] + ih[j] * gv[k][i];
		    H[i][j] = H[i][j] + ih[j] * ih[i];
		}
	    }			/* End lsq addins */
	}			/* End selected peaks */
    }				/* End first loop over spots */

    /* Now solve the least squares problem */
    /* inverse overwrites H with the inverse */
    k = inverse3x3(H);
    if (k == 0) {
	for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++)
		for (l = 0; l < 3; l++)
		    UB[i][j] += R[i][l] * H[l][j];
    }
    /* Now form ubi and copy to argument */
    if ((k == 0) && (inverse3x3(UB) == 0)) {
	for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++)
		ubi[i][j] = UB[i][j];
    } else {
	/* Determinant was zero - leave ubi as it is */
    }

    if (n > 0) {
	sumdrlv2 /= n;
    }

    /* return values */
    *n_arg = n;
    *sumdrlv2_arg = sumdrlv2;
}

int score_and_assign(vec * restrict ubi, vec * restrict gv, double tol,
		     double * restrict drlv2, int * restrict labels,
		     int label, int ng)
{

    double h0, h1, h2, t0, t1, t2, sumsq, tolsq;
    int k, n;
    tolsq = tol * tol;
    n = 0;
#pragma omp parallel for private(h0,h1,h2,t0,t1,t2,sumsq) reduction(+:n)
    for (k = 0; k < ng; k++) {
	h0 = ubi[0][0] * gv[k][0] + ubi[0][1] * gv[k][1] + ubi[0][2] * gv[k][2];
	h1 = ubi[1][0] * gv[k][0] + ubi[1][1] * gv[k][1] + ubi[1][2] * gv[k][2];
	h2 = ubi[2][0] * gv[k][0] + ubi[2][1] * gv[k][1] + ubi[2][2] * gv[k][2];
	t0 = h0 - conv_double_to_int_fast(h0);
	t1 = h1 - conv_double_to_int_fast(h1);
	t2 = h2 - conv_double_to_int_fast(h2);
	sumsq = t0 * t0 + t1 * t1 + t2 * t2;
	/* If this peak fits better than the one in drlv2 then we
	 * assign it:
	 */
	if ((sumsq < tolsq) && (sumsq < drlv2[k])) {
	    labels[k] = label;
	    drlv2[k] = sumsq;
	    n++;
	} else if (labels[k] == label) {
	    /* We thought it belonged but it does not */
	    labels[k] = -1;
	}
    }
    return n;
}

void refine_assigned(vec ubi[3], vec gv[], int labels[], int label,
		     int *npk, double *sumdrlv2, int ng)
{
    /* Skip the part about weights, not used */
    double sumsqtot, sumsq, h[3], t[3];
    double R[3][3], H[3][3], UB[3][3];
    int i, j, n, k, l, ih[3];
    n = 0;
    sumsqtot = 0;
    for (k = 0; k < ng; k++) {
	if (label != labels[k]) {
	    continue;
	}
	n++;
	for (j = 0; j < 3; j++) {
	    h[j] =
		ubi[j][0] * gv[k][0] + ubi[j][1] * gv[k][1] +
		ubi[j][2] * gv[k][2];
	    ih[j] = conv_double_to_int_fast(h[j]);
	    t[j] = h[j] - ih[j];
	}
	sumsq = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
	sumsqtot += sumsq;
	for (i = 0; i < 3; i++) {
	    for (j = 0; j < 3; j++) {
		R[i][j] = R[i][j] + ih[j] * gv[k][i];
		H[i][j] = H[i][j] + ih[j] * ih[i];
	    }
	}
    }
    /* outputs */
    *npk = n;
    if (n > 0) {
	*sumdrlv2 = sumsqtot / n;
    } else {
	*sumdrlv2 = 0.;
    }
    /* And the fitted matrix */
    k = inverse3x3(H);
    if (k == 0) {
	/* Form best fit UB */
	for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++)
		for (l = 0; l < 3; l++)
		    UB[i][j] = UB[i][j] + R[i][l] * H[l][j];
    }

    if (k == 0 && inverse3x3(UB) == 0) {
	/* Copy to output */
	for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++)
		ubi[i][j] = UB[i][j];
    }
}

/* weighted_refine was never used */

void put_incr(float data[], size_t ind[], float vals[], int boundscheck,
	      int n, int m)
{
    int k, ik;
    if (boundscheck == 0) {
	for (k = 0; k < n; k++)
	    data[ind[k]] += vals[k];
    } else {
	for (k = 0; k < n; k++) {
	    ik = ind[k];
	    if (ik < 0 || ik >= m) {
		printf("Array bounds error! k=%d ind[k]=%d\n", k, (int)ind[k]);
	    } else {
		data[ind[k]] += vals[k];
	    }
	}
    }
}


void cluster1d( double ar[], int n, int order[], double tol, // IN
		int *nclusters, int ids[], double avgs[]){   // OUT
  // Used in sandbox/friedel.py
  int i, ncl;
  double dv;
  // order is the order of the peaks to get them sorted
  avgs[0] = ar[order[0]];
  ncl = 1;               // number in this cluster
  ids[0] = 0;            // cluster assignments ( in order )
  for( i=1 ; i<n; i++){
    dv = ar[order[i]] - ar[order[i-1]];  // difference in values
    if( dv > tol ){    // make a new cluster
      if( ncl > 1 ) {  // make avg for the last one 
	avgs[ ids[i-1] ] = avgs[ids[i-1]] / ncl;
      }
      ids[i] = ids[i-1] + 1;         // increment id
      ncl = 1;                       // pks in this cluster
      avgs[ ids[i] ] = ar[order[i]]; // store value for avg
    } else { 
      ids[i] = ids[i-1];             // copy last id
      ncl = ncl + 1;
      avgs[ids[i]] = avgs[ids[i]] + ar[order[i]];  // sum on for avg
    }
  } // end for(i ...
  // make the last average if necessary
  if( ncl > 1 ) { 
    avgs[ ids[i-1] ] /= ncl;
  }
  *nclusters = ids[n-1]+1;
}



void score_gvec_z( vec ubi[3],       // in
		   vec ub[3],        // in
		   vec gv[],         // in
		   vec g0[],         // inout  normed(g)
		   vec g1[],         // inout  normed(axis x g)
		   vec g2[],         // inout  normed(axis x (axis x g))
		   vec e[],          // inout
		   int recompute,    // in
		   int n ){
  /*  Axis is z and we hard wire it here
   *     Compute errors in a co-ordinate system given by:
   *         parallel to gv  (gv*imodg)
   *         parallel to cross( axis, gv )
   *         parallel to cross( axis, cross( axis, gv ) )
   */
  int i;
  double t,txy;
  vec g, h, d;
  
#pragma omp parallel for private(i,t,txy,g,h,d)
  for(i=0; i<n; i++){
    g[0] = gv[i][0];
    g[1] = gv[i][1];
    g[2] = gv[i][2];
    // Test - is it faster to recompute or cache ?
    //        for many ubi ? Loop over peaks or ubis or ?
    if( recompute ){  // Fill in ub, modg, ax_x_gv, ax_ax_x_gv
      t = 1./sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
      g0[i][0] = g[0]*t;
      g0[i][1] = g[1]*t;
      g0[i][2] = g[2]*t;
      txy =  g[0]*g[0] + g[1]*g[1];
      t = 1./sqrt( txy );
      g1[i][0] = -g[1]*t;  // [-y,x,0]
      g1[i][1] =  g[0]*t;
      g1[i][2] =  0.;
      t = 1./sqrt(g[0]*g[0]*g[2]*g[2] + g[1]*g[1]*g[2]*g[2] + txy*txy);
      g2[i][0] = g[0]*g[2]*t;
      g2[i][1] = g[1]*g[2]*t;
      g2[i][2] = -(g[0]*g[0]+g[1]*g[1])*t; 
    } // end recompute
    
    // Find integer h,k,l
    h[0] = round( ubi[0][0]*g[0] + ubi[0][1]*g[1] + ubi[0][2]*g[2] );
    h[1] = round( ubi[1][0]*g[0] + ubi[1][1]*g[1] + ubi[1][2]*g[2] );
    h[2] = round( ubi[2][0]*g[0] + ubi[2][1]*g[1] + ubi[2][2]*g[2] );

    // Compute diff, the computed g-vector  - original
    d[0] = ub[0][0]*h[0] + ub[0][1]*h[1] + ub[0][2]*h[2]  - g[0];
    d[1] = ub[1][0]*h[0] + ub[1][1]*h[1] + ub[1][2]*h[2]  - g[1];
    d[2] = ub[2][0]*h[0] + ub[2][1]*h[1] + ub[2][2]*h[2]  - g[2];
    
    // Now dot this into the local co-ordinate system
    e[i][0] = g0[i][0]*d[0] + g0[i][1]*d[1] + g0[i][2]*d[2] ;
    e[i][1] = g1[i][0]*d[0] + g1[i][1]*d[1] + g1[i][2]*d[2] ;
    e[i][2] = g2[i][0]*d[0] + g2[i][1]*d[1] + g2[i][2]*d[2] ; 
  }
}
		 
