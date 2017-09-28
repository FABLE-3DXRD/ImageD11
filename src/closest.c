
#include <stdio.h>
#include <math.h>
typedef double vec[3];

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

static char moduledocs[] =
    "/* *******************************************************************\n"
    " * closest.c  Two useful functions for indexing grains \n"
    " *\n"
    " * 1. From python: closest(cosines_allowed, cosine_measured)\n"
    " *    Take 1 peak with hkl's assigned, compute the cosines to\n"
    " *    another powder ring (angles between this peak and peaks in\n"
    " *    that ring). This extension finds the best peak in the\n"
    " *    other ring to pair up with.\n"
    " *\n"
    " *    More generally closest(array, values) finds the closest\n"
    " *    item in values to one of the values in array.  Returns the\n"
    " *    difference and the index of the closest thing it finds.\n"
    " *\n"
    " *    Both arguments should be one dimensional Numeric arrays\n"
    " *    of type Numeric.Float\n"
    " *\n"
    " * 2. From python: score(ubi, gv, tol) where ubi is an orientation\n"
    " *    matrix and gv are an array of g-vectors. Returns the number\n"
    " *    of g-vectors which have integer hkl peaks within tolerance\n"
    " *    tol. Uses the conv_double_to_int_fast function in here for \n"
    " *    factor of !EIGHT! speed increase compared to rounding in \n"
    " *    C. In fact this gives the nearest even integer, instead\n"
    " *    of the nearest integer, but we don't care, as a peak having\n"
    " *    hkl of 0.5 is nowhere near being indexed anyway.\n"
    " *\n"
    " *    Returns and integer - number of peaks indexed\n"
    " *    UBI is a 3x3 Numeric.Float array (figure out the transposing yourself)\n"
    " *    GV is a nx3 Numeric.Float array, and you should try to make the 3 \n"
    " *       be the fast index for best performance\n"
    " *       \n"
    " * 3. score_and_refine(ubi, gv, tol) \n"
    " *    From python: same as score, I hope, but ubi is overwritten\n"
    " *    with refined matrix following paciorek algorithm which is \n"
    " *    in indexing.py\n"
    " *    \n"
    " * \n"
    " * 4. score_and_assign( ubi, gv, tol, drlv, labels, (int) label)\n"
    " *    as for score, but assignments are only made if drlv is lower than\n"
    " *    the previous. Fills in label with the new label\n"
    " *\n"
    " * 5. refine_assigned( ubi, gv, labels, label, weight) \n"
    " *    Use only the peaks where label matches labels in refinement \n"
    " *    ...following again the Paciorek algorithm \n"
    " *\n"
    " * 6. put_incr( data, indices, values)\n"
    " *    pretty much as numeric.put but as increments\n"
    " *\n"
    " * 7. weighted_refine(ubi, gv, tol, weights) \n"
    " *    From python: same as score, but with weights, and ubi is overwritten\n"
    " *    with refined matrix following paciorek algorithm which is \n"
    " *    in indexing.py\n"
    " *    \n"
    " * ****************************************************************** */ ";

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

int score_and_assign(vec * __restrict ubi, vec * __restrict gv, double tol,
		     double *__restrict drlv2, int *__restrict labels,
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
