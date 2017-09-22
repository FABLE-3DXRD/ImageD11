
#include <stdio.h>

#define DEBUG 0
/* turn to 1 to debug */
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


static char moduledocs[] = \
"/* *******************************************************************\n"\
" * closest.c  Two useful functions for indexing grains \n"\
" *\n"\
" * 1. From python: closest(cosines_allowed, cosine_measured)\n"\
" *    Take 1 peak with hkl's assigned, compute the cosines to\n"\
" *    another powder ring (angles between this peak and peaks in\n"\
" *    that ring). This extension finds the best peak in the\n"\
" *    other ring to pair up with.\n"\
" *\n"\
" *    More generally closest(array, values) finds the closest\n"\
" *    item in values to one of the values in array.  Returns the\n"\
" *    difference and the index of the closest thing it finds.\n"\
" *\n"\
" *    Both arguments should be one dimensional Numeric arrays\n"\
" *    of type Numeric.Float\n"\
" *\n"\
" * 2. From python: score(ubi, gv, tol) where ubi is an orientation\n"\
" *    matrix and gv are an array of g-vectors. Returns the number\n"\
" *    of g-vectors which have integer hkl peaks within tolerance\n"\
" *    tol. Uses the conv_double_to_int_fast function in here for \n"\
" *    factor of !EIGHT! speed increase compared to rounding in \n"\
" *    C. In fact this gives the nearest even integer, instead\n"\
" *    of the nearest integer, but we don't care, as a peak having\n"\
" *    hkl of 0.5 is nowhere near being indexed anyway.\n"\
" *\n"\
" *    Returns and integer - number of peaks indexed\n"\
" *    UBI is a 3x3 Numeric.Float array (figure out the transposing yourself)\n"\
" *    GV is a nx3 Numeric.Float array, and you should try to make the 3 \n"\
" *       be the fast index for best performance\n"\
" *       \n"\
" * 3. score_and_refine(ubi, gv, tol) \n"\
" *    From python: same as score, I hope, but ubi is overwritten\n"\
" *    with refined matrix following paciorek algorithm which is \n"\
" *    in indexing.py\n"\
" *    \n"\
" * \n"\
" * 4. score_and_assign( ubi, gv, tol, drlv, labels, (int) label)\n"\
" *    as for score, but assignments are only made if drlv is lower than\n"\
" *    the previous. Fills in label with the new label\n"\
" *\n"\
" * 5. refine_assigned( ubi, gv, labels, label, weight) \n"\
" *    Use only the peaks where label matches labels in refinement \n"\
" *    ...following again the Paciorek algorithm \n"\
" *\n"\
" * 6. put_incr( data, indices, values)\n"\
" *    pretty much as numeric.put but as increments\n"\
" *\n"\
" * 7. weighted_refine(ubi, gv, tol, weights) \n"\
" *    From python: same as score, but with weights, and ubi is overwritten\n"\
" *    with refined matrix following paciorek algorithm which is \n"\
" *    in indexing.py\n"\
" *    \n"\

" * ****************************************************************** */ ";

int conv_double_to_int_fast(double);

int conv_double_to_int_safe(double);

int inverse3x3(double A[3][3]);

/* Utils */


int conv_double_to_int_safe(double x){
   int a;
   a = (int) floor(x+0.5);
   return a;
}

typedef union {
	int i;
	double d;
	} a_union;
	
int conv_double_to_int_fast(double x){
	/*return conv_double_to_int_safe(x);*/
   /* This was benched as about eight times faster than the safe mode!! */
   /* Put in the reference for where this was found on the web TODO */
   const int p=52;
   const double c_p1 = (1L << (p/2));
   const double c_p2 = (1L << (p-p/2));
   const double c_mul = c_p1 * c_p2;
   const double cs = 1.5*c_mul;
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



void closest( double x[], double v[], int *ribest, double *rbest, int nx, int nv ){
   int i,j,ibest;
   double best;
   best=99.;
   ibest=0;
   for(i=0; i<nx;i++){
      for(j=0; j<nv;j++){
         if( fabs(x[i]-v[j]) < best){
            best=fabs(x[i]-v[j]);
            ibest=i;
         }
      }
   }
   *ribest = ibest;
   *rbest  = best;
}



int score(  vec ubi[3], vec gv[], double tol, int ng){
   double sumsq,h,t,atol;
   int n,k,j;
   n=0;
   atol=tol*tol;
   for( k=0 ; k<ng ; k++){ 
      sumsq = 0.;
      for( j = 0; j<3; j++){
          h  = ubi[j][0]*gv[k][0]+ubi[j][1]*gv[k][1]+ubi[j][2]*gv[k][2];
          h -= conv_double_to_int_fast( h );
          sumsq += h*h;
      }
      if (sumsq < atol){
          n=n+1;
      }
   }
   return n;
}


