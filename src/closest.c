
#define DEBUG 0
/* turn to 1 to debug */

#if defined( _MSC_VER ) && !defined( __cplusplus )
# define inline __inline
#endif // defined( _MSC_VER ) && !defined( __cplusplus )


/* 
# ImageD11_v1.0 Software for beamline ID11
# Copyright (C) 2005-2007  Jon Wright
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

#include <Python.h>                  /* To talk to python */
/* #include "Numeric/arrayobject.h"      Access to Numeric */
#include "numpy/arrayobject.h"     /*  upgrade to numpy */

int conv_double_to_int_fast(double);

int conv_double_to_int_safe(double);

int inverse3x3(double A[3][3]);

/* Utils */


int conv_double_to_int_safe(double x){
   int a;
   a = floor(x+0.5);
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



static PyObject *closest( PyObject *self, PyObject *args, PyObject *keywds){
   PyArrayObject *ar=NULL, *vals=NULL;
   double *x, *v;
   int ibest,i,j,nj;
   double best;
   if(!PyArg_ParseTuple(args,"O!O!",
                        &PyArray_Type, &ar,   /* array args */
                        &PyArray_Type, &vals))   /* array args */
      return NULL;

   if(ar->nd != 1 || ar->descr->type_num!=PyArray_DOUBLE){
      PyErr_SetString(PyExc_ValueError,
            "First array must be 1d and double");
      return NULL;
   }
   if(vals->nd != 1 || vals->descr->type_num!=PyArray_DOUBLE){
      PyErr_SetString(PyExc_ValueError,
            "Second array must be 1d and double");
      return NULL;
   }
      
   if(ar->strides[0] != sizeof(double) || vals->strides[0] != sizeof(double)){
      PyErr_SetString(PyExc_ValueError,
            "Arrays must be flat (strides == sizeof(double))");
      return NULL;
   }

   x =  (double*) ar->data;
   v =  (double*) vals->data;

   best=99.;
   ibest=0;
   nj=vals->dimensions[0];
   for(i=0; i<ar->dimensions[0];i++){
      for(j=0; j<nj;j++){
/*         printf("x[i]=%f v[i]=%f diff=%f\n",x[i],v[j],fabs(x[i]-v[j]));*/
         if( fabs(x[i]-v[j]) < best){
            best=fabs(x[i]-v[j]);
            ibest=i;
         }
      }
   }
/*   printf("Best %f ibest=%d\n",best,ibest); */
   return Py_BuildValue("if",ibest,best);
}

/* score starts here ==================================================== */

static PyObject *score( PyObject *self, PyObject *args, PyObject *keywds){
   PyArrayObject *ubi=NULL, *gv=NULL;
   double u00,u11,u22,u01,u02,u10,u12,u20,u21;
   double g0,g1,g2,h0,h1,h2,t0,t1,t2;
   double tol,sumsq;
   int n,k;

   
   if(!PyArg_ParseTuple(args,"O!O!d",
                        &PyArray_Type, &ubi,   /* array args */
                        &PyArray_Type, &gv,    /* array args */
                        &tol)) /* Tolerance */
      return NULL;

   if(ubi->nd != 2 || ubi->descr->type_num!=PyArray_DOUBLE){
      printf("first arg nd %d\n",ubi->nd);
      printf("first arg type %d\n",ubi->descr->type_num);
      PyErr_SetString(PyExc_ValueError,
            "First array must be 3x3 2d and double");
      return NULL;
   }
   if(gv->nd != 2 || gv->descr->type_num!=PyArray_DOUBLE){
      PyErr_SetString(PyExc_ValueError,
            "Second array must be 3xn 2d and double");
      return NULL;
   }
      
/*   if(gv->strides[0] != sizeof(double) || ubi->strides[0] != sizeof(double)){
      PyErr_SetString(PyExc_ValueError,
            "Arrays must be flat (strides == sizeof(double))");
      return NULL;
   } */

   if(ubi->dimensions[0] != 3 || ubi->dimensions[1] != 3){
      PyErr_SetString(PyExc_ValueError,
            "First array must be 3x3 2d, dimensions problem");
      return NULL;
   }

   if(gv->dimensions[1] != 3){
      PyErr_SetString(PyExc_ValueError,
            "Second array must be 3xn 2d and double");
      return NULL;
   }
/*   if(gv->strides[1] != sizeof(double)){
      PyErr_SetString(PyExc_ValueError,
            "Second array must be 3xn 2d and double, and you want it aligned!!!");
      return NULL;
   }
*/

   /*   for(i=0;i<3;i++)
      for(j=0;j<3;j++)
         u[i][j] = * (double *) (ubi->data + i*ubi->strides[0] + j*ubi->strides[1]); */

   u00=* (double *) (ubi->data + 0*ubi->strides[0] + 0*ubi->strides[1]);
   u01=* (double *) (ubi->data + 0*ubi->strides[0] + 1*ubi->strides[1]);
   u02=* (double *) (ubi->data + 0*ubi->strides[0] + 2*ubi->strides[1]);
   u10=* (double *) (ubi->data + 1*ubi->strides[0] + 0*ubi->strides[1]);
   u11=* (double *) (ubi->data + 1*ubi->strides[0] + 1*ubi->strides[1]);
   u12=* (double *) (ubi->data + 1*ubi->strides[0] + 2*ubi->strides[1]);
   u20=* (double *) (ubi->data + 2*ubi->strides[0] + 0*ubi->strides[1]);
   u21=* (double *) (ubi->data + 2*ubi->strides[0] + 1*ubi->strides[1]);
   u22=* (double *) (ubi->data + 2*ubi->strides[0] + 2*ubi->strides[1]);

   n=0;
   tol=tol*tol;

   for(k=0;k<gv->dimensions[0];k++){ /* Loop over observed peaks */
      /* Compute hkls of this peak as h = UBI g */
      g0 = * (double *) (gv->data + k*gv->strides[0] + 0*gv->strides[1]);
      g1 = * (double *) (gv->data + k*gv->strides[0] + 1*gv->strides[1]);
      g2 = * (double *) (gv->data + k*gv->strides[0] + 2*gv->strides[1]);
      h0 = u00*g0 + u01*g1 + u02*g2;
      h1 = u10*g0 + u11*g1 + u12*g2;
      h2 = u20*g0 + u21*g1 + u22*g2;
      t0=h0-conv_double_to_int_fast(h0);
      t1=h1-conv_double_to_int_fast(h1);
      t2=h2-conv_double_to_int_fast(h2);
      if ( (t0 > 0.501) ||
	       (t1 > 0.501) ||
	       (t2 > 0.501) ){
	        printf("Error in conv\nh0 = %f h1 = %f h2=%f \nt0 = %f t1 = %f t2=%f\n",
	        t0,t1,t2,h0,h1,h2);
	        return NULL;
	        }       
	        
      sumsq = t0*t0+t1*t1+t2*t2;
      if (sumsq < tol)
         n=n+1;
      /*
      printf("k=%d g=(",k);
      printf(" %6.3f ",g[0]);
      printf(" %6.3f ",g[1]);
      printf(" %6.3f )  h=(",g[2]);
      printf(" %6.3f ",h[0]);
      printf(" %6.3f ",h[1]);
      printf(" %6.3f ) ",h[2]);
      printf("int(h)=(%3d %3d %3d) sumsq=%g n=%d\n",conv_double_to_int_fast(h[0]),
            conv_double_to_int_fast(h[1]),conv_double_to_int_fast(h[2]),sumsq,n); 
      */
   }
   
   /*   printf("Exit score and assign\n");*/
   return Py_BuildValue("i",n);
}

/* score and refine  =========================================================== */

static PyObject *score_and_refine( PyObject *self, PyObject *args, PyObject *keywds){
   PyArrayObject *ubi=NULL, *gv=NULL;
   double u00,u11,u22,u01,u02,u10,u12,u20,u21;
   double g0,g1,g2,h0,h1,h2,t0,t1,t2;
   double tol,sumsq,tolsq,sumdrlv2;
   double R[3][3],H[3][3],ih[3],rh[3], UB[3][3];
   int n,k,i,j,l;

   for(i=0;i<3;i++){
     ih[i]=0.;
     rh[i]=0.;
     for(j=0;j<3;j++){
       R[i][j] = 0.;
       H[i][j] = 0.;
       UB[i][j] = 0.;
     }
   }

   
   if(!PyArg_ParseTuple(args,"O!O!d",
                        &PyArray_Type, &ubi,   /* array args */
                        &PyArray_Type, &gv,    /* array args */
                        &tol)) /* Tolerance */
     return NULL;
   
   if(ubi->nd != 2 || ubi->descr->type_num!=PyArray_DOUBLE){
     printf("first arg nd %d\n",ubi->nd);
     printf("first arg type %d\n",ubi->descr->type_num);
     PyErr_SetString(PyExc_ValueError,
		     "First array must be 3x3 2d and double");
     return NULL;
   }
   if(gv->nd != 2 || gv->descr->type_num!=PyArray_DOUBLE){
      PyErr_SetString(PyExc_ValueError,
		      "Second array must be 3xn 2d and double");
      return NULL;
   }
   
/*   if(gv->strides[0] != sizeof(double) || ubi->strides[0] != sizeof(double)){
      PyErr_SetString(PyExc_ValueError,
            "Arrays must be flat (strides == sizeof(double))");
      return NULL;
   } */

   if(ubi->dimensions[0] != 3 || ubi->dimensions[1] != 3){
      PyErr_SetString(PyExc_ValueError,
            "First array must be 3x3 2d, dimensions problem");
      return NULL;
   }

   if(gv->dimensions[1] != 3){
      PyErr_SetString(PyExc_ValueError,
            "Second array must be 3xn 2d and double");
      return NULL;
   }
/*   if(gv->strides[1] != sizeof(double)){
      PyErr_SetString(PyExc_ValueError,
            "Second array must be 3xn 2d and double, and you want it aligned!!!");
      return NULL;
   }
*/


   u00=* (double *) (ubi->data + 0*ubi->strides[0] + 0*ubi->strides[1]);
   u01=* (double *) (ubi->data + 0*ubi->strides[0] + 1*ubi->strides[1]);
   u02=* (double *) (ubi->data + 0*ubi->strides[0] + 2*ubi->strides[1]);
   u10=* (double *) (ubi->data + 1*ubi->strides[0] + 0*ubi->strides[1]);
   u11=* (double *) (ubi->data + 1*ubi->strides[0] + 1*ubi->strides[1]);
   u12=* (double *) (ubi->data + 1*ubi->strides[0] + 2*ubi->strides[1]);
   u20=* (double *) (ubi->data + 2*ubi->strides[0] + 0*ubi->strides[1]);
   u21=* (double *) (ubi->data + 2*ubi->strides[0] + 1*ubi->strides[1]);
   u22=* (double *) (ubi->data + 2*ubi->strides[0] + 2*ubi->strides[1]);

   n=0;
   sumdrlv2 = 0;
   tolsq=tol*tol;

   for(k=0;k<gv->dimensions[0];k++){ /* Loop over observed peaks */
      /* Compute hkls of this peak as h = UBI g */
      g0 = * (double *) (gv->data + k*gv->strides[0] + 0*gv->strides[1]);
      g1 = * (double *) (gv->data + k*gv->strides[0] + 1*gv->strides[1]);
      g2 = * (double *) (gv->data + k*gv->strides[0] + 2*gv->strides[1]);
      h0 = u00*g0 + u01*g1 + u02*g2;
      h1 = u10*g0 + u11*g1 + u12*g2;
      h2 = u20*g0 + u21*g1 + u22*g2;
      
      t0=h0-conv_double_to_int_fast(h0);
      t1=h1-conv_double_to_int_fast(h1);
      t2=h2-conv_double_to_int_fast(h2);
      sumsq = t0*t0+t1*t1+t2*t2;

      if ( (t0 > 0.501) ||
	   (t1 > 0.501) ||
	   (t2 > 0.501) ){
	printf("Error in conv\nh0 = %f h1 = %f h2=%f \nt0 = %f t1 = %f t2=%f\n",
	       t0,t1,t2,h0,h1,h2);
	return NULL;
      }       
	  
      
      if (sumsq < tolsq){
	n=n+1;
	sumdrlv2 += sumsq;
	/*   From Paciorek et al Acta A55 543 (1999)
	 *   UB = R H-1
	 *    where:
	 *    R = sum_n r_n h_n^t
	 *    H = sum_n h_n h_n^t
	 *    r = g-vectors
	 *    h = hkl indices
	 *    The hkl integer indices are: */
	ih[0] = conv_double_to_int_fast(h0);
	ih[1] = conv_double_to_int_fast(h1);
	ih[2] = conv_double_to_int_fast(h2);
	/* The g-vector was: */
	rh[0] = g0; rh[1] = g1 ; rh[2] = g2;
	for(i=0;i<3;i++){
	  for(j=0;j<3;j++){
	    /* Robust weight factor, fn(tol), would go here */
	    R[i][j] = R[i][j] + ih[j] * rh[i];
	    H[i][j] = H[i][j] + ih[j] * ih[i];
	  }
	}
      }
   }


   if (inverse3x3(H)==0){
     /* Form best fit UB */
     for(i=0;i<3;i++)
       for(j=0;j<3;j++)
	 for(l=0;l<3;l++)
	   UB[i][j] = UB[i][j] + R[i][l]*H[l][j];
     
     if (inverse3x3(UB)==0){
       /*printf("Got a new UB\n");
	 for(i=0;i<3;i++)
	 for(j=0;j<3;j++)
	 printf("UBi[%d][%d]=%f",i,j,UB[i][j]);
	 printf("\n");
       */

       /* Copy result into arg */
       for(i=0; i<3; i++)
	 for(j=0; j<3; j++)
	   *(double *)(ubi->data + 
		       i*ubi->strides[0] + 
		       j*ubi->strides[1]) = UB[i][j];

       
     }
     
   }
   
   if( n > 0){    
       return Py_BuildValue("id",n,sumdrlv2/n);
    }   else  {
       return Py_BuildValue("id",n,0);
    }
}


/* ================================================
 *
 * Add a new routine which takes a ubi matrix and integer label
 *
 * It assigns for each peak if this label gives a lower drlv2
 * and what is the new minimal drlv2 
 *
 *
 * This will allow grains to "fight over" the peaks
 * 4. score_and_assign( ubi, gv, tol, drlv, labels, (int) label)
 *    as for score, but assignments are only made if drlv is lower than
 *    the previous. Fills in label with the new label

 */

static PyObject *score_and_assign( PyObject *self, PyObject *args, PyObject *keywds){
  PyArrayObject *ubi=NULL, *gv=NULL, *labels=NULL, *drlv2=NULL;
  double u00,u11,u22,u01,u02,u10,u12,u20,u21;
  double g0,g1,g2, h0,h1,h2, t0,t1,t2;
  double tol,sumsq, *dp;
  npy_int32 n,k,label, *lp;
  

   
  if(!PyArg_ParseTuple(args,"O!O!dO!O!i",
		       &PyArray_Type, &ubi,   /* array args */
		       &PyArray_Type, &gv,    /* array args */
		       &tol,  /* Tolerance for assigning    */
		       &PyArray_Type, &drlv2, /* array args */
		       &PyArray_Type, &labels,/* array args */
		       &label))               /* integer */
    
    return NULL;
  if(DEBUG)printf("Welcome to score_and_assign\n");
  
  if(ubi->nd != 2 || ubi->descr->type_num!=PyArray_DOUBLE ||
     ubi->dimensions[0] != 3 || ubi->dimensions[1] != 3){
    printf("first arg nd %d\n",ubi->nd);
    printf("first arg type %d\n",ubi->descr->type_num);
    PyErr_SetString(PyExc_ValueError,
		    "First array must be 3x3 2d and double");
    return NULL;
  }
  if(DEBUG)printf("Got ubi\n");
  if(gv->nd != 2 || gv->descr->type_num!=PyArray_DOUBLE){
    PyErr_SetString(PyExc_ValueError,
		    "Second array must be 3xn 2d and double");
    return NULL;
  }
  if(DEBUG)printf("Got gv\n");  
  if(gv->dimensions[1] != 3){
    PyErr_SetString(PyExc_ValueError,
		    "Second array must be 3xn 2d and double");
    return NULL;
  }
  if(DEBUG)printf("Got gv-dim3\n");

  if(drlv2->nd != 1 || drlv2->descr->type_num!=PyArray_DOUBLE ||
     drlv2->dimensions[0] != gv->dimensions[0] ){
    PyErr_SetString(PyExc_ValueError,
		    "drlv2 array must be double and same length as gv");
    return NULL;
  }
  if(DEBUG)printf("Got drlv\n");

  if(labels->nd != 1 ||
    PyArray_EquivTypenums(labels->descr->type_num,NPY_INT32)!=NPY_TRUE ||
     labels->dimensions[0] != gv-> dimensions[0] ){
    printf("Problem with labels\n");
    printf("nd %d got %d want %d dims[0] %d  gv_dims[0] %d\n",labels->nd,
	   labels->descr->type_num,
       NPY_INT32,
	   (int) labels->dimensions[0],
	   (int) gv-> dimensions[0]) ;
    PyErr_SetString(PyExc_ValueError,
		    "Label array must be integer and same length as gv");
    return NULL;
  }
  if(DEBUG)printf("Got labels\n");
  if(DEBUG)printf("score_and_assign: copy UBI\n");
   /* for convenience */
   u00=* (double *) (ubi->data + 0*ubi->strides[0] + 0*ubi->strides[1]);
   u01=* (double *) (ubi->data + 0*ubi->strides[0] + 1*ubi->strides[1]);
   u02=* (double *) (ubi->data + 0*ubi->strides[0] + 2*ubi->strides[1]);
   u10=* (double *) (ubi->data + 1*ubi->strides[0] + 0*ubi->strides[1]);
   u11=* (double *) (ubi->data + 1*ubi->strides[0] + 1*ubi->strides[1]);
   u12=* (double *) (ubi->data + 1*ubi->strides[0] + 2*ubi->strides[1]);
   u20=* (double *) (ubi->data + 2*ubi->strides[0] + 0*ubi->strides[1]);
   u21=* (double *) (ubi->data + 2*ubi->strides[0] + 1*ubi->strides[1]);
   u22=* (double *) (ubi->data + 2*ubi->strides[0] + 2*ubi->strides[1]);

   n=0; /* This will be the number of peaks we assign */
   
   tol = tol*tol; /* To use drlv2 */

   if(DEBUG) printf("Into score and assign");

   for(k=0;k<gv->dimensions[0];k++){ /* Loop over observed peaks */
      /* Compute hkls of this peak as h = UBI g */
      g0 = * (double *) (gv->data + k*gv->strides[0] + 0*gv->strides[1]);
      g1 = * (double *) (gv->data + k*gv->strides[0] + 1*gv->strides[1]);
      g2 = * (double *) (gv->data + k*gv->strides[0] + 2*gv->strides[1]);
      h0 = u00*g0 + u01*g1 + u02*g2;
      h1 = u10*g0 + u11*g1 + u12*g2;
      h2 = u20*g0 + u21*g1 + u22*g2;
      t0=h0-conv_double_to_int_fast(h0);
      t1=h1-conv_double_to_int_fast(h1);
      t2=h2-conv_double_to_int_fast(h2);
      /* this error checking should be redundant nowadays (unit test?) */
      /*
      if ( (t0 > 0.501) ||
	   (t1 > 0.501) ||
	   (t2 > 0.501) ){
	printf("Error in conv\nh0 = %f h1 = %f h2=%f \nt0 = %f t1 = %f t2=%f\n",
	       t0,t1,t2,h0,h1,h2);
	return NULL;
	} */      
	        
      sumsq = t0*t0 + t1*t1 + t2*t2;

      /* Current label assignment (pointer to) */
      lp = (int *) (labels->data + k*labels->strides[0]);

      /* Current drlv assignment (pointer to) */
      dp = (double *) (drlv2->data + k*drlv2->strides[0]);

      if ( sumsq < *dp  && sumsq < tol ){
	/* This ubi matches this peak better than anything else so far
	 * ... we shall assign it 
	 */
	if(DEBUG>1)printf("Assigning label\n");
	*lp = label;
	*dp = sumsq;
	n++;
	continue; 
      } 

      /* if ( sumsq < tol ) {
       *   we don't care - the peak is taken by someone else already
       *  it only switches if the other guy fits worse than us
       */

      if ( *lp == label ){
	/* In case this peak thinks it belongs to this grain, but it doesn't
	 */ 
	*lp = -1;
      }
   }
   /* Nothing malloc'ed so nothing to free 
    *   We return the number of peaks "won" by this ubi
    */
   return Py_BuildValue("i",n);
}






/* ===================================================================== 
 *
 *
 * 5. refine_assigned( ubi, gv, labels, label, weight) 
 *    Use only the peaks where label matches labels in refinement 
 *    ...following again the Paciorek algorithm 
 */

static PyObject *refine_assigned( PyObject *self, 
				  PyObject *args, 
				  PyObject *keywds){
  PyArrayObject *ubi=NULL, *gv=NULL, *labels=NULL;
  int label;
  

  double u00,u11,u22,u01,u02,u10,u12,u20,u21;
  double g0,g1,g2,h0,h1,h2,t0,t1,t2;
  double sumsq,sumsqtot;
  double sigma, weight, wfac;
  double R[3][3],H[3][3],ih[3],rh[3], UB[3][3];
  int n,k,i,j,l;
  

  for(i=0;i<3;i++){
    ih[i]=0.;
    rh[i]=0.;
    for(j=0;j<3;j++){
      R[i][j] = 0.;
      H[i][j] = 0.;
      UB[i][j] = 0.;
    }
  }
  
   
  if(!PyArg_ParseTuple(args,"O!O!O!id",
		       &PyArray_Type, &ubi,   /* array args */
		       &PyArray_Type, &gv,    /* array args */
		       &PyArray_Type, &labels,    /* array args */
		       &label,
		       &sigma )) 
    return NULL;
  
  if (sigma > 0.0){
    weight = 0.5/sigma/sigma;
  }else {
    weight = -1;
  }


  if(ubi->nd != 2 || ubi->descr->type_num!=PyArray_DOUBLE ||
     ubi->dimensions[0] != 3 || ubi->dimensions[1] != 3){
    printf("first arg nd %d\n",ubi->nd);
    printf("first arg type %d\n",ubi->descr->type_num);
    PyErr_SetString(PyExc_ValueError,
		    "First array must be 3x3 2d and double");
    return NULL;
  }
  if(gv->nd != 2 || gv->descr->type_num!=PyArray_DOUBLE){
    PyErr_SetString(PyExc_ValueError,
		    "Second array must be 3xn 2d and double");
    return NULL;
  }
  if(gv->dimensions[1] != 3){
    PyErr_SetString(PyExc_ValueError,
		    "Second array must be 3xn 2d and double");
    return NULL;
  }
  if(labels->nd != 1 || labels->descr->type_num!=PyArray_INT32 ||
     labels->dimensions[0] != gv-> dimensions[0] ){
    PyErr_SetString(PyExc_ValueError,
		    "Label array must be integer and same length as gv");
    return NULL;
  }



   u00=* (double *) (ubi->data + 0*ubi->strides[0] + 0*ubi->strides[1]);
   u01=* (double *) (ubi->data + 0*ubi->strides[0] + 1*ubi->strides[1]);
   u02=* (double *) (ubi->data + 0*ubi->strides[0] + 2*ubi->strides[1]);
   u10=* (double *) (ubi->data + 1*ubi->strides[0] + 0*ubi->strides[1]);
   u11=* (double *) (ubi->data + 1*ubi->strides[0] + 1*ubi->strides[1]);
   u12=* (double *) (ubi->data + 1*ubi->strides[0] + 2*ubi->strides[1]);
   u20=* (double *) (ubi->data + 2*ubi->strides[0] + 0*ubi->strides[1]);
   u21=* (double *) (ubi->data + 2*ubi->strides[0] + 1*ubi->strides[1]);
   u22=* (double *) (ubi->data + 2*ubi->strides[0] + 2*ubi->strides[1]);

   sumsqtot = 0.0;
   n = 0;
   for(k=0;k<gv->dimensions[0];k++){ /* Loop over observed peaks */

     /* Are we using this peak? */
     if ( label != * (int *) (labels->data + k*labels->strides[0]) )
	  continue;

     n++; /* Number of peaks used */

     /* Compute hkls of this peak as h = UBI g */
     g0 = * (double *) (gv->data + k*gv->strides[0] + 0*gv->strides[1]);
     g1 = * (double *) (gv->data + k*gv->strides[0] + 1*gv->strides[1]);
     g2 = * (double *) (gv->data + k*gv->strides[0] + 2*gv->strides[1]);

     h0 = u00*g0 + u01*g1 + u02*g2;
     h1 = u10*g0 + u11*g1 + u12*g2;
     h2 = u20*g0 + u21*g1 + u22*g2;
      
     /*   From Paciorek et al Acta A55 543 (1999)
      *   UB = R H-1
      *    where:
      *    R = sum_n r_n h_n^t
      *    H = sum_n h_n h_n^t
      *    r = g-vectors
      *    h = hkl indices
      *    The hkl integer indices are: */
     ih[0] = conv_double_to_int_fast(h0);
     ih[1] = conv_double_to_int_fast(h1);
     ih[2] = conv_double_to_int_fast(h2);
     
     /* Error in fit */
     t0 = h0 - ih[0] ;
     t1 = h1 - ih[1] ;
     t2 = h2 - ih[2] ;

     sumsq = t0*t0 + t1*t1 + t2*t2 ;
     sumsqtot += sumsq;

     /* This is the soft cutoff factor */
     if ( weight > 0. ){
       wfac = exp( -sumsq*weight); /* see above weight == /2/sigma/sigma */ 
     } else {
       wfac = 1.0;
     }

     /* The g-vector was: */
     rh[0] = g0; 
     rh[1] = g1; 
     rh[2] = g2;
     
     for(i=0;i<3;i++){
       for(j=0;j<3;j++){
	 R[i][j] = R[i][j] + ih[j] * rh[i] * wfac;
	 H[i][j] = H[i][j] + ih[j] * ih[i] * wfac;
       }
     }
     
   }
   
   
   if (inverse3x3(H)==0){
     /* Form best fit UB */
     for(i=0;i<3;i++)
       for(j=0;j<3;j++)
	 for(l=0;l<3;l++)
	   UB[i][j] = UB[i][j] + R[i][l]*H[l][j];
     
     if (inverse3x3(UB)==0){
       /*printf("Got a new UB\n");
	 for(i=0;i<3;i++)
	 for(j=0;j<3;j++)
	 printf("UBi[%d][%d]=%f",i,j,UB[i][j]);
	 printf("\n");
       */
       /* Copy result into arg */
       for(i=0; i<3; i++)
	 for(j=0; j<3; j++)
	   *(double *)(ubi->data + 
		       i*ubi->strides[0] + 
		       j*ubi->strides[1]) = UB[i][j];

       
     }
     
   }
   
      
   return Py_BuildValue("id",n,sumsqtot/n);
}




/* weighted_refine  =========================================================== */

static PyObject *weighted_refine( PyObject *self, PyObject *args, PyObject *keywds){
   PyArrayObject *ubi=NULL, *gv=NULL, *w=NULL;
   double u00,u11,u22,u01,u02,u10,u12,u20,u21;
   double g0,g1,g2,h0,h1,h2,t0,t1,t2;
   double tol,sumsq,tolsq,sumdrlv2,wt;
   double R[3][3],H[3][3],ih[3],rh[3], UB[3][3];
   int n,k,i,j,l;

   for(i=0;i<3;i++){
     ih[i]=0.;
     rh[i]=0.;
     for(j=0;j<3;j++){
       R[i][j] = 0.;
       H[i][j] = 0.;
       UB[i][j] = 0.;
     }
   }

   
   if(!PyArg_ParseTuple(args,"O!O!dO!",
                        &PyArray_Type, &ubi,   /* array args */
                        &PyArray_Type, &gv,    /* array args */
                        &tol, /* Tolerance */
			&PyArray_Type, &w ))    /* array args */
     return NULL;
   
   if(ubi->nd != 2 || ubi->descr->type_num!=PyArray_DOUBLE){
     printf("first arg nd %d\n",ubi->nd);
     printf("first arg type %d\n",ubi->descr->type_num);
     PyErr_SetString(PyExc_ValueError,
		     "First array must be 3x3 2d and double");
     return NULL;
   }
   if(gv->nd != 2 || gv->descr->type_num!=PyArray_DOUBLE){
      PyErr_SetString(PyExc_ValueError,
		      "Second array must be 3xn 2d and double");
      return NULL;
   }
   
/*   if(gv->strides[0] != sizeof(double) || ubi->strides[0] != sizeof(double)){
      PyErr_SetString(PyExc_ValueError,
            "Arrays must be flat (strides == sizeof(double))");
      return NULL;
   } */

   if(ubi->dimensions[0] != 3 || ubi->dimensions[1] != 3){
      PyErr_SetString(PyExc_ValueError,
            "First array must be 3x3 2d, dimensions problem");
      return NULL;
   }

   if(gv->dimensions[1] != 3){
      PyErr_SetString(PyExc_ValueError,
            "Second array must be 3xn 2d and double");
      return NULL;
   }

   if(w->nd != 1 || w->descr->type_num!=PyArray_DOUBLE){
     PyErr_SetString(PyExc_ValueError,
            "Weight array must be 1D float");
      return NULL;
   }

   if(w->dimensions[0] != gv->dimensions[0]){
     PyErr_SetString(PyExc_ValueError,
            "Weight array must be dimension n where n matches gv");
      return NULL;
   }



/*   if(gv->strides[1] != sizeof(double)){
      PyErr_SetString(PyExc_ValueError,
            "Second array must be 3xn 2d and double, and you want it aligned!!!");
      return NULL;
   }
*/


   u00=* (double *) (ubi->data + 0*ubi->strides[0] + 0*ubi->strides[1]);
   u01=* (double *) (ubi->data + 0*ubi->strides[0] + 1*ubi->strides[1]);
   u02=* (double *) (ubi->data + 0*ubi->strides[0] + 2*ubi->strides[1]);
   u10=* (double *) (ubi->data + 1*ubi->strides[0] + 0*ubi->strides[1]);
   u11=* (double *) (ubi->data + 1*ubi->strides[0] + 1*ubi->strides[1]);
   u12=* (double *) (ubi->data + 1*ubi->strides[0] + 2*ubi->strides[1]);
   u20=* (double *) (ubi->data + 2*ubi->strides[0] + 0*ubi->strides[1]);
   u21=* (double *) (ubi->data + 2*ubi->strides[0] + 1*ubi->strides[1]);
   u22=* (double *) (ubi->data + 2*ubi->strides[0] + 2*ubi->strides[1]);

   n=0;
   sumdrlv2 = 0;
   tolsq=tol*tol;

   for(k=0;k<gv->dimensions[0];k++){ /* Loop over observed peaks */
      /* Compute hkls of this peak as h = UBI g */
      g0 = * (double *) (gv->data + k*gv->strides[0] + 0*gv->strides[1]);
      g1 = * (double *) (gv->data + k*gv->strides[0] + 1*gv->strides[1]);
      g2 = * (double *) (gv->data + k*gv->strides[0] + 2*gv->strides[1]);
      h0 = u00*g0 + u01*g1 + u02*g2;
      h1 = u10*g0 + u11*g1 + u12*g2;
      h2 = u20*g0 + u21*g1 + u22*g2;
      
      t0=h0-conv_double_to_int_fast(h0);
      t1=h1-conv_double_to_int_fast(h1);
      t2=h2-conv_double_to_int_fast(h2);
      sumsq = t0*t0+t1*t1+t2*t2;

      if ( (t0 > 0.501) ||
	   (t1 > 0.501) ||
	   (t2 > 0.501) ){
	printf("Error in conv\nh0 = %f h1 = %f h2=%f \nt0 = %f t1 = %f t2=%f\n",
	       t0,t1,t2,h0,h1,h2);
	return NULL;
      }       
	  
      
      if (sumsq < tolsq){
	n=n+1;
	sumdrlv2 += sumsq;
	/*   From Paciorek et al Acta A55 543 (1999)
	 *   UB = R H-1
	 *    where:
	 *    R = sum_n r_n h_n^t
	 *    H = sum_n h_n h_n^t
	 *    r = g-vectors
	 *    h = hkl indices
	 *    The hkl integer indices are: */
	ih[0] = conv_double_to_int_fast(h0);
	ih[1] = conv_double_to_int_fast(h1);
	ih[2] = conv_double_to_int_fast(h2);
	/* The g-vector was: */
	rh[0] = g0; rh[1] = g1 ; rh[2] = g2;
	wt = *( double*) (w->data  + k*w->strides[0]);
	for(i=0;i<3;i++){
	  for(j=0;j<3;j++){
	    /* Robust weight factor, fn(tol), would go here */
	    R[i][j] = R[i][j] + wt * ih[j] * rh[i];
	    H[i][j] = H[i][j] + wt * ih[j] * ih[i];
	  }
	}
      }
   }


   if (inverse3x3(H)==0){
     /* Form best fit UB */
     for(i=0;i<3;i++)
       for(j=0;j<3;j++)
	 for(l=0;l<3;l++)
	   UB[i][j] = UB[i][j] + R[i][l]*H[l][j];
     
     if (inverse3x3(UB)==0){
       /*printf("Got a new UB\n");
	 for(i=0;i<3;i++)
	 for(j=0;j<3;j++)
	 printf("UBi[%d][%d]=%f",i,j,UB[i][j]);
	 printf("\n");
       */

       /* Copy result into arg */
       for(i=0; i<3; i++)
	 for(j=0; j<3; j++)
	   *(double *)(ubi->data + 
		       i*ubi->strides[0] + 
		       j*ubi->strides[1]) = UB[i][j];

       
     }
     
   }
   
   if( n > 0){    
       return Py_BuildValue("id",n,sumdrlv2/n);
    }   else  {
       return Py_BuildValue("id",n,0);
    }
}

/* ================================================================ */



void ptype(int type){
   /* Print out the type of a Numeric array from the C point of view */
  printf("Your input type was ");
  switch (type){
  case    NPY_CHAR   : 
    printf("PyArray_CHAR *(char *)\n");
    break;
  case    NPY_SHORT  : 
    printf("PyArray_SHORT *(short *)\n");
    break;
  case    NPY_INT    : 
    printf("PyArray_INT *(int  *)\n");
    break;
  case    NPY_LONG   : 
    printf("PyArray_LONG *(long *)\n");
    break;
  case    NPY_FLOAT  : 
    printf("NPY_FLOAT *(float *)\n");
    break;
  case    NPY_DOUBLE : 
    printf("PyArray_DOUBLE *(double *)\n");
    break;
  case    NPY_UBYTE  : 
    printf("PyArray_UBYTE *(unsigned char *)\n");
    break;
  case    NPY_USHORT : 
    printf("PyArray_USHORT *(unsigned short*)\n");
    break;
  case    NPY_UINT   : 
    printf("PyArray_UINT *(unsigned int *)\n");
    break;
     }
}


static PyObject *put_incr( PyObject *self, 
			   PyObject *args, 
			   PyObject *keywds){
  PyArrayObject *data=NULL, *ind=NULL, *vals=NULL;
  int i;
  npy_intp j;
  float f;
  if(!PyArg_ParseTuple(args,"O!O!O!",
		       &PyArray_Type, &data,   /* array args */
		       &PyArray_Type, &ind,    /* array args */
		       &PyArray_Type, &vals     /* array args */
		       )) 
    return NULL;

   if(data->nd != 1 || data->descr->type_num!=NPY_FLOAT){
     PyErr_SetString(PyExc_ValueError,
		     "First array must be 1d and float32");
     return NULL;
   }
   if(ind->nd != 1 ){
     PyErr_SetString(PyExc_ValueError,
		     "Second array must be 1d");
     return NULL;
   }
   if( PyArray_EquivTypenums( NPY_INTP, data->descr->type_num  ) ){
     printf("ind->descr->type_num = %d , NPY_INTP %d, NPY_INT32 %d, NPY_INT %d\n",
             ind->descr->type_num, NPY_INTP, NPY_INT32, NPY_INT);
     ptype(ind->descr->type_num);       
     PyErr_SetString(PyExc_ValueError,
		     "Second array must be int pointer type");
     return NULL;
   }
   if(vals->nd != 1 || vals->descr->type_num!=NPY_FLOAT){
     PyErr_SetString(PyExc_ValueError,
		     "Third array must be 1d and float32");
     return NULL;
   }
   if(vals->dimensions[0] != ind->dimensions[0]){
     PyErr_SetString(PyExc_ValueError, 
		     "Second array must have same dims as first");
     return NULL;
   }
   /* Loop over points adding with increment */
   for(i = 0 ; i < ind->dimensions[0] ; i++){
     /* i applies to indices and vals */
     j = *(long *) ( ind->data + i* ind->strides[0]);
     
     if ( j < 0 || j > data->dimensions[0]){
       PyErr_SetString(PyExc_ValueError, 
		       "Array bounds exceeded");
       return NULL;
     }

     f =* (float *) (data->data + j*data->strides[0]) ; 
     * (float *) (data->data + j*data->strides[0]) = f + 
               * (float *) (vals->data + i*vals->strides[0]);
     /*     if( i < 20){
       printf("i %d  j %d f %f  tot %f\n",i,j,f,
	      * (float *) (vals->data + i*vals->strides[0]));
	      }*/
   }
   /* Nothing to return */
   Py_INCREF(Py_None);
   return Py_None;

}





int inverse3x3 ( double H[3][3]){
  double det, inverse[3][3];
  int i,j;
      /*
      # | a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
      # | a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
      # | a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |
      
      # DET=a11   (a33     a22    -a32     a23)-  
            a21   (a33      a12   -a32     a13)+      
            a31   (a23     a12    -a22     a13)
      */

      det=H[0][0]*(H[2][2]*H[1][1]-H[2][1]*H[1][2])-
          H[1][0]*(H[2][2]*H[0][1]-H[2][1]*H[0][2])+
 	  H[2][0]*(H[1][2]*H[0][1]-H[1][1]*H[0][2])  ;

      if (det != 0.) {


        inverse[0][0] =   (H[2][2]*H[1][1]-H[2][1]*H[1][2])/det;
	inverse[0][1] = -(H[2][2]*H[0][1]-H[2][1]*H[0][2])/det;
	inverse[0][2] =  (H[1][2]*H[0][1]-H[1][1]*H[0][2])/det; 
        inverse[1][0] =  -(H[2][2]*H[1][0]-H[2][0]*H[1][2])/det;
        inverse[1][1] =  (H[2][2]*H[0][0]-H[2][0]*H[0][2])/det;
	inverse[1][2] = -(H[1][2]*H[0][0]-H[1][0]*H[0][2])/det;
	inverse[2][0] = (H[2][1]*H[1][0]-H[2][0]*H[1][1])/det;
	inverse[2][1] = -(H[2][1]*H[0][0]-H[2][0]*H[0][1])/det;
        inverse[2][2] =  (H[1][1]*H[0][0]-H[1][0]*H[0][1])/det;

	for(i=0;i<3;i++)
	  for(j=0;j<3;j++)
	    H[i][j]=inverse[i][j];

	return 0;

      }
      else{
	return -1;
      }
}




/* TODO - Make a UBI matrix from two peaks and score it ?? */

   
   
static PyMethodDef closestMethods[] = {
   { "closest", (PyCFunction) closest, METH_VARARGS, 
     "int,float closest(x,v)\n"				\
     "    Find element in x closest to one in v\n"	\
     "    x,v 1D flat Float arrays"},
   { "score", (PyCFunction) score, METH_VARARGS, 
     "int score(ubi,gv,tol)\n"					\
     "    Find number of gv which are integer h within tol"},
   { "score_and_refine", (PyCFunction) score_and_refine, METH_VARARGS, 
     "int score_and_refine(ubi,gv,tol)\n"				\
     "    Find number of gv which are integer h within tol and overwrite UB with refined"},
   { "score_and_assign", (PyCFunction) score_and_assign, METH_VARARGS, 
     "int score_and_assign(ubi, gv, tol, drlv2, labels, label)\n"	\
     "    Find number of gv which have lower drlv2 than stored values in drlv2.\n" \
     "    Label them as label in labels array\n"},
   { "refine_assigned", (PyCFunction) refine_assigned, METH_VARARGS, 
     "int, float refine_assigned(ubi, gv, labels, label, sigma)\n"	\
     "    refines ubi to match gv where labels==label\n" \
     "    sigma is a soft cutoff on tolerance (drlv)\n" \
     "    returns npks, <drlv2> "},
   { "weighted_refine", (PyCFunction) weighted_refine, METH_VARARGS, 
     "int, float weighted_refine(ubi, gv, tol, weights )\n"	\
     "    refines ubi to match gv \n" \
     "    weights is the weighting to apply to peaks (eg, intensity)\n" \
     "    returns npks, <drlv2> "},
   { "put_incr", (PyCFunction) put_incr, METH_VARARGS,
     " None = put_incr( data, ind, vals )\n"\
     "    does the loop \n"\
     "    for i, v in ind, vals:\n"\
     "        data[ ind[i] ] += vals[i]\n"},
     
     



   { NULL, NULL, 0, NULL }  
};

void initclosest(void)
{
   PyObject *m, *d, *s;
   m=Py_InitModule("closest",closestMethods);
   import_array();
   d=PyModule_GetDict(m);
   s = PyString_FromString(moduledocs);
   PyDict_SetItemString(d,"__doc__",s);
   Py_DECREF(s);
   if (PyErr_Occurred())
      Py_FatalError("cant initialise closest module");
}
     


