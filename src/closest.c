

/* 
# ImageD11_v0.4 Software for beamline ID11
# Copyright (C) 2005  Jon Wright
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


/* *******************************************************************
 * closest.c  Two useful functions for indexing grains 
 *
 * 1. From python: closest(cosines_allowed, cosine_measured)
 *    Take 1 peak with hkl's assigned, compute the cosines to
 *    another powder ring (angles between this peak and peaks in
 *    that ring). This extension finds the best peak in the
 *    other ring to pair up with.
 *
 *    More generally closest(array, values) finds the closest
 *    item in values to one of the values in array.  Returns the
 *    difference and the index of the closest thing it finds.
 *
 *    Both arguments should be one dimensional Numeric arrays
 *    of type Numeric.Float
 *
 * 2. From python: score(ubi, gv, tol) where ubi is an orientation
 *    matrix and gv are an array of g-vectors. Returns the number
 *    of g-vectors which have integer hkl peaks within tolerance
 *    tol. Uses the conv_double_to_int_fast function in here for 
 *    factor of !EIGHT! speed increase compared to rounding in 
 *    C. In fact this gives the nearest even integer, instead
 *    of the nearest integer, but we don't care, as a peak having
 *    hkl of 0.5 is nowhere near being indexed anyway.
 *
 *    Returns and integer - number of peaks indexed
 *    UBI is a 3x3 Numeric.Float array (figure out the transposing yourself)
 *    GV is a nx3 Numeric.Float array, and you should try to make the 3 
 *       be the fast index for best performance
 *
 * ****************************************************************** */


#include <Python.h>                  /* To talk to python */
#include "Numeric/arrayobject.h"     /* Access to Numeric */

inline int conv_double_to_int_fast(double);
inline int conv_double_to_int_safe(double);
int inverse3x3(double A[3][3]);

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
/*    assert(t[0]<0.501);
      assert(t[1]<0.501);
      assert(t[2]<0.501);      */
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
   
   
   return Py_BuildValue("i",n);
}


static PyObject *score_and_refine( PyObject *self, PyObject *args, PyObject *keywds){
   PyArrayObject *ubi=NULL, *gv=NULL;
   double u00,u11,u22,u01,u02,u10,u12,u20,u21;
   double g0,g1,g2,h0,h1,h2,t0,t1,t2;
   double tol,sumsq;
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
      sumsq = t0*t0+t1*t1+t2*t2;

      if (sumsq < tol){
         n=n+1;
	 /* 	    From Paciorek et al Acta A55 543 (1999)
	    UB = R H-1
	    where:
	    R = sum_n r_n h_n^t
	    H = sum_n h_n h_n^t
	    r = g-vectors
	    h = hkl indices
	    The hkl integer indices are: */
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
     
	/* Return inverse in argument 
	   for(i=0;i<3;i++)
	   for(j=0;j<3;j++)
	   printf("UBi[%d][%d]=%f",i,j,UB[i][j]);
	*/
     if (inverse3x3(UB)==0){
       /*printf("Got a new UB\n");
       for(i=0;i<3;i++)
	  for(j=0;j<3;j++)
	    printf("UBi[%d][%d]=%f",i,j,UB[i][j]);
       printf("\n");
       */

       * (double *) (ubi->data + 0*ubi->strides[0] + 0*ubi->strides[1]) = UB[0][0];
       * (double *) (ubi->data + 0*ubi->strides[0] + 1*ubi->strides[1]) = UB[0][1];
       * (double *) (ubi->data + 0*ubi->strides[0] + 2*ubi->strides[1]) = UB[0][2];
       * (double *) (ubi->data + 1*ubi->strides[0] + 0*ubi->strides[1]) = UB[1][0];
       * (double *) (ubi->data + 1*ubi->strides[0] + 1*ubi->strides[1]) = UB[1][1];
       * (double *) (ubi->data + 1*ubi->strides[0] + 2*ubi->strides[1]) = UB[1][2];
       * (double *) (ubi->data + 2*ubi->strides[0] + 0*ubi->strides[1]) = UB[2][0];
       * (double *) (ubi->data + 2*ubi->strides[0] + 1*ubi->strides[1]) = UB[2][1];
       * (double *) (ubi->data + 2*ubi->strides[0] + 2*ubi->strides[1]) = UB[2][2];
       
     }
     
   }
   
      
   return Py_BuildValue("i",n);
}



inline int conv_double_to_int_safe(double x){
   int a;
   a = floor(x+0.5);
   return a;
}



inline int conv_double_to_int_fast(double x){
   /* This was benched as about eight times faster than the safe mode!! */

   /* Put in the reference for where this was found on the web TODO */
   const int p=52;
   const double c_p1 = (1L << (p/2));
   const double c_p2 = (1L << (p-p/2));
   const double c_mul = c_p1 * c_p2;
   const double cs = 1.5*c_mul;
   x+=cs;
   const int a = *(int *)(&x);
   return (a);
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
      "int,float closest(x,v)\n"\
      "    Find element in x closest to one in v\n"\
      "    x,v 1D flat Float arrays"},
   { "score", (PyCFunction) score, METH_VARARGS, 
      "int score(ubi,gv,tol)\n"\
      "    Find number of gv which are integer h within tol"},
   { "score_and_refine", (PyCFunction) score_and_refine, METH_VARARGS, 
      "int score_and_refine(ubi,gv,tol)\n"\
      "    Find number of gv which are integer h within tol and returned refined UB"},
   { NULL, NULL, 0, NULL }  
};

void initclosest(void)
{
   PyObject *m, *d;
   m=Py_InitModule("closest",closestMethods);
   import_array();
   d=PyModule_GetDict(m);
}
     

