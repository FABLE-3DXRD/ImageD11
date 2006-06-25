



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

static char moduledocs [] =\
"   _hist( data=Numeric.array(2D,UInt16) , hist=Numeric.array(1D,Int) , [Verbose=0] )\n"\
"   nothing returned (None)\n"\
"   Computes histogram (overwriting hist arg) of data \n"\
"   Assumes hist is dimensioned to hold all necessary bins (checks only with verbose)\n"\
"   \n"\
"   Good chance of dumping core, so use ImageD11/hist.py for interface ";
   

#include <Python.h>                  /* To talk to python */
#include "Numeric/arrayobject.h"     /* Access to Numeric */


/* Fill in a histogram of pixels */

static PyObject * _hist (PyObject *self, PyObject *args,  PyObject *keywds)
{
   PyArrayObject *dat=NULL; /* in (not modified) */
   PyArrayObject *hst=NULL; /* out ( modified) */
   static char *kwlist[] = {"data","hist","verbose", NULL};
   int verbose,i,j,k,f,s;
   verbose = 0;
   if(!PyArg_ParseTupleAndKeywords(args,keywds, "O!O!|i",kwlist,      
				   &PyArray_Type, &dat,   /* array arg */
				   &PyArray_Type, &hst,   /* array arg */
				   &verbose))        /* optional verbosity */
     return NULL;

   if(verbose!=0){
     printf("\n\nHello there from _hist, you wanted the verbose output...\n");
      }

   /* Check array is two dimensional and Ushort */
   if(dat->nd != 2 ){    
     PyErr_SetString(PyExc_ValueError,
		     "Data array must be 2d,first arg problem");
     return NULL;
   }
   if(dat->descr->type_num != PyArray_USHORT){
     PyErr_SetString(PyExc_ValueError,
		     "Data array must have type USHORT (sorry for now)");
     return NULL;
   }


   /* Check hst is one dimensional and Int */
   if(hst->nd != 1 && hst->descr->type_num != PyArray_INT){    
     PyErr_SetString(PyExc_ValueError,
		     "Hist array must be 1d, Int,first arg problem");
     return NULL;
   }
   /* Zero the histogram */
   for(i=0;i<hst->dimensions[0];i++)
     (*(int *)(hst->data + i*hst->strides[0])) = 0 ;

   /* Decide on fast/slow loop - inner versus outer */
   if(dat->strides[0] > dat->strides[1]) {
      f=1;  s=0;}
   else {
      f=0;  s=1;
   }
   if (verbose!=0){
        printf("Fast index is %d, slow index is %d, ",f,s);
        printf("strides[0]=%d, strides[1]=%d\n",dat->strides[0],dat->strides[1]);
   }

   if(verbose==0){
     printf("Starting\n");

     for( i = 0 ; i < dat->dimensions[s] ; i++ ){    /* i,j is looping along the indices data array */
       for( j = 0 ; j < dat->dimensions[f] ; j++ ){
	 k = (*(unsigned short *)(dat->data+ i*dat->strides[s] + j*dat->strides[f])); 
	 (*(int *)(hst->data + k*hst->strides[0])) = (*(int *)(hst->data + k*hst->strides[0])) + 1 ;
       }
     }

   }else{

     for( i = 0 ; i < dat->dimensions[s] ; i++ ){    /* i,j is looping along the indices data array */
       for( j = 0 ; j < dat->dimensions[f] ; j++ ){
	 k = (*(unsigned short *)(dat->data+ i*dat->strides[s] + j*dat->strides[f])); 
	 if (k>0 && k<hst->dimensions[0]){
	   (*(int *)(hst->data + k*hst->strides[0])) = (*(int *)(hst->data + k*hst->strides[0])) + 1 ;
	 }else{
	   PyErr_SetString(PyExc_ValueError,"Data does not fit in hist, check max/min values");
	   return NULL;
	 }
       }
     }

   }

   /* Arg modified so nothing to give back */
    Py_INCREF(Py_None);
    return Py_None;


}


static PyMethodDef _histMethods[] = {
   {"_hist", (PyCFunction) _hist, METH_VARARGS | METH_KEYWORDS,
    moduledocs},
   {NULL, NULL, 0, NULL} /* setinel */
};

void 
init_hist(void)
{
   PyObject *m, *d;

   m=Py_InitModule("_hist", _histMethods);
   import_array();
   d=PyModule_GetDict(m);
}
