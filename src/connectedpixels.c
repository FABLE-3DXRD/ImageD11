



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


#include <Python.h>                  /* To talk to python */
#include "Numeric/arrayobject.h"     /* Access to Numeric */
#include "dset.h"
#include <time.h>
#include <unistd.h>
#include <stdlib.h>
/*
int * dset_initialise(int size); 
int * dset_new(int ** S);
void dset_makeunion(int * S, int r1, int r2);
void dset_link(int *S, int r1, int r2);
int dset_find(int x, int * S);
*/

#undef PARANOID

int getval(char *p,int type);

int getval(char *p,int type){
  switch (type){
     case    PyArray_CHAR   : return *(char          *)p*1;
     case    PyArray_SBYTE  : return *(signed char   *)p*1;
     case    PyArray_SHORT  : return *(short         *)p*1;
     case    PyArray_INT    : return *(int           *)p*1;
     case    PyArray_LONG   : return *(long          *)p*1;
     case    PyArray_FLOAT  : return *(float         *)p*1;
     case    PyArray_DOUBLE : return *(double        *)p*1;
#ifdef PyArray_UNSIGNED_TYPES
     case    PyArray_UBYTE  : return *(unsigned char *)p*1;
     case    PyArray_USHORT : return *(unsigned short*)p*1;
     case    PyArray_UINT   : return *(unsigned int  *)p*1;
#endif
     }
   printf("Oh bugger in getval - unrecongnised numeric type\n");
   exit(1);
   return 0;
}

void ptype(int type){
  printf("Your input type was ");
  switch (type){
     case    PyArray_CHAR   : printf("PyArray_CHAR *(char *)\n");break;
     case    PyArray_SBYTE  : printf("PyArray_SBYTE *(signed char *)\n");break;
     case    PyArray_SHORT  : printf("PyArray_SHORT *(short *)\n");break;
     case    PyArray_INT    : printf("PyArray_INT *(int  *)\n");break;
     case    PyArray_LONG   : printf("PyArray_LONG *(long *)\n");break;
     case    PyArray_FLOAT  : printf("PyArray_FLOAT *(float *)\n");break;
     case    PyArray_DOUBLE : printf("PyArray_DOUBLE *(double *)\n");break;
#ifdef PyArray_UNSIGNED_TYPES
     case    PyArray_UBYTE  : printf("PyArray_UBYTE *(unsigned char *)\n");break;
     case    PyArray_USHORT : printf("PyArray_USHORT *(unsigned short*)\n");break;
     case    PyArray_UINT   : printf("PyArray_UINT *(unsigned int *)\n");break;
#endif
     }
}


/* make an image of peak assignement for pixels */
static PyObject * connectedpixels (PyObject *self, PyObject *args,  PyObject *keywds)
{
   PyArrayObject *dataarray=NULL,*results=NULL; /* in (not modified) and out */

   int i,j,k,l,ival,f,s,np;        
   int *T;                         /* for the disjoint set copy */
   int verbose=0,type;                /* whether to print stuff and the type of the input array */
   int percent,npover;                       /* for the progress indicator in verbose mode */
   static char *kwlist[] = {"data","results","threshold","verbose", NULL};
   clock_t tv1,tv1point5,tv2,tv3,tv4;
   double time;
   float threshold;
   if(!PyArg_ParseTupleAndKeywords(args,keywds, "O!O!f|i",kwlist,      
                        &PyArray_Type, &dataarray,   /* array args */
                        &PyArray_Type, &results,   /* array args */
                        &threshold, &verbose))        /* threshold and optional verbosity */
      return NULL;
   if(verbose!=0){
      printf("\n\nHello there from connected pixels, you wanted the verbose output...\n");
      }

   tv1=clock();  
   /* Check array is two dimensional and float */
   if(dataarray->nd != 2){    
      PyErr_SetString(PyExc_ValueError,
                       "Array must be 2d, first arg problem");
      return NULL;
      }
   type=dataarray->descr->type_num;
   if(verbose!=0)ptype(type);
   if(verbose!=0)printf("Thresholding at level %f\n",threshold);
   /* make an array to hold the integer peak assignments */
   /*   PyArray_FromDims(int n_dimensions, int dimensions[n_dimensions], int type_num) */
   /* Make a new array to hold the results... how to do optional arguments?? */
   /*   results=(PyArrayObject*)PyArray_FromDims(dataarray->nd, dataarray->dimensions, PyArray_INT); */
   /*   if(results==NULL){*/
   /*   printf("Could not make a results array\nGiving up\n");*/
   /*   exit(2);   } */
   /*   if(verbose!=0)printf("Made results array\n");*/

   if(results->nd != 2){    
      PyErr_SetString(PyExc_ValueError,
                       "Array must be 2d, second arg problem");
      return NULL;
      }
   if( results->dimensions[0] != dataarray->dimensions[0] ||
       results->dimensions[1] != dataarray->dimensions[1] ){    
      PyErr_SetString(PyExc_ValueError,
                       "Arrays must have same shape");
      return NULL;
      }

   S = dset_initialise(16384); /* Default number before reallocation */
   if(verbose!=0)printf("Initialised the disjoint set\n");
   npover=0; /* number of pixels over threshold */
   /* Decide on fast/slow loop - inner versus outer */
   if(dataarray->strides[0] > dataarray->strides[1]) {
      f=1;  s=0;}
   else {
      f=0;  s=1;
   }
   if (verbose!=0){
        printf("Fast index is %d, slow index is %d, ",f,s);
   printf("strides[0]=%d, strides[1]=%d\n",dataarray->strides[0],dataarray->strides[1]);
   }
   percent=(results->dimensions[s]-1)/80.0;
   if(percent < 1)percent=1;
   tv1point5=clock();   
   if(verbose!=0){
     time=(tv1point5-tv1)/CLOCKS_PER_SEC;
     printf("Ready to scan image, setup took %g and clock is only good to ~0.015 seconds\n",time);
     }

   if(verbose!=0)printf("Scanning image\n");
   for( i = 0 ; i < (results->dimensions[s]) ; i++ ){    /* i,j is looping along the indices data array */
      if(verbose!=0 && (i%percent == 0) )printf(".");


      for( j = 0 ; j < (results->dimensions[f]) ; j++ ){

	/* Set result for this pixel to zero */
	
	(*(int *)(results->data + i*results->strides[s] + j*results->strides[f])) = 0;

         ival= (* (unsigned short *) (dataarray->data + i*dataarray->strides[s] + j*dataarray->strides[f])) ;
	 /* val=getval((dataarray->data + i*dataarray->strides[s] + j*dataarray->strides[f]),type); */
         if( ival > threshold) {
             npover++;
             k=0;l=0;
                /* peak needs to be assigned */
	     /* i-1, j-1, assign same index */
        if(i!=0 && j!=0)
             if( ( k=*(int *)(results->data + (i-1)*results->strides[s] + (j-1)*results->strides[f] )) >0 ) {
                 (*(int *)(results->data + i*results->strides[s] + j*results->strides[f])) = k;
       l=k; /* l holds the assignment of the current pixel */
        }
           /* i-1,j */
        if(i!=0)
             if( ( k=*(int *)(results->data + (i-1)*results->strides[s] + j*results->strides[f] )) >0 ) {
            if(l==0){
                 (*(int *)(results->data + i*results->strides[s] + j*results->strides[f])) = k;
       l=k; /* l holds the assignment of the current pixel */
       } else {
       if(l!=k)dset_makeunion(S,k,l);
       }
        }
                /* i-1, j+1 */
        if(i!=0 && j!=(results->dimensions[f]-1))  
             if( ( k=*(int *)(results->data + (i-1)*results->strides[s] + (j+1)*results->strides[f] )) >0 ) {
            if(l==0){
                 (*(int *)(results->data + i*results->strides[s] + j*results->strides[f])) = k;
       l=k; /* l holds the assignment of the current pixel */
       } else {
       if(l!=k)dset_makeunion(S,k,l);
       }
        }
           /* i, j-1 */
        if(j!=0)
             if( ( k=*(int *)(results->data + i*results->strides[s] + (j-1)*results->strides[f] )) >0 ) {
            if(l==0){
                 (*(int *)(results->data + i*results->strides[s] + j*results->strides[f])) = k;
       l=k; /* l holds the assignment of the current pixel */
       } else {
       if(l!=k)dset_makeunion(S,k,l);
       }
        }
             if(l==0){ /* pixel has no neighbours thus far */
            S = dset_new(&S);
            (*(int *)(results->data + i*results->strides[s] + j*results->strides[f])) = S[S[0]-1];
        } 
    } /* active pixel */
    else { /* inactive pixel  - zero in results */
                 (*(int *)(results->data + i*results->strides[s] + j*results->strides[f])) = 0;
    }
      } /* j */
   } /* i */
   tv2=clock();   
   if(verbose!=0){
     printf("\nFinished scanning image, now make unique list and sums\n");
     printf("Number of pixels over threshold was : %d\n",npover);
     time=(tv2-tv1)/CLOCKS_PER_SEC;
     printf("That took %f seconds and %f in total being aware that the clock is only good to ~0.015 seconds\n",(tv2-tv1point5)/CLOCKS_PER_SEC,time);
     }
   /* First, count independent peaks */

   /* count peaks */
   np=S[S[0]-1];
   /* Now make each T[i] contain the unique ascending integer for the set */ 
   T=NULL;
   if(verbose!=0)printf("Now I want to malloc T=%p at size np=%d and sizeof(int)=%d\n",T,np,sizeof(int));
   if( ( T=(int *) ( malloc( (np+3)*sizeof(int) ) ) ) == NULL){
      printf("Memory allocation error in connected pixels\n");
      return 0;
   }
   if(verbose!=0)printf("Now T should be allocated, it is at %p\n",T);
   np=0;
   for(i=1;i<S[S[0]-1]+1;i++){
     if(S[i]==i){
      np++;
   T[i]=np;
   }
     else{  /* check */
        j=dset_find(i,S);
   T[i]=T[j];
   if(j>=i && verbose){
            printf("Oh dear - there was a problem compressing the disjoint set, j=%d, i=%d \n",j,i);
       }
   if(S[j]!=j && verbose!=0){
       printf("Oh dear - the disjoint set is squiff,S[j]=%d,j=%d\n",S[j],j);
       }
   }
     }
   if(verbose!=0)printf("\n");
   tv3=clock();
   if(verbose!=0){
     printf("Compressing list, time=%g, total time=%g\n",(tv3-tv2)/CLOCKS_PER_SEC,(tv3-tv1)/CLOCKS_PER_SEC);
     }

   if(verbose)printf("Found %d peaks\nNow assigning unique identifiers\n",np);


   /* Now loop through the results image again, assigning the corrected numbers to each pixel */
   for( i = 0 ; i < (results->dimensions[s]) ; i++ ){    /* i,j is looping along the indices data array */
      for( j = 0 ; j < (results->dimensions[f]) ; j++ ){
         ival=*(int *)(results->data + i*results->strides[s] + j*results->strides[f]);
         if(  ival > 0){
      if(T[ival]==0){
         printf("Something buggered up, ival=%d, T[ival]=%d, dset_find(ival,S)=%d\n",ival,T[ival],dset_find(ival,S));
         /* free(S); */
      exit(2);
      return NULL;
      }
      (*(int *)(results->data + i*results->strides[s] + j*results->strides[f]))=T[ival];
         }
      }
   }



 
   if(verbose!=0)printf("Freeing S=%p and T=%p\n",S,T);
   free(S);
   free(T);
   if(verbose!=0){
     printf("Freed S=%p and T=%p\n",S,T);
     printf("Finished assigning unique values and tidied up, now returning");
     tv4=clock();
     time=(tv4-tv3)/CLOCKS_PER_SEC;
     printf("\nThat took %g seconds, being aware that the clock is only good to ~0.015 seconds\n",time);
     printf("The total time in this routine was about %f\n",(tv4-tv1)/CLOCKS_PER_SEC); 
     }
   
/*   return PyArray_Return(results); */
/*   Py_DECREF(results);*/
     /*   Py_DECREF(dataarray);*/
   return Py_BuildValue("i", np); 
}


/* BLOB PROPERTIES BIT - things to measure are: */
/* Centre of mass (x,y) - optionally use x,y image */
/* Amount of intensity in blob */
/* std dev of intensity in blob */
/* ...? */


static PyObject * blobproperties (PyObject *self, PyObject *args,  PyObject *keywds)
{
   PyArrayObject *dataarray=NULL,*blobarray=NULL; /* in (not modified) */
   static char *kwlist[] = {"data","blob","npeaks","verbose", NULL};
   PyArrayObject *sum,*sumsq,*com0,*com1,*npix, *com00, *com11, *com01;
   double *asum,*asumsq,*acom0,*acom1, *acom00, *acom11, *acom01 ;
   double fval;
   int np,verbose=0,type,f,s,peak,bad;
   int i,j,safelyneed,*anpix, percent;
   if(!PyArg_ParseTupleAndKeywords(args,keywds, "O!O!i|i",kwlist,      
                        &PyArray_Type, &dataarray,   /* array args - data */
         &PyArray_Type, &blobarray,   /* blobs */
         &np,              /* Number of peaks to treat */
                        &verbose))        /* threshold and optional verbosity */
      return NULL;

   /* Check array is two dimensional */
   if(dataarray->nd != 2){     
      PyErr_SetString(PyExc_ValueError,
                       "data array must be 2d, first arg problem");
      return NULL;
      }
   if(verbose!=0)printf("Welcome to blobproperties\n");
   /* Check array is two dimensional and float */
   if(blobarray->nd != 2 && blobarray->descr->type_num != PyArray_INT){     
      PyErr_SetString(PyExc_ValueError,
                       "Blob array must be 2d and integer, second arg problem");
      return NULL;
      }
   type=dataarray->descr->type_num;
   /* Decide on fast/slow loop - inner versus outer */
   if(blobarray->strides[0] > blobarray->strides[1]) {
      f=1;  s=0;}
   else {
      f=0;  s=1;
   }
   if (verbose!=0){
        printf("Fast index is %d, slow index is %d, ",f,s);
   printf("strides[0]=%d, strides[1]=%d\n",dataarray->strides[0],dataarray->strides[1]);
   }
  /* results arrays */
  safelyneed=np+3;
  npix  = (PyArrayObject *)PyArray_FromDims(1,&safelyneed,PyArray_INT);
  sumsq = (PyArrayObject *)PyArray_FromDims(1,&safelyneed,PyArray_DOUBLE);
  sum   = (PyArrayObject *)PyArray_FromDims(1,&safelyneed,PyArray_DOUBLE);
  com0  = (PyArrayObject *)PyArray_FromDims(1,&safelyneed,PyArray_DOUBLE);
  com1  = (PyArrayObject *)PyArray_FromDims(1,&safelyneed,PyArray_DOUBLE);
  com00 = (PyArrayObject *)PyArray_FromDims(1,&safelyneed,PyArray_DOUBLE);
  com01 = (PyArrayObject *)PyArray_FromDims(1,&safelyneed,PyArray_DOUBLE);
  com11 = (PyArrayObject *)PyArray_FromDims(1,&safelyneed,PyArray_DOUBLE);

  if (npix == NULL || sumsq == NULL || sum == NULL || com0 == NULL || com1 == NULL) goto fail;
  anpix  = (int   *)  npix->data;
  asumsq = (double *) sumsq->data;
  asum   = (double *)   sum->data;
  acom0  = (double *)  com0->data;
  acom1  = (double *)  com1->data;
  acom00 = (double *)  com00->data;
  acom01 = (double *)  com01->data;
  acom11 = (double *)  com11->data;

  for ( i=0 ; i<safelyneed ; i++){
    anpix[i]  = 0;
    asumsq[i] = 0.;
    asum[i]   = 0.;
    acom0[i]  = 0.;
    acom1[i]  = 0.;
    acom00[i] = 0.;
    acom01[i] = 0.;
    acom11[i] = 0.;
  }   


   if(verbose!=0)printf("Got some arrays to put the results in\n");
   percent=(blobarray->dimensions[s]-1)/80.0;
   if(percent < 1)percent=1;
   bad=0;
   if(verbose!=0)printf("Scanning image\n");
   for( i = 0 ; i <= (blobarray->dimensions[s]-1) ; i++ ){    /* i,j is looping along the indices data array */
      if(verbose!=0 && (i%percent == 0) )printf(".");
      for( j = 0 ; j <= (blobarray->dimensions[f]-1) ; j++ ){
         peak=* (int *) (blobarray->data + i*blobarray->strides[s] + j*blobarray->strides[f]);
         if( peak > 0  && peak <=np ) {
             fval=(double)(getval((dataarray->data + i*dataarray->strides[s] + j*dataarray->strides[f]),type));
             asum[peak]  =asum[peak]+fval;
	     asumsq[peak]=asumsq[peak]+fval*fval;
	     acom0[peak] =acom0[peak]+i*fval;
	     acom1[peak] =acom1[peak]+j*fval;
	     acom00[peak] =acom00[peak]+i*i*fval;
	     acom01[peak] =acom01[peak]+i*j*fval;
	     acom11[peak] =acom11[peak]+j*j*fval;
	     anpix[peak] =anpix[peak]++;
#ifdef PARANOID
             if(verbose!=0){
	       printf("peak=%d i=%d j=%d val=%f acom00[peak]=%f 01=%f 11=%f\n",peak,i,j,fval,acom00[peak],
		      acom01[peak],acom11[peak]);}
#endif
	 }
     else{
       if(peak!=0){
          bad++;
          if(verbose!=0 && bad<10){printf("Found %d in your blob image at i=%d, j=%d\n",peak,i,j);}
          }
       }
      } /* j */
   } /* i */
   if(verbose){
     printf("\nFound %d bad pixels in the blob image\n",bad);
   }
   return Py_BuildValue("OOOOOOOO", PyArray_Return(npix),
			PyArray_Return(sum),PyArray_Return(sumsq),
			PyArray_Return(com0),PyArray_Return(com1),PyArray_Return(com00),
			PyArray_Return(com01),PyArray_Return(com11) ); 
   fail:
    Py_XDECREF(npix);
    Py_XDECREF(sum);
    Py_XDECREF(sumsq);
    Py_XDECREF(com0);
    Py_XDECREF(com1);
    Py_XDECREF(com00);
    Py_XDECREF(com01);
    Py_XDECREF(com11);
    return NULL;
}




static PyMethodDef connectedpixelsMethods[] = {
   {"connectedpixels", (PyCFunction) connectedpixels, METH_VARARGS | METH_KEYWORDS,
   "Assign connected pixels in image (first arg) above threshold (float second arg) to return image of peak number.Optional arguement verbose!=0 prints some info as it runs"},
   {"blobproperties", (PyCFunction) blobproperties, METH_VARARGS | METH_KEYWORDS,
   "Computes various properties of a blob image (created by connecetedpixels)"},
   
   {NULL, NULL, 0, NULL} /* setinel */
};

void 
initconnectedpixels(void)
{
   PyObject *m, *d;

   m=Py_InitModule("connectedpixels", connectedpixelsMethods);
   import_array();
   d=PyModule_GetDict(m);
}
