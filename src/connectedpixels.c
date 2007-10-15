



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


static char moduledocs[] =\
"   C extensions for image analysis, part of ImageD11\n"\
"   ";



#include <Python.h>                  /* To talk to python */
/* #include "Numeric/arrayobject.h"      Access to Numeric */
#include "numpy/arrayobject.h"    /* Upgrade numpy */

#include "blobs.h"   /* Disjoint sets thing for blob finding */

#include <time.h>  /* Benchmarking */
#include <unistd.h> /* ??? */
#include <stdlib.h> /* ??? */


#undef PARANOID

double getval(char *p,int type);
char* getptr(PyArrayObject *ar, int f, int s, int i, int j);



char* getptr(PyArrayObject *ar, int f, int s, int i, int j){
  return ar->data + i*ar->strides[s] + j*ar->strides[f] ;
}

double getval(char *p,int type){
   /* return the value of a python array element converted to the 
    * universal c type of double */

  switch (type){
     case    NPY_BYTE   : return *(signed char   *)p*1.;
     case    NPY_UBYTE  : return *(unsigned char   *)p*1.;
     case    NPY_SHORT  : return *(short         *)p*1.;
     case    NPY_USHORT  : return *(unsigned short         *)p*1.;
     case    NPY_INT    : return *(int           *)p*1.;
     case    NPY_LONG   : return *(long          *)p*1.;
     case    NPY_FLOAT  : return *(float         *)p*1.;
     case    NPY_DOUBLE : return *(double        *)p*1.;
     case    NPY_UINT   : return *(unsigned int  *)p*1.;
     }
   printf("Oh bugger in getval - unrecongnised numeric type\n");
   printf("type = %d ",type);
   exit(1);
   return 0;
}

void ptype(int type){
   /* Print out the type of a Numeric array from the C point of view */
  printf("Your input type was ");
  switch (type){
  case    PyArray_CHAR   : 
    printf("PyArray_CHAR *(char *)\n");
    break;
  case    PyArray_SHORT  : 
    printf("PyArray_SHORT *(short *)\n");
    break;
  case    PyArray_INT    : 
    printf("PyArray_INT *(int  *)\n");
    break;
  case    PyArray_LONG   : 
    printf("PyArray_LONG *(long *)\n");
    break;
  case    PyArray_FLOAT  : 
    printf("PyArray_FLOAT *(float *)\n");
    break;
  case    PyArray_DOUBLE : 
    printf("PyArray_DOUBLE *(double *)\n");
    break;
#ifdef PyArray_UNSIGNED_TYPES
  case    PyArray_UBYTE  : 
    printf("PyArray_UBYTE *(unsigned char *)\n");
    break;
  case    PyArray_USHORT : 
    printf("PyArray_USHORT *(unsigned short*)\n");
    break;
  case    PyArray_UINT   : 
    printf("PyArray_UINT *(unsigned int *)\n");
    break;
#endif
     }
}


/* ==== Begin roisum ========================================= */

static char roisum_doc[] =\
"  (float) roisum ( Numeric.array(2D), xl , xh , yl , yh , verbose=0 ) \n"\
"   sum( array[xl:xh , yl:yh] ) where x,y refer to SLOW/FAST indices   \n"\
"   ... NOT Numeric indices, but slow/fast \n"\
"   Unsure why this was written - can be done with Numeric anyway     ";
   

/* Fill in an image of peak assignments for pixels */

static PyObject * roisum (PyObject *self, PyObject *args,  PyObject *keywds)
{

   PyArrayObject *dat=NULL; /* in (not modified) */

   int i,j,f,s,np;
   int xl,xh,yl,yh;

   int verbose=0,type; /* whether to print stuff and 
			  the type of the input array */
   
   static char *kwlist[] = {"data","xl","xh","yl","yh","verbose", NULL};

   double sum;

   if(!PyArg_ParseTupleAndKeywords(args,keywds, "O!iiii|i",kwlist,      
				   &PyArray_Type, &dat,   /* array arg */
				   &xl,&xh,&yl,&yh,
				   &verbose))        /* optional verbosity */
     return NULL;

   if(verbose!=0){
     printf("\n\nHello there from roisum, you wanted the verbose output...\n");
   }

   /* Check array is two dimensional and Ushort */
   if(dat->nd != 2 ){    
     PyErr_SetString(PyExc_ValueError,
		     "Data array must be 2d,first arg problem");
     return NULL;
   }
   type=dat->descr->type_num;
   if(verbose!=0)ptype(type);

   /* Decide on fast/slow loop - inner versus outer */
   if(dat->strides[0] > dat->strides[1]) {
     f=1;  s=0;}
   else {
     f=0;  s=1;
   }
   if (verbose!=0){
     printf("Fast index is %d, slow index is %d, ",f,s);
     printf("strides[0]=%d, strides[1]=%d\n",
	    dat->strides[0],
	    dat->strides[1]);
   }

   /* Check ROI is sane */
   if(xl > xh || 
      yl > yh || 
      yl < 0  || 
      xl < 0  || 
      xh >dat->dimensions[s] || 
      yh >dat->dimensions[f]){
     PyErr_SetString(PyExc_ValueError,
		     "Problem with your ROI");
     return NULL;
   }

   if(verbose!=0)printf("Summing ROI\n");

   sum=0.;

   np=0;

   for( i = xl ; i < xh ; i++ ){    /* i,j is looping along the 
				       indices data array */
     for( j = yl ; j < yh ; j++ ){

       sum+=getval( getptr(dat,f,s,i,j),type); 
       np++;

       if(verbose>2){
	 printf("%d %d %f\n",
		i,j,
		getval(getptr(dat,f,s,i,j),
		       type)); 
       }
     }
   }

   if(verbose!=0){
	   printf("Sum %f np %d\n",sum,np);
   }

   return Py_BuildValue("d", sum/np); 
}

/* ==== End of roisum ======================================== */


/* ==== connectedpixels ======================================== */

static char connectedpixels_doc[] =\
"   nblobs = connectedpixels ( data=Numeric.array(data, 2D)  , \n"\
"                              results=Numeric.array(blob, 2D, Int)  ,\n"\
"                              threshold=float threshold ,\n"\
"                              verbose=Int verbose )\n"\
"   data is normally an image \n"\
"   blob is an array to receive pixel -> blob assignments\n"\
"   threshold is the value above which a pixel is considered to be in a blob\n"\
"   verbose flags printing on stdout";

                              
/* Fill in an image of peak assignments for pixels */

static PyObject * connectedpixels (PyObject *self, 
				   PyObject *args,  
				   PyObject *keywds)
{
  PyArrayObject *dataarray=NULL, *results=NULL; /* in (not modified) 
						  and out (modified) */

  int i, j,  /* looping over pixels */
    k,       /* used in looping neighouring pixels */
    *curr,   /* current blob assignment in loop */
    l,
    f, s, /* fast and slow directions in image */
    np,   /* number of pixels */
    ival; 
  int *T; /* for the disjoint set copy */

  int verbose=0,type;   /* whether to print stuff and the 
			   type of the input array */
  int percent,npover;    /* for the progress indicator in verbose mode */
  static char *kwlist[] = {
    "data",
    "results",
    "threshold",
    "verbose", 
    NULL};
  
  clock_t tv1,tv1point5,tv2,tv3,tv4;  /* for timing */
  
  double time;
  double val; /* data, flood and dark values */

  float threshold;
  if(!PyArg_ParseTupleAndKeywords(args,keywds, "O!O!f|i",kwlist,      
				  &PyArray_Type, &dataarray,   /* array args */
				  &PyArray_Type, &results,   /* array args */
				  &threshold,
				  &verbose))        /* threshold and 
						       optional verbosity */
    return NULL;
  if(verbose!=0){
    printf("\n\nHello there from connected pixels,"\
	   "you wanted the verbose output...\n");
  }

  tv1=clock();  
  /* Check array is two dimensional */
  if(dataarray->nd != 2){
    PyErr_SetString(PyExc_ValueError, "Array must be 2D!");
    return NULL;
  }
  /* Save type for getval */
  type = dataarray->descr->type_num;

  if(verbose!=0)ptype(type);
  if(verbose!=0)printf("Thresholding at level %f\n",threshold);
  
  /* Check results argument */
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
  /* Default number before reallocation, rather large 
   * FIXME - can we pass this in and out?
   */
  S = dset_initialise(16384); 
  if(verbose != 0)printf("Initialised the disjoint set\n");

  
  npover = 0; /* number of pixels over threshold */

  /* Decide on fast/slow loop - inner versus outer */
  if(dataarray->strides[0] > dataarray->strides[1]) {
    f=1;  s=0;}
  else {
    f=0;  s=1;
  }
  if (verbose!=0){
    printf("Fast index is %d, slow index is %d, ",f,s);
    printf("strides[0]=%d, strides[1]=%d\n",
	   dataarray->strides[0],
	   dataarray->strides[1]);
  }

  percent=(results->dimensions[s]-1)/80.0;
  
  if(percent < 1)percent=1;

  tv1point5=clock();   
  if(verbose!=0){
    time=(tv1point5-tv1)/CLOCKS_PER_SEC;
    printf("Ready to scan image, setup took %g and "\
	   "clock is only good to ~0.015 seconds\n",time);
  }

  if(verbose!=0)printf("Scanning image\n");

  /* === Mainloop === */
  for( i = 0 ; i < (results->dimensions[s]) ; i++ ){    /* i = slow */

    if(verbose!=0 && (i%percent == 0) )printf(".");
    for( j = 0 ; j < (results->dimensions[f]) ; j++ ){ /* j = fast*/

      /* Set result for this pixel to zero - Not needed here? */
      curr = (int *) getptr(results,f,s,i,j);
      *curr = 0;
      
      val = getval( getptr(dataarray,f,s,i,j), type);
	
      /* printf("curr %d val=%f",*curr,val); */
      if( val > threshold) {
	npover++;
	k = 0;
	l = 0;
	/* peak needs to be assigned */
	/* i-1, j-1, assign same index */
	if(i!=0 && j!=0)
	  if( (k = *(int *)(getptr(results,f,s,i-1,j-1))) > 0) {
	    *curr = k ;
	    l = k;
	  }
	/* i-1,j */
	if(i != 0)
	  if( (k = *(int *)(getptr(results,f,s,i-1,j))) > 0) {
	    if(l==0){
	      *curr = k;
	      l = k; /* l holds the assignment of the current pixel */
	    } else {
	      if(l != k)dset_makeunion(S,k,l);
	    }
	  }
	/* i-1, j+1 */
	if(i!=0 && j!=(results->dimensions[f]-1))  
	  if( (k = *(int *)(getptr(results,f,s,i-1,j+1))) > 0) {
	    if(l==0){
	      *curr = k;
	      l = k; /* l holds the assignment of the current pixel */
	    } else {
	      if(l != k)dset_makeunion(S,k,l);
	    }
	  }
	/* i, j-1 */
	if(j!=0)
	  if(  (k = *(int *)(getptr(results,f,s,i,j-1))) > 0) {
	    if(l==0){
	      *curr = k;
	      l = k; /* l holds the assignment of the current pixel */
	    } else {
	      if(l != k)dset_makeunion(S,k,l);
	    }
	  }
	if(l == 0){ /* pixel has no neighbours thus far */
	  S = dset_new(&S);
	  *curr = S[S[0]-1];
	} 
      } /* active pixel */
      else { /* inactive pixel  - zero in results */
	*curr = 0;		 
      }
    } /* j */
  } /* i */
  /* ===== End of loop over data assigning in disjoint set == */
  tv2=clock();   
  if(verbose!=0){
    printf("\nFinished scanning image, now make unique list and sums\n");
    printf("Number of pixels over threshold was : %d\n",npover);
    time=(tv2-tv1)/CLOCKS_PER_SEC;
    printf("That took %f seconds and %f in total being aware "\
	   "that the clock is only good to ~0.015 seconds\n",
	   1.*(tv2-tv1point5)/CLOCKS_PER_SEC,time*1.);
  }
  /* First, count independent peaks */

  /* count peaks */
  np = S[S[0]-1];

  /* Now make each T[i] contain the unique ascending integer for the set */ 
  T=NULL;
  if(verbose!=0)printf("Now I want to malloc T=%p at "\
		       "size np=%d and sizeof(int)=%d\n",
		       (void *) T ,np ,sizeof(int));

   if( ( T=(int *) ( malloc( (np+3)*sizeof(int) ) ) ) == NULL){
      printf("Memory allocation error in connected pixels\n");
      return 0;
   }
   if(verbose!=0)printf("Now T should be allocated, it is at %p\n",(void*)T);
   np = 0;
   for( i=1 ; i<S[S[0]-1]+1 ; i++){
     if(S[i] == i ){
       np++;
       T[i]=np;
     }
     else{  /* check */
       j=dset_find(i,S);
       T[i]=T[j];
       if(j>=i && verbose){
	 printf("Oh dear - there was a problem compressing"\
		" the disjoint set, j=%d, i=%d \n",j,i);
       }
       if(S[j]!=j && verbose!=0){
	 printf("Oh dear - the disjoint set is squiff,"\
		"S[j]=%d,j=%d\n",S[j],j);
       }
      }
   }
   if(verbose!=0)printf("\n");
   tv3=clock();
   if(verbose!=0){
     printf("Compressing list, time=%g, total time=%g\n",
	    1.*(tv3-tv2)/CLOCKS_PER_SEC,1.*(tv3-tv1)/CLOCKS_PER_SEC);
     }

   if(verbose)printf("Found %d peaks\nNow assigning unique identifiers\n",np);


   /* Now loop through the results image again,
      assigning the corrected numbers to each pixel, 
      and collecting the properties ??? */
   for( i = 0 ; i < (results->dimensions[s]) ; i++ ){    
     for( j = 0 ; j < (results->dimensions[f]) ; j++ ){
       curr = (int*)getptr(results,f,s,i,j);
       ival = *curr;
       if(  ival > 0){
	 if(T[ival]==0){
	   printf("Something buggered up, ival=%d, T[ival]=%d, "\
		  "dset_find(ival,S)=%d\n",
		  ival,T[ival],dset_find(ival,S));
	   free(S); 
	   exit(2);
	   return NULL;
	 }
	 *curr = T[ival];
       }
     }
   }

   if(verbose!=0)printf("Freeing S=%p and T=%p\n",(void*)S,(void*)T);
   free(S);
   free(T);
   if(verbose!=0){
     printf("Freed S=%p and T=%p\n",(void*)S,(void*)T);
     printf("Finished assigning unique values and tidied up, now returning");
     tv4=clock();
     time=(tv4-tv3)/CLOCKS_PER_SEC;
     printf("\nThat took %g seconds, being aware that the clock is "\
	    "only good to ~0.015 seconds\n",1.*time);
     printf("The total time in this routine was about %f\n",
	    1.*(tv4-tv1)/CLOCKS_PER_SEC); 
   }
   
   return Py_BuildValue("i", np);/* why the plus one?? */ 
}

/* ============================================================== */
/* BLOB PROPERTIES BIT - things to measure are: */
/* Centre of mass (x,y) - optionally use x,y image */
/* Amount of intensity in blob */
/* std dev of intensity in blob */
/* ...? */

static char blobproperties_doc[] =\
  "res = blobproperties ( Numeric.array(data, 2D)  , \n"	\
  "                          Numeric.array(blob, 2D, Int)  ,\n"	\
  "                          Int np , \n" \
  "                          Int verbose )\n"				\
  "\n"									\
  "   Computes various properties of a blob image "\
  "           (created by connectedpixels)\n" \
  "   data  = image data \n"						\
  "   blob  = integer peak assignments from connectedpixels \n"		\
  "   np    = number of peaks to treat \n"				\
  "   verbose  - flag about whether to print\n"				\
  "  \n";


   
static PyObject * blobproperties (PyObject *self, 
				  PyObject *args,  
				  PyObject *keywds)
{
   PyArrayObject *dataarray=NULL,*blobarray=NULL; /* in (not modified) */

   PyArrayObject *results=NULL;
   static char *kwlist[] = {
     "data",
     "blob",
     "npeaks",
     "omega",
     "verbose", 
     NULL};
   double *res;
   double fval;
   double omega=0;
   
   int np=0,verbose=0,type,f,s,peak,bad;

   int i,j, percent;
   int safelyneed[3];
   
   if(!PyArg_ParseTupleAndKeywords(args,keywds, "O!O!i|di",kwlist,      
				   &PyArray_Type, 
				   &dataarray,   /* array args - data */
				   &PyArray_Type, 
				   &blobarray,   /* blobs */
				   &np,  /* Number of peaks to treat */	   
				   &omega, /* omega angle - put 0 if unknown */
				   &verbose))     /* optional verbosity */
      return NULL;

   if(verbose)printf("omega = %f",omega);

   /* Check array is two dimensional */
   if(dataarray->nd != 2){     
     PyErr_SetString(PyExc_ValueError,
		     "data array must be 2d, first arg problem");
     return NULL;
   } 
   if(verbose!=0)printf("Welcome to blobproperties\n");
   /* Check array is two dimensional and int - 
      results from connectedpixels above */
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
     printf("strides[0]=%d, strides[1]=%d\n",
	    dataarray->strides[0],
	    dataarray->strides[1]);
   }
   /* results arrays */
   safelyneed[0]=np;
   safelyneed[1]=NPROPERTY;

   results = (PyArrayObject *) PyArray_FromDims(2, 
						safelyneed, 
						PyArray_DOUBLE);

   if( results == NULL){
     PyErr_SetString(PyExc_ValueError,
		     "Malloc failed fo results");
     return NULL;
   }

   if(verbose>0)printf("malloc'ed the results\n");

   res = (double*)results->data;

   /* Initialise the results */
   for ( i=0 ; i<safelyneed[0] ; i++) {
     for ( j=0 ; j<NPROPERTY; j++){
       res[i*NPROPERTY+j]=0.;
       res[i*NPROPERTY+bb_mn_f]=1.e9;
       res[i*NPROPERTY+bb_mn_s]=1.e9;
     }
   }
      
   if(verbose!=0)printf("Got some space to put the results in\n");

   percent=(blobarray->dimensions[s]-1)/80.0;
   if(percent < 1)percent=1;
   bad = 0;
   if(verbose!=0)printf("Scanning image\n");

   /* i,j is looping along the indices data array */
   for( i = 0 ; i <= (blobarray->dimensions[s]-1) ; i++ ){  
     if(verbose!=0 && (i%percent == 0) )printf(".");
     for( j = 0 ; j <= (blobarray->dimensions[f]-1) ; j++ ){

       peak =  *(int *) getptr(blobarray,f,s,i,j);
       
       if( peak > 0  && peak <=np ) {
	 fval = getval(getptr(dataarray,f,s,i,j),type);
	 /* printf("i,j,f,s,fval,peak %d %d %d %d %f %d\n",i,j,f,s,fval,peak); */
	 add_pixel( &res[NPROPERTY*(peak-1)] ,i , j ,fval , omega);
       }
       else{
	 if(peak!=0){
	   bad++;
	   if(verbose!=0 && bad<10){
	     printf("Found %d in your blob image at i=%d, j=%d\n",peak,i,j);}
	   /* Only complain 10 times - otherwise piles of crap go to screen */
	 }
       }
     } /* j */
   } /* i */
   if(verbose){
     printf("\nFound %d bad pixels in the blob image\n",bad);
   }

   
   return Py_BuildValue("O", PyArray_Return(results) ); 
}

/* ===  res2human     ==========================================================*/ 

static char blob_moments_doc [] =			\
  "   None = blob_moments(Numeric.array(peaks1))\n"	\
  "   \n"						\
  "   Loop over array filling out moments from sums\n";

static PyObject * blob_moments( PyObject *self,
				PyObject *args,
				PyObject *kwds){
  PyArrayObject *r = NULL;
  int verbose = 0 ;
  static char *kwlist[] = {
    "res",
    "verbose",
    NULL
  };
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!|i", kwlist,
				  &PyArray_Type, &r,
				  &verbose))
     return NULL;

  if (r->dimensions[1] != NPROPERTY ||
      r->descr->type_num != PyArray_DOUBLE ){
    PyErr_SetString(PyExc_ValueError,
		    "Results array looks corrupt\n");
    return NULL;
  }
  if(verbose){
    printf("Welcome to blob_moments\n");
    printf("%p r->data\n\n",r->data);
    printf("dim[0] = %d dim[1] = %d\n",r->dimensions[0], r->dimensions[1]);
  }
  compute_moments( (double *) r->data, r->dimensions[0] );

  Py_INCREF(Py_None);
  return Py_None;

}


/* ==============================================================================*/ 




/* ==============================================================================*/ 

static char bloboverlaps_doc[] =\
"   success = bloboverlaps (   Numeric.array(blob1, 2D, Int), n1  , res1 \n"\
"                              Numeric.array(blob2, 2D, Int), n2  , res2,\n"\
"                               verbose=0) \n"\
" \n"\
" merges the results from blob1/res1 into blob2/res2\n"\
" blob1 would be the previous image from the series, with its results\n"\
"   if it does not overlap it will stay untouched in res\n"\
"   if it overlaps with next image it is passed into that ones res\n";


 
   
static PyObject * bloboverlaps (PyObject *self, 
				PyObject *args,  
				PyObject *keywds)
{
  PyArrayObject *b1=NULL,*b2=NULL; /* in */
  PyArrayObject *r1=NULL, *r2=NULL; /* in/out */
  int i,j,s,f, verbose ; /* loop vars, flag */
  int p1,p2,n1,n2; /* peak num and npeaks in each */
  int *link , percent, safelyneed;
  int ipk, jpk;
  double *res1, *res2;
  static char *kwlist[] = {
    "blob1",
    "n1",
    "res1",
    "blob2",
    "n2",
    "res2",
    "verbose", NULL};                      /*   b nr b nr */
  if(!PyArg_ParseTupleAndKeywords(args,keywds, "O!iO!O!iO!|i",kwlist,      
				  &PyArray_Type, &b1,   /* blobs 1 */
				  &n1,                  /* npeaks 1 */
				  &PyArray_Type, &r1,    /* results 1 */
				  &PyArray_Type, &b2,   /* blobs 2 */
				  &n2,                  /* npeaks 2 */
				  &PyArray_Type, &r2,   /*results 2 */
				  &verbose))        /* optional verbosity */
    return NULL;

  if(verbose!=0)printf("Welcome to bloboverlaps\n");
  
  /* Check array is two dimensional and int - 
     results from connectedpixels above */
  
  if(b1->nd != 2 && b1->descr->type_num != PyArray_INT){     
    PyErr_SetString(PyExc_ValueError,
		    "Blob1 array must be 2d and integer, first arg problem");
    return NULL;
  }
  if(b2->nd != 2 && b2->descr->type_num != PyArray_INT){     
    PyErr_SetString(PyExc_ValueError,
		    "Blob2 array must be 2d and integer, second arg problem");
    return NULL;
  }
  if (b1->strides[0]!=b2->strides[0] ||
      b1->strides[1]!=b2->strides[1] ||
      b1->dimensions[0]!=b2->dimensions[0] ||
      b1->dimensions[1]!=b2->dimensions[1] ){
    PyErr_SetString(PyExc_ValueError,
		    "Blob1 and Blob2 array be similar (dims & strides)");
    return NULL;
  }       
  /* check results args */
  if (r1->dimensions[0] != n1 ||
      r1->dimensions[1] != NPROPERTY ||
      r2->dimensions[0] != n2 ||
      r2->dimensions[1] != NPROPERTY ||
      r1->descr->type_num != PyArray_DOUBLE ||
      r2->descr->type_num != PyArray_DOUBLE ){
    PyErr_SetString(PyExc_ValueError,
		    "Results arrays look corrupt\n");
    return NULL;
  }
  res1 = (double *)r1->data;
  res2 = (double *)r2->data;

  /* Decide on fast/slow loop - inner versus outer */
  if(b1->strides[0] > b1->strides[1]) {
    f=1;  s=0;}
  else {
    f=0;  s=1;
  }
  if (verbose!=0){
    printf("Fast index is %d, slow index is %d, ",f,s);
    printf("strides[0]=%d, strides[1]=%d\n",b1->strides[0],b2->strides[1]);
  }
  
  /* Initialise a disjoint set in link 
   * image 1 has peak[i]=i ; i=1->n1
   * image 2 has peak[i]=i+n1 ; i=1->n2
   * link to hold 0->n1-1 ; n1->n2+n2-1 */

  safelyneed=n1+n2+3;
  link =  (int *)malloc(safelyneed*sizeof(int));
  
  for(i=0;i<safelyneed;i++){link[i]=i;} /* first image */
  /* flag the start of image number 2 */
  link[n2+1]=-99999; /* ==n2=0 Should never be touched by anyone */
  
  /* results lists of pairs of numbers */
  if(verbose!=0){
    printf("Scanning image, n1=%d n2=%d\n",n1,n2);
  }
  percent = b1->dimensions[s] / 100.;
  for( i = 0 ; i <= (b1->dimensions[s]-1) ; i++ ){    /* i,j is looping 
							 along the indices 
							 of the data array */
    for( j = 0 ; j <= (b1->dimensions[f]-1) ; j++ ){
      p1 = * (int *) getptr(b1,f,s,i,j);
      if ( p1 == 0 ){ continue; } 
      p2 = * (int *) getptr(b2,f,s,i,j);
      if ( p2 == 0 ){ continue; } 
      /* Link contains the peak that this peak is */
      if(link[p2] < 0 || link[p1+n2+1]<0){
	printf("Whoops p1=%d p2=%d p2+n1=%d link[p1]=%d link[p2+n1]=%d\n",
	       p1,p2,p2+n1,link[p1],link[p2+n1+1]);
	return NULL;
      }
      if(verbose>2)printf("p1 %d p2 %d\n",p1,p2);
      dset_makeunion(link, p2, p1+n2+1);
      /* printf("link[60]=%d ",link[60]); */
      if(verbose>2)printf("link[p2=%d]=%d link[p1+n2=%d]=%d\n",
			  p2,link[p2],p2+n1+1,link[p1+n2+1]);
    } /* j */
  } /* i */
  if(verbose)printf("Finished scanning blob images\n");
  for(i=1; i < safelyneed; i++){
    if(link[i] != i && i != n2+1){
      j = dset_find(i, link);
      if(verbose>1){
	printf("i= %d link[i]= %d j= %d n1= %d n2=%d \n",i,link[i],j, n1,n2);
      }
      if(i > n2+1 && j<n2+1){
	/* linking between images */
	jpk = j-1;
	ipk = i-n2-2;
	if(jpk<0 || jpk>r2->dimensions[0])printf("bounds jpk %d",jpk);
	if(ipk<0 || ipk>r1->dimensions[0])printf("bounds ipk %d",ipk);
	merge( &res2[NPROPERTY*jpk], &res1[NPROPERTY*ipk] );
	if(verbose>2)printf("merged ipk %d jpk %d\n",ipk,jpk);
	continue;
      }
      if(i > n2+1 && j>n2+1){
	/* linking on same image */
	jpk = j-n2-2;
	ipk = i-n2-2;
	if(jpk<0 || jpk>r1->dimensions[0])printf("bounds jpk %d",jpk);
	if(ipk<0 || ipk>r1->dimensions[0])printf("bounds ipk %d",ipk);
	merge( &res1[NPROPERTY*jpk], &res1[NPROPERTY*ipk] );
	if(verbose>2)printf("merged ipk %d jpk %d\n",ipk,jpk);
	continue;
      }
      if(i < n2+1 && j < n2+1){
	/* linking on same image */
	jpk = j-1;
	ipk = i-1;
	if(jpk<0 || jpk>r2->dimensions[0])printf("bounds jpk %d",jpk);
	if(ipk<0 || ipk>r2->dimensions[0])printf("bounds ipk %d",ipk);
	merge( &res2[NPROPERTY*jpk], &res2[NPROPERTY*ipk] );
	if(verbose>2)printf("merged check ipk %d jpk %d\n",ipk,jpk);
	continue;
      }
      printf("bad logic in bloboverlaps\n");
      return NULL;
      
    }
  }



  if(verbose!=0){
    printf("\n");
    for(i=0;i<n1+n2;i++){
      if(i!=link[i]){
	printf("Link found!! %d %d\n",i,link[i]);
      }
    }
  }
  free(link);
  if(verbose)printf("returning from bloboverlaps\n");
  Py_INCREF(Py_None);
  return Py_None;
}


static char update_blobs_doc[] =\
  "   update_blobs (   Numeric.array(blob, 2D, Int), \n"		\
  "                    Numeric.array(set , 2D, Int) , verbose=0) \n"	\
  " \n"									\
  " updates blob image such that : \n"					\
  " if blob[i,j] > 0: \n"						\
  "     blob[i,j] = set[blob[i,j]]\n"					\
  "\n"									\
  " Used to update a blob image merging peaks which have overlapped due\n"
  " to the third dimension.\n";


static PyObject * update_blobs (PyObject *self, PyObject *args,  PyObject *keywds)
{
  PyArrayObject *bl=NULL,*set=NULL; /* in  */
  int i,j,s,f, verbose = 0, p1, v ; /* loop vars, flag */
  static char *kwlist[] = {"blobim","set","verbose", NULL};
  if(!PyArg_ParseTupleAndKeywords(args,keywds, "O!O!|i",kwlist,      
				  &PyArray_Type, &bl,   /* blobs */
				  &PyArray_Type, &set,   /* disjoint set array */
				  &verbose))        /* threshold and optional verbosity */
    return NULL;
  if(verbose!=0)printf("Welcome to update blobs\n");
  /* Check and validate args */
  if(bl->nd != 2 && bl->descr->type_num != PyArray_INT){     
    PyErr_SetString(PyExc_ValueError,
		    "Blob array must be 2d and integer, first arg problem");
    return NULL;
  }
  if(set->nd != 1 && bl->descr->type_num != PyArray_INT){     
    PyErr_SetString(PyExc_ValueError,
		    "Set array must be 1d and integer, second arg problem");
    return NULL;
  }
  /* Decide on fast/slow loop - inner versus outer */
  if(bl->strides[0] > bl->strides[1]) {
    f=1;  s=0;}
  else {
    f=0;  s=1;
  }
  if (verbose!=0){
    printf("Fast index is %d, slow index is %d, ",f,s);
    printf("strides[0]=%d, strides[1]=%d\n",bl->strides[0],bl->strides[1]);
  }
  for( i = 0 ; i <= (bl->dimensions[s]-1) ; i++ ){    /* i,j is looping along the indices data array */
    for( j = 0 ; j <= (bl->dimensions[f]-1) ; j++ ){
      p1 = * (int *) (bl->data + i*bl->strides[s] + j*bl->strides[f]);
      if (p1==0){ continue; }
      if (p1 < 0){
	PyErr_SetString(PyExc_ValueError,
			"Blob image contains negative number! Not allowed");
	return NULL;
      }
      if (p1 < set->dimensions[0]){
	/* Write in the value */
	v = *(int *)(set->data + p1*set->strides[0]);
	if (v > 0 && v < set->dimensions[0]){
	  (* (int *)(bl->data + i*bl->strides[s] + j*bl->strides[f])) =  p1;
	} else {
	  PyErr_SetString(PyExc_ValueError,
			  "Set contains a bad value");
	  return NULL;
	}         		
      }else{
	PyErr_SetString(PyExc_ValueError,
			"Blob image references overflows set array element");
	return NULL;
      }
    }
  }
  /* http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52309 */
  /* return None - side effect was on arg */
  Py_INCREF(Py_None);
  return Py_None;
}







static PyMethodDef connectedpixelsMethods[] = {
  {"connectedpixels", (PyCFunction) connectedpixels, 
   METH_VARARGS | METH_KEYWORDS,
   connectedpixels_doc},
  {"blobproperties", (PyCFunction) blobproperties, 
   METH_VARARGS | METH_KEYWORDS,
   blobproperties_doc},
  {"bloboverlaps", (PyCFunction) bloboverlaps,  
   METH_VARARGS | METH_KEYWORDS,
   bloboverlaps_doc},     
  {"update_blobs", (PyCFunction) update_blobs,  
   METH_VARARGS | METH_KEYWORDS,
   update_blobs_doc},     
  {"blob_moments", (PyCFunction) blob_moments,  
   METH_VARARGS | METH_KEYWORDS,
   blob_moments_doc},     
  {"roisum", (PyCFunction) roisum, 
   METH_VARARGS | METH_KEYWORDS,
   roisum_doc},
  {NULL, NULL, 0, NULL} /* setinel */
};

void initconnectedpixels(void) {
   PyObject *m, *d, *s;

   m=Py_InitModule("connectedpixels", connectedpixelsMethods);
   import_array();
   d=PyModule_GetDict(m);
   s=PyString_FromString(moduledocs);
   PyDict_SetItemString(d,"__doc__",s);
#define str(x) (#x)
   PyDict_SetItemString(d,str(s_1)    ,PyInt_FromLong(s_1));
   PyDict_SetItemString(d,str(s_I)    ,PyInt_FromLong(s_I));
   PyDict_SetItemString(d,str(s_I2)   ,PyInt_FromLong(s_I2));
   PyDict_SetItemString(d,str(s_fI)   ,PyInt_FromLong(s_fI));
   PyDict_SetItemString(d,str(s_ffI)  ,PyInt_FromLong(s_ffI));
   PyDict_SetItemString(d,str(s_sI)   ,PyInt_FromLong(s_sI));
   PyDict_SetItemString(d,str(s_ssI)  ,PyInt_FromLong(s_ssI));
   PyDict_SetItemString(d,str(s_sfI)  ,PyInt_FromLong(s_sfI));
   PyDict_SetItemString(d,str(s_oI)   ,PyInt_FromLong(s_oI));
   PyDict_SetItemString(d,str(s_ooI)  ,PyInt_FromLong(s_ooI));
   PyDict_SetItemString(d,str(s_soI)  ,PyInt_FromLong(s_soI));
   PyDict_SetItemString(d,str(s_foI)  ,PyInt_FromLong(s_foI));
   PyDict_SetItemString(d,str(mx_I)   ,PyInt_FromLong(mx_I));
   PyDict_SetItemString(d,str(mx_I_f) ,PyInt_FromLong(mx_I_f));
   PyDict_SetItemString(d,str(mx_I_s) ,PyInt_FromLong(mx_I_s));
   PyDict_SetItemString(d,str(mx_I_o) ,PyInt_FromLong(mx_I_o));
   PyDict_SetItemString(d,str(bb_mx_f),PyInt_FromLong(bb_mx_f));
   PyDict_SetItemString(d,str(bb_mx_s),PyInt_FromLong(bb_mx_s));
   PyDict_SetItemString(d,str(bb_mx_o),PyInt_FromLong(bb_mx_o));
   PyDict_SetItemString(d,str(bb_mn_f),PyInt_FromLong(bb_mn_f));
   PyDict_SetItemString(d,str(bb_mn_s),PyInt_FromLong(bb_mn_s));
   PyDict_SetItemString(d,str(bb_mn_o),PyInt_FromLong(bb_mn_o));

   PyDict_SetItemString(d,str(avg_i),PyInt_FromLong(avg_i));
   PyDict_SetItemString(d,str(f_raw),PyInt_FromLong(f_raw));
   PyDict_SetItemString(d,str(s_raw),PyInt_FromLong(s_raw));
   PyDict_SetItemString(d,str(f_cen),PyInt_FromLong(f_cen));
   PyDict_SetItemString(d,str(s_cen),PyInt_FromLong(s_cen));
   PyDict_SetItemString(d,str(o_raw),PyInt_FromLong(o_raw));
   PyDict_SetItemString(d,str(m_ff),PyInt_FromLong(m_ff));
   PyDict_SetItemString(d,str(m_ss),PyInt_FromLong(m_ss));
   PyDict_SetItemString(d,str(m_oo),PyInt_FromLong(m_oo));
   PyDict_SetItemString(d,str(m_sf),PyInt_FromLong(m_sf));
   PyDict_SetItemString(d,str(m_so),PyInt_FromLong(m_so));
   PyDict_SetItemString(d,str(m_fo),PyInt_FromLong(m_fo));
   PyDict_SetItemString(d,str(dety),PyInt_FromLong(dety));
   PyDict_SetItemString(d,str(detz),PyInt_FromLong(detz));
   
   PyDict_SetItemString(d,str(NPROPERTY),PyInt_FromLong(NPROPERTY));
   Py_DECREF(s);
   if(PyErr_Occurred())
     Py_FatalError("cant initialise connectedpixels");
}

