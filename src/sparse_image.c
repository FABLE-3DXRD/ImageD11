


#include <math.h>   /* sqrt */
#include <stdio.h>  /* printf */
#include <stdlib.h>  /* abs(int) */
#include <string.h> /* memset */

#include "blobs.h"




/**
 * Go through the data to compute sum and sum2 for mean and std 
 * Also find min and max. Uses openmp for reductions
 * 
 * @param img Input data array  (1D or 2D.ravel(), so sparse or dense)
 * @param npx length of data array (contiguous)
 * @param *minval minimum of the pixels
 * @param *maxval maximum of the pixels
 * @param *s1  Sum of all pixel
 * @param *s2  Sum of pixel^2
 */
void array_stats( float img[], int npx,
		  float *minval, float *maxval,
		  float *mean, float *var){
  int i;
  float mini, maxi;
  /* Use double to reduce rounding and subtraction errors */
  double t, s1, s2, y0; 
  
  s1 = 0.;
  s2 = 0.;
  mini = img[0];
  maxi = img[0];
  y0 = img[0];  
  
#pragma omp parallel for private(t) reduction(min:mini), reduction(max:maxi) reduction(+:s1,s2)
  for( i = 0; i < npx; i++ ){
    t = img[i] - y0;
    s1 = s1 + t;
    s2 = s2 + t * t;
    if( img[i] < mini ){
      mini = img[i];
    }
    if( img[i] > maxi ){
      maxi = img[i];
    }
  }
  /* results */
  *mean = s1 / npx + y0;
  *var  = (s2 - (s1*s1/npx))/npx;
  *minval = mini;
  *maxval = maxi;
}


/**
 * Go through the data to compute a histogram of the values
 * Previous call of array_stats would help to set up this call 
 *  compare to np.bincount - this does not gain much
 *  better implementations can be done 
 *
 * @param img[]  Input data array  (1D or 2D.ravel(), so sparse or dense)
 * @param npx    length of data array (contiguous)
 * @param low    Lower edge of first bin
 * @param high   Upper edge of last bin
 * @param hist[] Histogram to be output
 * @param nhist  Number of bins in the histogram
 */
void array_histogram( float img[],
		      int npx,
		      float low,
		      float high,
		      int32_t hist[],
		      int nhist){
  int i, ibin;
  float ostep;
  memset( hist, 0, nhist * sizeof( int32_t ) );
  /* Compute the multiplier to get the bin numbers */
  ostep = nhist/(high - low);
  for( i=0; i<npx; i++){
    ibin = (int) floor( (img[i] - low)*ostep );
    /* clip into range at ends */
    if(ibin < 0){ 
      ibin = 0;
    } 
    if(ibin >= nhist ){
      ibin=nhist-1;
    }
    hist[ibin] = hist[ibin] + 1; 
  }
}

/* 
 * https://fr.mathworks.com/help/matlab/ref/sparse.html
 * S = sparse(i,j,v,m,n,nz)
 *   i = array of index i
 *   j = array of index j
 *   v = array of values
 *   S( i(k), j(k) ) == v(k)
 *   m = dimension 0 of S
 *   n = dimension 1 of S
 *   nz = number of non-zero values
 *
 * https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
 *  scipy.sparse.coo_matrix(arg1, shape=None, dtype=None, copy=False)
 *   arg1 = ( v, i, j ) as above
 *   shape = (m, n)     as above
 *   dtype = float/double/int, etc
 */

/**
 * Takes a pair of i,j index arrays and checks if they are sorted
 *
 *  @param i, j index arrays
 *  @param nnz dimension of i, j
 *  returns 0 for all OK
 *          k for first non-sorted element
 *         -k for first duplicate
 */
int sparse_is_sorted( uint16_t i[],
		      uint16_t j[],
		      int nnz,
		      int verbose ){
  int k, i0, j0;
  i0 = i[0];
  j0 = j[0];
  if(verbose){printf("\n\nin sps : %d %d\n",nnz,verbose);}
  for( k=1 ; k<nnz ; k++){
    if(verbose) printf("k %d i0 %d j0 %d \n",k,i0,j0);
    if( i[k] < i0 ) {
      if(verbose) printf("k %d i[k] %d returning\n",k,i[k]);
      return k;  /* Not sorted */
    }
    if( i[k] > i0 ) { /* New row */
      i0 = i[k]; 
      j0 = j[k];
    } else {
      if( j[k] < j0 ) {
	if(verbose) printf("k %d j[k] %d returning\n",k,j[k]);
	return k;  /* Not sorted */
      }
      if( j[k] == j0 ){
	if(verbose) printf("k %d j[k] %d returning\n",k,j[k]);
	return -k; /* Duplicate */
      }
      j0 = j[k];
    }
  }
  return 0;
}


/** 
 * Connected pixels
 *
 */

int sparse_connectedpixels( float v[],
			    uint16_t i[],
			    uint16_t j[],
			    int nnz,
			    float threshold,
			    int32_t labels[]  /* nnz */
			    ){
  int k, p, pp, ir;
  int32_t *S, *T, np;
  /* Read k = kurrent
          p = prev */
  if(0){
    k = sparse_is_sorted( i, j, nnz, 0 );
    if( k!= 0){
      /* BAD */
      k = sparse_is_sorted( i, j, nnz, 1 );
      return k;
    } else {
      //printf( "OK, k=%d nnz=%d\n",k, nnz);
    }
  }
  pp=0;
  p=0;
  S = dset_initialise( 16384 );
  /* Main loop */
  for( k=0; k<nnz; k++){
    labels[k] = 0;
    if( v[k] <= threshold) {
      continue;
    }
    /* Decide on label for this one ... 
     *
     * 4 neighbors : k-1 is prev
     */
    p = k-1; /* previous pixel, same row */
    if( ( i[p]    == i[k])   &&
	((j[p]+1) == j[k])) {
      if ( labels[p] > 0 ) {
	match( &labels[k], &labels[p], S);
      }
    }
    ir = i[k]-1;
    /* pp should be on row above, on or after j-1 */
    while( ir > i[pp]  ) pp++;
    /* out if nothing on row above */
    if( i[pp] == i[k] ) {
      if( labels[k] == 0){
	S = dset_new( &S, &labels[k] );
	continue;
      }
    }
    /* Locate previous pixel on row above */ 
    while( ((j[k]-j[pp]) > 1)&&(i[pp]==ir) ) pp++;
    for( p = pp; p < (pp+3); p++ ){
      if( (j[p] - j[k]) > 1) break;
      if( (i[p] == ir)   &&   // same row
	  ( labels[p] > 0) ) {
	// Union p, k
	match( &labels[k], &labels[p], S);
      }
    }
    if( labels[k] == 0) dset_new( &S, &labels[k]);
  } // end loop over data
  T = dset_compress( &S, &np );
  // renumber labels
  for( k=0; k<nnz; k++){
    if( labels[k] > 0 ){
      /* if( T[labels[k]] == 0 ){
	printf("Error in sparse_connectedpixels\n");
	} */
      labels[k] = T[labels[k]];
    }
  }
  free(S);
  free(T);
  return np;
}


