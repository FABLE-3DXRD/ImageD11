


#include <math.h>   /* sqrt */
#include <stdio.h>  /* printf */
#include <stdlib.h>  /* abs(int) */
#include <string.h> /* memset */
#include <omp.h>
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
  /* Use double to reduce rounding and subtraction errors */  
  double t, s1, s2, y0; 
  float mini, maxi;
  mini = img[0];
  maxi = img[0];
  s1 = 0.;
  s2 = 0.;
  y0 = img[0];  
    /* Merge results - openmp 2.0 for windows has no min/max     */
#ifdef _MSC_VER
#pragma omp parallel for private(t) reduction(+:s1,s2)
    for( i = 0; i < npx; i++ ){
      t = img[i] - y0;
      s1 = s1 + t;
      s2 = s2 + t * t;
    }
    /* Just use serial - openmp 2.0 */
    for( i = 0; i < npx; i++ ){
      if( img[i] < mini) mini = img[i];
      if( img[i] > maxi) maxi = img[i];
    }
#else
#pragma omp parallel for private(t) reduction(+:s1,s2) reduction(min: mini) reduction(max: maxi)
    for( i = 0; i < npx; i++ ){
      t = img[i] - y0;
      s1 = s1 + t;
      s2 = s2 + t * t;
      if( img[i] < mini ) mini = img[i];
      if( img[i] > maxi ) maxi = img[i];
    }
#endif
    /* results */
  *mean = (float) (s1 / npx + y0);
  *var  = (float) ((s2 - (s1*s1/npx))/npx);
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
              int nnz ){
  int k, es, ed;
  es = nnz+1;
  ed = nnz+1;
#ifndef _MSC_VER
#pragma omp parallel for private(k) reduction(min: es, ed)
#endif
  for( k=1 ; k<nnz ; k++){
    if( i[k] < i[k-1] ){ /* bad, not sorted */
      es = (k < es) ? k : es;
      continue;
    }
    if (i[k] == i[k-1]) { /* Same row, j must be gt prev */
      if( j[k] < j[k-1] ){ /* bad */
    es = (k < es) ? k : es;
      } else if (j[k] == j[k-1] ){
    ed = (k < ed) ? k : ed;
      } else {
    continue;
      }
    }
  }
  if ((es == (nnz+1)) && ( ed == (nnz+1) ) ) return 0;
  if (es > ed) return -ed;
  else return es;
}



/** 
 * Connected pixels
 * Using sparse implementation
 */

#define NOISY 0
int sparse_connectedpixels( float * restrict v,
                            uint16_t * restrict i,
                            uint16_t * restrict j,
                            int nnz,
                            float threshold,
                            int32_t * restrict labels  /* nnz */
                            ){
  int k, p, pp, ir;
  int32_t *S, *T, np;
  /* Read k = kurrent
          p = prev */
  double start, mid, end;
  if(NOISY){
    start = my_get_time();
    k = sparse_is_sorted( i, j, nnz);
    if( k!= 0 ) return k;
  }
  pp=0;
  p=0;
  S = dset_initialise( 16384 );
  /* Main loop */
  if(NOISY)printf("ok to main loop\n");
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
    if( ((j[p]+1) == j[k])   &&
    ( i[p]    == i[k])   &&
    ( labels[p] > 0 )) {
      labels[k] = labels[p];
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
    for( p = pp; j[p] <= j[k]+1; p++ ){
      if( i[p] == ir ){
    if ( labels[p] > 0) {
      // Union p, k
      match( labels[k], labels[p], S);
    }
      } else {
    break; // not same row
      }
    }
    if( labels[k] == 0) dset_new( &S, &labels[k]);
  } // end loop over data
  if(NOISY) mid = my_get_time();
  T = dset_compress( &S, &np );
  // renumber labels
#pragma omp parallel for private(k) shared(labels)
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
  if(NOISY){
    end=my_get_time();
    printf("Time in sparse image %f ms %f ms\n", 1000*(end-mid),1000*(mid-start));
  }
  return np;
}



/** 
 * Connected pixels
 * Using internal dense implementation
 */

#define NOISY 0
int sparse_connectedpixels_splat( float * restrict v,
                  uint16_t * restrict i,
                  uint16_t * restrict j,
                  int nnz,
                  float threshold,
                  int32_t * restrict labels,  /* nnz */
                  int32_t * restrict Z, /* workspace, at least (imax+2)*(jmax+2) */
                  int imax,
                  int jmax
                  ){
  int k, p, pp, ir, jdim, ik, jk;
  int32_t *S, *T, np;
  /* Read k = kurrent
          p = prev */
  double start, mid;
  if(NOISY){
    start = my_get_time();
    k = sparse_is_sorted( i, j, nnz );
    if( k!= 0) return k;
  }
  if(NOISY){
    mid = my_get_time();
    printf("check sorted %.3f ms\n",(mid-start)*1000);
    start = my_get_time();
  }

  jdim = jmax + 2;
  /* This is not! delivered with zeros, we put a border in too 
   *  Z = (int32_t *) malloc(idim*jdim* sizeof(int32_t));
   * later we will write into Z as a scratch area for labels (filled at very end) */
  pp=0;
  p=0;
  S = dset_initialise( 16384 );
  if(NOISY){
    mid = my_get_time();
    printf("mallocs %.3f ms\n",(mid-start)*1000);
    start = my_get_time();
  }
  /* zero the parts of Z that we will read from (pixel neighbors) */
#pragma omp parallel for private(p, k, ik, jk) shared(Z)
  for( k=0; k<nnz; k++){
    ik = i[k]+1; /* the plus 1 is because we padded Z */
    jk = j[k]+1;
    p = ik*jdim+jk;
    Z[p] = 0;
    Z[p-1] = 0;
    Z[p-jdim-1] = 0;
    Z[p-jdim] = 0;
    Z[p-jdim+1] = 0;
  }
  if(NOISY){
    mid = my_get_time();
    printf("zeros %.3f ms\n",(mid-start)*1000);
    start = my_get_time();
  }
  
  /* Main loop */
  for( k=0; k<nnz; k++){
    if( v[k] <= threshold) {
      continue;
    }
    /* Decide on label for this one ... 
     *
     * 4 neighbors : k-1 is prev
     */
    ik = i[k]+1; /* the plus 1 is because we padded Z */
    jk = j[k]+1;
    p = ik*jdim+jk;
    /* previous pixel, same row */
    if( Z[p-1] > 0 ){
      Z[p] =  Z[p-1];
    }
    /* 3 pixels on previous row */
    ir = (ik-1)*jdim + jk;
    for( pp = ir-1; pp <= ir + 1; pp++ ){
      if( Z[pp] > 0 ){
        // Union p, k
        match( Z[p], Z[pp], S);
      }
    } 
    if( Z[p] == 0) dset_new( &S, &Z[p]);
  } // end loop over data
  if(NOISY){
    mid = my_get_time();
    printf("main loop %.3f ms\n",(mid-start)*1000);
    start = my_get_time();
  }
  T = dset_compress( &S, &np );
  // renumber labels
#pragma omp parallel for private(ik,jk,p) shared( labels )
  for( k=0; k<nnz; k++){
    ik = i[k]+1; /* the plus 1 is because we padded Z */
    jk = j[k]+1;
    p = ik*jdim+jk;
    if( Z[p] > 0 ){
      labels[k] = T[Z[p]];
    }
  }
  if(NOISY){
    mid=my_get_time();
    printf("Relabelling %f ms\n", 1000*(mid-start));
    start= my_get_time();
  }
  free(S);
  free(T);
  if(NOISY){
    mid=my_get_time();
    printf("Free %f ms\n", 1000*(mid-start));
  }
  return np;
}

