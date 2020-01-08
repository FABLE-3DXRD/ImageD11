/* stdint */

#include  <stdlib.h>
#include  <stdio.h>
#include  <string.h>

#include <float.h>
#include <math.h>

#include "cImageD11.h"

/* To compare to numpy (numba ought to do this anyway)

__m256i _mm256_cvtepu16_epi32 (__m128i a)
__m256 _mm256_cvtepi32_ps (__m256i a)

SSE: __m128 _mm_cvtpu16_ps (__m64 a)


AVX
__m256 _mm256_sub_ps (__m256 a, __m256 b)
__m256 _mm256_mul_ps (__m256 a, __m256 b)
  Architecture	Latency	Throughput (CPI)
  Skylake	4	0.5
  [divide is much slower]


  _m256_load_ps (float const * mem_addr)
  _m256_store_ps (float * mem_addr, __m256 a)
  aligned load/store (32 byte boundary)

  _mm256_loadu_ps (float const * mem_addr)
  _mm256_storeu_ps (float const * mem_addr)
  unaligned load/store (32 byte boundary)
*/


void uint16_to_float_darksub( float *restrict img,
                        const float *restrict drk,
                        const uint16_t *restrict data,
                        int npx ){
    int i;
#ifdef GOT_OMP_SIMD
#pragma omp parallel for simd
#else
#pragma omp parallel for
#endif
    for(i=0;i<npx;i++){
        img[i] = ((float)data[i]) - drk[i];
    }
}


void frelon_lines( float *img, int ns, int nf, float cut ){
    int i, j, p, npx;
    float rowsum, avg;
    avg = img[0];
    #pragma omp parallel for private(i,j,rowsum,npx,p) firstprivate(avg)
    for( i = 0; i < ns ; i++){
        rowsum = 0.;
        npx = 0;
        p = i*nf;
#ifdef GOT_OMP_SIMD
        #pragma omp simd
#endif
        for( j = 0; j<nf ; j++){
            if ( img[p+j] < cut ){
                rowsum += img[p+j];
                npx++;
            }
        }
        if( npx > 0 )
            avg = rowsum / npx;
#ifdef GOT_OMP_SIMD
        #pragma omp simd
#endif
        for( j = 0; j<nf ; j++)
           img[p+j] = img[p+j] - avg;
    }
}

void frelon_lines_sub( float * restrict img, float * restrict drk, int ns, int nf, float cut ){
    int i, j, p, npx;
    float rowsum, avg;
    avg = img[0];
    #pragma omp parallel for private(i,j,rowsum,npx,p) firstprivate(avg)
    for( i = 0; i < ns ; i++){
        rowsum = 0.;
        npx = 0;
        p = i*nf;
#ifdef GOT_OMP_SIMD
        #pragma omp simd
#endif
        for( j = 0; j<nf ; j++){
            img[p+j] = img[p+j] - drk[p+j];
            if ( img[p+j] < cut ){
                rowsum += img[p+j];
                npx++;
            }
        }
        if( npx > 0 )
            avg = rowsum / npx;
#ifdef GOT_OMP_SIMD
        #pragma omp simd
#endif
        for( j = 0; j<nf ; j++)
           img[p+j] = img[p+j] - avg;
    }
}

void array_mean_var_cut( float * restrict img, int npx, float *mean, float *std,
                         int n, float cut, int verbose ){
    int i, nactive;
    float t, s1, s2, wt, y0;
    y0 = img[0];
    s1 = 0;
    s2 = 0;
    if (verbose) printf("Args, img[0] %f npx %d n %d cut %f verbose %d\n", \
                           img[0], npx,  n, cut, verbose);
#ifdef GOT_OMP_SIMD
#pragma omp parallel for simd private(t) reduction(+:s1,s2)
#else
#pragma omp parallel for private(t) reduction(+:s1,s2)
#endif
    for (i = 0; i < npx; i++) {
         t = img[i] - y0;
         s1 = s1 + t;
         s2 = s2 + t * t;
    }
    /* mean and std */
    *mean = (float)(s1 / npx + y0);
    *std = sqrtf( (float)((s2 - (s1 * s1 / npx)) / npx) );
    if(verbose>0) printf("n=%d Mean %f, Std %f\n",n,*mean,*std);
    while(--n > 0){
        y0 = *mean;
        wt = y0 + cut*(*std);
        s1 = 0;
        s2 = 0;
        nactive = 0;
#ifdef GOT_OMP_SIMD
#pragma omp parallel for simd private(t) reduction(+:s1,s2,nactive)
#else
#pragma omp parallel for private(t) reduction(+:s1,s2,nactive)
#endif
        for (i = 0; i < npx; i++) {
            if( img[i] < wt ){
               t = img[i] - y0;
               s1 = s1 + t;
               s2 = s2 + t * t;
               nactive++;
            }
        }
        *mean = (float)(s1 / nactive + *mean);
        *std = sqrtf( ((s2 - (s1 * s1 / nactive)) / nactive) );
        if(verbose>0) printf("n=%d Mean %f, Std %f\n",n,*mean,*std);
    }
}


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
void array_stats(float img[], int npx,
		 float *minval, float *maxval, float *mean, float *var)
{
    int i;
    /* Use double to reduce rounding and subtraction errors */
    double t, s1, s2, y0,ts1, ts2;
    float mini, maxi, tmin, tmax ;
    mini = FLT_MAX;
    maxi = FLT_MIN;
    s1 = 0.; 
    s2 = 0.; 
    y0 = img[0];
    /* Merge results - openmp 2.0 for windows has no min/max     */
#pragma omp parallel private(i, t, ts1, ts2, tmin, tmax )
{
    tmin = FLT_MAX;
    tmax = FLT_MIN;
    ts1 = 0.;
    ts2 = 0.;
#ifdef GOT_OMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
    for (i = 0; i < npx; i++) {
	t = img[i] - y0;
	ts1 = ts1 + t;
	ts2 = ts2 + t * t;
	if (img[i] < tmin)
	    tmin = img[i];
	if (img[i] > tmax)
	    tmax = img[i];
    } // for
#pragma omp critical
{
    s1 += ts1;
    s2 += ts2;
    if( tmin < mini) mini = tmin;
    if( tmax > maxi) maxi = tmax;
}
}// parallel
    /* results */
    *mean = (float)(s1 / npx + y0);
    *var = (float)((s2 - (s1 * s1 / npx)) / npx);
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
void array_histogram(float img[],
		     int npx, float low, float high, int32_t hist[], int nhist)
{
    int i, ibin;
    float ostep;
    memset(hist, 0, nhist * sizeof(int32_t));
    /* Compute the multiplier to get the bin numbers */
    ostep = nhist / (high - low);
    for (i = 0; i < npx; i++) {
	    ibin = (int)floorf((img[i] - low) * ostep);
	    /* clip into range at ends */
	    if (ibin < 0) {
  	    ibin = 0;
	    }
	    if (ibin >= nhist) {
  	    ibin = nhist - 1;
	    }
	    hist[ibin] = hist[ibin] + 1;
    }
}


void reorder_u16_a32( uint16_t * restrict data,
		      uint32_t * restrict adr,
		      uint16_t * restrict out,
		      int N
		     ){
  int i;
  /*  printf("Hello, got N=%d\n",N);*/
#pragma omp parallel for
  for(i=0 ; i<N ; i++){
    out[ adr[i] ] = data[i];
  }
}

void reorderlut_u16_a32( uint16_t * restrict data,
			 uint32_t * restrict lut,
			 uint16_t * restrict out,
			 int N
		     ){
  int i;
  /*  printf("Hello, got N=%d\n",N);*/
#pragma omp parallel for
  for(i=0 ; i<N ; i++){
    out[ i ] = data[lut[i]];
  }
}


void reorder_u16_a32_a16( uint16_t * restrict data,
			  uint32_t * restrict a0,
			  int16_t  * restrict a1,
			  uint16_t * restrict out,
			  int ns, int nf
		     ){
  int i, j, p;
  /*  printf("Hello, got ns=%d nf=%d\n",ns, nf);*/
#pragma omp parallel for private(p, j) 
  for(i=0 ; i<ns ; i++){
    p = a0[i];
    for(j=0 ; j<nf ; j++){
      p += a1[i*nf + j];
      out[ p ] = data[i*nf+j];
    }
  }
}


