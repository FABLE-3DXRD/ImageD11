
#include <stdlib.h>
#include <stdio.h>
#include "cImageD11.h"

#include <omp.h>

// Most compilers should have this!
#ifdef _MSC_VER
 #include <intrin.h>
 #define cast128i(X)  (_mm_castps_si128((X)))
 #ifdef __AVX2__
  #define cast256i(X)  (_mm256_castps_si256(X))
 #else
  #warning "no avx2"
 #endif
#else				// gcc/clang
 #ifdef __SSE2__
  #include <x86intrin.h>
  #include <emmintrin.h>
  #define cast128i(X)  ((__m128i)(X))
  #ifdef __AVX2__
   #include <immintrin.h>
   #define cast256i(X)  ((__m256i)(X))
  #endif
 #else
  #warning "no sse2"
 #endif
#endif

/**
 * neighbormax assigns values in "l" to address highest pixel
 * of the 8 neighbors.
 *  l == 0 => this pixel is the max
 *  otherwise o[ l-1 ] give the offset
 */

void neighbormax(const float *restrict im,	// input
		 int32_t * restrict lout,	// output
		 uint8_t * restrict l,	// workspace temporary
		 int dim0,	// Image dimensions
		 int dim1,
		 int o[10]){
    int i, j, p, iq, k;
    float mx;    
    for (i = 0; i < dim0; i++) {	// first and last row here:
	lout[i] = 0;		// ends of rows are below
	l[i] = 0;
	lout[dim1 * (dim0 - 1) + i] = 0;
	l[dim1 * (dim0 - 1) + i] = 0;
    }
#pragma omp parallel for private( j, p, iq, k, mx ) 
    for (i = dim1; i < (dim0 - 1) * dim1; i = i + dim1) {	// skipping 1 pixel border
	lout[i] = 0;		// set edges to zero: pixel j=0:
	l[i] = 0;
	for (j = 1; j < dim1 - 1; j++) {
	    p = i + j;
	    mx = im[p + o[1]];
	    iq = 1;
	    for (k = 2; k < 10; k++) {
		if (im[p + o[k]] > mx) {
		    iq = k;
		    mx = im[p + o[k]];
		}
	    }
	    l[p] = iq;
	}
	// set edges to zero: pixel j=dim1:
	lout[i + dim1 - 1] = 0;
	l[i + dim1 - 1] = 0;
    }				// i
}





void neighbormax_sse2(const float *restrict im,	// input
		      int32_t * restrict lout,	// output
		      uint8_t * restrict l,	// workspace temporary
		      int dim0,	// Image dimensions
		      int dim1,
		      int o[10]){
  __m128i one, iqp, ik, imsk;
  __m128 mxp, mxq, msk;
  int i, j, p, iq, k;
  float mx;
  // Set edges to zero (background label)
  for (i = 0; i < dim0; i++) {	// first and last row here:
    lout[i] = 0;		// ends of rows are below
    l[i] = 0;
    lout[dim1 * (dim0 - 1) + i] = 0;
    l[dim1 * (dim0 - 1) + i] = 0;
  }
#pragma omp parallel for private( j, p, iq, k, mx, mxp, iqp, ik, one, msk, mxq, imsk)
  for (i = dim1; i < (dim0 - 1) * dim1; i = i + dim1) {	// skipping 1 pixel border
    lout[i] = 0;		// set edges to zero: pixel j=0:
    l[i] = 0;
    one = _mm_set1_epi32(1); // SSE2
    for (j = 1; j < dim1 - 1 - 4; j = j + 4) {	// do 4 pixels at a time with sse2
      p = i + j;
      mxp = _mm_loadu_ps(&im[p + o[1]]);	// SSE load 4 floats
      iqp = one;		// current selection
      ik = one;		// to increment
      for (k = 2; k < 10; k++) {
	ik = _mm_add_epi32(ik, one);	// SSE2 ik++ 
	mxq = _mm_loadu_ps(&im[p + o[k]]);	// load vector of floats
	msk = _mm_cmpgt_ps(mxq, mxp);	// SSE 1 for q > p
	// apply mask to max (mxq/mxp) and to index (iqp/ik)
	mxp = _mm_or_ps(_mm_and_ps(msk, mxq), _mm_andnot_ps(msk, mxp));
  imsk = cast128i( msk );	// SSE 1 for q > p
	iqp = _mm_or_si128(_mm_and_si128(imsk, ik),  _mm_andnot_si128(imsk, iqp));
      }			// k neighbors
      // Write results (note epi16 is sse2)
      l[p] = (uint8_t) _mm_extract_epi16(iqp, 0);
      l[p + 1] = (uint8_t) _mm_extract_epi16(iqp, 2);
      l[p + 2] = (uint8_t) _mm_extract_epi16(iqp, 4);
      l[p + 3] = (uint8_t) _mm_extract_epi16(iqp, 6);
      // Count peaks in here? ... was better in separate loop
    }
    for (; j < dim1 - 1; j++) {	// end of simd loop, continues on j from for
      p = i + j;
      mx = im[p + o[1]];
      iq = 1;
      for (k = 2; k < 10; k++) {
	if (im[p + o[k]] > mx) {
	  iq = k;
	  mx = im[p + o[k]];
	}
      }
      l[p] = iq;
    }
    // set edges to zero: pixel j=dim1:
    lout[i + dim1 - 1] = 0;
    l[i + dim1 - 1] = 0;
  }				// i
}

#ifdef __AVX2__
void neighbormax_avx2(const float *restrict im,	// input
		      int32_t * restrict lout,	// output
		      uint8_t * restrict l,	// workspace temporary
		      int dim0,	// Image dimensions
		      int dim1,
		      int o[10]){
  __m256i one, iqp, ik, imsk;
  __m256 mxp, mxq, msk;
  int i, j, p, iq, k;
  float mx;
  // Set edges to zero (background label)
  for (i = 0; i < dim0; i++) {	// first and last row here:
    lout[i] = 0;		// ends of rows are below
    l[i] = 0;
    lout[dim1 * (dim0 - 1) + i] = 0;
    l[dim1 * (dim0 - 1) + i] = 0;
  }
#pragma omp parallel for private( j, p, iq, k, mx, mxp, iqp, ik, one, msk, mxq, imsk)
  for (i = dim1; i < (dim0 - 1) * dim1; i = i + dim1) {	// skipping 1 pixel border
    lout[i] = 0;		// set edges to zero: pixel j=0:
    l[i] = 0;
    one = _mm256_set1_epi32(1); // AVX
    for (j = 1; j < dim1 - 1 - 8; j = j + 8) {	// do 8 pixels at a time with sse2
      p = i + j;
      mxp = _mm256_loadu_ps(&im[p + o[1]]);	// AVX load 8 floats
      iqp = one;		// current selection
      ik = one;		// to increment
      for (k = 2; k < 10; k++) {
	ik = _mm256_add_epi32(ik, one);	// AVX2+ 
	mxq = _mm256_loadu_ps(&im[p + o[k]]);	// load vector of floats
	msk = _mm256_cmp_ps(mxq, mxp, _CMP_GT_OQ);	// AVX for q > p
	// apply mask to max (mxq/mxp) and to index (iqp/ik)
	mxp = _mm256_or_ps(_mm256_and_ps(msk, mxq), _mm256_andnot_ps(msk, mxp));
	imsk = cast256i(msk);	// AVX for q > p
	iqp = _mm256_or_si256(_mm256_and_si256(imsk, ik), _mm256_andnot_si256(imsk, iqp));
      }			// k neighbors
      // Write results (note epi16 is sse2)
      l[p]     = (uint8_t) _mm256_extract_epi16(iqp, 0);
      l[p + 1] = (uint8_t) _mm256_extract_epi16(iqp, 2);
      l[p + 2] = (uint8_t) _mm256_extract_epi16(iqp, 4);
      l[p + 3] = (uint8_t) _mm256_extract_epi16(iqp, 6);
      l[p + 4] = (uint8_t) _mm256_extract_epi16(iqp, 8);
      l[p + 5] = (uint8_t) _mm256_extract_epi16(iqp,10);
      l[p + 6] = (uint8_t) _mm256_extract_epi16(iqp,12);
      l[p + 7] = (uint8_t) _mm256_extract_epi16(iqp,14);
      // Count peaks in here? ... was better in separate loop
    }
    for (; j < dim1 - 1; j++) {	// end of simd loop, continues on j from for
      p = i + j;
      mx = im[p + o[1]];
      iq = 1;
      for (k = 2; k < 10; k++) {
	if (im[p + o[k]] > mx) {
	  iq = k;
	  mx = im[p + o[k]];
	}
      }
      l[p] = iq;
    }
    // set edges to zero: pixel j=dim1:
    lout[i + dim1 - 1] = 0;
    l[i + dim1 - 1] = 0;
  }				// i
}
#else
void neighbormax_avx2(const float *restrict im,	// input
		      int32_t * restrict lout,	// output
		      uint8_t * restrict l,	// workspace temporary
		      int dim0,	// Image dimensions
		      int dim1,
		      int o[10]){
  printf("Not compiled with AVX2!\n");
}
#endif
//AVX2

int localmaxlabel(const float *restrict im,	// input
		  int32_t * restrict lout,	// output
		  uint8_t * restrict l,	// workspace temporary
		  int dim0,	// Image dimensions
		  int dim1,
		  int cpu
		  )
{
    // old msvc for python 2.7 requires ALL variables declared up here.  
  int i, j, p, k, q, npk, t, nt, tid, lo, hi;
    int *npka;
    int noisy=0;
    double tic, toc;
    int o[10] = { 0, // special case
		  -1 - dim1, -dim1, +1 - dim1,	// 1,2,3
		  -1, 0, +1,		// 4,5,6
		  -1 + dim1, +dim1, +1 + dim1
    };				// 7,8,9

    // Fill in l to point to the maximum neighbor
    // Fills edges in lout to zero too
    //   cpu 0 = vanilla, silent
    //   cpu 1 = sse
    //   cpu 2 = avx2
    //   
    if( cpu > 9){
      cpu -= 10;
      noisy=1;
      tic = my_get_time();
    }
    switch(cpu){
    case 2:
      if(noisy) printf("Using avx2 instructions\n");
      neighbormax_avx2( im, lout, l, dim0, dim1, o);
      break;
    case 1:
      if(noisy) printf("Using sse2 instructions\n");
      neighbormax_sse2( im, lout, l, dim0, dim1, o);
      break;
    case 0:
      if(noisy) printf("No intel intrinsics\n");
      neighbormax( im, lout, l, dim0, dim1, o);
      break;
    default:
      printf("cpu must be 0, 1, 2: %d\n",cpu);
    }
    if( noisy ){
      toc = my_get_time();
      printf("    neighbormax %.3f ms\n", 1000*(toc-tic));
      tic = toc;
    }
    // Now go through the byte array checking for max values
    // (row by row to be parallel)
    npka = (int *)malloc(sizeof(int) * dim0);
#pragma omp parallel for private( i, j, p, k ) shared( npka ) 
    for (i = 0; i < dim0; i++) {
      k = 0;
      for (j = 0; j < dim1; j++) {
	p = dim1 * i + j;
	if (l[p] == 5) {	// pointing to self : k+1==6
	  k ++;
	}
      }
#pragma omp critical 
      npka[i] = k;
    }
    npk = 0;
    if( noisy ){
      toc = my_get_time();
      printf("    Counting maxes %.3f ms\n", 1000*(toc-tic));
      tic = toc;
    }
    // Cumulating so that npka[i] holds cumulative sums 
    for (i = 0; i < dim0; i++) {
	t = npk;
	npk += npka[i];
	npka[i] = t;
	//    printf("%d %d\n",npk,npka[i]);
    }
    // Second pass with row offsets in place
#pragma omp parallel for private( t, j, p ) 
    for (i = 0; i < dim0; i++) {
        t = npka[i];
	for (j = 0; j < dim1; j++) {
	    p = dim1 * i + j;
	    if (l[p] == 5) {	// pointing to self : k+1==5
		t = t + 1;
		lout[p] = t;	// final label
		l[p] = 0;	// done, so tagged as zero
	    }
	}
    }
    free(npka);
    if( noisy ){
      toc = my_get_time();
      printf("    relabel %.3f ms\n", 1000*(toc-tic));
      tic = toc;
    }
    //
    // Now make all point to their max
    // If we re-write the paths we cannot run in parallel
    //  ... this was the slowest part, so try openmp
    // 105 ms to always walk to max
    // 40 ms if path relabelled
    //
    // This is the same for all versions (optimised or not)
    //  ... perhaps re-write to be a manual loop and fill in
    //  ... the steps that are thread local
    
#pragma omp parallel private( q, i, tid, nt, k, lo, hi )
{
  tid = omp_get_thread_num();
  nt = omp_get_num_threads();
  lo = dim0*dim1*tid/nt;
  hi = dim0*dim1*(tid+1)/nt;
  
  for (i = lo; i < hi ; i++) {
    if (l[i] == 0)
      continue;		// done
    // Now we need to walk to find a label
    k = 0;
    q = i + o[l[i]];
    while (l[q]) {
      q = q + o[l[q]];
      k++;
    }			// Now q addresses a max or a correct label
    lout[i] = lout[q];	// take label from max
    if( k > 0){ // relabel the path taken while we know the top label value
      q = i + o[l[i]];
      while (l[q]) {
	if( (q >= lo) && (q < hi)){
	  l[q] = 0;
	  lout[q] = lout[i];
	}
	q = q + o[l[q]];	
      }			
    }
    l[i] = 0;		// sharing problems??
  }
 }
    if( noisy ){
      toc = my_get_time();
      printf("    write %.3f ms\n", 1000*(toc-tic));
      tic = toc;
    } 
    return npk;
}
