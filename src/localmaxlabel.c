
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "cImageD11.h"

#include <omp.h>

// Most compilers should have this!
#ifdef _MSC_VER
 #include <intrin.h>
 #define cast128i(X)  (_mm_castps_si128((X)))
 #ifdef __AVX__
  #include <immintrin.h>
  #define cast256i(X)  (_mm256_castps_si256(X))

/* _mm256_extractf128_si256 */

 #endif
#else				// gcc/clang
 #ifdef __SSE2__
  #include <x86intrin.h>
  #include <emmintrin.h>
  #include <xmmintrin.h>
  #define cast128i(X)  ((__m128i)(X))
  #ifdef __AVX__
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

#define pick(A,B,I,J) \
  if( (A) > (B) ){    \
  (B)=(A);            \
  (I)=(J);            \
  }

void neighbormax(const float *restrict im,	// input
		 int32_t * restrict lout,	// output
		 uint8_t * restrict l,	// workspace temporary
		 int dim0,	// Image dimensions
		 int dim1,
		 int o[10]){
  int i, j, p, iq, k0, k1, k2;
  float mx0, mx1, mx2;
    assert( (o[1]+1) == o[4] );
    assert( (o[2]+1) == o[5] );
    assert( (o[3]+1) == o[6] );
    assert( (o[1]+2) == o[7] );
    assert( (o[2]+2) == o[8] );
    assert( (o[3]+2) == o[9] );
    
    /* The image borders */
    for (i = 0; i < dim0; i++) {	// first and last row here:
	lout[i] = 0;		// ends of rows are below
	l[i] = 0;
	lout[dim1 * (dim0 - 1) + i] = 0;
	l[dim1 * (dim0 - 1) + i] = 0;
    }  
#pragma omp parallel for private( j, p, iq, k0, k1, k2, mx0, mx1, mx2 )
    for (i = dim1; i < (dim0 - 1) * dim1; i = i + dim1) {	// skipping 1 pixel border
	lout[i] = 0;		// set edges to zero: pixel j=0:
	l[i] = 0;

	p=i+1;
	mx0 = im[ p + o[1] ];
	k0 = 1;
	pick( im[ p + o[2] ] , mx0 , k0, 2 );
	pick( im[ p + o[3] ] , mx0 , k0, 3 );
	mx1 = im[ p + o[4] ];
	pick( im[ p + o[5] ] , mx1 , k1, 5 );
	pick( im[ p + o[6] ] , mx1 , k1, 6 );
	for (j = 1; j < dim1 - 1; j++) {
	    p = i + j;
	    mx2 = im[p + o[7]];
	    k2 = 7;
	    pick( im[ p + o[8] ] , mx2 , k2, 8 );
	    pick( im[ p + o[9] ] , mx2 , k2, 9 );
	    pick( mx1, mx0, k0, k1);
	    pick( mx2, mx0, k0, k2);
	    l[p] = k0;
	    mx0 = mx1;
	    k0 = k1-3;
	    mx1 = mx2;
	    k1 = k2-3;
	}
	// set edges to zero: pixel j=dim1:
	lout[i + dim1 - 1] = 0;
	l[i + dim1 - 1] = 0;
    }				// i
}





void neighbormax_sse2_old(const float *restrict im,	// input
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


#define pick_sse2(A,B,I,J,M)	        \
  (M) = _mm_cmpgt_ps( (A), (B) );       \
  (B) = _mm_blendv_ps( (B), (A), (M) ); \
  (I) = _mm_blendv_ps( (I), (J), (M) ); \



void neighbormax_sse2(const float *restrict im,	// input
		      int32_t * restrict lout,	// output
		      uint8_t * restrict l,	// workspace temporary
		      int dim0,	// Image dimensions
		      int dim1,
		      int o[10]){
  __m128 mxq, msk, one, iqp, ik, mxp;
  __m128i smsk;
  int i, j, p, iq, k;
  float mx;
  // Set edges to zero (background label)
  for (i = 0; i < dim0; i++) {	// first and last row here:
    lout[i] = 0;		// ends of rows are below
    l[i] = 0;
    lout[dim1 * (dim0 - 1) + i] = 0;
    l[dim1 * (dim0 - 1) + i] = 0;
  }
  smsk = _mm_set_epi8( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 8, 4, 0 );
#pragma omp parallel for private( j, p, iq, k, mx, mxp, iqp, ik, one, msk, mxq)
  for (i = dim1; i < (dim0 - 1) * dim1; i = i + dim1) {	// skipping 1 pixel border
    lout[i] = 0;		// set edges to zero: pixel j=0:
    l[i] = 0;
    p = i+1;
    one = _mm_set1_ps( 1.0 );    
    for (j = 1; j < dim1 - 1 - 4; j = j + 4) {	// do 4 pixels at a time with sse2
      p = i + j;
      mxp = _mm_loadu_ps(&im[p + o[1]]);	// SSE load 4 floats
      iqp = one;
      ik  = one;
      for (k = 2; k < 10; k++) {
	ik  = _mm_add_ps(ik, one );
	mxq = _mm_loadu_ps(&im[p + o[k]]);
	pick_sse2( mxq,  mxp, iqp, ik, msk);
      }			// k neighbors
      // Write results (note epi16 is sse2)
      _mm_stream_si32( (int*) &l[p], _mm_cvtsi128_si32(_mm_shuffle_epi8( _mm_cvtps_epi32( iqp ), smsk)));
//      _mm_stream_si32( (int*) &l[p],  _mm_cvtsi64_si32( _mm_cvtps_pi8( iqp )));
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


#ifdef __AVX__
void neighbormax_avx(const float *restrict im,	// input
		      int32_t * restrict lout,	// output
		      uint8_t * restrict l,	// workspace temporary
		      int dim0,	// Image dimensions
		      int dim1,
		      int o[10]){
  __m256i iqpi;
  __m128i iqplo, iqphi, smsk;
  __m256 mxp, mxq, msk, one, iqp, ik;
  int i, j, p, iq, k;
  float mx;
    smsk = _mm_set_epi8( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 8, 4, 0 );

  // Set edges to zero (background label)
  for (i = 0; i < dim0; i++) {	// first and last row here:
    lout[i] = 0;		// ends of rows are below
    l[i] = 0;
    lout[dim1 * (dim0 - 1) + i] = 0;
    l[dim1 * (dim0 - 1) + i] = 0;
  }
#pragma omp parallel for private( j, p, iq, k, mx, mxp, iqp, ik, one, msk, mxq, iqpi, iqplo, iqphi)
  for (i = dim1; i < (dim0 - 1) * dim1; i = i + dim1) {	// skipping 1 pixel border
    lout[i] = 0;		// set edges to zero: pixel j=0:
    l[i] = 0;
//    one = _mm256_set1_epi32(1); // AVX
    one = _mm256_set1_ps(1.0); // AVX
    for (j = 1; j < dim1 - 1 - 8; j = j + 8) {	// do 8 pixels at a time with sse2
      p = i + j;
      mxp = _mm256_loadu_ps(&im[p + o[1]]);	// AVX load 8 floats
      iqp = one;		// current selection
      ik = one;		// to increment
      for (k = 2; k < 10; k++) {
	ik = _mm256_add_ps(ik, one);	// AVX 
	mxq = _mm256_loadu_ps(&im[p + o[k]]);	// load vector of floats
	msk = _mm256_cmp_ps(mxq, mxp, _CMP_GT_OQ);	// AVX for q > p
	// apply mask to max (mxq/mxp) and to index (iqp/ik)
	mxp = _mm256_or_ps(_mm256_and_ps(msk, mxq), _mm256_andnot_ps(msk, mxp));
	iqp = _mm256_or_ps(_mm256_and_ps(msk,  ik), _mm256_andnot_ps(msk, iqp));
      }			// k neighbors
      iqpi =  _mm256_cvtps_epi32(iqp);
      iqplo = _mm256_extractf128_si256( iqpi, 0);
      iqphi = _mm256_extractf128_si256( iqpi, 1);
      _mm_stream_si32( (int*) &l[p], _mm_cvtsi128_si32(_mm_shuffle_epi8( iqplo , smsk)));
      _mm_stream_si32( (int*) &l[p+4], _mm_cvtsi128_si32(_mm_shuffle_epi8( iqphi , smsk)));
/*
      l[p+0] = (uint8_t) _mm_extract_epi16(iqplo, 0);
      l[p+1] = (uint8_t) _mm_extract_epi16(iqplo, 2);
      l[p+2] = (uint8_t) _mm_extract_epi16(iqplo, 4);
      l[p+3] = (uint8_t) _mm_extract_epi16(iqplo, 6);
      iqphi = _mm256_extractf128_si256( iqpi, 1);
      l[p+4] = (uint8_t) _mm_extract_epi16(iqphi, 0);
      l[p+5] = (uint8_t) _mm_extract_epi16(iqphi, 2);
      l[p+6] = (uint8_t) _mm_extract_epi16(iqphi, 4);
      l[p+7] = (uint8_t) _mm_extract_epi16(iqphi, 6);*/
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
void neighbormax_avx(const float *restrict im,	// input
		      int32_t * restrict lout,	// output
		      uint8_t * restrict l,	// workspace temporary
		      int dim0,	// Image dimensions
		      int dim1,
		      int o[10]){
    neighbormax_sse2( im,	// input
		       lout,	// output
		       l,	// workspace temporary
		      dim0,	// Image dimensions
		       dim1,
		       o);
}
#endif
//AVX

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
    /*    int o[10] = { 0, // special case
		  -1 - dim1, -dim1, +1 - dim1,	// 1,2,3
		  -1, 0, +1,		// 4,5,6
		  -1 + dim1, +dim1, +1 + dim1
		  }; */				// 7,8,9
    int o[10] = { 0, // special case
		  -1 - dim1, -1, -1 + dim1,
		  -dim1, 0, +dim1,
		  +1 - dim1, +1, +1 +dim1
    };				// 7,8,9

    // Fill in l to point to the maximum neighbor
    // Fills edges in lout to zero too
    //   cpu 0 = vanilla, silent
    //   cpu 1 = sse
    //   cpu 2 = avx
    //   
    if( cpu > 9){
      cpu -= 10;
      noisy=1;
      tic = my_get_time();
    }
    switch(cpu){
    case 2:
      if(noisy) printf("Using avx instructions\n");
      neighbormax_avx( im, lout, l, dim0, dim1, o);
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
