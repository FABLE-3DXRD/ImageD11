
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>


#include <omp.h>

// Most compilers should have this!
#ifdef __SSE2__
#include <x86intrin.h>
#else
#warning "No SSE2 instructions for lmw!"
#endif

int localmaxlabel( const float* restrict im, // input
		   int32_t*  restrict lout,  // output
		   uint8_t*  restrict l,     // workspace temporary
		   int dim0,                 // Image dimensions
		   int dim1){
  // Using c99 to move private variables into parallel sections
  // ...so most things are declared where needed
  int  o[9]={ -1-dim1, -dim1, +1-dim1,  // 0,1,2
	      -1     ,     0,      +1,  // 3,4,5
	      -1+dim1, +dim1, +1+dim1 };// 6,7,8
  
  // Set edges to zero (background label)
  for( int i=0; i<dim0; i++){ // first and last row here:
    lout[i]=0;                // ends of rows are below
    l[i]=0;
    lout[dim1*(dim0-1)+i]=0;
    l[dim1*(dim0-1)+i]=0;
  }
#ifdef __SSE2__
  __m128i one = _mm_set1_epi32(1);
#endif
#pragma omp parallel for 
  for( int i=dim1; i<(dim0-1)*dim1; i=i+dim1){    // skipping 1 pixel border
    lout[i]=0;     // set edges to zero: pixel j=0:
    l[i]=0;
    int j;
#ifdef __SSE2__
    for( j=1; j<dim1-1-4; j=j+4){ // do 4 pixels at a time with sse2
      int p = i + j;
      __m128 mxp = _mm_loadu_ps(&im[p+o[0]]); // load 4 floats
      __m128i iqp = one;                      // current selection
      __m128i ik  = one;                      // to increment
      for(int k=1;k<9;k++){
	ik = _mm_add_epi32( ik, one );          // ik++ 
	__m128 mxq = _mm_loadu_ps(&im[p+o[k]]); // load vector of floats
	__m128 msk = _mm_cmpgt_ps( mxq, mxp );  // 1 for q > p
	// apply mask to max (mxq/mxp) and to index (iqp/ik)
	mxp = _mm_or_ps(_mm_and_ps( msk, mxq ),
			_mm_andnot_ps( msk, mxp ));
	iqp = (__m128i) _mm_or_ps(_mm_and_ps( msk, (__m128) ik ),
				  _mm_andnot_ps( msk, (__m128) iqp ));
      } // k neighbors
      // Write results (note epi16 is sse2)
      l[p]   = (uint8_t) _mm_extract_epi16( iqp, 0 );
      l[p+1] = (uint8_t) _mm_extract_epi16( iqp, 2 );
      l[p+2] = (uint8_t) _mm_extract_epi16( iqp, 4 );
      l[p+3] = (uint8_t) _mm_extract_epi16( iqp, 6 );
      // Count peaks in here? ... was better in separate loop
      }
#else 
      j = 1; // init next loop
#endif // SSE2
    for(  ; j<dim1-1; j++ ){ // end of simd loop, continues on j from for
      int p = i + j;
      float mx = im[p+o[0]];
      int iq = 1;
      for(int k=1;k<9;k++){
	if(im[p+o[k]]>mx){
	  iq = k+1;
	  mx = im[p+o[k]];
	}
      }
      l[p] = iq;
    }
    // set edges to zero: pixel j=dim1:
    lout[i+dim0-1]=0;
    l[i+dim0-1]=0;
  } // i
  // Now go through the byte array checking for max values
  // (row by row to be parallel)
  int npka[dim0];
#pragma omp parallel for
  for( int i=0; i<dim0; i++){
    npka[i] = 0;
    for( int j=0; j<dim1; j++){
      int p = dim1*i + j;
      if( l[p] == 5 ){ // pointing to self : k+1==5
	npka[i] = npka[i]+1;
      }
    }
  }
  int npk = 0;
  for(int i=0;i<dim0;i++){
    int t = npk;
    npk += npka[i];
    npka[i] = t;
    //    printf("%d %d\n",npk,npka[i]);
  }
  // Second pass with row offsets in place
#pragma omp parallel for
  for( int i=0; i<dim0; i++){
    for( int j=0; j<dim1; j++){
      int p = dim1*i + j;
      if( l[p] == 5 ){ // pointing to self : k+1==5
	npka[i] = npka[i]+1;
	lout[p] = npka[i]; // final label
	l[p] = 0;        // done tagged as zero
      }
    }
  }
  //
  // Now make all point to their max
  // If we re-write the paths we cannot run in parallel
  //  ... this was the slowest part, so try openmp
  // 105 ms to always walk to max
  // 40 ms if path relabelled
#pragma omp parallel for 
  for( int i=0; i<dim0*dim1; i++){
    if( l[i] == 0 ) continue; // done
    // Now we need to walk to find a label
    int q = i + o[l[i] - 1];
    while( l[q] ){
      q = q + o[l[q] - 1];
    }   // Now q addresses a max or a correct label
    lout[i] = lout[q]; // take label from max
    l[i] = 0; // sharing problem here??
  }
  return npk;
}

