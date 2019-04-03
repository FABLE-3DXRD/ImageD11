
#include <stdlib.h>
#include <stdio.h>
#include "cImageD11.h"

#include <omp.h>

// Most compilers should have this!
#ifdef _MSC_VER
#include <intrin.h>
#else				// gcc/clang
#ifdef __SSE2__
#include <x86intrin.h>
#else
#warning "no sse2"
#endif
#endif

int localmaxlabel(const float *restrict im,	// input
		  int32_t * restrict lout,	// output
		  uint8_t * restrict l,	// workspace temporary
		  int dim0,	// Image dimensions
		  int dim1)
{
    // Using c99 to move private variables into parallel sections
    // ...so most things are declared where needed
    int i;			// msvc requires openmp loop variable declared outside section
    // old msvc for python 2.7 requires ALL variables declared up here.  
    int j, p, iq, k, q, npk, t;
    int *npka;
    float mx;
    int o[9] = { -1 - dim1, -dim1, +1 - dim1,	// 0,1,2
	-1, 0, +1,		// 3,4,5
	-1 + dim1, +dim1, +1 + dim1
    };				// 6,7,8
#ifdef __SSE2__
    __m128i one, iqp, ik;
    __m128 mxp, mxq, msk;
#endif
    // Set edges to zero (background label)
    for (i = 0; i < dim0; i++) {	// first and last row here:
	lout[i] = 0;		// ends of rows are below
	l[i] = 0;
	lout[dim1 * (dim0 - 1) + i] = 0;
	l[dim1 * (dim0 - 1) + i] = 0;
    }
#ifdef __SSE2__
#pragma omp parallel for private( j, p, iq, k, mx, mxp, iqp, ik, one, msk, mxq)
#else
#pragma omp parallel for private( j, p, iq, k, mx )
#endif
    for (i = dim1; i < (dim0 - 1) * dim1; i = i + dim1) {	// skipping 1 pixel border
	lout[i] = 0;		// set edges to zero: pixel j=0:
	l[i] = 0;
	
#ifdef __SSE2__
	one = _mm_set1_epi32(1);
	//AVX __m256i _mm256_set1_epi32 (int a)
	for (j = 1; j < dim1 - 1 - 4; j = j + 4) {	// do 4 pixels at a time with sse2
	    p = i + j;
	    //AVX __m256 _mm256_loadu_ps (float const * mem_addr)
	    mxp = _mm_loadu_ps(&im[p + o[0]]);	// load 4 floats
	    iqp = one;		// current selection
	    ik = one;		// to increment
	    for (int k = 1; k < 9; k++) {
		ik = _mm_add_epi32(ik, one);	// ik++ 
		mxq = _mm_loadu_ps(&im[p + o[k]]);	// load vector of floats
		msk = _mm_cmpgt_ps(mxq, mxp);	// 1 for q > p
		//AVX __m256 _mm256_cmp_ps (__m256 a, __m256 b, const int imm8=_CMP_GT_OQ)
		// apply mask to max (mxq/mxp) and to index (iqp/ik)
		mxp = _mm_or_ps(_mm_and_ps(msk, mxq), _mm_andnot_ps(msk, mxp));
		//AVX __m256 _mm256_or_ps (__m256 a, __m256 b)
		//AVX __m256 _mm256_and_ps (__m256 a, __m256 b)
		//AVX __m256 _mm256_andnot_ps (__m256 a, __m256 b)
		iqp = (__m128i) _mm_or_ps(_mm_and_ps(msk, (__m128) ik),
					  _mm_andnot_ps(msk, (__m128) iqp));
	    }			// k neighbors
	    // Write results (note epi16 is sse2)
	    l[p] = (uint8_t) _mm_extract_epi16(iqp, 0);
	    //AVX int _mm256_extract_epi16 (__m256i a, const int index)
	    l[p + 1] = (uint8_t) _mm_extract_epi16(iqp, 2);
	    l[p + 2] = (uint8_t) _mm_extract_epi16(iqp, 4);
	    l[p + 3] = (uint8_t) _mm_extract_epi16(iqp, 6);
	    // Count peaks in here? ... was better in separate loop
	}
#else
	j = 1;			// init next loop
#endif				// SSE2
	for (; j < dim1 - 1; j++) {	// end of simd loop, continues on j from for
	    p = i + j;
	    mx = im[p + o[0]];
	    iq = 1;
	    for (k = 1; k < 9; k++) {
		if (im[p + o[k]] > mx) {
		    iq = k + 1;
		    mx = im[p + o[k]];
		}
	    }
	    l[p] = iq;
	}
	// set edges to zero: pixel j=dim1:
	lout[i + dim1 - 1] = 0;
	l[i + dim1 - 1] = 0;
    }				// i
    // Now go through the byte array checking for max values
    // (row by row to be parallel)
    npka = (int *)malloc(sizeof(int) * dim0);
#pragma omp parallel for private( j, p )
    for (i = 0; i < dim0; i++) {
	npka[i] = 0;
	for (j = 0; j < dim1; j++) {
	    p = dim1 * i + j;
	    if (l[p] == 5) {	// pointing to self : k+1==5
		npka[i] = npka[i] + 1;
	    }
	}
    }
    npk = 0;
    for (i = 0; i < dim0; i++) {
	t = npk;
	npk += npka[i];
	npka[i] = t;
	//    printf("%d %d\n",npk,npka[i]);
    }
    // Second pass with row offsets in place
#pragma omp parallel for private(j, p )
    for (i = 0; i < dim0; i++) {
	for (j = 0; j < dim1; j++) {
	    p = dim1 * i + j;
	    if (l[p] == 5) {	// pointing to self : k+1==5
		npka[i] = npka[i] + 1;
		lout[p] = npka[i];	// final label
		l[p] = 0;	// done tagged as zero
	    }
	}
    }
    free(npka);
    //
    // Now make all point to their max
    // If we re-write the paths we cannot run in parallel
    //  ... this was the slowest part, so try openmp
    // 105 ms to always walk to max
    // 40 ms if path relabelled
#pragma omp parallel for private( q )
    for (i = 0; i < dim0 * dim1; i++) {
	if (l[i] == 0)
	    continue;		// done
	// Now we need to walk to find a label
	q = i + o[l[i] - 1];
	while (l[q]) {
	    q = q + o[l[q] - 1];
	}			// Now q addresses a max or a correct label
	lout[i] = lout[q];	// take label from max
	l[i] = 0;		// sharing problem here??
    }
    return npk;
}

#define NOISY 0

#if defined(NOISY) && (NOISY == 1)
double tic, toc;
void clk(char *msg)
{
    toc = my_get_time();
    printf("%s %.6f ms\n", msg, 1e3 * (toc - tic));
    tic = my_get_time();
}
#else
#define clk(x)
#endif

int localmaxlabel_nosse(const float *restrict im,	// input
			int32_t * restrict lout,	// output
			uint8_t * restrict l,	// workspace temporary
			int dim0,	// Image dimensions
			int dim1)
{
    // Using c99 to move private variables into parallel sections
    // ...so most things are declared where needed
    int i;			// msvc requires openmp loop variable declared outside section
    // old msvc for python 2.7 requires ALL variables declared up here.  
    int j, p, iq, k, q, npk, t;
    int *npka;
    float mx;
    int o[9] = { -1 - dim1, -dim1, +1 - dim1,	// 0,1,2
	-1, 0, +1,		// 3,4,5
	-1 + dim1, +dim1, +1 + dim1
    };				// 6,7,8
    // Set edges to zero (background label)
    clk("localmaxlabel_nosse");
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
	    mx = im[p + o[0]];
	    iq = 1;
	    for (k = 1; k < 9; k++) {
		if (im[p + o[k]] > mx) {
		    iq = k + 1;
		    mx = im[p + o[k]];
		}
	    }
	    l[p] = iq;
	}
	// set edges to zero: pixel j=dim1:
	lout[i + dim1 - 1] = 0;
	l[i + dim1 - 1] = 0;
    }				// i
    // Now go through the byte array checking for max values
    // (row by row to be parallel)
    clk("First loop to neighbor");
    npka = (int *)malloc(sizeof(int) * dim0);
#pragma omp parallel for private( j, p )
    for (i = 0; i < dim0; i++) {
	npka[i] = 0;
	for (j = 0; j < dim1; j++) {
	    p = dim1 * i + j;
	    if (l[p] == 5) {	// pointing to self : k+1==5
		npka[i] = npka[i] + 1;
	    }
	}
    }
    clk("Second loop");
    npk = 0;
    for (i = 0; i < dim0; i++) {
	t = npk;
	npk += npka[i];
	npka[i] = t;
	//    printf("%d %d\n",npk,npka[i]);
    }
    clk("Reduction");
    // Second pass with row offsets in place
#pragma omp parallel for private(j, p )
    for (i = 0; i < dim0; i++) {
	for (j = 0; j < dim1; j++) {
	    p = dim1 * i + j;
	    if (l[p] == 5) {	// pointing to self : k+1==5
		npka[i] = npka[i] + 1;
		lout[p] = npka[i];	// final label
		l[p] = 0;	// done tagged as zero
	    }
	}
    }
    clk("Another pass");
    free(npka);
    //
    // Now make all point to their max
    // If we re-write the paths we cannot run in parallel
    //  ... this was the slowest part, so try openmp
    // 105 ms to always walk to max
    // 40 ms if path relabelled
#pragma omp parallel for private( q )
    for (i = 0; i < dim0 * dim1; i++) {
	if (l[i] == 0)
	    continue;		// done
	// Now we need to walk to find a label
	q = i + o[l[i] - 1];
	while (l[q]) {
	    q = q + o[l[q] - 1];
	}			// Now q addresses a max or a correct label
	lout[i] = lout[q];	// take label from max
	l[i] = 0;		// sharing problem here??
    }
    clk("Last pass");
    return npk;
}
