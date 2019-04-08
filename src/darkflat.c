

/* To compare to numpy

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

#include <immintrin.h>

#ifdef __FMA__
void darkflat_fma( float *img,
		   float *drk,
		   float *flm,
		   int n){
  /* in place dark / flat */
  // assert aligned
  __m256 i8, d8, f8, o8;
  int i;
#pragma omp parallel for private( i8, f8, d8, o8 )
  for( i = 0; i < n ; i = i+8 ){
    i8 = _mm256_load_ps( &img[i] );
    f8 = _mm256_load_ps( &flm[i] );
    d8 = _mm256_load_ps( &drk[i] );
    // a*b + c
    o8 = _mm256_fmadd_ps (i8, f8, d8);
    _mm256_store_ps( &img[i], o8 );
  }
}
#endif	       

#ifdef __AVX__
void darkflat_avx( float *img,
		   float *drk,
		   float *flm,
		   int n){
  /* in place dark / flat */
  // assert aligned
  __m256 i8, d8, f8, o8;
  int i;
#pragma omp parallel for private( i8, f8, d8, o8 )
  for( i = 0; i<n; i=i+8){
    i8 = _mm256_load_ps( &img[i] );
    f8 = _mm256_load_ps( &flm[i] );
    d8 = _mm256_load_ps( &drk[i] );
    // a*b + c
    o8 = _mm256_sub_ps( _mm256_mul_ps (i8, f8) , d8);
    _mm256_store_ps( &img[i], o8 );
  }
}
#endif


