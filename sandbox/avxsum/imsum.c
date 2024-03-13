
#include <immintrin.h>
#include <emmintrin.h>
#include <stdint.h> 

int64_t imsum_avx2_2( uint16_t im[], int64_t n ){
  __m256i s0,s1,s2,s3;
  __m128i l0, l1;
 int i;
 int64_t sum;
 s0 = _mm256_setzero_si256();
 s1 = _mm256_setzero_si256();
 s2 = _mm256_setzero_si256();
 s3 = _mm256_setzero_si256();
 for( i = 0; i < n; i = i+16){
     /* read 4*16=64 bytes, so not aligned 
      * goes to lower 4*32=128 in t
      */
   l0 = _mm_loadu_si128( (__m128i*) &im[i] ) ;
   s0 = _mm256_add_epi64( s0, _mm256_cvtepu16_epi64 ( l0 ) );
   s1 = _mm256_add_epi64( s1, _mm256_cvtepu16_epi64 ( _mm_bsrli_si128( l0,  8 )));
   l1 = _mm_loadu_si128( (__m128i*) &im[i+8] ) ;
   s2 = _mm256_add_epi64( s2, _mm256_cvtepu16_epi64 ( l1 ) );
   s3 = _mm256_add_epi64( s3, _mm256_cvtepu16_epi64 ( _mm_bsrli_si128( l1,  8 )));
 }
 sum = _mm256_extract_epi64 (s0, 0) + \
       _mm256_extract_epi64 (s0, 1) + \
       _mm256_extract_epi64 (s0, 2) + \
       _mm256_extract_epi64 (s0, 3) ;
 sum +=_mm256_extract_epi64 (s1, 0) +     \
       _mm256_extract_epi64 (s1, 1) + \
       _mm256_extract_epi64 (s1, 2) + \
       _mm256_extract_epi64 (s1, 3);
  sum += _mm256_extract_epi64 (s2, 0) + \
       _mm256_extract_epi64 (s2, 1) + \
       _mm256_extract_epi64 (s2, 2) + \
       _mm256_extract_epi64 (s2, 3);
  sum+=      _mm256_extract_epi64 (s3, 0) + \
       _mm256_extract_epi64 (s3, 1) + \
       _mm256_extract_epi64 (s3, 2) + \
       _mm256_extract_epi64 (s3, 3) ;
 return sum;
}


int64_t imsum_avx2_1( uint16_t im[], int64_t n ){
  __m256i s0,s1;
  __m128i l, lo, hi, z;
 int i;
 int64_t sum;
 s0 = _mm256_setzero_si256();
 s1 = _mm256_setzero_si256();
 for( i = 0; i < n; i = i+8){
     /* read 4*16=64 bytes, so not aligned 
      * goes to lower 4*32=128 in t
      */
   l = _mm_loadu_si128( (__m128i*) &im[i] ) ;
   s0 = _mm256_add_epi64( s0, _mm256_cvtepu16_epi64 ( l ) );
   s1 = _mm256_add_epi64( s1, _mm256_cvtepu16_epi64 ( _mm_bsrli_si128( l,  8 )));
 }
 sum = _mm256_extract_epi64 (s0, 0) + \
       _mm256_extract_epi64 (s0, 1) + \
       _mm256_extract_epi64 (s0, 2) + \
       _mm256_extract_epi64 (s0, 3) + \
       _mm256_extract_epi64 (s1, 0) +     \
       _mm256_extract_epi64 (s1, 1) + \
       _mm256_extract_epi64 (s1, 2) + \
       _mm256_extract_epi64 (s1, 3); 
 return sum;
}

int64_t imsum_avx2_0( uint16_t im[], int64_t n ){
  __m256i s;
  __m128i l;
 int i;
 int64_t sum;
 s = _mm256_setzero_si256();
 for( i = 0; i < n; i = i+4){
     /* read 4*16=64 bytes, so not aligned 
      * goes to lower 4*32=128 in t
      */
   l = _mm_loadu_si128( (__m128i*) &im[i] ) ;
   s = _mm256_add_epi64( s, _mm256_cvtepu16_epi64 ( l ) );
 }
 sum = _mm256_extract_epi64 (s, 0) + \
       _mm256_extract_epi64 (s, 1) + \
       _mm256_extract_epi64 (s, 2) + \
       _mm256_extract_epi64 (s, 3); 
 return sum;
}



int64_t imsum_sse2_0( uint16_t im[], int64_t n ){
  __m128i s0,s1,s2,s3,t,s32,l;
  int i;
  int64_t sum;
  s0 = _mm_setzero_si128();
  s1 = _mm_setzero_si128();
  s2 = _mm_setzero_si128();
  s3 = _mm_setzero_si128();
 for( i = 0; i < n; i = i+8){
     /* read 4*16=64 bytes, so not aligned 
      * goes to lower 4*32=128 in t
      */
   l = _mm_loadu_si128( (__m128i*) &im[i] );
   
   s0 = _mm_add_epi64( s0, _mm_cvtepu16_epi64( l ));
   s1 = _mm_add_epi64( s1, _mm_cvtepu16_epi64( _mm_bsrli_si128( l,  4 )));
   s2 = _mm_add_epi64( s2, _mm_cvtepu16_epi64( _mm_bsrli_si128( l,  8 )));
   s3 = _mm_add_epi64( s3, _mm_cvtepu16_epi64( _mm_bsrli_si128( l, 12 )));
 }
 sum = _mm_extract_epi64 (s0, 0) + \
       _mm_extract_epi64 (s0, 1) + \
    _mm_extract_epi64 (s1, 0) + \
       _mm_extract_epi64 (s1, 1) + \
    _mm_extract_epi64 (s2, 0) + \
       _mm_extract_epi64 (s2, 1) + \
    _mm_extract_epi64 (s3, 0) + \
   _mm_extract_epi64 (s3, 1) ;
 return sum;
}

int64_t imsum_sse2_1( uint16_t im[], int64_t n ){
  __m128i s1,s2,t,s32,l;
  int i;
  int64_t sum;
  s1 = _mm_setzero_si128();
  s2 = _mm_setzero_si128();
 for( i = 0; i < n; i = i+8){
     /* read 4*16=64 bytes, so not aligned 
      * goes to lower 4*32=128 in t
      */
   l = _mm_loadu_si128( (__m128i*) &im[i] );
   s32 = _mm_add_epi32( _mm_cvtepu16_epi32( l ),
			_mm_cvtepu16_epi32( _mm_bsrli_si128( l,  8 )));
   s1 = _mm_add_epi64( s1, _mm_cvtepi32_epi64( s32));
   s2 = _mm_add_epi64( s2, _mm_cvtepi32_epi64( _mm_bsrli_si128( s32,  8 )));
 }
 sum = _mm_extract_epi64 (s1, 0) + \
   _mm_extract_epi64 (s1, 1)+	   \
   _mm_extract_epi64 (s2, 0)+	   \
   _mm_extract_epi64 (s2, 1);
 return sum;
}



int64_t imsum_c( uint16_t im[], int64_t n ){
  int64_t sum = 0;
  for( int i = 0; i < n; i++)
    sum += im[i];
  return sum;
}

int64_t imsum_c_simd( uint16_t im[], int64_t n ){
  int64_t sum = 0;
#pragma omp simd reduction(+:sum)
  for( int i = 0; i < n; i++ )
    sum += im[i];
  return sum;
}

