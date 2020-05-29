/* Try to collect all the compiler workarounds in one place */

#ifndef _cImageD11_h
#define _cImageD11_h

#include <stdlib.h>

/* These did not expand properly anyway
 * #if _OPENMP >= 201307
 * #define SIMDCLAUSE simd
 * #define SIMDFOR omp simd
 * #else
 * #define SIMDCLAUSE
 * #define SIMDFOR ignore
 * #endif
 */

#ifdef _OPENMP
#include <omp.h>
#if _OPENMP >= 201307
#define GOT_OMP_SIMD
#endif
#endif

/* No really, thanks for standardising C99 20 years ago */
#define restrict __restrict

/* If we define functions as local they can be inlined at link time
 * in a shared library (e.g. not shared and overridden by LD_PRELOAD)
 * unlcear if/how this works currently across compilers...
 */
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#if __GNUC__ >= 4
#define DLL_PUBLIC __attribute__((visibility("default")))
#define DLL_LOCAL __attribute__((visibility("hidden")))
#else
#define DLL_PUBLIC
#define DLL_LOCAL
#endif
#else
#define DLL_PUBLIC
#define DLL_LOCAL
#endif

#ifdef _MSC_VER
typedef __int8 int8_t;
typedef __int16 int16_t;
typedef __int32 int32_t;
typedef __int64 int64_t;
typedef unsigned char uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;

#define inline __inline
#else
#include <stdint.h>
#endif

/* implemented in imaged11utils.c */
void cimaged11_omp_set_num_threads(int);
int cimaged11_omp_get_max_threads(void);
DLL_LOCAL
double my_get_time(void);

#endif
