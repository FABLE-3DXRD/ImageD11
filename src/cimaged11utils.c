
/* Various utility things that do not belong elsewhere */

#include "cImageD11.h"

#ifdef _OPENMP

#include <omp.h>

void cimaged11_omp_set_num_threads( int n){
  omp_set_num_threads(n);
}

int cimaged11_omp_get_max_threads( void ){
  return omp_get_max_threads();
}

#else

void cimaged11_omp_set_num_threads( int n){
}

int cimaged11_omp_get_max_threads( void ){
  return 0; /* e.g: zero parallelism via openmp */
}

#endif


#if defined(_MSC_VER) || defined(__MINGW32__)

#include <windows.h>
double my_get_time()
{
    LARGE_INTEGER t, f;
    QueryPerformanceCounter(&t);
    QueryPerformanceFrequency(&f);
    return (double)t.QuadPart / (double)f.QuadPart;
}

#else

#include <stdlib.h>
#include <sys/time.h>

double my_get_time()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec * 1e-6;
}
#endif
