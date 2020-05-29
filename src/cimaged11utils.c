
/* Various utility things that do not belong elsewhere */

#include "cImageD11.h"

#ifdef _OPENMP

#include <omp.h>

void cimaged11_omp_set_num_threads(int n) { omp_set_num_threads(n); }

int cimaged11_omp_get_max_threads(void) { return omp_get_max_threads(); }

#else

void cimaged11_omp_set_num_threads(int n) {}

int cimaged11_omp_get_max_threads(void) {
    return 0; /* e.g: zero parallelism via openmp */
}

#endif
/* F2PY_WRAPPER_START
    subroutine cimaged11_omp_set_num_threads(n)
!DOC cimaged11_omp_set_num_threads Sets the openmp number of
!DOC threads to use.
!DOC Change if you use multiprocessing or do not like
!DOC os.environ['OMP_NUM_THREADS']
!DOC see also: multiprocessing.cpu_count(),
!DOC see also os.cpu_count()
!DOC see docs/sphx/parallel.rst
        intent(c) cimaged11_omp_set_num_threads
        intent(c)
        integer, intent(in) :: n
    end subroutine cImaged11_set_omp_num_threads
F2PY_WRAPPER_END */

/* F2PY_WRAPPER_START
    function cimaged11_omp_get_max_threads()
!DOC cimaged11_omp_get_max_threads reads the openmp max number of threads
!DOC see docs/sphx/parallel.rst
        intent(c) cimaged11_omp_get_max_threads
        integer :: cimaged11_omp_get_max_threads
    end function cImaged11_omp_get_max_threads
F2PY_WRAPPER_END */

#if defined(_MSC_VER) || defined(__MINGW32__)

#include <windows.h>
double my_get_time() {
    LARGE_INTEGER t, f;
    QueryPerformanceCounter(&t);
    QueryPerformanceFrequency(&f);
    return (double)t.QuadPart / (double)f.QuadPart;
}

#else

#include <stdlib.h>
#include <sys/time.h>

double my_get_time() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec * 1e-6;
}
#endif
