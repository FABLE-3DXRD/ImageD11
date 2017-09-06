# Try some interfacing via cffi ...

from cffi import FFI


ffi = FFI()
ffi.set_source(
    "_ffi_diffraction",
    r"""
#include "omp.h"

double sum(double numbers[],int num_elements){
   int i;
   double s=0.0;
#pragma omp parallel for reduction(+: s)
   for (i=0; i<num_elements; i++)
   {
     s = s + numbers[i];
   }
   return s;
}
    """,
    Libraries = [] ,
    extra_compile_args = ["/openmp",], 
    extra_link_args = [], 
    )


ffi.cdef("""     
double sum(double [],int);
""")

if __name__ == "__main__":
    ffi.compile(verbose=True)

