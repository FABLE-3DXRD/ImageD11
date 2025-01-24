
/* Compile as
 * cl diffracCl.c /O2 /openmp /LD 
 * cl diffracCl.c /Ox /fp:fast /openmp /favor:AMD64 /Qfast_transcendentals /LD
 * mt -manifest diffracCl.dll.manifest -outputresource:diffracCl.dll;2
 * 0.21 s
 *
 * or
 * gcc diffracCl.c -O3 -fopenmp -march=native -ffast-math -static-libgcc -shared -o diffracCl.dll
 * 0.31 s
 */

/*
 * Center in fast, slow
 * size in fast, slow 
 * distance
 * tilts
 */


/* Positions in parameter array */

#define FC 0
#define FS 1
#define SC 2
#define SS 3
#define DISTANCE 4
#define R 5


#include <math.h>


#define PI 3.14159265358979323846f

#ifdef _OPENMP
#include <omp.h>
#endif


#ifdef _WIN32
/* This line is disgusting */
__declspec(dllexport)
#endif
void ttheta(  float * restrict tth,
              float * restrict eta,
              float * restrict p ,
              int ns,
              int nf)
{
    // data[i,j] => i == slow, j == fast
    int i, j, a,row;
    float s,f,x,y,z,th,et,OPI;
    OPI = 1.0/PI;
#pragma omp parallel for private(s,f,x,y,z,j,a,i)
    for(i=0; i<ns; i++) {
        s  = (i - p[SC]) * p[SS];
	a = i*nf;
        for(j=0; j<nf; j++, a++){
            f  = (j - p[FC]) * p[FS];
            /*  |R0  R1  R2| |0|
                |R3  R4  R5|.|f|
                |R6  R7  R8| |s| */
	    z  = p[R+7]*f + p[R+8]*s;
	    y  = p[R+4]*f + p[R+5]*s;
	    eta[a] = atan2f( z, y )*OPI;
	    x  = p[R+1]*f + p[R+2]*s + p[DISTANCE];
	    tth[a] = atan2f( sqrtf(y*y+z*z)  , x )*OPI;
        }
    }
}
    

