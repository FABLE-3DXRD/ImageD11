/* Try to collect all the compiler workarounds in one place */

#ifndef _cImageD11_h
#define _cImageD11_h

#ifdef __GNUC__
 #if __GNUC__ >= 4
   #define DLL_PUBLIC __attribute__ ((visibility ("default")))
   #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
 #else
   #define DLL_PUBLIC
   #define DLL_LOCAL 
 #endif
#else
   #define DLL_PUBLIC
   #define DLL_LOCAL 
#endif

DLL_LOCAL
double my_get_time(void);

#ifdef _MSC_VER
 typedef __int8 int8_t;
 typedef __int16 int16_t;
 typedef __int32 int32_t;
 typedef unsigned char uint8_t;
 typedef unsigned __int16 uint16_t;
 typedef unsigned __int32 uint32_t;
 #define restrict __restrict
 #define inline __inline
#else
#include <stdint.h>

/* If we define functions as local they can be inlined at link time
 * in a shared library (e.g. not shared and overridden by LD_PRELOAD)
 */

#endif

#endif
