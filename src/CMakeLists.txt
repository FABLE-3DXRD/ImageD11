
cmake_minimum_required(VERSION 3.0)


project(cImageD11lib)

set( COPTS -O2 -fPIC -fopenmp )
set( SSE2  -msse2 )
set( AVX2  -mavx2 )

if( MSVC )
  set( COPTS /Ox /WL /MD /openmp )
  set( SSE2  /arch:sse2 )
  set( AVX2  /arch:avx2 )
endif( MSVC )


include_directories( . )

set(CMAKE_C_FLAGS "-std=c99 -Wall")

add_library (cImageD11_sse2 STATIC blobs.c
  check_cpu_auto.c closest.c darkflat.c localmaxlabel.c
  cdiffraction.c connectedpixels.c deriv.c sparse_image.c)

target_compile_options(cImageD11_sse2 PRIVATE -msse2 ${COPTS})

add_library (cImageD11_avx2 STATIC blobs.c
  check_cpu_auto.c closest.c darkflat.c localmaxlabel.c
  cdiffraction.c connectedpixels.c deriv.c sparse_image.c)

target_compile_options(cImageD11_avx2 PRIVATE -mavx2 ${COPTS})



