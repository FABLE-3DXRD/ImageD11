

from __future__ import print_function, division

import ctypes, numpy, sys, time, os
if 0: 
   os.system("gcc -shared -fPIC -O3 -fopenmp -march=native imsum.c -o imsum.so")

timer = time.perf_counter

data = (numpy.random.random( 2048*2048*32 )*10000).astype(numpy.uint16)
if "win" in sys.platform:
    lib = ctypes.CDLL( "imsum.dll" )
    cfuns = [lib.imsum_sse2_0 ,lib.imsum_sse2_1 , lib.imsum_c_simd,
            lib.imsum_c,lib.imsum_avx2_0,lib.imsum_avx2_1,lib.imsum_avx2_2]
else:
    lib = ctypes.CDLL( "imsum.so" )
    cfuns = [lib.imsum_sse2_0 ,lib.imsum_sse2_1 , lib.imsum_c_simd,
            lib.imsum_c]
    if "avx2" in open("/proc/cpuinfo").read():
        print("You have avx2")
        cfuns += [lib.imsum_avx2_0,lib.imsum_avx2_1,lib.imsum_avx2_2]

lib.imsum_sse2_0.argtypes = [ ctypes.POINTER( ctypes.c_uint16 ), ctypes.c_int64 ]
lib.imsum_sse2_0.restype  = ctypes.c_int64
lib.imsum_sse2_1.argtypes = [ ctypes.POINTER( ctypes.c_uint16 ), ctypes.c_int64 ]
lib.imsum_sse2_1.restype  = ctypes.c_int64
lib.imsum_avx2_0.argtypes = [ ctypes.POINTER( ctypes.c_uint16 ), ctypes.c_int64 ]
lib.imsum_avx2_0.restype  = ctypes.c_int64
lib.imsum_avx2_1.argtypes = [ ctypes.POINTER( ctypes.c_uint16 ), ctypes.c_int64 ]
lib.imsum_avx2_1.restype  = ctypes.c_int64
lib.imsum_avx2_2.argtypes = [ ctypes.POINTER( ctypes.c_uint16 ), ctypes.c_int64 ]
lib.imsum_avx2_2.restype  = ctypes.c_int64
lib.imsum_c.argtypes = [ ctypes.POINTER( ctypes.c_uint16 ), ctypes.c_int64 ]
lib.imsum_c.restype  = ctypes.c_int64
lib.imsum_c_simd.argtypes = [ ctypes.POINTER( ctypes.c_uint16 ), ctypes.c_int64 ]
lib.imsum_c_simd.restype  = ctypes.c_int64




    
val = data.sum(dtype=numpy.int64)
    
N=100
t0=timer()
vals = [ data.sum(dtype=numpy.int64) for i in range(N)]
t1=timer()
avg = (t1-t0)/N
print("numpy",avg)
for cfunc in cfuns[::-1]:
    t0=timer()
    valc = [cfunc( data.ctypes.data_as( ctypes.POINTER( ctypes.c_uint16 ) ),
                   data.shape[0] )
            for i in range(N)]
    t1=timer()
    avc = (t1-t0)/N
    ok = (numpy.array(valc) == vals).all()
    print(cfunc.__name__,"%.6f"%(avc),"%.6f"%(avg/avc),ok)
    if not ok:
        print(valc[0],vals[0])

print(len(data)*5000)

"""
  gcc -shared -fPIC -Ofast -fopenmp -march=native -fprofile-generate=. imsum.c -o imsum.so
  python3 ./bench_imsum.py 
  gcc -shared -fPIC -Ofast -fopenmp -march=native -fprofile-use=. imsum.c -o imsum.so
  python3 ./bench_imsum.py 

  . /data/id11/jon/spack/
  spack load intel-oneapi-compilers@2021.3.0 
  icc -prof-gen -shared -fPIC -Ofast -fopenmp -xHost imsum.c -o imsum.so
  python3 ./bench_imsum.py 
  icc -prof-use -shared -fPIC -Ofast -fopenmp -xHost imsum.c -o imsum.so
  python3 ./bench_imsum.py 

  /usr/bin/clang -shared -fPIC -Ofast  -march=native imsum.c -o imsum.so -fprofile-instr-generate
  python3 bench_imsum.py 
  /usr/bin/llvm-profdata-10 merge default.profraw -o default.profdata
  /usr/bin/clang -shared -fPIC -Ofast  -march=native imsum.c -o imsum.so -fprofile-instr-use
  python3 bench_imsum.py 
 """