

# Try some interfacing via cffi ...

from _ffi_diffraction import ffi, lib
import numpy as np
import time
np.random.seed(42)
#       1K    1M     1G
size = 128 * 1024 * 1024 * 0.5
ar = np.random.random( int(size) ) 

print ar.shape


N=10
ar.sum()
start = time.time()
for i in range(N):
    sn = ar.sum()
print "Numpy",sn, (time.time() - start)/N
#raw_input()

N=10
lib.sum( ffi.cast('double *', ar.ctypes.data ),  len(ar) )
start = time.time()
for i in range(N):
    s = lib.sum( ffi.cast('double *', ar.ctypes.data ),  len(ar) )
print "FFI", s , (time.time()-start)/N
