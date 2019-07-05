

# check for avx2 and import cImageD11 module
import platform, sys, struct
from ImageD11 import cImageD11_sse2
if cImageD11_sse2.got_AVX:  
    from ImageD11 import cImageD11_avx
    print("Using AVX codes")
    from ImageD11.cImageD11_avx import *
else:
    from ImageD11.cImageD11_sse2 import *

# For 32 or 64 bits
nbyte = struct.calcsize("P") # 4 or 8

if nbyte == 8:
    put_incr = put_incr64

if nbyte == 4:
    put_incr = put_incr32
