

# check for avx2 and import cImageD11 module
import platform, sys
from ImageD11 import cImageD11_sse2

print(__file__, cImageD11_sse2.got_AVX2)
if cImageD11_sse2.got_AVX2:  
    from ImageD11 import cImageD11_avx2
    print("Using AVX2 codes")
    from ImageD11.cImageD11_avx2 import *
else:
    from ImageD11.cImageD11_sse2 import *

