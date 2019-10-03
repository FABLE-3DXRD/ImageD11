

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

# Add / fix the docstrings

from ImageD11 import cImageD11_docstrings

def fix_doc( oldstring, to_be_added ):
    if oldstring.find("Wrapper for"):
        head, tail = oldstring.split("Wrapper for")
        iname = tail.find("\n")
        return head + tail[:iname] + " " + to_be_added + tail[iname:]
    else:
        return oldstring + to_be_added


for name in cImageD11_docstrings.__all__:
    doc = getattr( cImageD11_docstrings, name )
    if name in globals():
        func = globals()[name]
        func.__doc__ = fix_doc( func.__doc__, doc )
        
