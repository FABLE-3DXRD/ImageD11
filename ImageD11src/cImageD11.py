

# check for avx2 and import cImageD11 module
import platform, sys, struct, logging
from ImageD11 import cImageD11_docstrings
if platform.machine() == 'x86_64':
    from ImageD11 import cImageD11_sse2
    if cImageD11_sse2.got_AVX:
        from ImageD11 import cImageD11_avx
        logging.info("Using AVX codes")
        from ImageD11.cImageD11_avx import *
    else:
        from ImageD11.cImageD11_sse2 import *
elif platform.machine() == 'ppc64le':
    from ImageD11.cImageD11_ppc64le import *

# For 32 or 64 bits
nbyte = struct.calcsize("P") # 4 or 8

if nbyte == 8:
    def put_incr(*a,**k):
        """ redirects to put_incr64 """
        return put_incr64(*a,**k)

if nbyte == 4:
    def put_incr(*a,**k):
        """ redirects to put_incr32 """
        return put_incr32(*a,**k)

# Add / fix the docstrings

def fix_doc( oldstring, to_be_added ):
    """ Adds a description of the function to the f2py string """
    if oldstring.find("Wrapper for"):
        try:
            head, tail = oldstring.split("Wrapper for")
        except:
            print( oldstring )
            print( to_be_added )
            raise
        iname = tail.find("\n")
        return head + tail[:iname] + " " + to_be_added + tail[iname:]
    else:
        return oldstring + to_be_added

def fill_in_docstrings():
    for name in cImageD11_docstrings.__all__:
        doc = getattr( cImageD11_docstrings, name )
        if name in globals():
            func = globals()[name]
            func.__doc__ = fix_doc( func.__doc__, doc )

fill_in_docstrings()



