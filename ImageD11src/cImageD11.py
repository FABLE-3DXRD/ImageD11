

# check for avx2 and import cImageD11 module
import platform, sys, struct, logging
from ImageD11 import cImageD11_docstrings


from ImageD11.cImageD11_safe import *
try:
    from ImageD11.cImageD11_fast import *
except ImportError:
    pass


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



