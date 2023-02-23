

# import cImageD11 compiled module and patch the docstrings

import sys
import struct
from ImageD11 import cImageD11_docstrings

try:
    from ImageD11._cImageD11 import *
except ImportError:
    print("Failed to import compiled C code for cImageD11 module")
    print("Are you running from a source directory?")
    print("Try something like:")
    print("   python -m pip install --editable .")
    print("or:")
    print("   python setup.py develop")
    print("or:")
    print("   python setup.py build_ext --inplace")
    raise

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

assert verify_rounding(20) == 0, "Problem with cImageD11 fast rounding code"

