# import cImageD11 compiled module and patch the docstrings

import struct
import warnings
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

# Check for the use of openmp interactions with os.fork and multiprocessing

def check_multiprocessing():
    """ You cannot safely use os.fork together with threads.
    But the cImageD11 codes uses threads via openmp.
    So please use forkserver or spawn for multiprocessing.
    
    https://discuss.python.org/t/concerns-regarding-deprecation-of-fork-with-alive-threads/33555
    https://github.com/FABLE-3DXRD/ImageD11/issues/177
    
    > Developers should respond by adjusting their use of multiprocessing or concurrent.futures
    > to explicitly specify either the "forkserver" or "spawn" start methods via a context.
    """
    import multiprocessing
    # Problem cases are:
    #     child processes -> we will set num threads to 1
    if multiprocessing.parent_process() is not None:
        cimaged11_omp_set_num_threads( 1 ) 
        # people wanting Nprocs * Mthreads need to reset after import
        # OMP_NUM_THREADS is not going to work for them
        # how are we going to remember this in the future??
    #
    # now check for the fork issue
    if ((multiprocessing.get_start_method(allow_none=False) == 'fork') and    # we have the problem
        (multiprocessing.get_start_method(allow_none=True) is None) and       # by accident
         ('forkserver' in multiprocessing.get_all_start_methods())):          # so fix it
            multiprocessing.set_start_method('forkserver')
    if ((multiprocessing.get_start_method(allow_none=False) == 'fork') and    #  we have the problem
        (parent is not None)):
        # Tell them about it.
        warnings.warn(__doc__)



if cimaged11_omp_get_max_threads() == 0:
    # The code was compiled without openmp
    OPENMP = False
else:
    # Openmp threading can be used
    OPENMP = True
    check_multiprocessing()

    
# For 32 or 64 bits
nbyte = struct.calcsize("P")  # 4 or 8

if nbyte == 8:

    def put_incr(*a, **k):
        """redirects to put_incr64"""
        return put_incr64(*a, **k)


if nbyte == 4:

    def put_incr(*a, **k):
        """redirects to put_incr32"""
        return put_incr32(*a, **k)


# Add / fix the docstrings


def fix_doc(oldstring, to_be_added):
    """Adds a description of the function to the f2py string"""
    if oldstring.find("Wrapper for"):
        try:
            head, tail = oldstring.split("Wrapper for")
        except:
            print(oldstring)
            print(to_be_added)
            raise
        iname = tail.find("\n")
        return head + tail[:iname] + " " + to_be_added + tail[iname:]
    else:
        return oldstring + to_be_added


def fill_in_docstrings():
    for name in cImageD11_docstrings.__all__:
        doc = getattr(cImageD11_docstrings, name)
        if name in globals():
            func = globals()[name]
            func.__doc__ = fix_doc(func.__doc__, doc)


fill_in_docstrings()

assert verify_rounding(20) == 0, "Problem with cImageD11 fast rounding code"
