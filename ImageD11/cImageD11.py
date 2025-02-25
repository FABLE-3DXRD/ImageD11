# import cImageD11 compiled module and patch the docstrings

import os
import struct
import warnings
from ImageD11 import cImageD11_docstrings
import numba, numpy as np

try:
    from ImageD11._cImageD11 import *
    import ImageD11._cImageD11
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


def check_multiprocessing(patch=False):
    """
    I tried to suppress the "fork" method for multiprocessing, but I could not.

    You cannot safely use os.fork together with threads.
    But the cImageD11 codes uses threads via openmp, and you are importing them.
    So please use forkserver or spawn for multiprocessing.

    https://discuss.python.org/t/concerns-regarding-deprecation-of-fork-with-alive-threads/33555
    https://github.com/FABLE-3DXRD/ImageD11/issues/177

    > Developers should respond by adjusting their use of multiprocessing or concurrent.futures
    > to explicitly specify either the "forkserver" or "spawn" start methods via a context.
    """
    # We only do this if os.fork exists
    if not hasattr(os, "fork"):
        return
    import multiprocessing

    #
    # now check for the fork issue
    if not hasattr(multiprocessing, "get_start_method"):
        # You are on python2.7. Give up?
        warnings.warn(
            "python2.7 with ImageD11.cImageD11: for multiprocessing use spawn\n"
        )
        return
    #
    # This has a side effect of fixing the start method if allow_none was not True
    method = multiprocessing.get_start_method(allow_none=True)
    if method == "fork":
        warnings.warn(check_multiprocessing.__doc__)
    # Problem cases are:
    #     child processes -> oversubscribe -> we will set num threads to 1
    #     your numpy (etc) are also going to be a problem if this was needed
    parent = None
    if hasattr(multiprocessing, "parent_process"):
        parent = multiprocessing.parent_process()
        # only for python 3.8 and up
        if parent is not None:
            if "OMP_NUM_THREADS" not in os.environ:
                # people wanting Nprocs * Mthreads may need to reset after import
                # or use OMP_NUM_THREADS as usual
                cimaged11_omp_set_num_threads(1)
            if method in ("fork", None):
                # warn if we are a child process
                warnings.warn(check_multiprocessing.__doc__)
    if method is None and patch:
        # we should change it away from fork
        # this is going to break code which relied on fork
        # needs one of those "warn people for a year or two and then break it"
        # messages. Perhaps wait on python3.14 when multiprocessing will get
        # fixed and break everyones code at that time?
        poss = multiprocessing.get_all_start_methods()
        if "forkserver" in poss:
            multiprocessing.set_start_method("forkserver")
        elif "spawn" in poss:
            multiprocessing.set_start_method("spawn")
        else:
            raise Exception("Could not set to forkserver or spawn")
    # It is doubtful that this function can ever please everyone.


def cores_available():
    """
    Return the number of CPU cores you can use
    First choice = os.sched_getaffinity
    Second choice = os.cpu_count (may be high)
    Third choice = os.environ['NUMBER_OF_PROCESSORS'] (windows+old)
    Guess at 1 when os.cpu_count someone returns None or 0
    """
    if hasattr(os, "sched_getaffinity"):
        ncpu = len(os.sched_getaffinity(os.getpid()))
    elif hasattr(os, "cpu_count"):
        ncpu = os.cpu_count()
    elif 'NUMBER_OF_PROCESSORS' in os.environ:
        ncpu = int(os.environ['NUMBER_OF_PROCESSORS'])
    if ncpu is None or ncpu < 1:
        return 1
    else:
        return ncpu


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



# dsb = np.floor( cf.ds * istep ).astype(int).clip(0,hbins-1)
@numba.njit(parallel=True)
def array_bin( ary, inverse_step, nbins ):
    """ Find integer bins for array reals from 0 to nbins and 1/step
    returns int( floor( ary * inverse_step ) ).clip(0,nbins-1)
    """
    out = np.zeros( ary.shape[0], dtype=np.dtype(np.intp) )
    for i in numba.prange( len(ary) ):
        j = int(np.floor(ary[i] * inverse_step)) # find bin
        j = min(max(0,j), nbins-1) # clip to range
        out[i] = j
    return out


@numba.njit(parallel=True)
def array_lt( ary, cut, out):
    """ parallel ary < cut """
    for i in numba.prange(ary.size):
        out[i] = ary[i] < cut
