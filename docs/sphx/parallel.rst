================================
Parallel code in ImageD11
================================

The main reason to write this code was to have something which is as fast,
or faster, than the data collection at the beamline. Performance is 
supposed to be a feature. If you cannot see what you are doing during
an experiment then you can waste a lot of beamtime. Over time this
is easier as the detectors get older.

Now, in 2020, ID11 is about to get a new detector and so speed is a 
bit of problem again.

Compiled C codes
----------------

In the beginning there were hand written wrappers which compiled 
against the Numeric python module. These were never ported
to numarray but following the great numpy reunification they
were ported to numpy. For a very long time we resisted
the move to python3 because it broke these wrappers. Accepting
that python3 is coming to end of life, we looked into a range
of options here : https://github.com/jonwright/cython_example

Many of the experts suggested cython (and now numba). There 
is also ctypes and cffi as plausible options. At the time we
had some Fortran code in place and this was being built by f2py.
Fortran has the huge advantage of being a standardised language.
The meaning of the code does not change when someone decides to 
deprecate a feature you are using. For long term maintenance
this is quite appealing. Sadly no one has Fortran comilers and
so it is painful to install. Also there are no unsigned or
small integer types.

The C code is nominally portable/ansi conforming and can run without
openmp if needed. There are a few nasty things which rely on ieee 
and an optimiser which is not too agressive. These should have 
testcases.

Looking to the future: it would be handy to find a nice parallel
language which is easier than opencl. The main drawback
of the current system is that you need to compile and optimise
the C code for the machine you are running on. If you run on a
cluster where each CPU is different then you miss having a JIT.


Performance measurements
------------------------

To check the C code is actually performing OK we plan try to set something
up using "airspeed velocity" : https://github.com/airspeed-velocity/asv

Over the years there have been a few different cases of introducing parallel
programming into ImageD11. They have all introduced regrettable bugs or 
other problems to distribute the code. It should be worth the cost.

To avoid exploding the git repository, this is planned for a seperate github 
project http://github.com/jonwright/ImageD11_benchmarks . A bunch 
of things that are currently in testcases (like examples and so on)
will probably go there. Eventually.

OpenMP Threading
----------------

Many of the compiled C codes have openmp loops inside. In the vast
majority of cases this is a simple for-loop decorator that should
not be too hard to debug. The main problem is when you are calling
several C codes at the same time from the main python process.
Or if you have a lot of python processes running at the same 
time, then you overload the CPU's. To avoid that you can:

- set OMP_NUM_THREADS to the number per process
- set dynamically via: cImageD11.cimaged11_omp_set_num_threads

Python Threading + Queue modules
--------------------------------

Used in ImageD11.peaksearcher. There is a thread for each of the following
tasks:
- reading frames (blocking IO)
- correcting frames for dark / flat
- one thread per threshold of a multi-threshold peaksearch

Within each of these tasks there may be underlying openmp threads.

Python Multiprocessing module
-----------------------------

Used in grid_index_parallel.py and also refine_em.py. Convenient
but remember to set OMP_NUM_THREADS to something sensible you
you get access to a big CPU. Problems show up as Ncpu^2 if
openmp and multiprocessing co-exist.

SIMD parallelism
----------------

Many CPU's have an instruction level parallelism that can 
give large speedups and also do what you would want. See, e.g.: 
https://software.intel.com/sites/landingpage/IntrinsicsGuide/

Take some examples from avx2::

  __m256i _mm256_subs_epu16 (__m256i a, __m256i b)
  Subtract packed unsigned 16-bit integers in b from packed unsigned 
  16-bit integers in a using saturation, and store the results in dst.

This is what you might want to subtract a dark image
and it does 16 pixels in a single instruction. Notably it 
sounds like it deals with the overflow problems correctly instead of 
wrapping. Another example, if you want to normalise a vector you have::

  __m256 _mm256_rsqrt_ps (__m256 a)
  Compute the approximate reciprocal square root of packed 
  single-precision (32-bit) floating-point elements in a, and store 
  the results in dst. The maximum relative error for this 
  approximation is less than 1.5*2^-12.
  
This one is harder to justify. It depends on how often you need to
repeat the computation as to whether you need it.

For various reasons of portability, these are NOT currently used 
in ImageD11 source code. If there were to be a nice portable library
that gave access to these across different CPUs with runtime
instruction selection then it could have a big impact.

The current hope is that the C compilers will use these things
if they are appropriate. Then we re-compile to code for each
machine and put it in a different conda or virtualenv. This
is controlled via CFLAGS when compiling.