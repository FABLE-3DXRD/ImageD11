
from __future__ import print_function
import multiprocessing, sys, timeit, os

timer = timeit.default_timer
ncpu = multiprocessing.cpu_count()

r = {}

for ib in range(1,ncpu+1):
    for io in  range(1,ncpu+1):
        cmd0 = "OPENBLAS_NUM_THREADS=%d OMP_NUM_THREADS=%d "%(ib,io)
        cmd1 = "python idx.py > /dev/null"
        start = timer()
        x = os.system( cmd0+cmd1 )
        end = timer()
        print( cmd0, "%.2f/s "%(end-start))
        if x!=0:
            break

"""
jon@brainsIntel:~/Documents/ImageD11/test/simul_1000_grains$ python bencher.py 
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1  5.79/s 
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2  5.02/s 
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=3  5.17/s 
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4  5.51/s 
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=1  6.01/s 
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2  6.08/s 
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=3  5.92/s 
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=4  56.00/s 
OPENBLAS_NUM_THREADS=3 OMP_NUM_THREADS=1  6.02/s 
OPENBLAS_NUM_THREADS=3 OMP_NUM_THREADS=2  6.46/s 
OPENBLAS_NUM_THREADS=3 OMP_NUM_THREADS=3  23.10/s 
OPENBLAS_NUM_THREADS=3 OMP_NUM_THREADS=4  58.63/s 
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=1  11.46/s 
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=2  11.29/s 
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=3  20.34/s 
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4  60.67/s 
"""
