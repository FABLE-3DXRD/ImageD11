from __future__ import print_function, division
from ImageD11 import cImageD11, indexing


try:
    from time import perf_counter_ns
except:
    from time import time
    def perf_counter_ns():
        return time()*1e9

import numpy as np
i = indexing.indexer()
i.readgvfile("eu3.gve")

gve = i.gv.copy().astype(float)
#gve = np.random.random( size=(100000,3) )-0.5
rgba = np.zeros( (512, 512, 4), np.uint8 )
u = np.eye(3, dtype=float).ravel()
t0 = perf_counter_ns()
cImageD11.splat( rgba, gve, u, 1 )
t1 = perf_counter_ns()
for j in range(5):
    cts = []
    for i in range(10):
        t2 = perf_counter_ns()
        cImageD11.splat( rgba, gve, u, j )
        t3 = perf_counter_ns()
        cts.append(t3-t2)
    print("Ran that for npx=",j,"s, fps:", 1e9/np.mean(cts),1e9/np.max(cts),1e9/np.min(cts))
import pylab as pl
pl.imshow(rgba)
pl.show()
