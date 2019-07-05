from __future__ import division

import numpy as np
cimport numpy as np

cdef extern from "cImageD11.h":
     void compute_xlylzl( float s[], float f[], double p[4],
		   	    double r[9], double dist[3],
		  	       float xlylzl[][3], int n)


def compute_xlylzl_wrap(
    np.ndarray[np.float32_t, ndim=1] s,
    np.ndarray[np.float32_t, ndim=1] f,
    np.ndarray[np.float_t,   ndim=1] p,
    np.ndarray[np.float_t,   ndim=1] r,
    np.ndarray[np.float_t,   ndim=1] dist,
    np.ndarray[np.float32_t, ndim=2] xlylzl,
    ):
    assert xlylzl.shape[1] == 3
    assert dist.shape[0] == 3
    assert r.shape[0] == 9
    assert p.shape[0] == 4
    assert xlylzl.shape[0] == s.shape[0] == f.shape[0]
    n = xlylzl.shape[0]
    compute_xlylzl( <float *> s.data, <float *> f.data, 
                    <double *> p.data, <double *> r.data, 
		    <double *> dist.data ,  <float (*)[3]> xlylzl.data, n )