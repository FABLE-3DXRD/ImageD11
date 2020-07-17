

from __future__ import print_function, division

import numpy as np
import timeit

timer = timeit.default_timer

# spots at x,y,z positions in 3D space.
# computed tth/eta in the past
# ... tth = angle with X axis,            atan2( sqrt( Y*Y + Z*Z), X ) 
# ... eta = angle with Z axis, sense from atan2( -Y, Z )

# later we use:
#  -sin(theta)
#  -cos(theta) sin(eta)   [ note this is also minus due to the earlier -Y ]
#   cos(theta) cos(eta)

# Given X, Y, Z
def k_from_xyz(xyz, wvln):
    """
    x,y,z = 3D peak positions, origin as 0,0,0
    wvln  = wavelength
    returns the k-vectors

    this was already buried inside src/cdiffraction.c
    """
    w1 = 1.0 / wvln
    M  = w1 /  np.linalg.norm( xyz, axis=0 ) 
    k  = xyz * M
    k[0] -=  w1
    return k 

from ImageD11.transform import compute_k_vectors, compute_tth_eta_from_xyz


def test1():
    np.random.seed(42)
    w = 40 / 12.
    nv = 1024*1024
    xyz = 1e5 * (np.random.random( (nv, 3) ) + ( 0.2, -0.5, -0.5 )).T
    old = [] ; new = []
    for i in range(5):
        start = timer()
        t, e = compute_tth_eta_from_xyz( xyz, None )
        k = compute_k_vectors( t, e, w )
        old.append( timer()-start )
        knew = k_from_xyz(xyz, w )
        start = timer()
        ok = np.allclose( k, knew )
        new.append( timer() - start )
    if ok:
        print("Matches", nv)
    else:
        print("Fails")
    print("old perf", np.min(old), np.mean(old), np.std(old))
    print("new perf", np.min(new), np.mean(new), np.std(new))
    print("speedup", np.min(old)/np.min(new))    


if __name__ == "__main__":
    test1()


