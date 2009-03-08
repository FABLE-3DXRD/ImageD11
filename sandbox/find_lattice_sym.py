

"""
Routines for finding the lattice symmetry from a ubi file
"""
import numpy as np

from ImageD11.lattice_reduction import mod as lattice_mod

def cosangle_rr_vec( ubi, v ):
    """
    Returns the cosine of the angle between the vector, v, in
    real and reciprocal space.

    In monoclinic-b symmetry the real space b and reciprocal space
    b* axes are parallel.
    This generalises, so all 2-fold lattice symmetries are parallel
    See:
    J. Appl. Cryst. (1982). 15, 255-259    [ doi:10.1107/S0021889882011959 ]
    "The derivation of the axes of the conventional unit cell from the
    dimensions of the Buerger-reduced cell"
    """
    real = np.dot( ubi.T, v )
    reci = np.dot( np.linalg.inv(ubi) , v )
    return np.dot( real, reci )/np.sqrt(
        np.dot(real, real) * np.dot(reci, reci) )

def find_2_folds( ubi , tol = 2.0, debug = False):
    """
    Scans over different directions to find all the potential
    2-fold axes.
    Default tolerance is 2 degrees.
    """
    hr = range(-2,3)
    tol_c = abs(1.0 - np.cos( tol * np.pi / 180.0))
    axes = []
    for h in hr:
        for k in hr:
            for l in hr:
                if h==0 and k==0 and l==0:
                    continue
                c = cosangle_rr_vec( ubi, [h,k,l] )
                if debug:
                    print h, k, l, c, np.arccos(c)*180/np.pi
                if abs(1 - c) < tol_c:
                    # print h, k, l, c, np.arccos(c)*180/np.pi
                    axes.append( (h,k,l) )
    return axes

def r2sweep( vl ):
    """
    A sweep of subtracting vectors from each other
    """
    vn = np.asarray( vl ).copy()
    nv = len(vl)
    for i in range(nv):
        for j in range(i+1, i+nv):
            k = j%nv
            assert i != k
            vn[k] = lattice_mod( vn[k], vn[i] )
    return vn


def reduce_two_folds( vecs ):
    """
    Try to make a uniq list of 2-fold axes (eg, recognise
    2x and -1x etc)
    """
    vn = r2sweep( vecs )
    uniq = []
    for v in vn:
        if np.dot( v, v) > 0:
            if v.sum( ) > 0:
                uniq.append( v )
            else:
                uniq.append( -v )
    return uniq

import xfab.tools

def test_monoclinic( ):
    ubi = xfab.tools.form_a_mat( [ 1, 2, 3, 89.8, 101, 90.02 ] )
    twofolds = find_2_folds( ubi )
    uniq = reduce_two_folds( twofolds )
    print "Monoclinic", uniq

def test_monoclinic_c( ):
    """ a=c , alpha=gamma """
    ubi = xfab.tools.form_a_mat( [ 4., 2., 4., 81., 95., 81. ] )
    twofolds = find_2_folds( ubi, debug = True )
    uniq = reduce_two_folds( twofolds )
    print "Monoclinic_c", uniq


def test_orthorhombic( ):
    ubi = xfab.tools.form_a_mat( [ 1, 2, 3, 89.8, 90.1, 90.1 ] )
    twofolds = find_2_folds( ubi )
    uniq = reduce_two_folds( twofolds )
    print "Orthorhombic", uniq

    

def test():
    test_monoclinic()
    test_monoclinic_c()
    test_orthorhombic()



    
if __name__=="__main__":
    test()
