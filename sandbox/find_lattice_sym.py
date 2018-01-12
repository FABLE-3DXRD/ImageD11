
from __future__ import print_function


"""
Routines for finding the lattice symmetry from a ubi file
"""
import numpy as np



def cosangle_rr_vec( ubi, v1, v2 ):
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
    real = np.dot( ubi.T, v1 )
    reci = np.dot( np.linalg.inv(ubi) , v2 )
    return np.dot( real, reci )/np.sqrt(
        np.dot(real, real) * np.dot(reci, reci) )

# All the possible indices of lattice rows and planes
rows_planes = [ (h, k, l) for h in range(-2,3)
                          for k in range(-2,3)
                          for l in range(-2,3) ]
rows_planes.remove( (0,0,0) )
rows_planes = np.array( rows_planes )

from ImageD11.rc_array import rc_array

def find_2_folds( ubi , tol = 1.0, debug = False):
    """
    Scans over different directions to find all the potential
    2-fold axes.
    Default tolerance is 0.5 degrees.
    """
    tol_c = abs(1.0 - np.cos( tol * np.pi / 180.0))
    ub = np.linalg.inv( ubi )

    real_vecs = np.dot( ubi, rows_planes.T ).T
    reci_vecs = np.dot( ub.T , rows_planes.T ).T
    
    real_lens = np.sqrt((real_vecs*real_vecs).sum(axis = 1))
    reci_lens = np.sqrt((reci_vecs*reci_vecs).sum(axis = 1))
    nvec = len(rows_planes)
    pairs = []
    for v1, l1, i in zip(real_vecs, real_lens, list(range(nvec))):
        print(v1,l1,i)
        for v2, l2, j in zip(reci_vecs, reci_lens, list(range(nvec))):
            #if j < i:
            #    continue
            c = np.dot(v1, v2)/ l1/ l2
            if abs( 1 - c ) < tol_c:
                print(i, j, c, rows_planes[i],rows_planes[j])
                pairs.append( (i,j) )
    print(pairs)
    return pairs

    


from ImageD11.lattice_reduction import mod as lattice_mod


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


mono = ( 1, 2, 3, 89.8, 107, 90.02 )
ortho = ( 1.5, 2, 3, 90., 90., 90. )
def transform( cell, matrix):
    return xfab.tools.ubi_to_cell(
        np.dot(  np.linalg.inv(matrix), 
                 xfab.tools.form_a_mat( cell )))

testcases = {
    
    'monoclinic'   : [ mono , (0,1,0) ],
    'monoclinic_c'   : [ transform( mono, [[1,1,0],[-1,1,0],[0,0,1]]),
                         (1,-1,0) ],
    
    'orthorhombic' : [ ortho ,
                       (1,0,0),  (0,1,0),  (0,0,1) ],
    
    'orthorhombic_c' : [transform( ortho, [[1,1,0],[-1,1,0],[0,0,1]]) ,
                        (1,1,0),  (1,-1,0),  (0,0,1) ],
    
    'orthorhombic_i' : [ transform( ortho, [[0,1,1],[1,0,1],[1,1,0]]),
                          (1,1,-1),  (1,-1,1),  (-1,1,1) ],
    
    'orthorhombic_f' : [ transform( ortho, [[-1,1,1],[1,-1,1],[1,1,-1]] ),
                         (1,0,1),  (0,1,1),  (1,1,0) ],
    
    
    }


def test_reduce( name, cell , expected ):
    print(("%6.3f  "*6)%tuple(cell), end=' ')
    ubi = xfab.tools.form_a_mat( cell ).T
    twofolds = find_2_folds( ubi )
    lines  = [ rows_planes[t[0]] for t in twofolds]
    planes = [ rows_planes[t[1]] for t in twofolds]
    uniql = reduce_two_folds( planes )
    uniqp = reduce_two_folds( lines ) 
    # print name, uniq
    if len(expected) != len( uniql ):
        print(name + " wrong number of axes ")
        print(uniql)
        print(uniqp)
    for ul, up in zip(uniql, uniqp):
        ok = (tuple(ul) in expected) or (tuple(-ul) in expected)
        if not ok:
            print("PROBLEM",name, ul, up, expected)
            
def runtests():
    np.seterr(all='raise')
    names = list(testcases.keys())
    names.sort()
    for name in names:
        cell = testcases[name][0]
        expected = testcases[name][1:]

        print("Testing",name, end=' ')
        test_reduce( name, cell, expected)
        print("passed")
        

    
if __name__=="__main__":
    runtests()
