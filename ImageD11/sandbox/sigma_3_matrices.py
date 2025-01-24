from __future__ import print_function


import ImageD11.grain
import ImageD11.grid_index_parallel
import scipy.spatial.transform
import numpy as np


def make_sigma3_mats( ):
    """ Rotations of 60 degrees on 111 """
    v = 60/np.sqrt( 3 )
    S3_matrices = [ scipy.spatial.transform.Rotation.from_rotvec( (x,y,z), degrees=True ).as_matrix()
            for x in (-v,v)
            for y in (-v,v)
            for z in (-v,) ]  # only need 1 z as -60 == 60
    return S3_matrices


def applytwinmat( g, mat ):
    """ U.B -> U.M'.B """
    ubinew = np.linalg.inv( g.U.dot( mat.T ).dot( g.B ) )
    return ImageD11.grain.grain( ubinew, g.translation )


def twin_check( ind, matrices, symmetry='cubic', poserr=1000, angerr=0.5 ):
    gl = ImageD11.grid_index_parallel.uniq_grain_list( symmetry, poserr, angerr )
    gl.add( [ ImageD11.grain.grain( ubi, np.array((0,0,0)) ) for ubi in ind.ubis ] )
    for ubi in ind.ubis:
        g = ImageD11.grain.grain( ubi )
        sg1 = ImageD11.cImageD11.score( ubi, ind.gvflat, ind.hkl_tol )
        print("%4d   "%(sg1), end='')
        for mat in matrices:
            gnew = applytwinmat( g, mat )
            s1 = ImageD11.cImageD11.score( gnew.ubi, ind.gvflat, ind.hkl_tol )
            if s1 > ind.minpks:  # should be accepted assuming uniq
                b4 = len(gl.uniqgrains)
            gl.add( [gnew,] )
            if len(gl.uniqgrains) == b4:
                e = 'r'
            else:
                e = '!'
            print("%4d "%(s1), end=e)
        print() 
