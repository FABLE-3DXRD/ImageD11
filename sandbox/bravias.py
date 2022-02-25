
from __future__ import print_function, division
"""
Determine Bravias lattice for UBI

Return transformation matrix

Due to:
ferraris_and_ivaldi

    Acta Cryst (1983) A39 595-596

    Determine the multiplicities of vectors in a computed powder pattern
    Match these to multiplicities according to lattice types
    ... appeals to powder person
    Tolerance is for grouping vector lengths.
    

"""

import sys
import numpy as np, pylab as pl


from ImageD11.unitcell import unitcell
from ImageD11.grain import read_grain_file
from ImageD11.indexing import ubitocellpars



def fi_triclinic(ubi):
    """ Ensure we have the 3 shortest non-coplanar vectors """
    return ubi

def fi_monoclinic(ubi, cell, mult):
    """ 
    Choose vectors having mult == 2
    These are either along b or in ac plane
    If centred then one mult 4 vector is shorter and cell par
    """
    possibles = []
    for i in range(len(mult)):
        if mult[i] == 2:
            possibles.append( cell.ringhkls(cell.ringds[i])[0] )
    possibles = [cell.ring]
    

def ferraris_and_ivaldi(ubi, tol = 0.001):
    """
    Acta Cryst (1983) A39 595-596

    Determine the multiplicities of vectors in a computed powder pattern
    Match these to multiplicities according to lattice types
    ... appeals to powder person
    Tolerance is for grouping vector lengths.
    """
    cpar = ubitocellpars(ubi)
    cell = unitcell( cpar )
    dlim = cell.ds( np.array( [7,7,7] ))
    cell.makerings( dlim , tol)
    mult = [ len(cell.ringhkls[d]) for d in cell.ringds ]
    mx = max(mult)
    mn = min(mult)
    cases = {
        (2,  2) : fi_triclinic,
        (4,  2) : fi_monoclinic,
#        (8,  2) : fi_orthorhombic,
#        (12, 2) : fi_rhombohedral,
#        (16, 2) : fi_tetragonal,
#        (24, 2) : fi_hexagonal,
#        (48, 2) : fi_cubic
    }

def vecang( v1, v2 ):
    cosa = np.dot(v1,v2)/np.sqrt( np.dot( v1, v1 ) * np.dot( v2, v2 ))
    return np.degrees( np.arccos( cosa ))

def lattice_powder( ubi, tol=0.01):
    a,b,c = np.asarray( ubi )
    nmax = 5
    vecs = []
    lens = []
    for ai in range(-nmax,nmax):
        for bi in range(-nmax,nmax):
            for ci in range(-nmax,nmax):
                vec = ai*a + bi*b + ci*c
                modv = np.sqrt( np.dot(vec, vec) )
                lens.append( modv )
                vecs.append( vec )
    order = np.argsort( lens )
    lens  = np.take( lens, order )
    vecs  = np.take( vecs, order ,axis=0 )
    steps = (lens[1:]-lens[:-1])/lens[1:]   # diff divided by upper
    gaps  = steps > tol
    rings = np.cumsum(gaps)
    mults, bins = np.histogram( rings, bins=np.arange(-0.5, rings[-1]+0.5, 1 ) )
    print(rings)
    print(mults)
    nr = min( 10, rings.max() )
    for i in range( nr, 0 , -1 ):
        # this ring:
        M = mults[i]
        print("Powder ring",i,"multiplicity", M )
        print("   Angles between vectors:")
        vl = vecs[ 1: ][ rings==i ]
        vll = lens[ 1: ][ rings==i ]
        for j in range(M):
            print ("%2d:"%(j),vl[j],vll[j])
        print("  ",end="")
        for j in range(M-1,0,-1):
            print("%7d"%(j),end=" ")
        print()
        for j in range(M-1):
            print(" %2d: "%(j), end="")
            for k in range(M-1, j, -1):
                print( "%7.3f"%( vecang( vl[j],vl[k] ) ) , end=" ")
            print()
    print( "Cell pars before",ubitocellpars(ubi))

        
def katayama(ubi):
    """
    J. Appl. Cryst (1986) 19 69-72
    
    Take four vectors a,b,c and d=-(a+b+c)
    Compute the 6 pairwise dot products (ab, ac, ad, bc, bd, cd)
    Apply transformations (up to 12) and recompute
    
    etc
    """
    pass


if __name__=="__main__":
    for grain in read_grain_file( sys.argv[1] ):
        print ("ubi:\n",grain.ubi)
        lattice_powder(grain.ubi)
