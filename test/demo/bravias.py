

"""
Determine Bravias lattice for UBI

Return transformation matrix


"""


from ImageD11.unitcell import unitcell
from ImageD11.indexing import ubitocellpars
from numpy import array


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
    cell = unitcell( cp )
    dlim = cell.ds( array( [7,7,7] ))
    u.makerings( dlim , tol)
    mult = [ len(cell.ringhkls[d]) for d in cell.ringds ]
    mx = max(mult)
    mn = min(mult)
    cases = {
        (2,  2) : fi_triclinic,
        (4,  2) : fi_monoclinic,
        (8,  2) : fi_orthorhombic,
        (12, 2) : fi_rhombohedral,
        (16, 2) : fi_tetragonal,
        (24, 2) : fi_hexagonal,
        (48, 2) : fi_cubic }
        

def katayama(ubi):
    """
    J. Appl. Cryst (1986) 19 69-72
    
    Take four vectors a,b,c and d=-(a+b+c)
    Compute the 6 pairwise dot products (ab, ac, ad, bc, bd, cd)
    Apply transformations (up to 12) and recompute
    
    etc
    """
    pass

