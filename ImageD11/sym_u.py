
from __future__ import print_function

## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py




import numpy as np
import logging

DEBUG = False

def m_from_string(s):
    """
    Creates a symmetry operator from a string
    """
    m = []
    t = np.array(eval("lambda x,y,z: ( %s )"%(s))(0,0,0))
    for v1,v2,v3 in [ [ 1,0,0] , [ 0,1,0], [0,0,1] ]:
        r = eval("lambda x,y,z: ( %s )"%(s))(v1,v2,v3)
        m.append(np.array(r)-t)
    return np.array(m)

def fmt(c):
    if c == 1:
        return "+"
    if c == -1:
        return "-"
    else:
        return "%f"%(c)

def m_to_string(m):
    """
    Creates a symmetry operator string from a matrix
    """
    st = []
    for i in range(3):
        for v,s in zip([ [ 1,0,0] , [ 0,1,0], [0,0,1] ],
                            "xyz"):
            c = np.dot(v,m.T[i])
            if abs(c)>0:
                st.append( "%s%s"%(fmt(c),s))
        if i<2:
            st.append(",")
    return "".join(st)







class group:
    """ An abstract mathematical finite(?) point rotation groups """
    def __init__(self, tol=1e-5):
        """
        Basic group is identity
        tol is for numerical comparison of group membership
        """
        self.group = [ np.identity(3, float) ]
        self.tol = 1e-5
    def op(self, x, y):
        """
        Normally multiplication ?
        Means of generating new thing from two others
        """
        m = np.dot(x, y)
        #
        # d = np.linalg.det(m)
        # Only appears to make sense for pure rotation matrices
        # assert abs(d-1)<1e-6, (str((d,m,x,y)))
        return m
    def comp(self, x, y):
        """
        Compare two things for equality
        """
        return np.allclose(x, y, rtol = self.tol, atol=self.tol)
    def isMember(self, item):
        """
        Decide if item is already in the group
        """
        for g in self.group:
            if self.comp(g, item):
                return True
        return False
    def additem(self, item):
        """
        add a new member
        """
        item = np.asarray(item)
        if not self.isMember(item):
            self.group.append(item)
        #else:
        #    logging.warning(str(item)+" is already a group member")
        self.makegroup()
    def makegroup(self):
        """
        ensure all items = op(x,y) are in group
        """
        #flake8: global DEBUG
        if DEBUG:
            print("making new group")
        new = True
        while new:
            for a in self.group:
                for b in self.group:
                    c = self.op(a,b)
                    new = True
                    if self.isMember(c):
                        new=False
                    if new:
                        if DEBUG: print("adding",c,"to group")
                        self.group.append(c)

symcache = {}

def generate_group(*args):
    #flake8: global symcache
    if args in symcache:
        return symcache[args]
    g=group()
    for a in args:
        g.additem(m_from_string(a))
    symcache[args]=g
    return g


# From International tables the "Laue classes" are:
#
#  -1  Triclinic
# 2/m  Monoclinic
# mmm  Orthorhombic
# 4/m  or 4/mmm Tetragonal
# -3 or -3/m Trigonal
# 6/m or 6/mmm Hexagonal
# m-3 or m-3m Cubic
#
#  Note: the Laue class is the symmetry of the reflection intensities assuming friedel's
#        law is obeyed. e.g. these are the 32 point groups combined with inversion.
#        Needed for peak intensities
#
# These groups are used for lattice permutations, they should preserve the handedness
#  and so they DO NOT CONTAIN MIRROR PLANES
#
crystallographic_point_groups = """
# index Crystal_system  Hermann-Mauguin	Order
1   Triclinic	1	1
2   Triclinic	-1	2
3   Monoclinic	2	2
4   Monoclinic	m	2
5   Monoclinic	2/m	4
6   Orthorhombic	222	4
7   Orthorhombic    mm2	4
8   Orthorhombic    mmm 8
9   Tetragonal	4   4
10  Tetragonal	-4   4
11  Tetragonal	4/m   8
12  Tetragonal	422   8
13  Tetragonal	4mm   8
14  Tetragonal	-42m   8
15  Tetragonal	4/mmm   16
16  Trigonal    3	3
17  Trigonal    -3	3
18  Trigonal    32	6
19  Trigonal    3m  6
20  Trigonal    -3m 12
21  Hexagonal	6   6
22  Hexagonal	-6	6
23  Hexagonal	6/m 12
24  Hexagonal	622 12
25  Hexagonal	6mm 12
26  Hexagonal	-62m 12
27  Hexagonal	6/mmm   24
28  Cubic	23  12
29  Cubic   m-3  24
30  Cubic   432   24
31  Cubic   -43m 24
32  Cubic   m3m  48"""

proper_point_groups = """1 2 3 4 6 222 322 422 622 23 432""".split()
assert len(proper_point_groups) == 11
improper_centrosymmetrical_point_groups = """-1 2/m -3 4/m 6/m mmm -3m 4/mmm 6/mmm m-3 m3m""".split()
assert len(improper_centrosymmetrical_point_groups) == 11
improper_acentric_point_groups = """m -4 -6 mm2 3mm 4mm -42m 6mm -62m -43m""".split()
assert len(improper_acentric_point_groups) == 10

def cubic():
    """ P32 """
    return generate_group( "z,x,y",  "-y,x,z" )

def hexagonal():
    #""" P6 168 """
    #return generate_group ( "-y,x-y,z", "-x,-y,z" )
    """ P6/mmm 191"""
    return generate_group ( "-y,x-y,z", "-x,-y,z", "y,x,-z" )

def trigonal():
    """ P321 150 """
    return generate_group ( "y,-x-y,z", "y,x,-z" )

def rhombohedralP():
    """ R3 primitive """
    return generate_group("z,x,y", "-z,-y,-x")

trigonalP = rhombohedralP

def tetragonal():
    """ P4 75"""
    return generate_group ( "-y,x,z", "-x,y,-z" )

def orthorhombic():
    """ P222 16 """
    return generate_group( "-x,-y,z", "-x,y,-z" )

def monoclinic_c():
    return generate_group("-x,-y,z" )

def monoclinic_a():
    return generate_group("x,-y,-z" )

def monoclinic_b():
    return generate_group("-x,y,-z" )

def triclinic():
    return generate_group("x, y, z" )


def find_uniq_u(u, grp, debug=0, func=np.trace):
    uniq = u
    tmax = func(uniq)
    for o in grp.group:
        cand = grp.op(o, u)
        t = func(cand)
        if debug: print(t)
        if func(cand) > tmax:
            uniq = cand
            tmax = t
    return np.array(uniq)


def hklmax(h, hmax=1000):
    # Assumes |h| < hmax
    return (h[0]*hmax + h[1])*hmax + h[2]

def find_uniq_hkls( hkls, grp, func=hklmax):
    assert hkls.shape[0] == 3, 'hkls must be 3xn array'
    uniq = hkls.copy()
    tmax = func( hkls )
    for o in grp.group:
        cand = grp.op( o, hkls)
        t = func(cand)
        msk = t > tmax
        for i in range(3):
            uniq[i] = np.where( msk , cand[i], uniq[i] )
        tmax = np.where( msk , t, tmax )
    return uniq


class trans_group(group):
    """
    Translation group (eg crystal lattice)

    FIXME - this is mostly done in lattice_reduction.py instead now
    """
    def __init__(self, tol = 1e-5):
        """
        Identity is to not move at all
        """
        group.__init__(self,tol)
        
    def op(self, x, y):
        """
        Means of generating new thing from two others
        In this case add them and mod by group members
        """
        return self.reduce(x + y)
    def reduce(self, v):
        """
        Perform lattice reduction
        """
        vc = np.array(v).copy() # copies
        for o in self.group:
            vc = self.mod(vc, o)
        # if DEBUG: print "reduced",v,vc
        return vc
    def additem(self, x):
        """ Do lattice reduction before adding as infinite group"""
        t = self.reduce(x)
        group.additem(self, self.reduce(x))
        # Now try to remove anything which is spare??
        return
    def mod(self, x, y):
        """
        Remove y from x to give smallest possible result
        Find component of x || to y and remove it
        """
        ly2 = np.dot(y, y)
        if ly2 > 1e-9:
            ny = np.dot(x,y)/ly2
            parl = ny * y
            ints = np.round_(ny)
            return x - ints * y
        else:
            return x
    def isMember(self, x):
        return group.isMember(self, self.reduce(x))

def test():
    assert np.allclose( m_from_string( "x,y,z" ), np.identity(3))
    assert np.allclose( m_from_string( "-y,x,z" ), np.array([ [ 0,1,0],
                                                          [-1,0,0],
                                                          [ 0,0,1]] ))
    assert np.allclose( m_from_string( "-y,y-x,z" ), np.array([[ 0,-1,0],
                                                             [ -1, 1,0],
                                                             [ 0, 0,1]] ))
    print("testing1")
    for op in [ "x,y,z", "-y,x-y,z", "-y,x,z"]:
        d = np.linalg.det(m_from_string(op))
        assert d == 1.0, "Determinant = %f %s"%(d,op)
    print("testing2")
    assert len(cubic().group) == 24, "not 24 ops found for cubic !"
    assert len(hexagonal().group) == 12 ,"not 12 ops found for hexagonal !"
    assert len(trigonal().group) == 6 ,"not 7 ops found for trigonal !"+\
        str(trigonal().group)
    assert len(tetragonal().group) == 8 ,"not 8 ops found for tetragonal !"
    assert len(orthorhombic().group) == 4 ,"not 4 ops found for orthorhombic !"
    print("testing3")
    for f in [ monoclinic_a, monoclinic_b, monoclinic_c]:
        r = f().group
        assert len(r) == 2, " not 2 ops for monoclinic "

    assert np.allclose(
        find_uniq_u( np.array(
                [[0,1,0],[-1,0,0],[0,0,1]]),cubic()),
                     np.identity(3) ), "Should easily get this unique choice"

    # translational groups
    g1 = trans_group()
    g2 = trans_group()
    ops = [ np.array( [ 1,0,0], float) ,
            np.array( [ 0,1,0], float) ,
            np.array( [ 0,0,1], float) ]
    for op in ops:
        g1.additem(op)
        g2.additem(op)
    g2.additem( np.array( [ 5,6,7], float)  )
    for op2 in g2.group:
        found = False
        for op1 in g1.group:
            if ( op1 == op2).all():
                found = True
        if not found:
            raise Exception ("Bad translation groups")
    assert not g2.isMember([0.1,0.5,10]), "Not a member"
    assert g2.isMember([99,-1e5,4e7]), "Not a member"
    global DEBUG
    DEBUG = True
    g2.additem([0.1, 0.45, 10])
    print(g2.group)
    DEBUG = False


def getgroup(s):
    """
    convert a user supplied string to a group
    ... a little vague still
    """
    groups =['cubic', 'hexagonal','trigonal','tetragonal',
             'orthorhombic','monoclinic_c','monoclinic_a',
             'monoclinic_b','triclinic','rhombohedralP'] 
    if s in groups:
        return globals()[s]
    else:
        raise Exception( s,"not in", groups )



if __name__=="__main__":
    test()

    u = np.array([[ 0.71850787 , 0.69517833,  0.02176059],
                 [-0.62925889 , 0.66306714, -0.40543213],
                 [-0.29627636 , 0.27761313 , 0.91386611]])

    find_uniq_u(u,cubic(),debug=0)




