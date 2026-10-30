
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

    @property
    def _group_stack(self):
        if getattr(self, "_gs", None) is None or len(self._gs) != len(self.group):
            self._gs = np.array(self.group)
        return self._gs

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
    return generate_group ( "-y,x-y,z", "y,x,-z" )

trigonalH = trigonal

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

# Generators of the rotations of each Laue class, in the ITA standard setting.
# A Laue class is these rotations and their products with the inversion.
# Rotations only: they act on orientations (op.ubi stays right handed).
laue_generators = {
    "-1"    : ("x, y, z",),                          # triclinic()
    "2/m"   : ("-x,y,-z",),                          # monoclinic_b(), unique axis b
    "mmm"   : ("-x,-y,z", "-x,y,-z"),                # orthorhombic()
    "4/m"   : ("-y,x,z",),
    "4/mmm" : ("-y,x,z", "-x,y,-z"),                 # tetragonal()
    "-3"    : ("-y,x-y,z",),                         # hexagonal axes
    "-3m1"  : ("-y,x-y,z", "y,x,-z"),                # trigonal(), 2-folds along a
    "-31m"  : ("-y,x-y,z", "-y,-x,-z"),              # 2-folds along a*
    "6/m"   : ("-y,x-y,z", "-x,-y,z"),
    "6/mmm" : ("-y,x-y,z", "-x,-y,z", "y,x,-z"),     # hexagonal()
    "m-3"   : ("z,x,y", "-x,-y,z"),
    "m-3m"  : ("z,x,y", "-y,x,z"),                   # cubic()
    "-3R"   : ("z,x,y",),                            # rhombohedral axes
    "-3mR"  : ("z,x,y", "-z,-y,-x"),                 # rhombohedralP()
}

def laue_class(spacegroup, cell_choice="standard"):
    """
    The Laue class of a space group (number 1-230), a key of laue_generators.
    cell_choice "rhombohedral" puts R space groups on rhombohedral axes.
    """
    from xfab import sg
    spg = sg.sg(sgno=spacegroup)
    if cell_choice == "rhombohedral" and spg.name.startswith("R"):
        return {"-3": "-3R", "-3m1": "-3mR"}[spg.Laue]
    return spg.Laue

def laue_rotations(spacegroup, cell_choice="standard"):
    """
    Rotations of the Laue class of a space group (number 1-230).
    These are the orientations diffraction cannot tell apart (Friedel's law).
    cell_choice is "standard" or "rhombohedral" (R groups on rhombohedral axes)
    Returns a group acting on ubi as op.ubi, like cubic() etc
    """
    return generate_group(*laue_generators[laue_class(spacegroup, cell_choice)])

def cartesian(grp, B):
    """
    The operations of grp as rotations S of the crystal frame of B:
    for ubi = (U.B)^-1, (op.ubi)^-1 = U.S.B
    Returns an array (n, 3, 3) in the order of grp.group
    """
    Bi = np.linalg.inv(B)
    return np.array([np.dot(B, np.dot(np.linalg.inv(o), Bi)) for o in grp.group])

def _axis_order(r):
    """ Unit axis and order of a proper rotation, (None, 1) for the identity """
    angle = np.arccos(np.clip((np.trace(r) - 1) / 2, -1, 1))
    if angle < 1e-6:
        return None, 1
    w, v = np.linalg.eig(r)
    axis = np.real(v[:, np.argmin(np.abs(w - 1))])
    return axis / np.linalg.norm(axis), int(round(2 * np.pi / angle))

def _same_set(a, b, tol=1e-6):
    return len(a) == len(b) and all(min(np.abs(x - y).max() for y in b) < tol for x in a)

def _frame(u, v):
    """ Orthonormal frame (columns) with u first and v in the first plane """
    e1 = u / np.linalg.norm(u)
    e2 = v - np.dot(v, e1) * e1
    e2 = e2 / np.linalg.norm(e2)
    return np.column_stack((e1, e2, np.cross(e1, e2)))

def align(rotations, target):
    """
    A rotation C with { C.r.C^T for r in rotations } == target, for two
    lists of Cartesian rotations (e.g. from cartesian). The identity when
    they are already the same. Raises ValueError if there is no such C.
    """
    if _same_set(rotations, target):
        return np.eye(3)
    axr = [ak for ak in map(_axis_order, rotations) if ak[0] is not None]
    axt = [ak for ak in map(_axis_order, target) if ak[0] is not None]
    if len(axr) == len(axt) > 0:
        # take the highest order axis of target to one of rotations ...
        n = max(k for a, k in axt)
        ut = [a for a, k in axt if k == n][0]
        # ... and a second axis of target, if there is one, to another
        vts = [a for a, k in axt if abs(abs(np.dot(a, ut)) - 1) > 1e-6]
        for ur in [s * a for a, k in axr if k == n for s in (1, -1)]:
            if vts:
                vt = vts[0]
                vrs = [s * a for a, k in axr for s in (1, -1)
                       if abs(np.dot(ur, s * a) - np.dot(ut, vt)) < 1e-6]
            else:   # one axis : any perpendicular will do
                vt = np.eye(3)[np.argmin(np.abs(ut))]
                vrs = [np.eye(3)[np.argmin(np.abs(ur))]]
            for vr in vrs:
                C = np.dot(_frame(ut, vt), _frame(ur, vr).T)
                if _same_set([np.dot(np.dot(C, r), C.T) for r in rotations], target):
                    return C
    raise ValueError("No rotation takes these rotations to the target")

def find_uniq_u(u, grp, debug=0, func=np.trace):
    if func is np.trace and not debug:
        # trace(o.u) = sum_ij o_ij u_ji, so all the candidate traces at once
        ops = grp._group_stack                      # (n, 3, 3), cached below
        traces = ops.reshape(len(ops), 9).dot(np.asarray(u).T.ravel())
        k = int(traces.argmax())
        if traces[k] > np.trace(u):
            return np.array(grp.op(ops[k], u))
        return np.array(u)
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
            ints = np.round(ny)
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
             'monoclinic_b','triclinic','rhombohedralP',
             'trigonalH','trigonalP']
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




