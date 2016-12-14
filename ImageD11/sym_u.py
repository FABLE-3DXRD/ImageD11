
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
        needplus = 0
        for v,s in zip([ [ 1,0,0] , [ 0,1,0], [0,0,1] ],
                            "xyz"):
            c = np.dot(v,m.T[i])
            if abs(c)>0:
                st.append( "%s%s"%(fmt(c),s))
                needplus = 1
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
        self.group = [ np.identity(3, np.float) ]
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
        global DEBUG
        if DEBUG:
            print "making new group"
        new = True            
        while new:
            for a in self.group:
                for b in self.group:
                    c = self.op(a,b)
                    new = True
                    if self.isMember(c):
                        new=False
                    if new:
                        if DEBUG: print "adding",c,"to group"
                        self.group.append(c)

def generate_group(*args):
    g=group()
    for a in args:
        g.additem(m_from_string(a))
    return g



def cubic():
    return generate_group( "z,x,y",  "-y,x,z" )

def hexagonal():
    #""" P6 168 """
    #return generate_group ( "-y,x-y,z", "-x,-y,z" )
    """ P6/mmm 191"""
    return generate_group ( "-y,x-y,z", "-x,-y,z", "y,x,-z" )

def trigonal():
    """ P3 143 """
    return generate_group ( "y,-x-y,z" )

def rhombohedralP():
    """ R3 primitive """
    return generate_group("z,x,y", "-z,-y,-x")

def tetragonal():
    """ P4 75"""
    return generate_group ( "-y,x,z", "-x,-y,z" )
    
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
    return generate_group("-x,-y,-z" )


def find_uniq_u(u, grp, debug=0, func=np.trace):
    uniq = u
    tmax = func(uniq)
    for o in grp.group:
        cand = grp.op(o, u)
        t = func(cand)
        if debug: print t
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
        self.group = [ np.zeros(3, np.float) ]
        self.tol = tol
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
    print "testing1"
    for op in [ "x,y,z", "-y,x-y,z", "-y,x,z"]:
        d = np.linalg.det(m_from_string(op)) 
        assert d == 1.0, "Determinant = %f %s"%(d,op)
    print "testing2"    
    assert len(cubic().group) == 24, "not 24 ops found for cubic !"
    print  len(hexagonal().group)
    assert len(hexagonal().group) == 12 ,"not 6 ops found for hexagonal !"
    assert len(trigonal().group) == 3 ,"not 3 ops found for trigonal !"+\
        str(trigonal().group)
    assert len(tetragonal().group) == 4 ,"not 8 ops found for tetragonal !"
    assert len(orthorhombic().group) == 4 ,"not 4 ops found for orthorhombic !"
    print "testing3"
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
    ops = [ np.array( [ 1,0,0], np.float) ,
            np.array( [ 0,1,0], np.float) ,
            np.array( [ 0,0,1], np.float) ]
    for op in ops:
        g1.additem(op)
        g2.additem(op)
    g2.additem( np.array( [ 5,6,7], np.float)  )
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
    print g2.group
    DEBUG = False
    

def getgroup(s):
    """
    convert a user supplied string to a group
    ... a little vague still
    """
    if s in ['cubic', 'hexagonal','trigonal','tetragonal',
             'orthorhombic','monoclinic_c','monoclinic_a',
             'monoclinic_b','triclinic','rhombohedralP']:
        import ImageD11.sym_u
        return getattr(ImageD11.sym_u, s)

    

if __name__=="__main__":
    test()
    

    u = np.array([[ 0.71850787 , 0.69517833,  0.02176059],
                 [-0.62925889 , 0.66306714, -0.40543213],
                 [-0.29627636 , 0.27761313 , 0.91386611]])

    find_uniq_u(u,cubic(),debug=0)

    

        
