## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py



from Numeric import *
from numpy.oldnumeric.linear_algebra import determinant

def m_from_string(s):
    m = []
    t = array(eval("lambda x,y,z: ( %s )"%(s))(0,0,0))
    for v1,v2,v3 in [ [ 1,0,0] , [ 0,1,0], [0,0,1] ]:
        r = eval("lambda x,y,z: ( %s )"%(s))(v1,v2,v3)
        m.append(array(r)-t)
    return array(m)



class group:
    def __init__(self):
        """
        Basic group is identity
        """
        self.group = [ identity(3,Float) ]
    def op(self, x, y):
        """
        Means of generating new thing from two others
        """
        m = matrixmultiply(x,y)
        d = determinant(m)
        assert abs(d-1)<1e-6, (str((d,m,x,y)))
        return m
    def comp(self,x,y):
        """
        Compare two things for equality
        """
        return allclose(x,y)
    def isMember(self,item):
        """
        Decide if item is already in the group
        """
        for g in self.group:
            if self.comp(g,item):
                return True
    def additem(self,item):
        """
        add a new member
        """
        if not self.isMember(item):
            self.group.append(item)
        self.makegroup()
    def makegroup(self):
        """
        ensure all items = op(x,y) are in group
        """
        new = True
        while new:
            for a in self.group:
                for b in self.group:
                    c = self.op(a,b)
                    new = True
                    if self.isMember(c):
                        new=False
                    if new:
                        self.group.append(c)
    
def generate_group(*args):
    g=group()
    for a in args:
        g.additem(m_from_string(a))
    return g

def cubic():
    return generate_group( "z,x,y",  "-y,x,z" )

def hexagonal():
    return generate_group ( "x+y,-x,z", "x,-x-y,-z")

def trigonal():
    return generate_group ("z,x,y","-y,-x,-z")

def tetragonal():
    return generate_group ("-y,x,z", "x,-y,-z")
    
def orthorhombic():
    return generate_group("-x,-y,z","x,-y,-z")

def monoclinic_c():
    return generate_group("-x,-y,z")

def monoclinic_a():
    return generate_group("x,-y,-z")

def monoclinic_b():
    return generate_group("-x,y,-z")

def triclinic():
    return generate_group("-x,-y,-z")


def find_uniq_u(u,grp,debug=0):
    uniq = u
    tmax = trace(uniq)
    for o in grp.group:
        cand = grp.op(o,u)
        t = trace(cand)
        if debug: print t
        if trace(cand) > tmax:
            uniq = cand
            tmax = t
    return array(uniq)


def test():
    assert allclose( m_from_string( "x,y,z" ), identity(3))
    assert allclose( m_from_string( "-y,x,z" ), array([ [ 0,1,0],
                                                        [-1,0,0],
                                                        [ 0,0,1]] ))
    assert allclose( m_from_string( "-y,y-x,z" ), array([[ 0,-1,0],
                                                         [ -1, 1,0],
                                                         [ 0, 0,1]] ))
    
    for op in [ "x,y,z", "-y,x-y,z", "-y,x,z"]:
        d = determinant(m_from_string(op)) 
        assert d == 1.0, "Determinant = %f %s"%(d,op)
    
    assert len(cubic().group) == 24, "not 24 ops found for cubic !"
    assert len(hexagonal().group) == 12 ,"not 12 ops found for hexagonal !"
    assert len(trigonal().group) == 6 ,"not 6 ops found for hexagonal !"
    assert len(tetragonal().group) == 8 ,"not 8 ops found for tetragonal !"
    assert len(orthorhombic().group) == 4 ,"not 4 ops found for orthorhombic !"
    for f in [ monoclinic_a, monoclinic_b, monoclinic_c]:
        r = f().group
        assert len(r) == 2, " not 2 ops for monoclinic "

    assert allclose( find_uniq_u( array( [[0,1,0],[-1,0,0],[0,0,1]]),cubic()),
                     identity(3) ), "Should easily get this unique choice"
    

if __name__=="__main__":
    test()
    

    u = array([[ 0.71850787 , 0.69517833,  0.02176059],
               [-0.62925889 , 0.66306714, -0.40543213],
               [-0.29627636 , 0.27761313 , 0.91386611]])

    find_uniq_u(u,cubic(),debug=1)

    

        
