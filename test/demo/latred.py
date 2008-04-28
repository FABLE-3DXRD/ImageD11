
from numpy import dot, round_, array, float, allclose, asarray, fabs,\
        argmax, sqrt
from numpy.linalg import inv

def fparl(x,y):
    """ fraction of x parallel to y"""
    ly2=dot(y,y)
    if ly2 > 1e-9:
        return dot(x,y)/ly2
    return 0

def mod(x,y):
    n = round_(fparl(x,y))
    return x - n*y

MIN_VEC2 = 1e-9*1e-9

class lattice(object):
    def __init__(self, v1, v2, v3):
        """ Make the lattice 
        Currently v1, v2, v3 are reciprocal space vectors
        TODO : apply to pattersons      
        """
        vl=[v1, v2, v3]
        for i in range(3):
            assert dot(vl[i], vl[i]) > MIN_VEC2, vl[i]
            for j in range(3):
                if i==j: continue
                rem = mod(vl[i], vl[j])
                if dot(rem, rem) > MIN_VEC2:
                    vl[i] = rem
        self.v = array(vl)
        self.vi = inv(self.v.T)
    def hkls(self, x):
        """ Give back the hkl indices """
        #f1 = [ fparl( x , v ) for v in self.v ]
        #h1 = round_(f1)
        f2 = dot(self.vi, x)
        h2 = round_(f2)
        # Should fail for non-orthogonal lattice
        #assert (h2 == h1).all(), [str(z) for z in [x, h1, f1, h2, f2 ] ]
        return h2
    def rem(self, x):
        """ Remainder of x after removing closest lattice point """
        h = self.hkls(x)
        r1 = x - dot(self.v.T, h)
        #r2 = x.copy()
        #for v in self.v:
        #    r2 = mod(r2, v)
        # Should fail for non-orthogonal lattice
        #assert allclose(r1, r2), [str(z) for z in [ x, h, r1, r2] ]
        return r1
    def withvec(self, x):
        """ 
        Try to fit x into the lattice 
        Make the remainder together with current vectors
        Index it as hkl indices
        whichever vector has the biggest projection is replaced
        remake the lattice with these 3 vectors
        """
        r = self.rem(x)
        vd = argmax(fabs(dot(self.vi, x)))
        v = list(self.v)
        v[vd] = r
        return lattice( v[0], v[1], v[2] )
        
     

def test_eu():
    from ImageD11.indexing import indexer
    #0 = (0, 1,  1)
    #1 = (0, 1, -1)
    #6 = (1, 0, -1)
    o = indexer()
    o.readgvfile("eu3.gve")
    v1 = o.gv[0]
    v2 = o.gv[1]
    v3 = o.gv[6]
    print v1, v2, v3
    l = lattice ( v1, v2, v3)
    print l.v
    print l.vi
    for v in o.gv:
        print ("%8.5f "*3+"%4d "*3)%tuple( list(v)+list(l.hkls(v))),
        err = l.rem(v)
        print "%.5f"%(sqrt(dot(err,err)))
        
def test2():
    """ Adding a fourth vector """
    o = lattice( array( [ 1, 0, 0] , float) ,
                 array( [ 0, 1, 0] , float) ,
                 array( [ 0, 0, 2] , float) )
    # ... how to do it?
    r  = o.rem( array( [ 0, 0, 1] , float) ) 
    assert ( r == array( [ 0, 0, 1] , float) ).all()
    o = o.withvec( r )
    assert allclose(o.v[0], array( [ 1, 0, 0] , float) )
    assert allclose(o.v[1], array( [ 0, 1, 0] , float) )
    assert allclose(o.v[2], array( [ 0, 0, 1] , float) )
    
def test1():
    """ Make a lattice from 3 vectors"""
    o = lattice( array( [ 1, 0, 0] , float) ,
                 array( [ 2, 1, 0] , float) ,
                 array( [ 3, 4, 1] , float) )
    assert allclose(o.v[0], array( [ 1, 0, 0] , float) )
    assert allclose(o.v[1], array( [ 0, 1, 0] , float) )
    assert allclose(o.v[2], array( [ 0, 0, 1] , float) )
    """ Assign hkl indices in a strange way (no matrix inversing) """
    for x in [ [ 1,3,4 ], [ -9, 10, 13] ]:
        h = o.hkls(x)
        assert  (asarray(h) == asarray(x)).all(), x
    """ 3 more difficult vectors """
    o = lattice( array( [ 1, 0, 1] , float) ,
                 array( [ 0, 1, 0] , float) ,
                 array( [ 0, 0, 1] , float) )
    assert allclose(o.v[0], array( [ 1, 0, 0] , float) )
    assert allclose(o.v[1], array( [ 0, 1, 0] , float) )
    assert allclose(o.v[2], array( [ 0, 0, 1] , float) )

        

if __name__=="__main__":
    test1()
    test2()
    test_eu()
