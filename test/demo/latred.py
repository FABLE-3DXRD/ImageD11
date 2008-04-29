
from numpy import dot, round_, array, float, allclose, asarray, fabs,\
        argmax, sqrt, argsort, take
from numpy.linalg import inv

DEBUG = False

# Some sort of round off
MIN_VEC2 = 1e-9*1e-9


def fparl(x, y):
    """ fraction of x parallel to y """
    ly2 = dot(y ,y)
    if ly2 > 1e-9:
        return dot(x, y) / ly2
    return 0

def mod(x,y):
    """ 
    Returns the part of x with integer multiples of y removed
    assert that ||mod(x,y)|| <= ||x|| for all y 
    """
    if __debug__: 
        b4  = dot(x,x)
    n   = round_(fparl(x,y))
    ret = x - n * y
    if __debug__:
        af  = dot(ret,ret)
        assert b4 >= af , "Bad mod "+str(x) + " " + str(y)
    return ret


def sortvec_len( vl ):
    """
    Sorts according to length (d-s-u)
    """
    ul = [ (dot(v,v),tuple(v)) for v in vl ]
    # Choose x,y,z
    ul.sort()
    return [ asarray(v[1]) for v in ul ]

def sortvec_xyz( vl ):
    """
    For each choose between x, -x
    Sorts according to x then y then z components
    """
    ul = []
    # Choose x/-x
    for v in vl:
        # Wonder how to do this with numpy?
        ul.append( max( tuple(v),tuple(-asarray(v)) ) )
    # Choose x,y,z
    ul.sort()
    return [ asarray(v) for v in ul[::-1] ]



class lattice(object):
    """
    Represents a 3D crystal lattice built from 3 vectors
    """

    def __init__(self, v1, v2, v3, space='reciprocal', min_vec2=MIN_VEC2):
        """ Make the lattice 
        Currently v1, v2, v3 are vectors - which gives 3D
        space [ 'real' | 'reciprocal' ]
        ... means they are real space lattice vectors or measured scattering
            vectors from crystallography
        It will attempt to find a primitive basis from the supplied vectors
        """
        # 3 supplied vectors (these are lattice vectors)
        vl=[v1, v2, v3] 
        # 9 difference vectors (these are also lattice vectors)
        # Assuming : We can always add or subtract integer amounts of any
        #            lattice vector to get another lattice vector.
        #            We look for the 3 shortest which are non-coplanar
        #
        # We could do this recursively adding and subtracting
        # ... We *need* to do it only until the solution is stable 
        # one would hope this is enough, but more testcases should be added?
        # ... currently the algorith is vague
        dsu = sortvec_len( vl + [  vi - vj for vi in vl for vj in vl ] +
                                [  vi + vj for vi in vl for vj in vl ]  )
        # Sort by length
        # Take shortest as first vector
        vl = []
        for i in range(len(dsu)):
            t = dsu[i]
            for v in vl: # remove previous vectors
                t = mod( t, v )
            if len(vl) == 2: # And try to catch hexagonal planes too
                t = mod( t, vl[0] + vl[1] )
                t = mod( t, vl[0] - vl[1] )
            if dot(t,t) < min_vec2: # Skip if zeroed out
                continue
            vl.append( t )
            if len(vl) == 3: break
            # Remove this from future vectors ???
            dsu = [ mod( d , t ) for d in dsu]
        assert len(vl) == 3, "Have not got 3 vectors"
        # print "vl:",vl
        vl = sortvec_xyz(vl)
        # print "vl",vl
        if space == 'reciprocal':  
            self.v  = array(vl).T   # We got column vectors
            self.vi = inv(self.v)
        elif space == 'real':
            self.vi = array(vl)     # We got row vectors
            self.v  = inv(self.vi) 
        else:
            raise Exception("Space must be real or reciprocal "+str(space))
    def ijks(self, u):
        """
        Give back the nearest lattice point indices
        ... only present for symmetry with hkls
        """
        return round_(dot(self.v, u))
    def rem_ijks(self, u):
        """
        Remainder of u (distance from a lattice point)
        """
        return g - dot(self.vi, self.hkls(g))
    def hkls(self, g):
        """ 
        Give back the hkl indices 
        ... these are the nearest lattice point
        """
        return round_(dot(self.vi, g))
    def rem_hkls(self, g):
        """ 
        Remainder of x after removing closest lattice point in reciprocal space
        """
        return g - dot(self.v, self.hkls(g))
    def withvec(self, x, space="reciprocal"):
        """ 
        Try to fit x into the lattice 
        Make the remainder together with current vectors
        Index it as hkl indices
        whichever vector has the biggest projection is replaced
        remake the lattice with these 3 vectors
        """
        x = asarray(x)
        if space == 'reciprocal':
            r = self.rem_hkls(x)
            vd = argmax(fabs(dot(self.vi, x)))
            v = list(self.v.T)  # reciprocal space is columns
            v[vd] = r
        if space == 'real':
            r = self.rem_ijks(x)
            vd = argmax(fabs(dot(self.vi, x)))
            v = list(self.vi)
            v[vd] = r
        return lattice( v[0], v[1], v[2] , space=space )
        

def get_eu_gv():
    from ImageD11.indexing import indexer
    o = indexer()
    o.readgvfile("eu3.gve" , quiet = True)
    return o.gv

def test_fft():
    gv = get_eu_gv()
    from ImageD11.fft_index_refac import grid
    g = grid( np = 128,
              mr = 1.0,
              nsig = 20,
              minlen = 3. )
    g.gv_to_grid_new(gv)
    g.fft()
    g.props()
    g.peaksearch(open("eu.patterson_pks","w"))
    g.read_peaks("eu.patterson_pks")
    vecs = g.UBIALL
    order = argsort( g.colfile.sum_intensity )[::-1]
    try:
        for i in order:
            print vecs[i], g.colfile.sum_intensity[i]
            for j in order[i:]:
                for k in order[j,:]:
                    try:
                        l = lattice( vec[i], vecs[j], vecs[k] )
                    except:
                        pass
                    raise Exception("Really bad control structure")
    except:
        print "Got something"

def test_eu():
    #Conventional cell gives:
    #0 = (0, 1,  1)
    #1 = (0, 1, -1)
    #6 = (1, 0, -1)
    gv = get_eu_gv()
    v1 = gv[0]
    v2 = gv[1]
    v3 = gv[6]
    l = lattice ( v1, v2, v3)
    esum = 0.0
    for v in gv:
        # print ("%8.5f "*3+"%4d "*3)%tuple( list(v)+list(l.hkls(v))),
        err   = l.rem_hkls(v)
        esum += sqrt(dot(err,err))
        # print "%.5f"%(sqrt(dot(err,err)))
    # print esum / len(o.gv)
    assert esum / len(gv) < 0.0102, "Did not fit"

def test2():
    """ Adding a fourth vector """
    o = lattice( array( [ 1, 0, 0] , float) ,
                 array( [ 0, 1, 0] , float) ,
                 array( [ 0, 0, 2] , float) )
    # ... how to do it?
    r  = o.rem_hkls( array( [ 0, 0, 1] , float) ) 
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
    assert allclose(o.v[:,0], array( [ 1, 0, 0] , float) )
    assert allclose(o.v[:,1], array( [ 0, 1, 0] , float) )
    assert allclose(o.v[:,2], array( [ 0, 0, 1] , float) )
    """ Assign hkl indices in a strange way (no matrix inversing) """
    for x in [ [ 1,3,4 ], [ -9, 10, 13] ]:
        h = o.hkls(x)
        assert  (asarray(h) == asarray(x)).all(), x
    """ 3 more difficult vectors """
    o = lattice( array( [ 1, 0, 1] , float) ,
                 array( [ 0, 1, 0] , float) ,
                 array( [ 0, 0, 1] , float) )
    assert allclose(o.v[:,0], array( [ 1, 0, 0] , float) )
    assert allclose(o.v[:,1], array( [ 0, 1, 0] , float) )
    assert allclose(o.v[:,2], array( [ 0, 0, 1] , float) )

    """ 3 really hard vectors """
    global DEBUG
    DEBUG = True
    o = lattice( array( [ 991, 990, 990] , float) ,
                 array( [ 990, 991, 990] , float) ,
                 array( [ 990, 990, 991] , float) )
    # print "o.v",o.v
    o = o.withvec( [990,990,990] ) # Remove long
    assert allclose(o.v[:,0], array( [ 1, 0, 0] , float) )
    assert allclose(o.v[:,1], array( [ 0, 1, 0] , float) )
    assert allclose(o.v[:,2], array( [ 0, 0, 1] , float) )
    DEBUG = False

        

if __name__=="__main__":
    import time
    start = time.time()
    test1()
    test2()
    test_eu()
    test_fft()
    print time.time()-start
    
