

from numpy import dot, round_, array, float, allclose, asarray, fabs,\
    argmin, argmax, sqrt, argsort, take, sum, where, ndarray, eye,\
    zeros, cross
from numpy.linalg import inv, LinAlgError

# Confirm that dot'ting a 3x3 matrix with a 3x10 gives a 3x10
assert dot(eye(3), zeros( (3, 10) ) ).shape == (3, 10), \
    "Numpy dot insanity problem"
           
# It is unclear why it is not a 10x3 result (row/col vectors)

try:
    dot(eye(3), zeros( (10, 3) ) )
    raise Exception("Numpy dot insanity problem")
except ValueError: 
    pass
except:
    print "Unexpected exception when checking numpy behaviour"
    raise

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

class BadVectors(Exception):
    """ Raised by lattice class when vectors are coplanar or colinear"""
    pass

def checkvol(v1, v2, v3):
    """ Check whether vectors are singular """
    v = dot(v1, cross(v2,v3))
    assert abs(v)>0., "Volume problem"
    return True


# TODO move this rc stuff to it's own module

#  From http://www.scipy.org/Subclasses
class rc_array(ndarray):
    def __new__(subtype, data, direction=None, dtype=None, copy=False):
        subarr = array(data, dtype=dtype, copy=copy)
        subarr = subarr.view(subtype)
        if direction is not None:
            subarr.direction = direction
            # print "set it   ",
        elif hasattr(data, 'info'):
            subarr.direction = data.direction
            # print "who knows",
        # print "__new__ received direction",direction
        return subarr
    
    def __array_finalize__(self, obj):
        # 3rd argument is default 
        self.direction = getattr(obj, 'direction', 'row' )
        # print "called __array_finalize__",self.shape
        # assert check(self)

    def __str__(self):
        desc = \
"""{%(data)s in %(direction)s direction}"""
        return desc % { 'data': super.__str__(self), 
                        'direction' : self.direction }

    def __iter__(self):
        """
        Iterate depending on rows columns
        """
        # print "iter called"
        if self.direction == 'row':
            return ndarray.__iter__(self)
        elif self.direction == 'col':
            return ndarray.__iter__(self.T)
        else:
            raise Exception("rc_array with direction not in row|col")

def this_direction_axis(direction):
    """ The axis which has the 3 on it """
    return ['col', 'row'].index(direction)

def other_direction(direction):
    """ The axis having the 'n' on it """
    return ['col', 'row'][1 - this_direction_axis(direction)]

def rc_norm2(v):
    """ len(v*v) for row or col """
    # print type(self)
    return sum( v*v, axis=this_direction_axis(v.direction))


def rsweep( vl ):
    """ 
    One sweep subtracting each from other two 
    This idea comes from Knuth TAOCP sec 3.3.4C 
    """
    vn = asarray(vl).copy()
    for i in range(3):
        for j in range(i+1,i+3):
            k = j%3
            assert i!=k
            vn[k] = mod(vn[k], vn[i])
    return vn


def reduce(v1, v2, v3, min_vec2):
    assert checkvol(v1,v2,v3)
    vl = array([v1, v2, v3])
    vn = rsweep(vl)
    i = 0
    while not allclose(vn ,vl) :
        vl = [ v.copy() for v in vn ]
        vn = rsweep( vl )
        i += 1
        if i>10:
            raise Exception("Algorithmic flaw")
    # choose the "bigger" compared to -v
    for i in range(3):
        vn[i] = sortvec_xyz( [vn[i], -vn[i]] )[0]
    return vn


def check(v):
    """ 
    Ensure we have an rc_array which is well behaved 
    Pattern assert(check(v)) should disappear in optimiser
    """
    assert hasattr(v, 'direction')
    assert v.direction in ['row', 'col']
    if len(v.shape) == 1:
        assert v.shape[0] == 3
    elif len(v.shape) == 2:
        if v.direction == 'row':
            assert v.shape[1] == 3
        if v.direction == 'col':
            assert v.shape[0] == 3
    else:
        raise Exception("Only 1D or 2D rc_arrays allowed so far")
    return True



class lattice(object):
    """
    Represents a 3D crystal lattice built from 3 vectors
    """

    def __init__(self, v1, v2, v3, direction=None, min_vec2=MIN_VEC2):
        """ Make the lattice 
        Currently v1, v2, v3 are vectors - which gives 3D
        direction [ 'row' | 'col' ] - if none read from v1, v2, v3
        ... means they are row direction lattice vectors or measured scattering
            vectors from crystallography
        It will attempt to find a primitive basis from the supplied vectors

        These are put together in the following way:

           gv = [ [ g0x , g0y, g0z ],
                  [ g1x , g1y, g1z ],     <-- g is a row vector
                  [ g2x , g2y, g2z ],
                  [ g3x , g3y, g3z ],
                  [ g4x , g4y, g4z ] ]

        [ h, k, l ]  = [ [ col0x  col0y  col0z ] ] [ gx ]
                       |   -------------------   | [    ]
                       | [ col1x  col1y  col1z ] | [ gy ]
                       |   -------------------   | [    ]
                       [ [ col2x  col2y  col2z ] ] [ gz ]

        This matrix of column vectors is often called 

        or ...  h.T  =  dot(r2c, g.T)
           ...  since h is column vector - appears here as a row vecot
           ...  Note that numpy.dot transposes its result
                eg:   s = (3,4) ; 
                assert dot(zeros( (3,3) ), zeros( s )).shape == s
                This is a double head screwing up error. 
                            
        [ gx  gy  gz ]  = [  -------   -------   -------      [ hx ]
                          [ [ row0x ] [ row1y ] [ row2z ] ]   [    ]
                          [ | row0y | [ row1y ] [ row2z ] ]   [ hy ]
                          [ [ row0z ] [ row1y ] [ row2z ] ]   [    ]
                          [  -------   -------   -------      [ hz ]

        ... due to the double head screwing (dot transposes its result)
            we then have to transpose this too
        
        

        """
        if direction is None:
            assert hasattr(v1, 'direction')
            assert hasattr(v2, 'direction')
            assert hasattr(v3, 'direction')
            assert v1.direction == v2.direction
            assert v1.direction == v3.direction
            direction = v1.direction
        # Direction is irrelevant for reduction to shortest 3
        vl =    reduce( v1   ,    v2,    v3, min_vec2 )
        again = reduce( vl[0], vl[1], vl[2], min_vec2 )  
        # Check reduction is stable
        assert allclose( array(vl),array(again) ), "Bad reduction %s %s"%(
                str(vl), str(again))
        try:
            if direction == 'col':  
                # print "Supplied col direction vectors"
                self.r2c  = array(vl)
                self.c2r = inv(self.r2c)
            elif direction == 'row':
                # Supplied with g-vectors
                # print "Supplied row direction vectors"
                self.c2r = array(vl).T
                self.r2c  = inv(self.c2r) 
            else:
                raise Exception("Direction must be row or col "+str(direction))
        except LinAlgError:
            print "problem with vectors"
            print v1,v2,v3
            print "Reduced to"
            print vl
            raise
        assert self.c2r.shape == (3, 3)
        assert self.r2c.shape == (3, 3)
        


    def flip(self, v):
        """
        See also __init__.__doc__
        """
        assert check(v)
        # Case of many vectors, shape == (n, 3) 
        if v.direction == 'row':
            # print "flipping",v,'row in so using r2c'
            # print "r2c",self.r2c
            ra = dot(self.r2c, v.T)
            # Note that ra here has the same shape as v.T
            # This appears to be a nonsensical behaviour of dot
        elif v.direction == 'col':
            # print "flipping",v,'recip in so using c2r'
            # print "c2r",self.c2r
            ra = dot(self.c2r, v).T
            # The .T is because dot doesn't when it should (mathematically)
        # This implies transposing
        # single vector is also OK
        ret = rc_array( ra, direction = other_direction(v.direction))
        assert check(ret)
        assert ret.shape == v.shape[::-1], "Shape mismatch, %s %s %s %s"%(
            str(v.shape[::-1]),str(v.shape), str(ra.shape), v.direction)
        return ret




    def nearest(self, vecs):
        """ Give back the nearest lattice point indices, 
        in the same direction """
        new_vecs = self.flip( vecs )
        int_vecs = round_(new_vecs)
        return self.flip( int_vecs )
        
    def remainders(self, vecs):
        """ Difference between vecs and closest lattice points """
        check(vecs)
        return vecs - self.nearest(vecs)

    def withvec(self, x, direction="col"):
        """ 
        Try to fit x into the lattice 
        Make the remainder together with current vectors
        Index it as hkl indices
        whichever vector has the biggest projection is replaced
        remake the lattice with these 3 vectors
        """
        #print self.r2c
        assert hasattr(x, 'direction')
        #print "x",x
        r = self.remainders( x )
        #print "r",r
        worst = argmax(fabs(r))
        if r.direction == 'col':
            # r is g-vector errors
            v = list(self.r2c.T)
        if r.direction == 'row':
            # r is hkl errors
            v = list(self.c2r)
        #print 'v',v
        v[worst]=r
        #print 'worst',worst
        l_new = lattice( v[0], v[1], v[2] , direction=r.direction )
        return l_new

    def score(self, vecs, tol=0.1):
        """
        How many peaks have rem less than tol
        """
        assert check(vecs)
        # print "In score, tol=",tol
        diffs = self.remainders(vecs)
        r2 = rc_norm2( diffs )
        #for i in range(10):
        #    print vecs[i], diffs[i], r2[i]
        #print r2.shape, diffs.shape
        s = sum( where( r2 < tol * tol, 1, 0) )
        return s


def iter3d(n):
    """
    Generate all possible unordered combinations of vectors i,j,k
    for i,j,k < n
    """
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                yield i,j,k 


def find_lattice(vecs, 
                 min_vec2=1,
                 n_try=None, 
                 test_vecs = None, 
                 tol = 0.1,
                 fraction_indexed=0.9):
    """
    vecs - vectors to use to generate the lattice
    min_vec2 - squared length of min_vec (units as vec)
    n_try - max number of vecs to test (clips vecs for generating)
    test_vecs - vectors to index with the unit cell (else use vecs)
    fraction_indexed  - whether or not to return cells
    """
    assert hasattr( vecs, 'direction' )
    if n_try is None:
        n_try = len(vecs)
    if test_vecs is None:
        test_vecs = vecs
    assert hasattr( test_vecs, 'direction' )
    gen_dir = vecs[0].direction
    for i,j,k in iter3d(n_try):
        try:
            if gen_dir == 'row':
                l = lattice(vecs[i], vecs[j], vecs[k], direction=gen_dir)
            elif gen_dir == 'col':
                l = lattice(vecs[:,i], vecs[:,j], vecs[:,k], direction=gen_dir)
            else:
                raise Exception("Logical impossibility")
            scor = l.score( test_vecs, tol )
            frac = 1.0 * scor / len(test_vecs)
            if frac > fraction_indexed:
                return l
        except BadVectors:
            pass
    return None

def get_eu_gv():
    from ImageD11.indexing import indexer, ubitocellpars
    o = indexer()
    o.readgvfile("eu3.gve" , quiet = True)
    return o.gv

def test_fft():
    gv = get_eu_gv()
    from ImageD11.fft_index_refac import grid
    from ImageD11.indexing import ubitocellpars, write_ubi_file,\
        refine
    g = grid( np = 128,
              mr = 1.0,
              nsig = 20,
              minlen = 3. )
    g.gv_to_grid_new(gv)
    g.fft()
    g.props()
    g.peaksearch(open("eu.patterson_pks","w"))
    g.read_peaks("eu.patterson_pks")
    vecs = rc_array(g.UBIALL.T , direction='col')
    assert vecs.shape == (3, g.colfile.nrows)
    order = argsort( g.colfile.sum_intensity )[::-1]
    vecs = take( vecs, order, axis = 1)
    min_pks = 300
    ntry = 0

    assert g.gv.shape[1] == 3
    tv = rc_array( g.gv, direction='row' )
    print "Finding lattice l1 from patterson"
    l1 = find_lattice( vecs,
                       min_vec2 = 1,
                       n_try = 20 )
    print "r2c == ubi matrix"
    print l1.r2c 
    print "scores", l1.score( tv )
    print "cell",ubitocellpars(l1.r2c)
    l1.r2c = refine( l1.r2c, g.gv, tol=0.1)
    l1.c2r = inv(l1.r2c)
    print "With refine",l1.score( tv )
    print "cell",ubitocellpars(l1.r2c)
    print "Finding lattice l2 with gvectors to test"
    l2 = find_lattice( vecs,
                       min_vec2 = 1,
                       n_try = 20,
                       test_vecs = tv)
    print "r2c == ubi matrix"
    print l2.r2c       
    print "scores", l2.score (tv)
    print "cell",ubitocellpars(l2.r2c)
    l2.r2c = refine( l2.r2c, g.gv, tol=0.1)
    l2.c2r = inv(l2.r2c)
    print "With refine",l2.score( tv )
    print "cell",ubitocellpars(l2.r2c)
    return

def test_eu():
    #Conventional cell gives:
    #0 = (0, 1,  1)
    #1 = (0, 1, -1)
    #6 = (1, 0, -1)
    gv = get_eu_gv()
    v1 = gv[0]
    v2 = gv[1]
    v3 = gv[6]
    assert len(v1) == 3
    # This means that the g-vectors are in row direction
    l = lattice ( v1, v2, v3, direction='row')
    esum = 0.0
    gv = rc_array( gv , direction='row' )
    # print v1, v2, v3
    print "ubi/r2c",l.r2c
    print "ub /c2r",l.c2r
    print "dot(l.r2c, gv[0])",dot(l.r2c.T, gv[0])
    for v in gv:
        #print ("%8.5f "*3+"%8.5f "*3)%tuple( list(v)+list(l.flip(v))),
        assert len(v) == 3
        err   = l.remainders( v )
        esum += sqrt(dot(err,err))
        #print "%.5f"%(sqrt(dot(err,err)))
    #import sys
    #sys.exit()
    # print 
    assert esum / len(gv) < 0.0102, "Did not fit"
    s = l.score(gv, tol = 0.1)
    assert  s == 602, "Expecting to index 602 peaks, got %s"%(s) 
    print "Indexing of eu3.gve, Average",esum / len(gv),"for",s,"peaks"
    print "UBI is",l.r2c


def test2():
    """ Adding a fourth vector """
    o = lattice( array( [ 1, 0, 0] , float) ,
                 array( [ 0, 1, 0] , float) ,
                 array( [ 0, 0, 2] , float) , direction = 'row')
    # ... how to do it?

    u = rc_array( [ 0, 0, 1] , dtype=float, direction = 'row') 
    r  = o.remainders( u )
    assert ( r == array( [ 0, 0, 1] , float) ).all()
    o = o.withvec( r )
    assert allclose(o.c2r[0], array( [ 1, 0, 0] , float) )
    assert allclose(o.c2r[1], array( [ 0, 1, 0] , float) )
    assert allclose(o.c2r[2], array( [ 0, 0, 1] , float) )
    assert allclose(o.c2r, o.r2c)
    
def test1():
    """ Make a lattice from 3 vectors"""
    o = lattice( array( [ 1, 0, 0] , float) ,
                 array( [ 2, 1, 0] , float) ,
                 array( [ 3, 4, 1] , float) , direction = 'col')
    assert allclose(o.c2r[:,0], array( [ 1, 0, 0] , float) )
    assert allclose(o.c2r[:,1], array( [ 0, 1, 0] , float) )
    assert allclose(o.c2r[:,2], array( [ 0, 0, 1] , float) )
    # Both are identity
    assert allclose(o.c2r, eye(3))
    assert allclose(o.r2c, eye(3))

    """ 3 more difficult vectors """
    o = lattice( array( [ 1, 0, 1] , float) ,
                 array( [ 0, 1, 0] , float) ,
                 array( [ 0, 0, 1] , float) , direction = 'col')
    assert allclose(o.c2r[:,0], array( [ 1, 0, 0] , float) ), str(o.c2r)
    assert allclose(o.c2r[:,1], array( [ 0, 1, 0] , float) ), str(o.c2r)
    assert allclose(o.c2r[:,2], array( [ 0, 0, 1] , float) ), str(o.c2r)
    assert allclose(o.c2r, eye(3))
    assert allclose(o.r2c, eye(3))
    

    """ 3 shorter vectors """
    # print '3 short'
    o = lattice( array( [1e-2,  0,  0] , float) ,
                 array( [  0,1e-2,  0] , float) ,
                 array( [  0,  0,1e-2] , float) , direction='row' )
    # Consider those as 3 g-vectors
    # supplied rows go into c2r matrix as columns
    # The r2c matrix (==ubi) is then the inverse of this
    assert allclose(o.r2c[:,0], array( [ 100, 0, 0] , float) ), str(o.r2c)
    assert allclose(o.r2c[:,1], array( [ 0, 100, 0] , float) ), str(o.r2c)
    assert allclose(o.r2c[:,2], array( [ 0, 0, 100] , float) ), str(o.r2c)

    # g-vectors written as rows here - note transpose
    g = array( [[1e-2,   0.,0.], 
                [1e-2, 1e-2,0.], 
                [1e-2,   0.,1e-2],
                [  0.,   0.,1e-2]])
    #print g, g.shape
    g = rc_array( g , direction='row')
    assert g.shape == ( 4, 3)
    check(g)
    #print g.shape
    #print g
    hkl = o.flip( g ) 
    assert hkl.direction == 'col', hkl.direction
    # print g,"\nFlips to give\n",hkl
    assert allclose(hkl, 
                    rc_array(array([[1.,0.,0.],
                                    [1.,1.,0.],
                                    [1.,0.,1.],
                                    [0.,0.,1.]]).T,
                             direction='col') ), str(hkl)


    # print '3 long'
    """ 3 longer vectors """
    o = lattice( array( [ 10,  0,  0] , float) ,
                 array( [  0, 10,  0] , float) ,
                 array( [  0,  0, 10] , float) , direction='col' )
    assert allclose(o.r2c[:,0], array( [ 10, 0, 0] , float) ), str(o.r2c)
    assert allclose(o.r2c[:,1], array( [ 0, 10, 0] , float) ), str(o.r2c)
    assert allclose(o.r2c[:,2], array( [ 0, 0, 10] , float) ), str(o.r2c)
    test_col_vec = rc_array([10.,10.,10.], direction='col') 
    flipped = o.flip( test_col_vec )
    assert flipped.direction == 'row'
    assert allclose( flipped,
                     rc_array([1.,1.,1.], direction='col') ), str(flipped)

    """ 3 shorter vectors row"""
    # print '3 short'
    o = lattice( array( [1e-2,  0,  0] , float) ,
                 array( [  0,1e-2,  0] , float) ,
                 array( [  0,  0,1e-2] , float) , direction='row' )
    # Consider those as 3 g-vectors
    assert allclose(o.r2c[:,0], array( [ 100, 0, 0] , float) ), str(o.r2c)
    assert allclose(o.r2c[:,1], array( [ 0, 100, 0] , float) ), str(o.r2c)
    assert allclose(o.r2c[:,2], array( [ 0, 0, 100] , float) ), str(o.r2c)
    g = rc_array([[1e-2,   0.,0.], 
                  [1e-2, 1e-2,0.], 
                  [1e-2,   0.,1e-2],
                  [  0.,   0.,1e-2]],
                 direction='row')
    hkl = o.flip( g ) 
    # print g,"\nFlips to give\n",hkl

    assert allclose(hkl, 
                    rc_array([[1.,0.,0.],
                              [1.,1.,0.],
                              [1.,0.,1.],
                              [0.,0.,1.]],
                             direction='col').T ), str(hkl)
    assert hkl.direction == 'col'


    # print '3 long'
    """ 3 longer vectors """
    o = lattice( array( [ 10,  0,  0] , float) ,
                 array( [  0, 10,  0] , float) ,
                 array( [  0,  0, 10] , float) , direction='col' )
    assert allclose(o.r2c[:,0], array( [ 10, 0, 0] , float) ), str(o.r2c)
    assert allclose(o.r2c[:,1], array( [ 0, 10, 0] , float) ), str(o.r2c)
    assert allclose(o.r2c[:,2], array( [ 0, 0, 10] , float) ), str(o.r2c)
    hkl = o.flip( rc_array([10.,10.,10.], direction='col') )
    assert allclose(hkl,
                    rc_array([1.,1.,1.], direction='row') ), str(hkl)

    
    # FIXME - check is transpose of row is OK or not 
    # this break the row / col symmetry?
    
    # TODO - think on whether h / g are really col pairs
    #        Find the Sands book.



    # print hkl
    
    """ 3 really hard vectors """
    global DEBUG
    DEBUG = True
    o = lattice( array( [ 991, 990, 990] , float) ,
                 array( [ 990, 991, 990] , float) ,
                 array( [ 990, 990, 991] , float) ,direction='row')
    # Eg, these were long g-vectors far from origin but making a basis
    # Anticipate them making a ubi matrix with 990
    assert allclose( o.c2r , array( [[ 991. ,   1. ,   1.],
                                     [ 990. ,  -1. ,   0.],
                                     [ 990. ,   0. ,  -1.],] )), str(o.c2r)
    # OK
    u = rc_array([990.,990.,990.], direction='row') 
    assert (o.nearest( u ) == rc_array([991.,990.,990], 
                                       direction='row')).all(),\
        str(o.nearest( u ))+" "+str( u )


    o = o.withvec( u ) # Remove long
    # print o.c2r
    print o.r2c
    assert allclose(o.c2r[:,0], array( [ 1, 0, 0] , float) ), o.c2r
    assert allclose(o.c2r[:,1], array( [ 0, 1, 0] , float) ), o.c2r
    assert allclose(o.c2r[:,2], array( [ 0, 0, 1] , float) ), o.c2r
    DEBUG = False

        

if __name__=="__main__":
    import time
    start = time.time()
    test1()
    test2()
    test_eu()
    test_fft()
    print time.time()-start
    
