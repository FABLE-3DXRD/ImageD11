
from __future__ import print_function

from .rc_array import rc_array

from numpy import dot, array, allclose, asarray, fabs, round, \
    argmin, argmax, sqrt, argsort, take, sum, where, ndarray, eye,\
    zeros, cross, pi, arccos, floor
from numpy.linalg import inv, LinAlgError, det

import logging

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
    print("Unexpected exception when checking numpy behaviour")
    raise

DEBUG = False

# Some sort of round off
MIN_VEC2 = 1e-9*1e-9



def fparl(x, y):
    """ fraction of x parallel to y """
    ly2 = dot(1.0*y ,y)
    if ly2 > 1e-9:
        return dot(1.0*x, y) / ly2
    return 0

def mod(x,y):
    """ 
    Returns the part of x with integer multiples of y removed
    assert that ||mod(x,y)|| <= ||x|| for all y 
    """
    if __debug__: 
        b4  = dot(x,x)
    n   = round(fparl(x,y))
    ret = x - n * y
    if __debug__:
        af  = dot(ret,ret)
    if b4 < af and n != 0 : 
        print("Bad mod "+str(x) + " " + str(y))
        print(ret, b4, af, n)
        raise Exception("problem in mod")
    return ret

def sortvec_len( vl ):
    """
    Sorts according to length (d-s-u)
    """
    # Here v is ALWAYS a shape==(3,) vector
    ul = [ ( dot(v,v), tuple(v)) for v in vl ]
    ul.sort()
    return asarray([ v[1] for v in ul[::-1] ])

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

def checkvol(v1, v2, v3, min_vec2=MIN_VEC2):
    """ Check whether vectors are singular """
    v = dot(v1, cross(v2,v3))
    assert abs(v)> pow(min_vec2, 1.5), "Volume problem"
    return True

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


def reduce(v1, v2, v3, min_vec2=MIN_VEC2):
    """
    Try to find a reduced lattice basis of v1, v2, v3
    """
    assert checkvol(v1,v2,v3,min_vec2)
    vl = array([v1, v2, v3])
    vn = rsweep(vl)
    i = 0
    while not allclose(vn ,vl) :
        vl = [ v.copy() for v in vn ]
        vn = rsweep( vl )
        i += 1
        if i > 10:
            raise Exception("Algorithmic flaw")
    # choose the "bigger" compared to -v
    for i in range(3):
        vn[i] = sortvec_xyz( [vn[i], -vn[i]] )[0]
    return vn


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
        try:
            vl =    reduce( v1   ,    v2,    v3, min_vec2 )
            again = reduce( vl[0], vl[1], vl[2], min_vec2 )

        except:
            raise BadVectors()
        # Check reduction is stable
        # print vl
        assert allclose( array(vl),array(again) ), "Bad reduction %s %s"%(
                str(vl), str(again))
        # This cause a problem - why?
        # vl = sortvec_len( vl )
        try:
            if direction == 'col':  
                # print "Supplied col direction vectors"
                self.r2c  = array(vl)
                if det( self.r2c ) < 0:
                    self.r2c = array( [vl[0], vl[2], vl[1] ] )
                self.c2r = inv(self.r2c)
            elif direction == 'row':
                # Supplied with g-vectors
                # print "Supplied row direction vectors"
                self.c2r = array(vl).T
                if det( self.c2r ) < 0:
                    self.c2r = array( [vl[0], vl[2], vl[1] ] ).T
                self.r2c  = inv(self.c2r) 
            else:
                raise Exception("Direction must be row or col "+str(direction))
        except LinAlgError:
            print("problem with vectors")
            print(v1,v2,v3)
            print("Reduced to")
            print(vl)
            raise
        assert self.c2r.shape == (3, 3)
        assert self.r2c.shape == (3, 3)
        
    def flip(self, v):
        """
        See also __init__.__doc__
        """
        assert v.check()
        ret = v.flip(self.matrix(v.direction))
        assert isinstance( ret, rc_array )
        assert ret.check()
        assert ret.shape == v.shape[::-1], "Shape mismatch, %s %s %s %s"%(
            str(v.shape[::-1]),str(v.shape), str(ret.shape), v.direction)
        return ret

    def matrix(self, direction):
        if direction == 'row': return self.r2c
        if direction == 'col': return self.c2r
        raise Exception("direction not in row|col")


    def nearest(self, vecs):
        """ Give back the nearest lattice point indices, 
        in the same direction """
        new_vecs = vecs.flip(self.matrix(vecs.direction))
        int_vecs = round(new_vecs)
        return int_vecs.flip(self.matrix(int_vecs.direction))
        
    def remainders(self, vecs):
        """ Difference between vecs and closest lattice points """
        vecs.check()
        return vecs - self.nearest(vecs)

    def withvec(self, x, direction="col"):
        """ 
        Try to fit x into the lattice 
        Make the remainder together with current vectors
        Index it as hkl indices
        whichever vector has the biggest projection is replaced
        remake the lattice with these 3 vectors
        """
        assert hasattr(x, 'direction')
        r = self.remainders( x )
        worst = argmax(fabs(r))
        if r.direction == 'col':
            # supplied vector is from a patterson
            v = list(self.r2c)
        if r.direction == 'row':
            # supplied vector is a g-vector
            v = list(self.c2r.T)
        v[worst]=r
        l_new = lattice( v[0], v[1], v[2] , direction=r.direction )
        return l_new

    def score(self, vecs, tol=0.1, debug=False):
        """
        How many peaks have integer error less than tol?
        """
        assert vecs.check()
        # These are in units of g-vector
        diffs = self.remainders(vecs)
        # Put into other space to compare to tol
        # ... works for g-vectors as hkl
        int_err = diffs.flip( self.matrix( diffs.direction ) )
        r2 = int_err.norm2()
        if debug:
            print(vecs.shape, r2.shape, tol*tol)
            print(r2[:10])
        s = sum( where( r2 < tol * tol, 1, 0) )
        return s


def iter3d_old(n):
    """
    Generate all possible unordered combinations of vectors i,j,k
    for i,j,k < n
    """
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                yield i,j,k 

def iter3d(n):
    """
    Generate all possible unordered combinations of vectors i,j,k
    for i,j,k < n

    This looping was rewritten thanks to:
    TAOCP V4 fascicle 3, section 7.2.1.3.
    
    It gives much nicer ordering than previously, as is gradually
    expands down the list, instead of hammering the start
    """
    for k in range(2,n):
        for j in range(1,k):
            for i in range(j):
                yield i,j,k

#t1 = [ l for l in iter3d_old(10) ]
#t2 = [ l for l in iter3d(10) ]
#print t1
#print t2
#assert t1 == t2


def find_lattice(vecs, 
                 min_vec2=1,
                 n_try=None, 
                 test_vecs = None, 
                 tol = 0.1,
                 fraction_indexed=0.9,
                 noisy = False ):
    """
    vecs - vectors to use to generate the lattice
    min_vec2 - squared length of min_vec (units as vec)
    n_try - max number of vecs to test (clips vecs for generating)
    test_vecs - vectors to index with the unit cell (else use vecs)
    fraction_indexed  - whether or not to return cells
    """
    assert isinstance(vecs, rc_array)
    if n_try is None:
        n_try = len(vecs)
    else:
        if n_try > vecs.nvectors():
            n_try = vecs.nvectors()
            logging.warning("Adjusting number of trial vectors to %d"%(n_try))
    if test_vecs is None:
        test_vecs = vecs
    assert isinstance(test_vecs, rc_array)
    gen_dir = vecs[0].direction
    if noisy:
        print("Finding with dir",gen_dir)
        for i,v in enumerate(vecs):
            print(i, v)
            if i > n_try:
                break
        print("min_vec2",min_vec2)
    for i,j,k in iter3d(n_try):
        # if (i,j,k) == (0,1,6):
        #    print vecs[i],vecs[j],vecs[k]
        #    print gen_dir, min_vec2
        if noisy:
            print("Try",i,j,k, end=' ')
        try:
            if gen_dir == 'row':
                if dot(vecs[i], vecs[i]) < min_vec2: continue
                if dot(vecs[j], vecs[j]) < min_vec2: continue
                if dot(vecs[k], vecs[k]) < min_vec2: continue
                print(i,j,k, end=' ')
                l = lattice(vecs[i], vecs[j], vecs[k], 
                            direction = gen_dir,
                            min_vec2 = min_vec2)
            elif gen_dir == 'col':
                try:
                    if dot(vecs[:,i], vecs[:,i]) < min_vec2: continue
                    if dot(vecs[:,j], vecs[:,j]) < min_vec2: continue
                    if dot(vecs[:,k], vecs[:,k]) < min_vec2: continue
                    # print i,j,k,dot(vecs[:,i], vecs[:,i]),dot(vecs[:,j], vecs[:,j]),dot(vecs[:,k], vecs[:,k])
                    l = lattice(vecs[:,i], vecs[:,j], vecs[:,k], 
                            direction = gen_dir,
                            min_vec2 = min_vec2)
                except IndexError:
                    print(i,j,k,n_try,vecs.shape)
                    raise
                    
            else:
                raise Exception("Logical impossibility")
            # First test on vecs
            
            scor = l.score( vecs, tol )
            frac = 1.0 * scor / vecs.nvectors()
            if noisy:
                print("Score on vecs",scor,frac, end=' ')
            
            scor = l.score( test_vecs, tol )
            frac = 1.0 * scor / test_vecs.nvectors()
            if noisy:
                print("score on test_vecs",scor,frac)
            if frac > fraction_indexed:
                if noisy:
                    print("Returning")
                return l
        except BadVectors:
            pass
    return None


def cosangle_vec( ubi, v ):
    """
    Angle between v in real and reciprocal space
    eg, is a* parallel to a or not?
    """
    real = dot( ubi.T, v )
    reci = dot( inv(ubi) , v )
    return dot( real, reci )/sqrt(
        dot(real, real) * dot(reci, reci) )


def search_2folds( ubi ):
    """
    Inspired by the Yvon Lepage's method for finding lattice symmetry
    Check for 2 fold axes by measuring the directions between real
    and reciprocal vectors with the same indices. In the case of 2 fold
    axes they should be parallel
    """
    hr = list(range(-2,3))
    for h in hr:
        for k in hr:
            for l in hr:
                if h==0 and k==0 and l==0:
                    continue
                c = cosangle_vec( ubi, [h,k,l] )
                if abs(c - floor( c + 0.5)) < 0.001:
                    print(h, k, l, c, arccos(c)*180/pi)



def get_options(parser):
    parser.add_argument('-v', '--min_vec2',
                      action='store',
                      type=float,
                      dest="min_vec2",
                      help='Minimum axis length ^2, (angstrom^2) [1.5]',
                      default = 1.5)
    parser.add_argument('-m', '--n_try',
                      action='store',
                      type=int,
                      dest="n_try",
                      default=None,
                      help='Number of vectors to test in finding lattice [all]')
    parser.add_argument('-f', '--fraction_indexed',
                      action='store',
                      type=float,
                      dest="fraction_indexed",
                      default=0.9,
                      help='Fraction of peaks to be indexed')
    parser.add_argument('-t','--tol',
                      action='store',
                      type=float,
                      default = 0.1,
                      dest="tol",
                      help='tolerance in hkl error for indexing')
    return parser
                      
                      

