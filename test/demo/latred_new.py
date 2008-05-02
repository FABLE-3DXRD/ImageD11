
from numpy import dot, round_, array, float, allclose, asarray, fabs,\
        argmax, sqrt, argsort, take, sum, where, ndarray
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

class BadVectors(Exception):
    pass


def this_space_axis(space):
    return ['real', 'reciprocal'].index(space)

def other_space(space):
    return ['real', 'reciprocal'][1 - this_space_axis(space)]

#  From http://www.scipy.org/Subclasses
class rr_array(ndarray):
    def __new__(subtype, data, space=None, dtype=None, copy=False):
        subarr = array(data, dtype=dtype, copy=copy)
        subarr = subarr.view(subtype)
        if space is not None:
            subarr.space = space
            # print "set it   ",
        elif hasattr(data, 'info'):
            subarr.space = data.space
            # print "who knows",
        # print "__new__ received space",space
        return subarr
    
    def __array_finalize__(self, obj):
       self.space = getattr(obj, 'space', 'real' )

    def __repr__(self):
        desc = """\
array(data=%(data)s, 
  in %(space)s space)"""
        return desc % { 'data': str(self), 
                        'space' : self.space }

    def length_of_vectors(self):
        # print type(self)
        return sum( self*self, axis=this_space_axis(self.space))


    


class flipper:
    def __init__(self, ub=None, ubi=None):
        assert ub.shape == (3, 3)
        assert ubi.shape == (3, 3)
        self.ub = ub
        self.ubi = ubi
    def flip(self, v):
        assert hasattr(v, 'space')
        if v.space == 'real':
            print "flipping",v,'real in so using ubi'
            print "ubi",self.ubi
            ra = dot(self.ubi, v.T)
        elif v.space == 'reciprocal':
            print "flipping",v,'recip in so using ub'
            print "ub",self.ub
            ra = dot(self.ub, v)
        else:
            raise Exception("Bad space")
        return rr_array( ra, space=other_space(v.space))


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
        if len(vl) != 3:
            raise BadVectors()
        # print "vl:",vl
        vl = sortvec_xyz(vl)
        # print "vl",vl
        if space == 'reciprocal':  
            # We got column vectors
            self.ub  = array(vl).T 
            self.ubi = inv(self.ub)
        elif space == 'real':
            # We got row vectors
            self.ubi = array(vl)
            self.ub  = inv(self.ubi) 
        else:
            raise Exception("Space must be real or reciprocal "+str(space))
        self.flipper = flipper( ub=self.ub, ubi=self.ubi )

    def nearest(self, vecs):
        """ Give back the nearest lattice point indices, in the same space """
        assert hasattr(vecs, 'space')
        print 'vecs',vecs
        print "ubi",self.ubi
        print "ub",self.ub
        print "Going to flip"
        new_vecs = self.flipper.flip( vecs )
        assert hasattr(new_vecs, 'space')
        print "new_vecs",new_vecs
        int_vecs = round_(new_vecs)
        print 'int_vecs',int_vecs
        int_vecs.space = new_vecs.space

        return self.flipper.flip( int_vecs )
        
    def remainders(self, vecs):
        """ Difference between vecs and closest lattice points """
        assert hasattr(vecs, 'space')
        int_in_other_space = self.nearest(vecs)
        # print "ints",int_in_other_space
        computed_in_this_space = self.flipper.flip( int_in_other_space )
        err = vecs - computed_in_this_space
        return err

    def withvec(self, x, space="reciprocal"):
        """ 
        Try to fit x into the lattice 
        Make the remainder together with current vectors
        Index it as hkl indices
        whichever vector has the biggest projection is replaced
        remake the lattice with these 3 vectors
        """
        #print self.ubi
        assert hasattr(x, 'space')
        #print "x",x
        r = self.remainders( x )
        #print "r",r
        worst = argmax(r)
        if r.space == 'reciprocal':
            # r is g-vector errors
            v = list(self.ub.T)
        if r.space == 'real':
            # r is hkl errors
            v = list(self.ubi)
        v[worst]=r
        print v
        l_new = lattice( v[0], v[1], v[2] , space=r.space )
        print l_new.ubi
        print l_new.ub
        return l_new

    def score(self, vecs, tol=0.1):
        """
        How many peaks have rem less than tol
        """
        assert hasattr(vecs, 'space')
        r2 = self.remainders( vecs, tol ).length_of_vectors2()
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
                 space='reciprocal', 
                 min_vec2=1,
                 n_try=None, 
                 test_vecs = None, 
                 testspace=None,
                 tol= 0.1,
                 fraction_indexed=0.9):
    """
    vecs - vectors to use to generate the lattice
    space - whether vecs are real space or reciprocal
    min_vec2 - squared length of min_vec (units as vec)
    n_try - max number of vecs to test (clips vecs for generating)
    test_on - vectors to index with the unit cell
    testspace - real/recip for test set
    fraction_indexed  - whether or not to return cells
    """
    if n_try is None:
        n_try = len(vecs)
    if test_vecs is None:
        assert testspace is None, "testspace arg implies test_on arg"
        test_vecs = vecs
    if testspace is None:
        testspace = space
    assert testspace in ['real' , 'reciprocal' ], 'space must be real or reciprocal'
    for i,j,k in iter3d(len(vecs)):
        try:
            l = lattice(vecs[i], vecs[j], vecs[k], space=space)
            if testspace == 'real':
                scor = l.score_real( test_vecs, tol )
            elif testspace == 'reciprocal':
                scor = l.score_reciprocal( test_vecs, tol )
            else:
                raise Exception("Logic error")
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

    vecs = rr_array(g.UBIALL, space='real')
    assert vecs.shape == (g.colfile.nrows, 3)
    order = argsort( g.colfile.sum_intensity )[::-1]
    vecs = take( vecs, order, axis = 1)
    min_pks = 300
    ntry = 0
    print ""
    l1 = find_lattice( vecs,
                       space='real',
                       min_vec2 = 1,
                       n_try = 20 )
                  
    l2 = find_lattice( vecs,
                       space='real',
                       min_vec2 = 1,
                       n_try = 20 )
                  
                  


    for i,j,k in iter3d( min(len(vecs),50)):
        ntry += 1
        # if ntry> 10: break
        try:
            l = lattice( vecs[order[i]], vecs[order[j]], vecs[order[k]],
                         space='real', min_vec2=1)
            t = 0.1
            ref = refine( l.vi, gv, t)
            l = lattice( ref[0], ref[1], ref[2], space='real')
            print "i,j,k",i,j,k
            print "UBI:"
            print l.vi
            print ubitocellpars(l.vi)
            s = l.score_recip( gv , t )
            print "Indexes",s ,"at tolerance",t
        except BadVectors:
            continue
        except:
            import traceback
            print traceback.print_exc()
            raise Exception("error")
        if s > min_pks:
            write_ubi_file( "trial.ubi", [ l.vi ] )
            break


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
    assert allclose(o.ub[0], array( [ 1, 0, 0] , float) )
    assert allclose(o.ub[1], array( [ 0, 1, 0] , float) )
    assert allclose(o.ub[2], array( [ 0, 0, 1] , float) )
    
def test1():
    """ Make a lattice from 3 vectors"""
    o = lattice( array( [ 1, 0, 0] , float) ,
                 array( [ 2, 1, 0] , float) ,
                 array( [ 3, 4, 1] , float) )
    assert allclose(o.ub[:,0], array( [ 1, 0, 0] , float) )
    assert allclose(o.ub[:,1], array( [ 0, 1, 0] , float) )
    assert allclose(o.ub[:,2], array( [ 0, 0, 1] , float) )
    """ 3 more difficult vectors """
    o = lattice( array( [ 1, 0, 1] , float) ,
                 array( [ 0, 1, 0] , float) ,
                 array( [ 0, 0, 1] , float) )
    assert allclose(o.ub[:,0], array( [ 1, 0, 0] , float) )
    assert allclose(o.ub[:,1], array( [ 0, 1, 0] , float) )
    assert allclose(o.ub[:,2], array( [ 0, 0, 1] , float) )

    """ 3 really hard vectors """
    global DEBUG
    DEBUG = True
    o = lattice( array( [ 991, 990, 990] , float) ,
                 array( [ 990, 991, 990] , float) ,
                 array( [ 990, 990, 991] , float) ,space='reciprocal')
    assert ( o.ub  == array( [[ 990. ,   1. ,   1.],
                              [ 990. ,  -0. ,  -1.],
                              [ 991. ,  -1. ,  -0.],] )).all(), str(o.ub)

    u = rr_array([990.,990.,990.], space='reciprocal') 
    assert (o.nearest( u ) == rr_array([991.,990.,000.], 
                                       space='reciprocal')).all(),\
        str(o.nearest( u ))+" "+str( u )
    
    o = o.withvec( u ) # Remove long
    assert allclose(o.ub[:,0], array( [ 1, 0, 0] , float) )
    assert allclose(o.ub[:,1], array( [ 0, 1, 0] , float) )
    assert allclose(o.ub[:,2], array( [ 0, 0, 1] , float) )
    DEBUG = False

        

if __name__=="__main__":
    import time
    start = time.time()
    test1()
    test2()
    test_eu()
    test_fft()
    print time.time()-start
    
