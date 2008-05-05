
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


def this_direction_axis(direction):
    return ['row', 'col'].index(direction)

def other_direction(direction):
    return ['row', 'col'][1 - this_direction_axis(direction)]

#  From http://www.scipy.org/Subclasses
class rr_array(ndarray):
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

    

    def __str__(self):
        desc = \
"""{%(data)s in %(direction)s direction}"""
        return desc % { 'data': super.__str__(self), 
                        'direction' : self.direction }


    def length_of_vectors2(self):
        # print type(self)
        return sum( self*self, axis=this_direction_axis(self.direction))


    


class flipper:
    def __init__(self, ub=None, ubi=None):
        assert ub.shape == (3, 3)
        assert ubi.shape == (3, 3)
        self.ub = ub
        self.ubi = ubi
    def flip(self, v):
        assert hasattr(v, 'direction')
        try:
            if v.direction == 'row':
                # print "flipping",v,'row in so using ubi'
                # print "ubi",self.ubi
                ra = dot(self.ub, v).T
            elif v.direction == 'col':
                # print "flipping",v,'recip in so using ub'
                # print "ub",self.ub
                ra = dot(self.ubi, v.T)
            else:
                raise Exception("Bad direction")
        except:
            print v.shape, v.direction
            raise
        assert ra.shape == v.shape[::-1], "Shape mismatch, %s %s %s"%(
            str(v.shape), str(ra.shape), v.direction)
        return rr_array( ra, direction=other_direction(v.direction))


class lattice(object):
    """
    Represents a 3D crystal lattice built from 3 vectors
    """

    def __init__(self, v1, v2, v3, direction='col', min_vec2=MIN_VEC2):
        """ Make the lattice 
        Currently v1, v2, v3 are vectors - which gives 3D
        direction [ 'row' | 'col' ]
        ... means they are row direction lattice vectors or measured scattering
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
        if direction == 'col':  
            # We got g-vectors
            # print "Supplied col direction vectors"
            self.ub  = array(vl).T 
            self.ubi = inv(self.ub)
        elif direction == 'row':
            # We were supplied with peaks from patterson peaksearch
            # print "Supplied row direction vectors"
            self.ubi = array(vl)
            self.ub  = inv(self.ubi) 
        else:
            raise Exception("Direction must be row or col "+str(direction))
        self.flipper = flipper( ub=self.ub, ubi=self.ubi )
        self.flip = self.flipper.flip

    def nearest(self, vecs):
        """ Give back the nearest lattice point indices, 
        in the same direction """
        assert hasattr(vecs, 'direction')
        new_vecs = self.flip( vecs )
        int_vecs = round_(new_vecs)
        return self.flip( int_vecs )
        
    def remainders(self, vecs):
        """ Difference between vecs and closest lattice points """
        assert hasattr(vecs, 'direction')
        
        return vecs - self.nearest(vecs)

    def withvec(self, x, direction="col"):
        """ 
        Try to fit x into the lattice 
        Make the remainder together with current vectors
        Index it as hkl indices
        whichever vector has the biggest projection is replaced
        remake the lattice with these 3 vectors
        """
        #print self.ubi
        assert hasattr(x, 'direction')
        #print "x",x
        r = self.remainders( x )
        #print "r",r
        worst = argmax(r)
        if r.direction == 'col':
            # r is g-vector errors
            v = list(self.ub.T)
        if r.direction == 'row':
            # r is hkl errors
            v = list(self.ubi)
        v[worst]=r
        l_new = lattice( v[0], v[1], v[2] , direction=r.direction )
        return l_new

    def score(self, vecs, tol=0.1):
        """
        How many peaks have rem less than tol
        """
        assert hasattr(vecs, 'direction')
        r2 = self.remainders( vecs ).length_of_vectors2()
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
            l = lattice(vecs[i], vecs[j], vecs[k], direction=vecs.direction)
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
    vecs = rr_array(g.UBIALL.T , direction='row')
    assert vecs.shape == (3, g.colfile.nrows)
    order = argsort( g.colfile.sum_intensity )[::-1]
    vecs = take( vecs, order, axis = 1)
    min_pks = 300
    ntry = 0
    print ""
    l1 = find_lattice( vecs,
                       min_vec2 = 1,
                       n_try = 20 )
                  
    l2 = find_lattice( vecs,
                       min_vec2 = 1,
                       n_try = 20 )
                  
                  


    for i,j,k in iter3d( min(len(vecs),50)):
        ntry += 1
        # if ntry> 10: break
        try:
            l = lattice( vecs[order[i]], vecs[order[j]], vecs[order[k]],
                         direction='row', min_vec2=1)
            t = 0.1
            ref = refine( l.vi, gv, t)
            l = lattice( ref[0], ref[1], ref[2], direction='row')
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
    # G-vectors are in row direction
    l = lattice ( v1, v2, v3, direction='col')
    esum = 0.0
    gv = rr_array( gv, direction='col' )
    #print v1, v2, v3
    #print l.ubi
    #print l.ub
    for v in gv:
        # print ("%8.5f "*3+"%8.5f "*3)%tuple( list(v)+list(l.flip(v))),
        err   = l.remainders( v )
        esum += sqrt(dot(err,err))
        # print "%.5f"%(sqrt(dot(err,err)))
    #import sys
    #sys.exit()
    #print "Average",esum / len(gv)
    # print 
    assert esum / len(gv) < 0.0102, "Did not fit"
    assert l.score(gv, tol = 0.1) == 602, "Expecting to index 602 peaks" 

def test2():
    """ Adding a fourth vector """
    o = lattice( array( [ 1, 0, 0] , float) ,
                 array( [ 0, 1, 0] , float) ,
                 array( [ 0, 0, 2] , float) , direction = 'row')
    # ... how to do it?
    u = rr_array( [ 0, 0, 1] , dtype=float, direction = 'row') 
    r  = o.remainders( u )
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
    assert allclose(o.ub[:,0], array( [ 1, 0, 0] , float) ), str(o.ub)
    assert allclose(o.ub[:,1], array( [ 0, 1, 0] , float) ), str(o.ub)
    assert allclose(o.ub[:,2], array( [ 0, 0, 1] , float) ), str(o.ub)


    """ 3 shorter vectors """
    # print '3 short'
    o = lattice( array( [1e-2,  0,  0] , float) ,
                 array( [  0,1e-2,  0] , float) ,
                 array( [  0,  0,1e-2] , float) , direction='col' )
    # Consider those as 3 g-vectors
    assert allclose(o.ubi[:,0], array( [ 100, 0, 0] , float) ), str(o.ubi)
    assert allclose(o.ubi[:,1], array( [ 0, 100, 0] , float) ), str(o.ubi)
    assert allclose(o.ubi[:,2], array( [ 0, 0, 100] , float) ), str(o.ubi)
    g = rr_array([1e-2,0.,0.], direction='col') 
    hkl = o.flip( g ) 
    assert hkl.direction == 'row', hkl.direction
    # print g,"\nFlips to give\n",hkl
    assert allclose(hkl, 
                    rr_array([1.,0.,0.], direction='row') ), str(hkl)


    # print '3 long'
    """ 3 longer vectors """
    o = lattice( array( [ 10,  0,  0] , float) ,
                 array( [  0, 10,  0] , float) ,
                 array( [  0,  0, 10] , float) , direction='row' )
    assert allclose(o.ubi[:,0], array( [ 10, 0, 0] , float) ), str(o.ubi)
    assert allclose(o.ubi[:,1], array( [ 0, 10, 0] , float) ), str(o.ubi)
    assert allclose(o.ubi[:,2], array( [ 0, 0, 10] , float) ), str(o.ubi)
    hkl = o.flip( rr_array([10.,10.,10.], direction='row') )
    assert hkl.direction == 'col'
    assert allclose(hkl,
                    rr_array([1.,1.,1.], direction='col') ), str(hkl)

    """ 3 shorter vectors row"""
    # print '3 short'
    o = lattice( array( [1e-2,  0,  0] , float) ,
                 array( [  0,1e-2,  0] , float) ,
                 array( [  0,  0,1e-2] , float) , direction='row' )
    # Consider those as 3 g-vectors
    assert allclose(o.ub[:,0], array( [ 100, 0, 0] , float) ), str(o.ub)
    assert allclose(o.ub[:,1], array( [ 0, 100, 0] , float) ), str(o.ub)
    assert allclose(o.ub[:,2], array( [ 0, 0, 100] , float) ), str(o.ub)
    g = rr_array([1e-2,0.,0.], direction='row') 
    hkl = o.flip( g ) 
    # print g,"\nFlips to give\n",hkl
    assert allclose(hkl, 
                    rr_array([1.,0.,0.], direction='row') ), str(hkl)
    assert hkl.direction == 'col'


    # print '3 long'
    """ 3 longer vectors """
    o = lattice( array( [ 10,  0,  0] , float) ,
                 array( [  0, 10,  0] , float) ,
                 array( [  0,  0, 10] , float) , direction='col' )
    assert allclose(o.ub[:,0], array( [ 10, 0, 0] , float) ), str(o.ub)
    assert allclose(o.ub[:,1], array( [ 0, 10, 0] , float) ), str(o.ub)
    assert allclose(o.ub[:,2], array( [ 0, 0, 10] , float) ), str(o.ub)
    hkl = o.flip( rr_array([10.,10.,10.], direction='col') )
    assert allclose(hkl,
                    rr_array([1.,1.,1.], direction='row') ), str(hkl)

    
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
                 array( [ 990, 990, 991] , float) ,direction='col')
    # Eg, these were long g-vectors far from origin but making a basis
    # Anticipate them making a ubi matrix with 990
    assert ( o.ub  == array( [[ 990. ,   1. ,   1.],
                              [ 990. ,  -0. ,  -1.],
                              [ 991. ,  -1. ,  -0.],] )).all(), str(o.ubi)
    # OK
    u = rr_array([992.,990.,990.], direction='col') 
    assert (o.nearest( u ) == rr_array([991.,990.,990], 
                                       direction='col')).all(),\
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
    
