

from numpy import dot, round_, array, float, allclose, asarray, fabs,\
    argmin, argmax, sqrt, argsort, take, sum, where, ndarray, eye,\
    zeros, cross
from numpy.linalg import inv, LinAlgError

from ImageD11.lattice_reduction import lattice, find_lattice
from ImageD11.rc_array import rc_array

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
              nsig = 20)
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
                       min_vec2 = 9,
                       test_vecs = tv,
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
                       min_vec2 = 9,
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
    print v1,v2,v3
    # This means that the g-vectors are in row direction
    l = lattice ( v1, v2, v3, direction='row')
    esum = 0.0
    gv = rc_array( gv , direction='row' )
    # print v1, v2, v3
    print "ubi/r2c",l.r2c
    print "ub /c2r",l.c2r
    print "dot(l.r2c, gv[0])",dot(l.r2c, gv[0])
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
    assert  s == 583, "Expecting to index 582 peaks, got %s"%(s) 
    print "Indexing of eu3.gve, Average",esum / len(gv),"for",s,"peaks"
    print "UBI is",l.r2c

def test_eu_find():
    gv = get_eu_gv()
    vecs = rc_array(gv, direction="row")
    l = find_lattice( vecs,
                      # This next bit is important - tolerance for gve 
                      # is very different to patterson
                      min_vec2=1./81.,
                      n_try=20,
                      test_vecs=vecs)
    s = l.score(vecs, tol=0.1)
    print "Eu3 using find_lattice scores",s

def test2():
    """ Adding a fourth vector """
    o = lattice( array( [ 1, 0, 0] , float) ,
                 array( [ 0, 1, 0] , float) ,
                 array( [ 0, 0, 2] , float) , direction = 'row')
    # ... how to do it?
    #print o.r2c
    #print o.c2r
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
    g.check()
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
    #assert allclose( o.c2r , array( [[ 991. ,   1. ,   1.],
    #                                 [ 990. ,   0. ,  -1.],
    #                                 [ 990. ,  -1. ,   0.],] )), str(o.c2r)
    # OK
    u = rc_array([990.,990.,990.], direction='row') 
    assert (o.nearest( u ) == rc_array([991.,990.,990], 
                                       direction='row')).all(),\
        str(o.nearest( u ))+" "+str( u )


    o = o.withvec( u ) # Remove long
    # print o.c2r
    # print o.r2c
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
    test_eu_find()
    test_fft()

    print time.time()-start
    
