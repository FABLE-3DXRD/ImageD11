

"""
Some simplistic ART algorithms
"""
import numpy as np
from ImageD11 import closest

def update_nearest( recon, proj, angle, debug=True ):
    """
    recon is the current reconstruction
    proj is a projection at angle (degrees)

    the middle pixel of proj matches the middle pixel of recon
    TODO: angle + offset as more general geometric term
    """
    assert len(recon.shape) == 2
    assert len(proj.shape) == 1
    rads = angle * np.pi / 180.0
    cth = np.cos( rads )
    sth = np.sin( rads )
    # For each pixel in recon, compute the array index in proj
    nx = recon.shape[0]
    ny = recon.shape[1]
    x , y = np.mgrid[ -nx/2:nx/2, -ny/2:ny/2 ]
    # pos is a 2D array of positions in reconstruction
    pos = cth * x + sth * y
    # Nearest pixel:
    idx = pos.round().astype(np.int)
    # idx will go from -xxx to xxx
    mn = idx.ravel().min()
    mx = idx.ravel().max()
    calc_proj = np.zeros( (mx-mn), np.float32 )
    closest.put_incr( calc_proj, (idx-mn).ravel(), recon.ravel() )
    if debug:
        print "angle",angle,"recon.shape",recon.shape
        print "idx.shape",idx.shape,"min_i, max_i",mn,mx
        print "calc_proj.shape",calc_proj.shape, proj.shape
        import pylab as pl
        pl.figure()
        pl.hist( (idx-mn).ravel(), bins = 100)
        pl.figure( )
        pl.imshow(recon)
        pl.figure()
        pl.plot( calc_proj )
    return calc_proj


def update_wtd( recon, proj, angle, dbg=True ):
    """
    recon is the current reconstruction
    proj is a projection at angle (degrees)

    the middle pixel of proj matches the middle pixel of recon
    TODO: angle + offset as more general geometric term
    """
    assert len(recon.shape) == 2
    assert len(proj.shape) == 1
    rads = angle * np.pi / 180.0
    cth = np.cos( rads ).astype(np.float32)
    sth = np.sin( rads ).astype(np.float32)
    # For each pixel in recon, compute the array index in proj
    nx = recon.shape[0]
    ny = recon.shape[1]
    x , y = np.mgrid[ -nx/2:nx/2, -ny/2:ny/2 ]
    # pos is a 2D array of positions in reconstruction
    pos = np.zeros( x.shape, np.float32).ravel()
    np.add( x.ravel()*cth, y.ravel()*sth, pos) 
    # Lower & upper pixel:
    idx_lo = np.floor(pos).astype(np.int)
    idx_hi = np.ceil( pos).astype(np.int)
    # idx will go from -xxx to xxx
    mn = idx_lo.min()
    mx = idx_hi.max()+1
    # Weights:
    wt_lo = np.subtract( idx_hi, pos, np.zeros( pos.shape, np.float32) )
    wt_hi = np.subtract( pos, idx_lo, np.zeros( pos.shape, np.float32) )
    #
    wt_lo = np.where( idx_lo == idx_hi, 1, wt_lo)
    calc_proj = np.zeros( (mx-mn), np.float32 )
    r = recon.ravel()
    np.subtract( idx_lo, mn, idx_lo) # Move these to offsets in calc_proj
    np.subtract( idx_hi, mn, idx_hi)
    closest.put_incr( 
            calc_proj, 
            idx_lo, 
            r*wt_lo
            )
    closest.put_incr( 
            calc_proj, 
            idx_hi, 
            r*wt_hi
            )
    if dbg:
        assert wt_lo.ravel().min() >= 0
        assert wt_hi.ravel().min() >= 0
        assert wt_lo.ravel().max() <= 1
        assert wt_hi.ravel().max() <= 1
        print "Avg total weight", (wt_lo + wt_hi).mean()
        import pylab as pl
        print "angle",angle,"recon.shape",recon.shape
        print "calc_proj.shape",calc_proj.shape, proj.shape
        pl.figure( )
        pl.imshow(recon)
        pl.figure()
        pl.plot( calc_proj )
    # Now make calc_proj align with the supplied proj.
    # Convention will be to make the central pixel of our image
    # match the central pixel of proj. Logically this is the central
    # pixel of calc_proj too...
    error = np.zeros( calc_proj.shape, np.float32 )
    start = (len(calc_proj)- len(proj))/2
    error[ start:start+len(proj) ] = proj
    error = error - calc_proj
    # Now take this error off of recon according to the weights
    wt_all = np.zeros( calc_proj.shape, np.float32 )
    closest.put_incr( 
            wt_all,
            idx_lo,
            wt_lo*wt_lo )
    closest.put_incr( 
            wt_all,
            idx_hi,
            wt_hi*wt_hi )
    werr = error/np.clip(wt_all, 1, error.max()  )
    update = np.zeros( r.shape, np.float32)
    if dbg:
        pl.figure()
        pl.title( 'werr' )
        pl.plot( werr )
    print "lo, err.shape",wt_lo.shape, werr.shape
    update = wt_lo*np.take( werr, idx_lo ) + wt_hi*np.take( werr, idx_hi )
    if dbg:
        pl.figure()
        pl.plot( update )
        pl.plot( wt_lo + wt_hi )
    update.shape = recon.shape
    print "Sum of update",update.sum()
    print "Sum of error ",error.sum()
    print "Ratio",update.sum()/error.sum()
    if dbg:
        pl.figure()
        pl.imshow( np.reshape(update, recon.shape) )
        pl.colorbar()
    return update

def make_data():
    import rad, Image
    im = Image.open("phantom3.gif")
    a = np.asarray(im).astype(np.float32)
    theta = range(45,1800,107)
    r = rad.radon( a, theta )
    return a, r, theta

def test_u_n():
    ideal, projections, theta = make_data()
    recon = np.zeros(ideal.shape,np.float32)
    import time, pylab as pl
    sumtime = ncalc = 0
    dbg = False
    x, y = np.mgrid[ 0:recon.shape[0], 0:recon.shape[1] ]
    x = x-recon.shape[0]/2
    y = y-recon.shape[1]/2
    n = min(recon.shape)/2
    msk = np.where( x*x + y*y < n*n , 1, 0)
    #    pl.imshow(msk)
    #pl.show()
    for j in range(10):
        for proj, angle in zip(projections.T, theta):
            print j, angle
            start = time.clock()
            # print proj.shape
            update = update_wtd( recon, proj, angle, dbg=False)
            print "Recon",recon.min(),recon.max(),recon.mean()
            recon = recon + update
            # clip to zero - constructing positive intensity
            recon = np.where( recon > 0, recon, 0)
            np.multiply( recon , msk, recon )
            err = (ideal - recon).ravel()
            print "Error, mean, min, max",err.mean(), err.min(), err.max(), time.clock()-start
            if True:
                pl.figure(1)
                pl.clf()
                pl.subplot(121)
                pl.imshow(recon)
                pl.subplot(122)
                pl.imshow(recon-ideal)
                pl.colorbar()
                pl.show()
            end = time.clock()
            sumtime += end-start
            ncalc +=1
        pl.show()
        raw_input()
        # recon = recon.clip( 0, 10)
    print sumtime/ncalc

if __name__=="__main__":
    test_u_n()




