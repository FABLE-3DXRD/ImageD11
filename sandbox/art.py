
from __future__ import print_function, division


"""
Some simplistic ART algorithms
"""
import numpy as np
from ImageD11 import cImageD11
import pylab as pl


def update_wtd( recon, proj, angle, msk, dbg=True ):
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
    idx_lo = np.floor(pos).astype(int)
    idx_hi = idx_lo+1
    # idx will go from -xxx to xxx
    mn = idx_lo.min()
    mx = idx_hi.max()+1
    # Weights:
    wt_lo = np.subtract( idx_hi, pos, np.zeros( pos.shape, np.float32) )
    wt_hi = np.subtract( pos, idx_lo, np.zeros( pos.shape, np.float32) )
    wt_lo = np.where(msk.ravel(), wt_lo, 0)
    wt_hi = np.where(msk.ravel(), wt_hi, 0)
    #
    wtc = wt_lo + wt_hi
    calc_proj = np.zeros( (mx-mn), np.float32 )    
    r = recon.ravel()
    np.subtract( idx_lo, mn, idx_lo) # Move these to offsets in calc_proj
    np.subtract( idx_hi, mn, idx_hi)
    cImageD11.put_incr( calc_proj, idx_lo, r*wt_lo )
    cImageD11.put_incr( calc_proj, idx_hi, r*wt_hi )
    error = np.zeros( calc_proj.shape, np.float32 )
    start = int((len(calc_proj)- len(proj))/2)
    error[ start:start+len(proj) ] = proj
    error = error - calc_proj
    npcalc_proj = np.zeros( (mx-mn), np.float32 )
    cImageD11.put_incr( npcalc_proj, idx_hi, wt_hi )
    cImageD11.put_incr( npcalc_proj, idx_lo, wt_lo )
    npcalc_proj = np.where(npcalc_proj >0 ,npcalc_proj, 1)
    # update =  np.zeros( r.shape, np.float32)
    update =  wt_hi*np.take( error/npcalc_proj, idx_hi ) 
    update += wt_lo*np.take( error/npcalc_proj, idx_lo )    

    # Now make calc_proj align with the supplied proj.
    # Convention will be to make the central pixel of our image
    # match the central pixel of proj. Logically this is the central
    # pixel of calc_proj too...
    # Now take this error off of recon according to the weights

    #tot    =  np.zeros( (mx-mn), np.float32 )
    #cImageD11.put_incr( tot, idx_hi, wt_hi)
    #cImageD11.put_incr( tot, idx_lo, wt_lo)
    
    
    
    
    #pl.imshow(np.reshape(update,recon.shape))
    #pl.colorbar()
    #pl.show()

    # update = update 
    update.shape = recon.shape
    print(update.sum(),error.sum(), end=' ')
    if dbg:
        pl.figure()
        pl.imshow( np.reshape(update, recon.shape) )
        pl.colorbar()
        # pl.show()
    return update



def make_data():
    import rad
    from PIL import Image
    im = Image.open("phantom3.gif")
    a = np.asarray(im).astype(np.float32)
    print(a.min(),a.max())
    # 70 random projection
    np.random.seed(42)
    theta = np.random.random(70)*180.
    r = []
    for t in theta:
        r.append( np.asarray( im.rotate( -t, resample=Image.BILINEAR) ).astype(np.float32).sum(axis=1) )
    return a, np.array(r), theta

def test_u_n():
    import time, pylab as pl

    ideal, projections, theta = make_data()

    rshape = (ideal.shape[0]*2, ideal.shape[1]*2)
    recon = np.zeros(rshape,np.float32)

    sumtime = ncalc = 0
    dbg = False
    x, y = np.mgrid[ 0:recon.shape[0], 0:recon.shape[1] ]
    x = x-recon.shape[0]/2
    y = y-recon.shape[1]/2
    n = min(recon.shape)/2
    msk = np.where( x*x + y*y < n*n , 1, 0)
    realstart = time.time()
    pl.ion()
    for j in range(3):
        for proj, angle in zip(projections, theta):
            print(j, angle, end=' ')
            start = time.time()
            update = update_wtd( recon, proj, angle, msk, dbg=False)
            recon = recon + update*0.9
            # clip to zero - constructing positive intensity
            recon = np.where( recon > 0, recon, 0)
            # err = (ideal - recon).ravel()
            if False:
                pl.figure(1)
                pl.clf()
                pl.imshow(recon)
                pl.colorbar()
            #    pl.show()
            end = time.time()
            sumtime += end-start
            ncalc +=1
            # print err.sum()
            print()
        # print "Waiting"
        # raw_input()
        # recon = recon.clip( 0, 10)
    print(sumtime/ncalc)
    np.save("recon.npy",recon)
    end = time.time()
    print("Time take",realstart - end)
    pl.figure(1)
    pl.clf()
    pl.subplot(221)
    pl.title("recon")
    pl.imshow(recon)
    pl.colorbar()
    pl.subplot(222)
    pl.title("ideal")
    pl.imshow( ideal)
    pl.colorbar()
#    pl.subplot(223)
#    pl.title("difference")
#    pl.imshow(ideal-recon)
#    pl.colorbar()
    print(len(theta),"projections")
    pl.show()


if __name__=="__main__":
    test_u_n()




