
import numpy as np, fabio

M = 2048   # Overall image size
N = 64     # tile size
SIG = 3.   # sigma for peaks (width)

# tile co-ordinates
i,j = np.mgrid[0:N,0:N]-N/2.

def pk(x,y,s):
    # 2D peak shape (with tails)
    dx = i-x
    dy = j-y
    gg= np.exp( -(dx*dx+dy*dy)/s/s )
    ll=s*s/(s*s+dx*dx+dy*dy)
    return ll+gg

# Image for testing
im = np.zeros( (M,M) )

# Make a grid of points (tiles)
for ii in range(0,M,N):
    for jj in range(0,M,N):
        # offset of centroids
        oi = 3*SIG*ii/M
        oj = SIG*ii/M
        # Ratio of heights
        sc = jj*1.0/M
        # put in two peaks, offset increases one direction, scale in the other
        im[ii:ii+N,jj:jj+N] = pk( -oi, -oj, SIG ) + pk( oi, oj, SIG)*sc

if 1:
    import pylab as pl
    pl.imshow( pl.log(im ))
    pl.show()

s = 30000/im.max()
int_im = (im*s).astype( np.uint16)
    
fabio.edfimage.edfimage( int_im ).write( "testoverlaps0000.edf" )
    

