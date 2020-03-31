

import numba, numpy as np
import timeit
timer = timeit.default_timer

def nbsplat( 
    rgba, gve, M, 
    npx = 2,
    zsort=True, 
    fg = (255,255,255,255),
    bg = (0,0,0,255),
    colors = None ):
    """ Type checking wrapper - does the z sort and clip """
    astart = start = timer()
    h,w,d = rgba.shape
    rve = np.dot( M, gve )
    msk = (abs(rve[0]) < (w//2 - npx)) & (abs(rve[1]) < (h//2 - npx) )
    if zsort:
        start = timer()
        order = np.argsort( rve[2] )[msk]
    else:
        order = np.arange( len(rve[2]), dtype(np.int) )[msk]
    start = timer()
    i = (rve[1]+h//2).astype(np.int16)
    j = (rve[0]+w//2).astype(np.int16)
    if bg is not None:
        # surprising but quite a lot quicker to use fill
        rgba.view( np.uint32 ).fill( np.array( bg, np.uint8).view(np.uint32)[0] )
    if colors is None:
        fgpx = np.array( fg, np.uint8 )
        nbsplat_onecolor( rgba, i, j, fgpx, order, npx )
    else:
        nbsplat_manycolor( rgba, i, j, colors, order, npx )
    

@numba.njit(boundscheck = False )
def nbsplat_onecolor( rgba, i, j, fg, order, npx):
    for k in order:
        rgba[ i[k]-npx:i[k]+npx+1 , j[k]-npx:j[k]+npx+1  ] = fg
    

@numba.njit(boundscheck = False )
def nbsplat_manycolor( rgba, i, j, c, order, npx):
    for k in order:
        rgba[ i[k]-npx:i[k]+npx+1 , j[k]-npx:j[k]+npx+1  ] = c[k]

