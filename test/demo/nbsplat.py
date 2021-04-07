

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
    order = np.empty( len(rve[0]), np.int )
    tmp = np.empty( len(rve[0]), np.int )
    npk = do_zsort_clip( rve, tmp, w//2-npx, h//2-npx, zsort, order )
    i = (rve[1]+h//2).astype(int)
    j = (rve[0]+w//2).astype(int)
    if bg is not None:
        # surprising but quite a lot quicker to use fill
        rgba.view( np.uint32 ).fill( np.array( bg, np.uint8).view(np.uint32)[0] )
    if colors is None:
        fgpx = np.array( fg, np.uint8 )
        nbsplat_onecolor( rgba, i, j, fgpx, order[:npk], npx )
    else:
        nbsplat_manycolor( rgba, i, j, colors, order[:npk], npx )

@numba.njit(boundscheck=False)
def do_zsort_clip( ar, tmp, xlim, ylim, zsort, order ):
    k = 0
    zmin =  1e20
    zmax = -1e20
    if zsort:
        for i in range(len(ar[0])):
            if abs( ar[0][i] ) < xlim and abs( ar[1][i] ) < ylim:
                tmp[k] = i
                k += 1
                zmin = min(zmin, ar[2][i])
                zmax = max(zmax, ar[2][i])
            npk = k
        s = 255.0/(zmax - zmin)
        ha = np.zeros( 257, numba.int32 )
        v8 = np.empty( ar[2].shape, numba.int32 )
        h = ha[1:]
        for i in tmp[:npk]:
            v8[i] = int(np.round((ar[2][i] - zmin)*s))
            h[v8[i]] = h[v8[i]] + 1
        for i in range(1,len(ha)):
            ha[i] = ha[i] + ha[i-1]
        for i in tmp[:npk]:
            order[ha[v8[i]]] = i
            ha[v8[i]] += 1
    else:
        for i in range(len(ar[0])):
            if abs( ar[0][i] ) < xlim and abs( ar[1][i] ) < ylim:
                order[k] = i
                k += 1
                zmin = min(zmin, ar[2][i])
                zmax = max(zmax, ar[2][i])
        npk = k
    return npk
    

@numba.njit(boundscheck = False )
def nbsplat_onecolor( rgba, i, j, fg, order, npx):
    for k in order:
        rgba[ i[k]-npx:i[k]+npx+1 , j[k]-npx:j[k]+npx+1  ] = fg
    

@numba.njit(boundscheck = False )
def nbsplat_manycolor( rgba, i, j, c, order, npx):
    for k in order:
        rgba[ i[k]-npx:i[k]+npx+1 , j[k]-npx:j[k]+npx+1  ] = c[k]

