
import numpy as np

def localmax( im, out=None ):
    neighbours = [ (-1,-1),(-1,0), (-1,1),
                   (0,-1),         (0,1),
                   (1,-1), (1,0),  (1,1) ]
    out = np.arange( im.shape[0]*im.shape[1], dtype=np.int )
    out.shape = im.shape
    out0 = out.copy()
    for i in range(1,im.shape[0]-1):
#        print i,
        for j in range(1,im.shape[1]-1):
            mx = im[i,j]
            for k,l in neighbours:
                if im[i+k,j+l] > mx:
                    out[i,j] = out0[i+k,j+l]
                    mx = im[i+k,j+l]
    return out

import ctypes
l=ctypes.CDLL("lmw.so")
_lmw = l.lmw
_lmw.argtypes = [ ctypes.c_void_p,
                 ctypes.c_void_p,
                  ctypes.c_float,
                 ctypes.c_int,
                 ctypes.c_int ]
_lmw.restype = ctypes.c_int

def lmw( im, threshold, out=None):
    im = im.astype(np.float32)
    out = np.zeros( im.shape, np.int32 )
    n = _lmw( im.ctypes.data, out.ctypes.data,
              threshold,
              im.shape[0], im.shape[1])
    return out,n


if __name__=="__main__":
    import fabio, pylab as pl, time
    import scipy.ndimage
    im = fabio.open("/data/id11/jon/1607/spottysucr/spottysucr0000.edf").data.astype(np.float32).copy()
    im = scipy.ndimage.gaussian_filter(im,2)
    print(im.shape)
    start=time.time()
    l,n = lmw( im, 0 )
    print(time.time()-start,n)
#    t = localmax( im )
    r = np.reshape(np.arange( im.shape[0]*im.shape[1], dtype = np.int),
                   im.shape)
    pl.imshow((l*199)%255)
    pl.colorbar()
#    pl.figure()
#    pl.imshow(t-r)
#    pl.figure()
#    pl.imshow(l-t)
#    pl.colorbar()
    pl.show()
    
