
"""
fast radial regrouping method proposed by Peter Boesecke
"""
import numpy as np, pylab as pl
import timeit
from ImageD11 import cImageD11

timer = timeit.default_timer


try:
    from numba import jit, prange
    GOTNUMBA = True
    print("Got numba")
except:
    print("No numba")
    GOTNUMBA = False
    def jit(pyfunc=None, **kwds):
        def wrap(func):
            return func
        if pyfunc is not None:
            return wrap(pyfunc)
        else:
            return wrap
    def prange(*a):
        return range(*a)

# 2018 code:
#   compute radius / azimuth
#        integer bin for radius
#        sort on radial bin, then azimuth (or vice versa)
#
#   Re-order an input array to write output
#    no splits or sums n->n mapping
#
@jit(nopython=True, parallel=True)
def compute_r_chi( dims, pars ):
    """
    Replace this with pyFAI / ImageD11 / whatever
      dims = array dimension
      pars = center in slow/fast
    returns r/chi (float)
    """
    r = np.empty( dims, np.float32 )
    c = np.empty( dims, np.float32 )
    p0 = pars[0]
    p1 = pars[1]
    ns = dims[0]
    nf = dims[1]
    for i in prange( ns ):
      di = i - p0
      for j in range( nf ):
        dj = j - p1
        r[i,j] = np.sqrt( di*di + dj*dj )
        c[i,j] = np.arctan2( di, dj )
    return r, c




N = 2048
M = 2048
dims = ( N, M )
pars = 1022.1, 1018.5, 1.25
t0 = timer()
rad, chi = compute_r_chi( dims, pars )
t1 = timer()
rad, chi = compute_r_chi( dims, pars )
t2 = timer()
print(t1-t0,t2-t1)

arg = np.round( rad.ravel() / pars[2] ) + chi.ravel()/2/np.pi
start=timer()
lut = np.argsort(arg).astype(np.uint32)
adrout = np.argsort(lut).astype(np.uint32)
# required to save
#   length of radial rows, r, IR and chi value for each
# TODO: optimise the order to apply lut
#      assert len(set(lut)) == len(lut)

gdata = (np.random.random(  (N,M)  )*256).astype( np.uint16 )
for i in range(0,N,256):
    gdata[i:i+120] += 256

p = np.empty( (N,M), np.float32 )
pk = 50
while pk < rad.max():
    dr = ((rad - pk)**2)/20.
    p += 500*np.exp( -dr )
    pk = pk + 200 + pk/10.
gdata += p.astype(np.uint16)

def measure( func, name, lut, adrout ):
    data = gdata.copy()
    out  = np.zeros_like(data)
    t = [timer(),]
    w, r = func( lut, adrout, data, out )
    t.append( timer() )
    for i in range(50):
        func ( lut, adrout, data, out )
        t.append(timer())
        if(t[-1]-t[1]) > 0.5:
            break
    t = np.array(t)
    dt = t[1:] - t[:-1]
    # print(dt)
    t = np.min(dt)
    te = dt[1:].std()
    s = name + "% 6.3f ms (std %.3f), write % 7.0f MB/s, read % 7.0f MB/s %g MB"%(
        t*1e3, te*1e3, w/t/1e6, r/t/1e6,r/1e6)
    print(s)
    return out



def fazit_numpy_write( lut, adrout, data, out ):
    out.ravel()[:] = data.ravel()[lut]
    writes = out.itemsize * out.size
    reads =  data.itemsize*data.size + lut.itemsize*lut.size
    return writes, reads

trans = measure( fazit_numpy_write, "Numpy  write ", lut, adrout  )
assert (trans.ravel() == gdata.flat[lut]).all()



@jit(nopython=True, parallel=True, nogil=True)
def fazit_sortedwrite( lut, adrout, data, out ):
    npx = lut.size
    for i in prange(npx):
        out.flat[i] = data.flat[lut.flat[i]]
    writes = out.itemsize * out.size
    reads =  data.itemsize*data.size + lut.itemsize*lut.size
    return writes, reads
    
trans = measure( fazit_sortedwrite, "Sorted write  ", lut, adrout  )
assert (trans.ravel() == gdata.flat[lut]).all()

@jit(nopython=True, parallel=True, nogil=True)
def fazit_sortedread( lut, adrout, data, out ):
    npx = lut.size
    for i in prange(npx):
        out.flat[adrout.flat[i]] = data.flat[i]
    writes = out.itemsize * out.size
    reads =  data.itemsize*data.size + adrout.itemsize*adrout.size
    return writes, reads

trans = measure( fazit_sortedread,  "Sorted read   ", lut, adrout  )
assert (trans.ravel() == gdata.flat[lut]).all()


@jit(nopython=True, parallel=True, nogil=True)
def fazit_both( lut, adrout, data, out ):
    npx = lut.size
    for i in prange(npx):
        out.flat[adrout.flat[i]] = data.flat[lut.flat[i]]
    writes = out.itemsize * out.size
    reads =  data.itemsize*data.size + adrout.itemsize*adrout.size + lut.itemsize*lut.size
    return writes, reads

ao = np.arange( N*M, dtype=np.uint32 )
trans = measure( fazit_both,        "sort out both ", lut, ao  )
assert (trans.ravel() == gdata.flat[lut]).all()

trans = measure( fazit_both,        "sort in  both ", ao, adrout  )
assert (trans.ravel() == gdata.flat[lut]).all()

@jit(nopython=True, parallel=True, nogil=True)
def fazit_u32_i16( u32, i16, data, out ):
    ns = i16.shape[0]
    nf = i16.shape[1]
    for i in prange( ns ):
        a0 = u32[i]
        i0 = i*nf
        for j in range( nf ):
            a0 = a0 + i16.flat[i0 + j]
            out.flat[a0] = data.flat[i0 + j]
    writes = out.itemsize * out.size
    reads =  data.itemsize*data.size + u32.itemsize*u32.size + i16.itemsize*i16.size
    return writes, reads

trans = measure( fazit_sortedread,  "Sorted read   ", lut, adrout  )
assert (trans.ravel() == gdata.flat[lut]).all()


def fazit_f2py( lut, adrout, data, out):
    cImageD11.reorder_u16_a32( data.ravel(), adrout, out.ravel())
    writes = out.nbytes
    reads = data.nbytes + adrout.nbytes
    return writes, reads

trans= measure( fazit_f2py, "f2py sort read", lut, adrout )
assert (trans.ravel() == gdata.flat[lut]).all()

def fazit_f2py_lut( lut, adrout, data, out):
    cImageD11.reorderlut_u16_a32( data.ravel(), lut.ravel(), out.ravel())
    writes = out.nbytes
    reads = data.nbytes + lut.nbytes
    return writes, reads

trans= measure( fazit_f2py_lut, "f2py sortwrite", lut, adrout )
assert (trans.ravel() == gdata.flat[lut]).all()


def fazit_f2py_i16( u32, i16, data, out):
    cImageD11.reorder_u16_a32_a16( data.reshape(i16.shape),
                                   u32,
                                   i16,
                                   out.reshape(i16.shape) )
    writes = out.nbytes
    reads = data.nbytes + u32.nbytes + i16.nbytes
    return writes, reads


k = 1
while k<2048:
    try:
        a =  adrout.reshape(N*k,M//k)
    except:
        k *=2
        continue
    u32 = a[:,0].copy().astype(np.int32)
    i16 = np.zeros(a.shape, np.int16)
    da = a[:,1:].astype(int) - a[:,:-1]
    i16[:,1:] = da
    assert np.alltrue( i16[:,1:] == da ) , "clipping issue"
    if k == 1:
        print("highest jump was",abs(da).max())
    trans = measure( fazit_u32_i16,        "u32 i16  %4d "%(k), u32, i16  )
    assert (trans.ravel() == gdata.flat[lut]).all()

    trans = measure( fazit_f2py_i16,        "f2pyi16  %4d "%(k), u32, i16  )
    assert (trans.ravel() == gdata.flat[lut]).all()
  
    k *= 2


if 0:
    pl.subplot(221)
    pl.imshow(arg.reshape(N,N))
    pl.title("radius")
    pl.colorbar()
    pl.subplot(223)
    pl.imshow(chi)
    pl.title("azimuth")
    pl.colorbar()
    pl.subplot(222)
    pl.imshow(lut.reshape(dims),aspect='auto',
           interpolation='nearest')
    pl.title("Source address to come from")
    pl.subplot(224)
    pl.imshow(adrout.reshape(dims),aspect='auto',
           interpolation='nearest')
    pl.title("Output address to go to")
    pl.figure()
    pl.imshow(gdata,interpolation='nearest' ,aspect='auto')
    pl.figure()
    pl.imshow( trans.reshape(N,N).T, interpolation='nearest', aspect='auto' )
    pl.show()
