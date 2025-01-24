
from __future__ import print_function, division
"""
fast radial regrouping method proposed by Peter Boesecke
"""
import numpy as np, pylab as pl
import timeit, sys
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

if sys.version_info[0] == 2:
    prange = range

# 2018 code:
#   compute radius / azimuth
#        integer bin for radius
#        sort on radial bin, then azimuth (or vice versa)
#
#   Re-order an input array to write output
#    no splits or sums n->n mapping
#
@jit(nopython=True, parallel=sys.version_info[0]>2)
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
print("Setting up, r_chi",t1-t0,t2-t1)

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

def measure( func, name, lut, adrout, f32=False ):
    NFRAME = 500//8 # overfill the cache (500 MB)
    if f32:
        data = [gdata.copy().astype(np.float32) for i in range(NFRAME)]
    else:
        data = [gdata.copy() for i in range(NFRAME)]
    out  = np.zeros_like(data[0])
    w, r = func( lut, adrout, data[NFRAME-1], out )
    dt = []
    for i in range(NFRAME):
        d = data[i]
        t0 = timer()
        func ( lut, adrout, d, out )
        dt.append(timer()-t0)
        if np.sum(dt) > 15:
            break
    dt = np.array(dt)
    t =  np.median(dt)#.mean()
    te = dt.std()
    s = name + "% 6.3f ms (std % 6.2f), write % 7.0f MB/s, read % 7.0f MB/s %g MB"%(
        t*1e3, te*1e3, w/t/1e6, r/t/1e6,r/1e6)
    print(s)
    if f32:
        return out.astype( gdata.dtype )
    else:
        return out

@jit(nopython=True, parallel=sys.version_info[0]>2, nogil=True)
def memcp(a,b):
    for i in prange(a.size):
	    b.flat[i] = a.flat[i]

dest = np.empty_like(gdata)
memcp(gdata,dest)
t = [timer(),]
for i in range(10):
   memcp(gdata, dest)
   t.append(timer())
t = np.array(t)
t = t[1:]-t[:-1]
print("memcpy %.2f s %.2f MB/s"%( t.min(), gdata.nbytes/t.min()/1e6 ))

ao = np.arange( N*M, dtype=np.uint32 )
# lut is sorted for reads. Use cachelen blocks instead
a32 = lut.copy()
a32.shape = a32.size//32 , 32
maxread = a32.max( axis=1 )
run_order = np.argsort( maxread )
print(run_order.shape)
l32 = a32[ run_order ].ravel()
print(l32.shape)
i32 = ao.reshape( a32.shape )[ run_order ].ravel()
assert l32.shape == lut.shape

# lut is sorted for reads. Use cachelen blocks instead
d32 = adrout.copy()
d32.shape = d32.size//32 , 32
maxread = d32.max( axis=1 )
run_order = np.argsort( maxread )
print(run_order.shape)
m32 = d32[ run_order ].ravel()
print(l32.shape)
j32 = ao.reshape( d32.shape )[ run_order ].ravel()
assert j32.shape == lut.shape
assert m32.shape == lut.shape



def fazit_numpy_write( lut, adrout, data, out ):
    out.ravel()[:] = data.ravel()[lut]
    writes = out.itemsize * out.size
    reads =  data.itemsize*data.size + lut.itemsize*lut.size
    return writes, reads

trans = measure( fazit_numpy_write, "Numpy  write ", lut, adrout  )
assert (trans.ravel() == gdata.flat[lut]).all()



@jit(nopython=True, parallel=sys.version_info[0]>2, nogil=True)
def fazit_sortedwrite( lut, adrout, data, out ):
    npx = lut.size
    for i in prange(npx):
        out.flat[i] = data.flat[lut.flat[i]]
    writes = out.itemsize * out.size
    reads =  data.itemsize*data.size + lut.itemsize*lut.size
    return writes, reads
    
trans = measure( fazit_sortedwrite, "Sorted write  ", lut, adrout  )
assert (trans.ravel() == gdata.flat[lut]).all()

@jit(nopython=True, parallel=sys.version_info[0]>2, nogil=True)
def fazit_sortedread( lut, adrout, data, out ):
    npx = lut.size
    for i in prange(npx):
        out.flat[adrout.flat[i]] = data.flat[i]
    writes = out.itemsize * out.size
    reads =  data.itemsize*data.size + adrout.itemsize*adrout.size
    return writes, reads

trans = measure( fazit_sortedread,  "Sorted read   ", lut, adrout  )
assert (trans.ravel() == gdata.flat[lut]).all()


@jit(nopython=True, parallel=sys.version_info[0]>2, nogil=True)
def fazit_both( lut, adrout, data, out ):
    npx = lut.size
    for i in prange(npx):
        out.flat[adrout.flat[i]] = data.flat[lut.flat[i]]
    writes = out.itemsize * out.size
    reads =  data.itemsize*data.size + adrout.itemsize*adrout.size + lut.itemsize*lut.size
    return writes, reads

@jit(nopython=True, parallel=sys.version_info[0]>2, nogil=True)
def fazit_b32out( lut, a32, data, out ):
    npx = lut.size
    for i in prange(0,npx//32):
        p = a32[i]
        for j in range(32):
            out.flat[p+j] = data.flat[lut.flat[i*32+j]]
    writes = out.itemsize * out.size
    reads =  data.itemsize*data.size + a32.itemsize*a32.size + lut.itemsize*lut.size
    return writes, reads





trans = measure( fazit_both,        "sort write 2Id", lut, ao  )
assert (trans.ravel() == gdata.flat[lut]).all()

trans = measure( fazit_both,        "sort reads 2Id", ao, adrout  )
assert (trans.ravel() == gdata.flat[lut]).all()

trans = measure( fazit_both,        "sort blk32r2Id", j32 , m32)
assert (trans.ravel() == gdata.flat[lut]).all()

trans = measure( fazit_both,        "sort blk32o2Id", l32 , i32)
assert (trans.ravel() == gdata.flat[lut]).all()

io32 = i32.reshape(i32.size//32,32)[:,0].copy()
trans = measure( fazit_b32out ,     "sort b32o   Id", l32 , io32)
assert (trans.ravel() == gdata.flat[lut]).all()



@jit(nopython=True, parallel=sys.version_info[0]>2, nogil=True)
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

othr = cImageD11.cimaged11_omp_get_max_threads()
trans= measure( fazit_f2py_lut, "f sortwrite %d"%(othr), lut, adrout )
assert (trans.ravel() == gdata.flat[lut]).all()


def fazit_f2py32( lut, adrout, data, out ):
    cImageD11.reorder_f32_a32( data.ravel(), adrout, out.ravel())
    writes = out.nbytes
    reads = data.nbytes + adrout.nbytes
    return writes, reads

trans= measure( fazit_f2py32, "f32 sort read", lut, adrout, f32=True )
assert (trans.ravel() == gdata.flat[lut]).all()

def fazit_f2py_lut32( lut, adrout, data, out):
    cImageD11.reorderlut_f32_a32( data.ravel(), lut.ravel(), out.ravel())
    writes = out.nbytes
    reads = data.nbytes + lut.nbytes
    return writes, reads

othr = cImageD11.cimaged11_omp_get_max_threads()
trans= measure( fazit_f2py_lut32, "f32 sortwrt %d"%(othr), lut, adrout, f32=True )
assert (trans.ravel() == gdata.flat[lut]).all()

    

for nthr in range(7):
    cImageD11.cimaged11_omp_set_num_threads(2**nthr)
    gthr = cImageD11.cimaged11_omp_get_max_threads()
    trans= measure( fazit_f2py_lut, "f sortwrite %d"%(gthr), lut, adrout )
    assert (trans.ravel() == gdata.flat[lut]).all()

cImageD11.cimaged11_omp_set_num_threads(othr)
gthr = cImageD11.cimaged11_omp_get_max_threads()
trans= measure( fazit_f2py_lut, "f sortwrite %d"%(gthr), lut, adrout )
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
    #    if k  == 1:
    print("highest jump was",abs(da).max())
    trans = measure( fazit_u32_i16,        "u32 i16  %4d "%(k), u32, i16  )
    assert (trans.ravel() == gdata.flat[lut]).all()
    for nthr in range(7):
        cImageD11.cimaged11_omp_set_num_threads(2**nthr)
        gthr = cImageD11.cimaged11_omp_get_max_threads()
        trans = measure( fazit_f2py_i16,        "f2pyi16 %2d%4d"%(k, gthr), u32, i16  )
        assert (trans.ravel() == gdata.flat[lut]).all()

    cImageD11.cimaged11_omp_set_num_threads(othr)
    gthr = cImageD11.cimaged11_omp_get_max_threads()
    trans = measure( fazit_f2py_i16,        "f2pyi16 %2d%4d"%(k, gthr), u32, i16  )
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
