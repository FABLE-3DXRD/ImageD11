import ctypes, fabio, os, time, numpy as np

def compile(x):
    if os.path.exists("%s.so"%(x)):
        os.remove("%s.so"%(x))
    os.system("gcc -shared -O3 -fPIC %s.c -o %s.so"%(x,x))
    
x = "../src/sigma_cut"
compile(x)
sc=ctypes.CDLL("%s.so"%(x))
sc.threshold_image.restype = None
sc.threshold_image.argtypes = [
    np.ctypeslib.ndpointer( dtype=np.float32,
                            ndim=1, flags="CONTIGUOUS"),
    ctypes.c_int,
    ctypes.c_float,
    np.ctypeslib.ndpointer( dtype=np.int8,
                            ndim=1, flags="CONTIGUOUS") ]


im=fabio.open("testoverlaps0000.edf").data.astype(np.float32)
mask = np.zeros( im.shape, np.int8 )
start = time.time()
sc.threshold_image( im.ravel(), len(im.ravel()), 3., mask.ravel())
end = time.time()
print(end-start)
