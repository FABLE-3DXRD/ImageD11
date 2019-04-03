import ctypes, fabio, os, time, numpy as np, pylab as pl

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

sc.histogram_image.restype = None
sc.histogram_image.argtypes = [
    np.ctypeslib.ndpointer( dtype=np.float32,
                            ndim=1, flags="CONTIGUOUS"),
    ctypes.c_int,
    ctypes.c_float,
    ctypes.c_float,
    np.ctypeslib.ndpointer( dtype=np.int32,
                            ndim=1, flags="CONTIGUOUS"),
    ctypes.c_int,]



if 1:
    fname = "testoverlaps0000.edf"
    im=fabio.open(fname).data.astype(np.float32)
    hh = np.zeros(256, np.int32)
    low = 24
    step = 2
    

else:
    fname = "/data/id11/nanoscope/Commissioning/2017Feb/gold/Au6_s0_048_b/Au6_s0_048_b0432.edf.gz"
    bname = "/data/id11/nanoscope/Commissioning/2017Feb/gold/pks/bkg.edf"
    im=fabio.open(fname).data.astype(np.float32)
    im = im - fabio.open(bname).data.astype(np.float32)
    h,b = np.histogram( im , bins=np.arange(-128,128))
    bc = (b[1:]+b[:-1])/2
    pl.plot( bc, h, "-")
    hh = np.zeros(256, np.int32)
    low = -50
    step = 1

sc.histogram_image( im.ravel(), len(im.ravel()), low, step, hh, len(hh))    
# Notes:
#
#  The cutoff in the iterative approach depends if we come from below or above
#  ... the process was not converging properly for the testoverlaps image
#
#  Alternative idea:
#    Compute the histogram of the image
#    Use the histogram to determine the cutoff (e.g. fit low side of "peak")
#  TODO : use the histogram...
    

pl.plot( np.arange(low,low+len(hh)*step,step), hh, "-")
pl.semilogy()
pl.show()
mask = np.zeros( im.shape, np.int8 )
start = time.time()
sc.threshold_image( im.ravel(), len(im.ravel()), 4., mask.ravel())
end = time.time()
print(end-start)
