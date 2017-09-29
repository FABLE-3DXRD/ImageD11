
import pyopencl as cl
import numpy
import time, sys, os
import ctypes


if "win" in sys.platform:
    timer = time.clock
    dll = ctypes.CDLL("diffracCl.dll")
    print("# time.clock", sys.platform, sys.version)
else:
    timer = time.time
    dll = ctypes.CDLL("./diffracCl.so")
    print("# time.time", sys.platform, sys.version)

import pyopencl.version
print("# pyopencl:",pyopencl.version.VERSION)

class myCL:
    def __init__(self, npx):
        self.ctx = cl.create_some_context()
        for d in self.ctx.devices:
            print("#",d.platform.name)
            print("#",d.vendor)
            print("#",d.name)
        self.npx = npx
        self.queue = cl.CommandQueue(self.ctx)
        self.pars = numpy.zeros(14, dtype=numpy.float32)

    def loadProgram(self, filename):
        #read in the OpenCL source file as a string
        f = open(filename, 'r')
        fstr = "".join(f.readlines())
        #create the program
        self.program = cl.Program(self.ctx, fstr).build(
               [ '-cl-fast-relaxed-math' ])
#                '-cl-mad-enable',
#                '-cl-no-signed-zeros',
#                '-w'] )
                

    def popCorn(self):
        mf = cl.mem_flags
        nb = self.npx*self.npx*4
        #create OpenCL buffers
        self.par_buf = cl.Buffer(self.ctx, 
                mf.READ_ONLY | mf.COPY_HOST_PTR,
                hostbuf = self.pars )

        self.tth_buf = cl.Buffer(self.ctx, mf.WRITE_ONLY, nb )
        self.eta_buf = cl.Buffer(self.ctx, mf.WRITE_ONLY, nb )



    def execute(self):
        # start = timer()
        evtcompute = self.program.tthetaf4(self.queue, 
                (self.npx/4, self.npx),
                None, 
                self.tth_buf, 
                self.eta_buf,
                self.par_buf)
        #evtcompute.wait()
        #print timer()-start
        self.tthl,evtt = cl.enqueue_map_buffer( self.queue,
                    self.tth_buf,
                    cl.map_flags.READ,
                    0,
                    (self.npx,self.npx),
                    numpy.float32,
                    'C',
                    wait_for = [evtcompute],
                    is_blocking=False)
        self.etal,evte = cl.enqueue_map_buffer( self.queue,
                    self.eta_buf,
                    cl.map_flags.READ,
                    0,
                    (self.npx,self.npx),
                    numpy.float32,
                    'C',
                    wait_for = [evtcompute],
                    is_blocking=False)
        evtcompute.wait()
        evtt.wait()
        evte.wait()
        return self.tthl, self.etal

    def setpars(self, pars ):
        self.pars = pars
        # print pars
        cl.enqueue_write_buffer( self.queue,
                         self.par_buf,   
                         self.pars).wait()




class ctp:
    def __init__(self, npx):
        self.npx = npx
        if "win" in sys.platform:
            fname = "diffracCl.dll"
        else:
            fname = "./diffracCl.so"
        self.dll = ctypes.CDLL( fname )
        self.dll.ttheta.argtypes=[ctypes.POINTER(ctypes.c_float),
                             ctypes.POINTER(ctypes.c_float),
                             ctypes.POINTER(ctypes.c_float),
                             ctypes.c_int,
                             ctypes.c_int]

    def compute( self, tth, eta, p):
        t = tth.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        e = eta.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        p = pars.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        times = []
        # warmup
        self.dll.ttheta(t , e , p, self.npx, self.npx )
        for i in range(10):
            start = timer()
            self.dll.ttheta(t , e , p, self.npx, self.npx )
            times.append( timer() - start )
        return times

def make_pars( parfile ):
    from ImageD11.parameters import parameters
    from ImageD11 import transform
    p = parameters( )
    p.loadparameters( parfile )
    rmat = transform.detector_rotation_matrix(
            float(p.parameters["tilt_x"]),
            float(p.parameters["tilt_y"]),
            float(p.parameters["tilt_z"]) )
    fmat = numpy.array( [[ 1 , 0 , 0],
          [ 0 , float(p.parameters['o11']), float(p.parameters['o12']) ],
          [ 0 , float(p.parameters['o21']), float(p.parameters['o22']) ]],
          numpy.float32)
    pars = numpy.array( [ float(p.parameters["y_center"]) ,
                          float(p.parameters["y_size"]) ,
                          float(p.parameters["z_center"]) ,
                          float(p.parameters["z_size"]) ,
                          float(p.parameters["distance"]) ] +
                          list( numpy.dot( rmat, fmat).ravel() ) ,
                                                      numpy.float32)
    return pars



if __name__ == "__main__":
    start = timer()
    npx = 4096
    example = myCL(npx)
    example.loadProgram("diffracCl.cl")
    example.popCorn()
    pars = make_pars(sys.argv[1]) 
    example.setpars( pars )
    print("# Init", timer()-start)
    times = []
    # Warmup
    tth, eta = example.execute()
    for i in range(100):
        start = timer()
        tth, eta = example.execute()
        times.append(timer()-start)
    times = numpy.array(times)
    print("# mean    min    max     std")
    print("%.4f  %.4f  %.4f  %.4f"%( times.mean(), times.min(),
            times.max(), times.std()))
    t = numpy.median(times)
    print("%.1f ms, %.1f fps,"%(1e3*t,1.0/t), end=' ')
    print(tth.max(),tth.min())
    eta_cl = eta.copy()
    eta_ct = eta.copy()
    tth_cl = tth.copy()
    tth_ct = tth.copy()
    o = ctp( npx )
    times = numpy.array( o.compute( tth_ct, eta_ct, pars ) )

    print("# ctypes module, hopefully with openmp")
    print("# mean    min    max     std")
    print("%.4f  %.4f  %.4f  %.4f"%( times.mean(), times.min(),
            times.max(), times.std()))
    t = numpy.median(times)
    print("%.1f ms, %.1f fps,"%(1e3*t,1.0/t), end=' ')
    print(tth.max(),tth.min())

    # Check same ness
    eta_err = (abs(eta_cl - eta_ct)).mean()
    tth_err = (abs(tth_cl - tth_ct)).mean()
    print("Mean diff tth,eta",tth_err,eta_err)
    if len(sys.argv)>2:
     from matplotlib.pylab import imshow, figure, show, colorbar, title
     figure(1)
     title("OpenCL")
     imshow(eta_cl)
     colorbar()
     figure(2)
     title("Ctypes")
     imshow(eta_ct)
     colorbar()
     show()



