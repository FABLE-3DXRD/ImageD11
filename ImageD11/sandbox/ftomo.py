
from __future__ import print_function
import math, numpy, time
from ImageD11 import cImageD11
from fabio.openimage import openimage
print("Using class version")




class fourier_radial(object):
    """ Cache results for re-use where possible on next layer """
    def __init__(self, dims, theta=None):
        self.dims = dims
        self.theta = theta
        if self.theta is not None:
            assert len(self.theta) == dims[1]
            self.theta = numpy.array(theta)*numpy.pi/180.0
            self.make_indices()

    def set_theta(self, theta):
        """
        th is the list of angles in degrees of the projections
        assumed 1 degree steps otherwise
        """
        self.theta = theta
        assert len(self.theta) == self.dims[1]
        self.theta = numpy.array(theta)*numpy.pi/180.0
        self.make_indices()

    def make_indices(self):
        arshape = self.dims[0]/2+1, self.dims[1]
        nv = (arshape[0]-1)*2
        nh = arshape[0]
        print("NV,NH",nv,nh)
        self.ftimshape = (nv, nh)
        self.ftimlen = nv*nh
        n1 = (self.dims[0]/2+1)*self.dims[1]
        xv   = numpy.arange(0, self.dims[0]/2+1, 1, 
                            dtype=numpy.float32 )
        # dimensions?
        cth = numpy.cos( self.theta ) # 1D
        sth = numpy.sin( self.theta ) # 1D
        ia = numpy.round(numpy.outer( cth, xv )).astype(numpy.int)
        ja = numpy.round(numpy.outer( sth, xv )).astype(numpy.int)
        on = numpy.array([1.0],numpy.float32)
        jm = numpy.where(ja < 0, -on, on)
        numpy.multiply( ia, jm, ia ) #  if j<0: i=-i
        numpy.multiply( ja, jm, ja ) #  if j<0: j=-j
        #  if j<0: f=f.conj()
        ia = numpy.where( ia < 0, nv+ia, ia)
        inds = (ia*nh + ja).ravel()
        self.conjer = jm
        self.inds = inds
        nim  = numpy.zeros( ( nv* nh), numpy.float32 )
        wons = numpy.ones( (len(inds)), dtype=numpy.float32 )
        # This is now more dense - bincount?
        cImageD11.put_incr( nim , inds, wons )
        nim = nim.astype(numpy.int)
        self.nim_div = nim + (nim==0)
        

    def process_sinogram( self,
                          sinogram, 
                          do_interpolation=False):

        """
        sinogram is from the data
            dimensions [npixels, nangles]
            do_interp - tries to fill in some of the missing data in
            fourier space
        returns the radon transform
        """
        assert sinogram.shape == self.dims
        ar = numpy.fft.rfft(sinogram, axis=0)
        faprojr = (ar.T.real.astype(numpy.float32))
        faprojc = (ar.T.imag.astype(numpy.float32))
        numpy.multiply( faprojc, self.conjer, faprojc)
        fimr = numpy.zeros( self.ftimlen , numpy.float32 )
        fimc = numpy.zeros( self.ftimlen , numpy.float32 )
        cImageD11.put_incr( fimr, self.inds, faprojr.ravel())
        cImageD11.put_incr( fimc, self.inds, faprojc.ravel())
        fim = fimr + fimc*1j
        fim = numpy.divide( fimr + fimc*1j, self.nim_div)
        fim.shape = self.ftimshape 
        return fim

    def sino2im(self, sinogram, centrepixel ):
        # Take out high frequency in mean (some ring artifacts)
        s = sinogram
        cp = centrepixel
        d  = numpy.concatenate( (  s[cp:,:], s[:cp,:], ), axis=0)
        im = self.process_sinogram( d , centrepixel )    
        # figure out what the scale factor really is
        ret = numpy.fft.irfft2( im ) * im.shape[0] * im.shape[1]
        ret = numpy.fft.fftshift( ret )
        return ret


if __name__=="__main__":
    import sys
    if len(sys.argv) != 5:
        print("Usage: sinogram startangle step centrepixel")
        sys.exit()
    
    fname = sys.argv[1]
    star = time.time()
    sino = openimage( fname )
    na,nx = sino.data.shape
    start = float(sys.argv[2])
    step  = float(sys.argv[3])
    
    centrepixel = int( sys.argv[4] )
    end   = na*step + start
    print("start, step, end",start, step, end)
    angles = numpy.arange(start, end, step)
    assert len(angles) == na,"%d %d  ... %d"%(nx,na,len(angles))
    print("%.2f setup"%(time.time()-star))
    d =  sino.data.T[:1200]
    o = fourier_radial( d.shape, angles )
    start = time.time()
    im = o.sino2im( d, centrepixel )
    sino.data = im
    sino.write(fname+"_r", force_type=numpy.float32)
    import pylab
    pylab.imshow(im, interpolation='nearest',aspect='auto')
    pylab.show()
    print("per image",time.time()-start)

