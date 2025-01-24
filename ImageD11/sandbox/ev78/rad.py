

import math, numpy, time
import h5py, sys


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
        ia = numpy.round(numpy.outer( cth, xv )).astype(numpy.int).ravel()
        ja = numpy.round(numpy.outer( sth, xv )).astype(numpy.int).ravel()
        on = numpy.ones(ja.shape,numpy.float32)
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
        # closest.put_incr( nim , inds, wons )
        nim[inds] += wons
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
        faprojr = (ar.T.real.astype(numpy.float32)).ravel()
        faprojc = (ar.T.imag.astype(numpy.float32)).ravel()
        numpy.multiply( faprojc, self.conjer, faprojc)
        fimr = numpy.zeros( self.ftimlen , numpy.float32 )
        fimc = numpy.zeros( self.ftimlen , numpy.float32 )
        #closest.put_incr( fimr, self.inds, faprojr)
        #closest.put_incr( fimc, self.inds, faprojc)
        fimr[self.inds] += faprojr
        fimc[self.inds] += faprojc
        fim = fimr + fimc*1j
        fim = numpy.divide( fimr + fimc*1j, self.nim_div)
        fim.shape = self.ftimshape 
        return fim

    def sino2im(self, sinogram, centrepixel ):
        # Take out high frequency in mean (some ring artifacts)
        d = numpy.zeros(self.dims,numpy.float32)
        n = sinogram.shape[0]
        cp = centrepixel
        d[:n-cp] = sinogram[cp:] # end
        d[-cp:] = sinogram[:cp] # beginning
        im = self.process_sinogram( d , centrepixel )    
        # figure out what the scale factor really is
        ret = numpy.fft.irfft2( im ) * im.shape[0] * im.shape[1]
        ret = numpy.fft.fftshift( ret )
        return ret




def find_offset(a1, a2):
    assert len(a1) == len(a2)
    c = numpy.convolve(a1-a1.mean(),a2-a2.mean())
    if 0:
        import pylab
        pylab.plot(c)
        pylab.figure()
        pylab.plot(a1)
        pylab.plot(a2[::-1])
        pylab.show()
        print(numpy.argmax(c)/2.0)
        1/0
    return numpy.argmax(c)/2.0

def savearray(ar,name,grp):
    if name not in grp:
        dataset = grp.create_dataset(
            name=name,
            shape=ar.shape,
            dtype="float32")
    else:
        dataset = grp[name]
    dataset[:] = ar[:]

def difftomo(N=256):
    h = h5py.File(sys.argv[1])
    # sinograms
    s = h["DiffTomo/NXdata/sinogram"][:]
    # angles
    a = h["DiffTomo/NXdata/xaxis"][:]
    # use a 256 grid (significantly bigger than the data)
    o = fourier_radial( (N, len(a)), a)
    # Look for a strong peak
    ipk = s[s.shape[0]/2,s.shape[1]/2].argmax()
    offset = find_offset( s[:,0,ipk], s[:,-1,ipk] )
    ioff = int(numpy.round(offset))
    print("Offset (centre pixel) seems to be",offset,ioff)
    recon = numpy.zeros((N,N,s.shape[2]),numpy.float32)
    i0 = numpy.zeros(s.shape[2],numpy.float32)
    for i in range(s.shape[2]):
        if i%10 == 0:
            sys.stdout.write("%4d\r"%(i))
            sys.stdout.flush()
        i0[i] = numpy.median( s[-1,:,i]  )
        recon[:,:,i] = o.sino2im( s[:,:,i] - i0[i] , ioff )
        if 0:
            import pylab
            pylab.ion()
            pylab.figure(1)
            pylab.clf()
            pylab.subplot(121)
            pylab.title(str(i))
            pylab.imshow(recon[:,:,i])
            pylab.colorbar()
            pylab.subplot(122)
            pylab.imshow(s[:,:,i] - i0[i])
            pylab.colorbar()
            input("next?")
    grp = h["DiffTomo/NXdata"]
    savearray( i0, "recon_bg", grp )
    savearray( recon, "recon", grp)
    h.close()
    
        
    

if __name__=="__main__":
    difftomo()
