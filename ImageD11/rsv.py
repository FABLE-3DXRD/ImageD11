
"""
Reciprocal Space Volume

This class is just for holding a volume and loading/saving it to disk
Another class (rsv_mapper) will add images into it.
"""


# ImageD11_v1.0 Software for beamline ID11
# Copyright (C) 2005-2007  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  0211-1307  USA


import logging, numpy, h5py

class rsv(object):
    """
    A reciprocal space volume
    """
    def __init__(self, dimensions, bounds, np, **kwds ):
        """
        dimensions = NX*NY*NZ grid for the space
        uspace = a 3x3 matrix describing the grid
         eg:
            pixel at vol[i,j,k] comes from:
               [i,j,k] = (uspace).gvec
          gvec is the scattering vector in reciprocal space
          so uspace are vectors in real space 
        uorigin = a 3 vector giving the position of the [0,0,0] pixel
        """
        assert len(dimensions) == 3
        self.SIG = None # signal
        self.MON = None # monitor
        self.NR =  dimensions # dimensions      
        self.NORMED = None
        self.bounds = bounds  # boundary in reciprocal space
        self.np = np          # px per hkl
        self.metadata = kwds
        self.allocate_vol()

    def allocate_vol( self ):
        """ 
        Allocates memory for a volume data
        """
        if self.NR is None:
            raise Exception("Cannot allocate rsv")
        total = self.NR[0]*self.NR[1]*self.NR[2]
        print "rsv: memory used = %.2f MB"%(total*8.0/1024/1024)
        self.SIG = numpy.zeros( total, numpy.float32 )
        self.MON = numpy.zeros( total, numpy.float32 )
        
        
    def normalise( self ):
        """
        Return the normalised but avoid divide by zero
        """
        self.NORMED = numpy.where( self.MON > 0.1, 
                              self.SIG/(self.MON+1e-32), 
                              0.0)

    plnames = {
        0 :0, "h":0, "H":0,
        1 :1, "k":1, "K":1,
        2 :2, "l":2, "L":2,
        }
        

    def slice(self, plane, num):
        """
        return signal on plane index num 
        """
        if not self.plnames.has_key(plane):
            raise Exception("Plane should be one of %s"%(
                        str(self.plnames)))
        p = self.plnames[plane]
        # floor(x+0.5) is nearest integer
        ind = numpy.floor( num * self.np + 0.5) - self.bounds[p][0]
        # convert this back to num
        testnum = 1.0*(self.bounds[p][0] + ind )/self.np
        if abs(testnum - num)>1e-6:
            logging.info("Nearest plane to %f is %f"%(num, testnum))
        if ind < 0 or ind >= self.NR[p]:
            print ind,num,self.np,self.bounds
            raise Exception("slice is out of volume bounds")
        if self.NORMED is None:
            self.normalise()
        if p==0:
            return self.NORMED.reshape(self.NR)[ind, :, :]
        if p==1:
            return self.NORMED.reshape(self.NR)[:, ind, :]
        if p==2:
            return self.NORMED.reshape(self.NR)[:, :, ind]


        
        


def writevol(vol, filename):
    """
    Write volume in vol to filename
    
    Compress -1 is for the zeros, which there might be a lot of
    """
    if not isinstance( vol, rsv ):
        raise Exception("First arg to writevol should be an rsv object")
    if None in [vol.NR, vol.SIG, vol.MON]:
        raise Exception("Cannot save rsv, has not data in it")
    volout = h5py.File( filename,"w")
    if vol.SIG.dtype != numpy.float32:
        logging.warning("rsv SIG was not float32, converting")
        vol.SIG = vol.SIG.astype(numpy.float32)
    volout.create_dataset( "signal",
                           (vol.NR[0],vol.NR[1],vol.NR[2]),
                           vol.SIG.dtype,
                           data = vol.SIG,
                           compression = 'gzip', 
                           compression_opts = 1)
    if vol.MON.dtype != numpy.float32:
        logging.warning("rsv MON was not float32, converting")
        vol.MON = vol.MON.astype(numpy.float32)
    volout.create_dataset( "monitor",
                           (vol.NR[0],vol.NR[1],vol.NR[2]),
                           vol.MON.dtype,
                           data = vol.MON,
                           compression = 'gzip',
                           compression_opts = 1)
    volout.attrs['bounds'] = vol.bounds
    volout.attrs['np'] = vol.np
    for key, value in vol.metadata.iteritems():
        volout.attrs[key]=value
    volout.flush()
    volout.close()


def readvol(filename):
    """
    Read volume from a file
    returns an rsv object
    """
    volf = h5py.File(filename)
    try:
        mon = volf['monitor']
        sig = volf['signal']
    except:
        raise Exception("Your file %s is not an rsv"%(filename))
    mona= mon[:,:,:]
    siga= sig[:,:,:]
    assert mona.shape == siga.shape
    assert len(mona.shape)==3
    bounds = volf.attrs['bounds']
    np = volf.attrs['np']
    vol = rsv( mona.shape, bounds, np )
    vol.MON = mona
    vol.SIG = siga
    for name, value in volf.attrs.iteritems():
        vol.metadata[name] = value
    volf.close()
    return vol




