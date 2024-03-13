#!/usr/bin/env python

from __future__ import print_function

"""
Reciprocal space volume mapper
Transfers images into reciprocal space by pixel mapping
"""

import numpy, logging
from ImageD11 import parameters, transform, indexing, \
    cImageD11, blobcorrector, rsv, ImageD11options


class rsv_mapper(object):
    """
    Maps images into a reciprocal space volume

    Basic idea is for a fixed detector, scattering vectors in lab
    frame are fixed. When rotating sample, just rotate these into
    crystal frame and use as array indices into volume
    """
    def __init__(self, dims, pars, ubi, 
                 splinefile = None,
                 np=16, 
                 border = 10, 
                 omegarange = list(range(360)),
                 maxpix = None,
                 mask = None):
        """
        Create a new mapper intance. It will transform images into 
        reciprocal space (has its own rsv object holding the space)
        dims - image dimensions
        par  - ImageD11 parameter filename for experiment
        ubi  - Orientation matrix (ImageD11 style)
        np   - Number of pixels per hkl index [16]
        border - amount to add around edge of images [10]
        omegarange - omega values to be mapped (0->360)
        maxpix - value for saturated pixels to be ignored
        mask - fit2d style mask for removing bad pixels / border
        """
        if len(dims)!=2: raise Exception("For 2D dims!")
        self.dims = dims
        print(dims)
        # Experiment parameters
        if not isinstance( pars, parameters.parameters):
            raise Exception("Pars should be an ImageD11 parameters object")
        for key in ["distance", "wavelength"]: #etc
            assert key in pars.parameters
        self.pars = pars

        # Orientation matrix
        self.ubi = ubi

        # Saturation
        self.maxpix = maxpix

        # Mask
        self.mask = mask
        if self.mask is not None:
            assert self.mask.shape == self.dims, "Mask dimensions mush match image"
        
        # spatial
        if splinefile is None:
            self.spatial = blobcorrector.perfect()
        else:
            self.spatial = blobcorrector.correctorclass(splinefile)

        # npixels
        self.np  = np

        self.uspace = np*ubi

        self.find_vol( border = border, omegarange = omegarange )

        self.rsv.metadata['ubi'] = ubi
        self.rsv.metadata['uspace'] = self.uspace
        # Make and cache the k vectors
        self.make_k_vecs()
        
        
    def find_vol( self, border , omegarange ):
        """
        find limiting volume
        The four image corners over 360 degrees

        np is the number of pixels per hkl

        returns (RSV, NR)
           RSV min/max in reciprocal space (np*ubi).gv
           NR number of points == 1+RSV[i][1]-RSV[i][0]
        """
        # Note that ImageD11 peaks are [slow, fast]
        #  1   2   3
        #  4   5   6    
        #  7   8   9
        p1 = [ -border, self.dims[0]/2, self.dims[0]+border ,
                self.dims[0]/2, self.dims[0]+border ,-border,
                self.dims[0]+border, -border, self.dims[0]/2 ]
        p2 = [ -border, self.dims[1]/2, self.dims[1]+border ,
                -border, self.dims[1]/2, self.dims[1]+border ,
                -border, self.dims[1]/2, self.dims[1]+border ]
        for i in range(9):
            p1[i], p2[i] = self.spatial.correct( p1[i], p2[i] )
        peaks   = [p1*len(omegarange), p2*len(omegarange)]
        om = numpy.array(list(omegarange)*9, numpy.float32)
        tth, eta = transform.compute_tth_eta( peaks,
                                              **self.pars.get_parameters() )
        #print "tth",tth.min(),tth.max()
        #print "eta",eta.min(),eta.max()
        assert om.shape == tth.shape
        gv = transform.compute_g_vectors( tth, eta, om,
                                          self.pars.get('wavelength'),
                                          float(self.pars.get('wedge')),
                                          float(self.pars.get('chi')))
                                          
    
        # Rotate g-vectors into target volume
        hkls = numpy.dot( self.uspace, gv )

        # print "Ranges for RSV"
        bounds=numpy.zeros((3,2))
        npv =numpy.zeros(3)
        for i in range(3):
            bounds[i][0] = numpy.floor(hkls[i].min())
            bounds[i][1] = numpy.ceil(hkls[i].max())
            npv[i] = (bounds[i][1]-bounds[i][0]) + 1
        self.bounds = bounds
        self.rsv = rsv.rsv( npv , bounds=bounds, np=self.np )
        # Cross your fingers and....
        self.rsv.allocate_vol()

                            

    def make_k_vecs( self ):
        """
        Generate the k vectors from the experiment parameters
        given in constructor
        """
        xim, yim = self.spatial.make_pixel_lut( self.dims )
        peaks = [ numpy.ravel(xim), numpy.ravel(yim) ]
        # First, x, is the slow pixel direction, should not change
        # when raveled
        #for i in range(10):
        #    print "slow, fast"
        #    print peaks[0][i],peaks[1][i]
        assert abs(peaks[0][10]-peaks[0][0]) < 3
        # Second, y, is the fast, should change by ~ 1 pixel per pixel
        assert abs(peaks[1][10]-peaks[1][0]-10) < 3
        tth, eta = transform.compute_tth_eta( peaks, 
                                              **self.pars.get_parameters() )
        self.k = transform.compute_k_vectors(tth, eta,
                                             self.pars.get('wavelength'))
                                             
        # FIXME
        # This should be something like domega/dk where
        #      dk [ k(omega=0) - k(omega=1) ]
        
        self.lorfac = numpy.ones( self.dims[0]*self.dims[1],
                                  numpy.float32)
                        


    def add_image( self, om, data ):
        """
        RSV = bounds of reciprocal space vol
        NR = dims of RSV
        k = scattering vector for image at om == 0
        data = 2D image (dims == k.dims)
        SIG = signal
        MON = monitor
        """
        dat = numpy.ravel(data).astype(numpy.float32)
        assert len(dat) == len(self.k[0]), "dimensioning issue"
    
        # hkl = ubi.( gtok. k )
        gvm = transform.compute_g_from_k(numpy.eye(3),
                                         # this transform module sucks
                                         om*float(self.pars.get('omegasign')),
                                         wedge = float(self.pars.get('wedge')),
                                         chi = float(self.pars.get('chi')))
        tmat = numpy.dot( self.uspace , gvm )
        hkls = numpy.dot( tmat , self.k )

        # Find a way to test if we are doing the transform OK
    
        # Bounds checks
#        for i in range(3):
#            assert hkls[i].min() > self.bounds[i][0], \
#                "%d %s %s"%(i, str(hkls[i].min()),str( self.bounds[i][0]))
#            assert hkls[i].max() < self.bounds[i][1], \
#                "%d %s %s"%(i, str(hkls[i].max()),str( self.bounds[i][1]))

        NR = self.rsv.NR        
        # hkls[0] is the slowest index. integer steps of NR[1]*NR[2]
        ind = numpy.floor(hkls[0]+0.5-self.bounds[0][0]).astype(numpy.intp)
        numpy.multiply( ind, NR[1]*NR[2] , ind )
        assert ind.dtype == numpy.intp
        # hkls[1] is faster. Steps by NR[2] only
        numpy.add( ind, NR[2]*numpy.floor(
                hkls[1] + 0.5 - self.bounds[1][0]).astype(numpy.intp),
                   ind )
        numpy.add( ind, numpy.floor(
                hkls[2] + 0.5 - self.bounds[2][0]).astype(numpy.intp),
                   ind )
        #
        #
        #
        if self.maxpix is not None:
            msk =  numpy.where( dat > self.maxpix, 0, 1).astype(numpy.uint8)
        else:
            msk = None

        if self.mask is not None:
            # This excludes saturated pixels
            if msk is None:
                msk = self.mask
            else:
                numpy.multiply(msk, numpy.ravel(self.mask), msk)
                
        # cases:
        #    maxpix only     == msk
        #    mask only       == msk
        #    maxpix and mask == msk
        #    neither      


        if msk is not None:
            numpy.multiply(dat, msk, dat)
            cImageD11.put_incr( self.rsv.SIG,
                              ind,
                              dat )
            cImageD11.put_incr( self.rsv.MON,
                              ind,
                              self.lorfac * msk)
        else:
            cImageD11.put_incr( self.rsv.SIG,
                              ind,
                              dat )
            cImageD11.put_incr( self.rsv.MON,
                              ind,
                              self.lorfac)
        return

    def writevol(self, filename):
        """
        Save the volume in a hdf file
        """
        rsv.writevol( self.rsv, filename )



from ImageD11 import ImageD11_file_series

def get_options(parser):
    """
    Command line interface for making a mapping
    Add our options to a parser object
    """
    parser = ImageD11_file_series.get_options( parser )

    parser.add_argument("-p", "--pars", action="store",
                      dest = "pars", default = None,
                      type=ImageD11options.ParameterFileType(mode='r'),
                      help = "ImageD11 parameter file for experiment")
    
    parser.add_argument("-o", "--output", action="store",
                      dest = "output", default = None,
                      type=ImageD11options.HdfFileType(mode='r'),
                      help = "Name of hdf5 output file")

    parser.add_argument("-s", "--splinefile", action="store", 
                      dest = "spline", default = None,
                      type=ImageD11options.SplineFileType(mode='r'),
                      help = "Name of fit2d spline file for spatial dist")

    parser.add_argument("-u", "--ubifile", action="store", 
                      dest = "ubifile", default = None,
                      type = ImageD11options.UbiFileType(mode='r'),
                      help = "Name of ubi file (first matrix is used)")

    parser.add_argument("-x", "--npixels", action="store", type=int,
                      dest = "npixels", default = 16,
      help = "Number of pixels in reciprocal space map per integer hkl [16]")

    parser.add_argument("-i", "--images", action="store", type=int,
                      dest = "images", default = None,
                      help = "Number of images to process [all]")

    parser.add_argument("-b", "--border", action="store", type=int,
                       dest = "border", default = 10,
                       help = "Border around images to allocate space, px [10]")
    parser.add_argument("-t", "--saturation", action="store", type=float,
                      dest = "maxpix", default = None,
                      help = "Saturation value for excluding pixels")


    #parser.add_argument("-t", "--testcolfile", action="store", type="string",
    #                  dest = "testcolfile", default=None,
    #                  help = "A columnfile to test geometry")

    parser.add_argument("-c", "--subslice", action="store", type=int,
                      dest = "subslice", default=1,
                      help = "Number of omega subslices to repeat images")

    parser.add_argument("--maskfilename", action="store", type=str,
                      dest = "maskfilename", default=None,
                      help = "Mask image (fit2d style)" )
    
    return parser
                      
    
                    

def main():
    """
    A user interface
    """
    import sys, time, os, logging
    start = time.time()

    root = logging.getLogger('')
    root.setLevel(logging.WARNING)

    try:
        from argparse import ArgumentParser
        parser = ArgumentParser()
        parser = get_options( parser )
        options, args = parser.parse_known_args()
    except SystemExit:
        raise
    except:
        parser.print_help()
        print("\nProblem with your options:")
        raise
    
    if options.output is None:
        print("You must supply an output file (-o vol.h5)")
        sys.exit()

    if os.path.exists( options.output ):
        print("I would overwrite your output file",options.output)
        print("If you really want that then delete it first and re-run")
        sys.exit()
    
    try:
        if options.pars is None:
            print("You must supply a parameter file, -p file.pars")
            sys.exit()
        pars = parameters.parameters()
        pars.loadparameters(options.pars)
        print("Got parameters from",options.pars)
        pd = pars.get_parameters()
        names = list(pd.keys())
        names.sort()
        for name in names:
            print("%30s   %s"%(name, pd[name]))
    except:
        print("Problem with parameters:",options.pars)
        raise

    try:
        if options.ubifile is None:
            print("You must supply an input ubifile")
        ubi = indexing.readubis(options.ubifile)[0]
        print("UBI:\n",ubi)
        print("Cell parameters:")
        print("%.5f %.5f %.5f %.4f %.4f %.4f" % \
              indexing.ubitocellpars(ubi))
    except:
        print("Problem with ubi file:",options.ubifile)
        raise

    if options.maskfilename is not None:
        from fabio.openimage import openimage
        try:
            mask = ( openimage( options.maskfilename ).data == 0 )
        except:
            print("Problem with your mask image",options.maskfilename)
            raise
        print("Using a mask from",options.maskfilename)
        print("percent of image used %.3f"%(100.0*mask.sum()/mask.shape[0]/mask.shape[1]))
    else:
        mask = None
        
    
    first_image = True
    nimage = 0

    imagefiles = ImageD11_file_series.get_series_from_options( options, args )

    print("Subslicing by",options.subslice)

    try:
        for fim in imagefiles:
            
            if first_image: # allocate volume, compute k etc
                
                first_image = False
                
                mapper = rsv_mapper( fim.data.shape, 
                                     pars , ubi, 
                                     options.spline,
                                     np = options.npixels,
                                     border = options.border,
                                     maxpix = options.maxpix,
                                     mask = mask
                                     # FIXME omegarange
                                     )
                                
                logging.info( "Setting up time %.4f s"%(time.time()-start))

                
            ltp = time.time()
            om = float(fim.header['Omega'])
            oms = float(fim.header['OmegaStep'])
            for i in range(options.subslice):
                print(".", end=' ')
                omv = om + i*oms/options.subslice
                # ==1 : 0*s/1
                # ==2 : 0*s/2 , 1*s/2
                # ==3 : 0*s/3 , 1*s/3, 2*s/3 etc
                mapper.add_image( omv, fim.data )
                    
                    
            nimage = nimage + 1
            
            print("  %d %.3f %.4f s, %.4f s"%(nimage, om, 
                                            time.time()-ltp,
                                            time.time()-start))
            if options.images is not None:
                if nimage >= options.images:
                    break

    except KeyboardInterrupt:
        print("\nCaught a control-c")
        if nimage > 0:
            print("Problem, trying to save volume so far to:",options.output)
            mapper.writevol( options.output )
            print("Saved what I had")
        sys.exit()
    except:
        print("\nAn error occured")
        if nimage > 0:
            print("Problem, trying to save volume so far to:",options.output)
            mapper.writevol( options.output )
            print("Saved what I had")
        raise

    if nimage > 0:
        mapper.writevol( options.output )
        
            


        
if __name__=="__main__":
    main()
