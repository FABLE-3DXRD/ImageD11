
from __future__ import print_function


"""
Try to compute the radius/arc images required for the
fast radial regrouping method proposed by Peter Boesecke

We wish to have the x/y displacement images for the data
image to put it into an output array as if it were a spatial
distortion

(Uses the approach of the rsv_mapper method which does 3D
regrouping)

Jon Wright 18/11/2009
"""

from ImageD11 import parameters, blobcorrector, transform
from fabio import edfimage

import numpy, math


class xydisp:
    required_pars = ["wavelength", "distance", # etc
                     ]
                     
                     
    def __init__(self, splinefile = None,
                 parfile = None):
        """
        splinefile = fit2d spline file, or None for images
        that are already corrected

        parfile = ImageD11 parameter file. Can be fitted
        using old ImageD11_gui.py or newer fable.transform
        plugin.
        """
        self.splinefile = splinefile
        
        if self.splinefile is None:
            self.spatial = blobcorrector.perfect()
        else:
            self.spatial = blobcorrector.correctorclass( splinefile )

        self.parfile = parfile
        self.pars = parameters.parameters()
        self.pars.loadparameters( parfile )
        for key in self.required_pars:
            
            if key not in self.pars.parameters:
                raise Exception("Missing parameter "+str(key))

    def compute_tth_eta(self, dims):
        """
        Find the twotheta and azimuth images
        """
        assert len(dims) == 2
        xim, yim = self.spatial.make_pixel_lut( dims )
        self.dims = dims
        peaks =  [ numpy.ravel(xim), numpy.ravel(yim) ]
        tth, eta =  transform.compute_tth_eta( peaks, 
                                               **self.pars.get_parameters() )
        assert len(tth) == dims[0]*dims[1]
        assert len(eta) == dims[0]*dims[1]

        # Now we have the twotheta and azimuth images in memory
        # they are in degrees
        
        self.tth = numpy.reshape( tth, dims )
#        self.eta = numpy.mod(numpy.reshape( eta, dims ), 360)-180
        self.eta = numpy.reshape( eta, dims )-eta.mean()
        self.compute_rad_arc()
        
    def compute_rad_arc(self):
        """
        This part needs more work - how to properly define the output
        after spd is patched
        For now we aim to make a triangle filling the same array size

        # FIXME - this remains unclear to me, what are the output
        units supposed to be??
        """
        tth_rad = self.tth * math.pi / 180.0
        eta_rad = self.eta * math.pi / 180.0
        arclength = tth_rad * eta_rad

        # x-axis, eg [0], is tth
        tthmax =  numpy.max( self.tth )
        tthmin =  numpy.min( self.tth )
        tthstep = (tthmax - tthmin)/(self.dims[0] - 1)
        self.tthbin = numpy.floor( (self.tth - tthmin)/tthstep )
        self.tthvals = numpy.arange(tthmin,tthmax+tthstep*0.5,tthstep)
        # Ideally we want the arc bins to vary with tth?
        #arcmax = numpy.max( arclength )
        #arcmin = numpy.min( arclength )
        # 4 corners of image
        arcmin = arclength.min()
        arcmax = arclength.max()
        #from matplotlib.pylab import imshow,show, colorbar
        #imshow(arclength)
        #colorbar()
        #show()
        arcstep = (arcmax - arcmin)/(self.dims[1] - 1)
        arcmid  = 0.5*(arcmax+arcmin)
        # Make integer pixel id images

        self.arcbin = numpy.floor((arclength - arcmid)/arcstep )+self.dims[1]/2

        assert self.tthbin.min() >= 0
        assert self.tthbin.max() < self.dims[0], self.tthbin.max()
        assert self.arcbin.min() >= 0, self.arcbin.min()
        assert self.arcbin.max() < self.dims[1]

        # Now convert these to displacements compared to input image
        # Use the same code as for the spline case to get the "x/y" images
        ideal = blobcorrector.perfect()
        idealx, idealy = ideal.make_pixel_lut( self.dims )
        self.dx = self.tthbin - idealx
        self.dy = self.arcbin - idealy

    def write(self, stemname):
        """
        save the dx, dy images
        """
        im = edfimage.edfimage()
        im.data = self.dx
        im.write( "%s_dx.edf"%(stemname), force_type = numpy.float32)
        im = edfimage.edfimage()
        im.data = self.dy
        im.write( "%s_dy.edf"%(stemname), force_type = numpy.float32)
        numpy.save("%s_tth.npy"%(stemname),self.tthvals)

def get_options(parser):

    parser.add_option("-p", "--pars", action="store",type="string",
                      dest = "pars", default = None,
                      help = "ImageD11 parameter file for experiment")
    
    parser.add_option("-o", "--output", action="store", type="string",
                      dest = "output", default = None,
                      help = "stem name for output x/y edf images")

    parser.add_option("-s", "--splinefile", action="store", type="string",
                      dest = "spline", default = None,
                      help = "Name of fit2d spline file for spatial dist")

    parser.add_option("--nf", action="store", type="int",
                      dest = "nf", default = 2048,
                      help = "Number of pixels in fast direction, eg 2048")
    
    parser.add_option("--ns", action="store", type="int",
                      dest = "ns", default = 2048,
                      help = "Number of pixels in slow direction, eg 2048")
    
    return parser


def main():
    """
    A CLI interface
    """
    import sys, time, os, logging
    start = time.time()
    root = logging.getLogger('')
    root.setLevel( logging.WARNING )
    try:
        from optparse import OptionParser
        parser = OptionParser()
        parser = get_options( parser )
        options, args = parser.parse_args()
    except SystemExit:
        raise
    except:
        parser.print_help()
        print("\nSorry, there was a problem interpreting your command line")
        raise

        
    if options.pars is None:
        print("Failed: You must supply a parameters file, -p option")
        sys.exit()
    if not os.path.exists( options.pars ):
        print("Cannot find your file",options.pars)
        sys.exit()

    worker = xydisp(
        splinefile = options.spline,
        parfile = options.pars
        )
    
    worker.compute_tth_eta( (options.nf, options.ns) )
    worker.write( options.output )
        

        
if __name__ == "__main__":
    main()
