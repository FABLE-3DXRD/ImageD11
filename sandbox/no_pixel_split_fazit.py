
from __future__ import print_function


import sys, numpy, time
from ImageD11 import cImageD11, blobcorrector, ImageD11_file_series
from fabio.openimage import openimage
from fabio.edfimage import edfimage

def get_options(parser):
    """
    Command line interface for making a mapping
    Add our options to a parser object
    """
    parser = ImageD11_file_series.get_options( parser )
    parser.add_option("--lookup", action="store",type="string",
                      dest="lookup",default="None",
                      help="Lookup tables stem name")
    parser.add_option("--outdir", action="store",type="string",
                      dest="outdir",default=None,
                      help="Output directory")
    parser.add_option("--mask", action="store",type="string",
                      dest="mask",default=None,
                      help="fit2d mask file")

    return parser


def main():
    """
    A CLI user interface
    """
    import sys, time, os, logging
    start = time.time()

    root = logging.getLogger('')
    root.setLevel(logging.WARNING)

    try:
        from optparse import OptionParser
        parser = OptionParser()
        parser = get_options( parser )
        options, args = parser.parse_args()
    except SystemExit:
        raise
    except:
        parser.print_help()
        print("\nProblem with your options:")
        raise

    if options.mask is not None:
        fit2dmask = (1-openimage(options.mask).data).ravel()
    else:
        fit2dmask = 1.0

    first_image = True

    imagefiles = ImageD11_file_series.get_series_from_options( options, args )

    tthvals = numpy.load(options.lookup+"_tth.npy")

    try:
        for fim in imagefiles:
            dataim = fim.data
            print(fim.filename)
            if first_image: # allocate volume, compute k etc

                first_image = False


                dxim = openimage(options.lookup+"_dx.edf").data
                dyim = openimage(options.lookup+"_dy.edf").data

                outsum = numpy.ravel( numpy.zeros( dataim.shape,
                                                   numpy.float32 ) )
                outnp = numpy.ravel( numpy.zeros( dataim.shape,
                                                  numpy.float32 ) )

                e = edfimage()


                # C code from rsv_mapper (not intended to be obfuscated)
                o = blobcorrector.perfect( )
                idealx, idealy = o.make_pixel_lut( dataim.shape )

                destx, desty = idealx + dxim, idealy + dyim
    
                assert destx.min() >= 0
                assert destx.max() < dataim.shape[1]
                assert desty.min() >= 0
                assert desty.max() < dataim.shape[1]

                imageshape = dataim.shape
                
                indices = numpy.ravel( destx ).astype( numpy.intp )
                numpy.multiply( indices, dataim.shape[1], indices )
                numpy.add( indices,
                           numpy.ravel( desty ).astype( numpy.intp ),
                           indices )

                assert indices.min() >= 0
                assert indices.max() < dataim.shape[0]*dataim.shape[1]
                on = numpy.ones(len(outnp), numpy.float32)
                if fit2dmask is not None:
                    on = on * fit2dmask
                    
                # Number of pixels and mask are constant
                cImageD11.put_incr( outnp,
                                  indices,
                                  on )

                mask = outnp < 0.1
                scalar = ( 1.0 - mask ) / ( outnp + mask )
                    
                flatshape = outsum.shape

                # 

                arsorted = mask.copy()
                outmask = mask.copy()
                outmask = outmask * 1e6
                outmask.shape = imageshape
                arsorted.shape = imageshape
                arsorted.sort(axis=1)
                minds = numpy.array([ l.searchsorted(0.5) for l in arsorted ])

                

            # ENDIF firstimage


            
            start = time.time()

            numpy.multiply( outsum, 0, outsum )
            outsum.shape = flatshape

            dm = (dataim.ravel()*fit2dmask).astype(numpy.float32)
            
            cImageD11.put_incr( outsum,
                              indices,
                              dm )

            # outsum = outsum.reshape( dataim.shape )
            # outnp = outnp.reshape( dataim.shape ).astype(numpy.int32)
            # Normalise
            numpy.multiply( outsum, scalar, outsum )
            print(dataim.max(),dataim.min(), end=' ')
            print(scalar.max(),scalar.min(),outsum.min(), outsum.max( ))

            outsum.shape = imageshape
            # saving edf
            e.data=outsum
            e.write( "r_"+fim.filename  , force_type=numpy.float32)
            
            print(time.time()-start)

            


    except:
        raise

#f = open("ar.bin","w")
#f.write(outsum.tostring())
#f.write(outnp.tostring())


    
        
if __name__ == "__main__":
    main()

