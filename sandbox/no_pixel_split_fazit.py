

import sys, numpy, time
from ImageD11 import closest, blobcorrector
from fabio.openimage import openimage
from fabio.edfimage import edfimage

def get_options(parser):
    """
    Command line interface for making a mapping
    Add our options to a parser object
    """
    parser = ImageD11_file_series.get_options( parser )
    


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
        print "\nProblem with your options:"
        raise

    first_image = True

    imagefiles = ImageD11_file_series.get_series_from_options( options, args )

    try:
        for fim in imagefiles:

            if first_image: # allocate volume, compute k etc

                first_image = False

                dataim = fim.data
                dxim = openimage(options.dxim).data
                dyim = openimage(options.dyim).data

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
                
                # Number of pixels and mask are constant
                closest.put_incr( outnp,
                                  indices,
                                  on )

                mask = outnp == 0
                scalar = ( 1.0 - mask ) / ( outnp + mask )

                flatshape = outsum.shape

                # 
            # ENDIF firstimage

            start = time.time()

            numpy.multiply( outsum, 0, outsum )
            outsum.shape = flatshape
            
            closest.put_incr( outsum,
                              indices,
                              numpy.ravel(dataim).astype(numpy.float32) )
            
            outsum.shape = imageshape

            print "The actual rebin took around",time.time()-start,"seconds"

            # outsum = outsum.reshape( dataim.shape )
            # outnp = outnp.reshape( dataim.shape ).astype(numpy.int32)


            # Normalise
            numpy.multiply( outsum, scalar, outsum )

            e.data=outsum

            e.write( , force_type=numpy.float32)



#f = open("ar.bin","w")
#f.write(outsum.tostring())
#f.write(outnp.tostring())


    

