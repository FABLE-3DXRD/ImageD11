

import sys, numpy
from ImageD11 import closest, blobcorrector
from fabio.openimage import openimage
from fabio.edfimage import edfimage


dataim = openimage(sys.argv[1]).data
dxim = openimage(sys.argv[2]).data
dyim = openimage(sys.argv[3]).data




if False:
    # Super slow
    outsum = numpy.zeros( dataim.shape, numpy.float32 )
    outnp = numpy.zeros( dataim.shape, numpy.int32 )
    for i in range(dataim.shape[0]):
        print i,
        sys.stdout.flush()
        for j in range(dataim.shape[1]):
            dest_i = i + dxim[i,j]
            dest_j = j + dyim[i,j]
            outsum[ dest_i, dest_j ] += dataim[i,j]
            outnp[ dest_i, dest_j ]  += 1
else:

    outsum = numpy.ravel( numpy.zeros( dataim.shape, numpy.float32 ) )
    outnp = numpy.ravel( numpy.zeros( dataim.shape, numpy.float32 ) )

    # C code from rsv_mapper (not intended to be obfuscated)
    o = blobcorrector.perfect( )
    idealx, idealy = o.make_pixel_lut( dataim.shape )

    destx, desty = idealx + dxim, idealy + dyim
    
    assert destx.min() >= 0
    assert destx.max() < dataim.shape[1]
    assert desty.min() >= 0
    assert desty.max() < dataim.shape[1]


    indices = numpy.ravel( destx ).astype( numpy.intp )
    numpy.multiply( indices, dataim.shape[1], indices )
    numpy.add( indices,
               numpy.ravel( desty ).astype( numpy.intp ),
               indices )

    assert indices.min() >= 0
    assert indices.max() < dataim.shape[0]*dataim.shape[1]

    import time
    start = time.time()

    closest.put_incr( outsum,
                      indices,
                      numpy.ravel(dataim).astype(numpy.float32) )

    on = numpy.ones(len(outnp), numpy.float32) 
    closest.put_incr( outnp,
                      indices,
                      on )

    print "The actual rebin took around",time.time()-start,"seconds"

    outsum = outsum.reshape( dataim.shape )
    outnp = outnp.reshape( dataim.shape ).astype(numpy.int32)


# Normalise
mask = outnp == 0
signal = outsum / ( outnp + mask )
signal = signal*(1-mask) 

e = edfimage()
e.data=signal
e.write( sys.argv[4], force_type=numpy.float32)



#f = open("ar.bin","w")
#f.write(outsum.tostring())
#f.write(outnp.tostring())


    

