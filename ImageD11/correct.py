

import numpy
# From PIL
import ImageFilter

# These don't work
filternames = [ "BLUR", "CONTOUR", "DETAIL", "EDGE_ENHANCE",
                "EDGE_ENHANCE_MORE", "EMBOSS", "FIND_EDGES", "SMOOTH",
                "SMOOTH_MORE", "SHARPEN"]
filters = {}
for f in filternames:
    filters[f] = getattr( ImageFilter, f )

filternames.append("MedianFilter(3)")
filters["MedianFilter(3)"] = ImageFilter.MedianFilter(3)
filternames.append("MedianFilter(5)")
filters["MedianFilter(5)"] = ImageFilter.MedianFilter(5) 

# fixme - subtracting median filtered
# coarser medians - eg rebinned too

def correct(data_object, dark = None, flood = None, do_median = False,
            filterlist = [] ):
    """
    Does the dark and flood corrections
    Also PIL filters
    """
    picture = data_object.data.astype(numpy.float32)
    if dark != None:
        # This is meant to be quicker
        picture = numpy.subtract( picture , dark, picture )
        data_object.data = picture
    if flood != None:
        picture = numpy.divide( picture, flood, picture )
        data_object.data = picture
    if do_median:
        # We do this after corrections
        # The expectation is that this is a median on the azimuth
        # direction of a previously radially transformed image
        # Gives the liquid contribution 
        med = numpy.median( picture )
        if True: # Suboption - save the median or not?
            obj = fabio.deconstruct_filename( data_object.header['filename'] )
            obj.extension = ".bkm"
            medfilename = obj.tostring()
            med.tofile( medfilename , sep = "\n")
        picture = numpy.subtract( picture , med, picture )
    # Apply series of PIL filters
    if len( filterlist ) > 0:
        pim = data_object.toPIL16()
        print "Applied",
        for item in filterlist:
            if item in filternames:
                try:
                     pim = pim.filter( filters[ item ])
                except:
                    raise
                print item,
        data_object.data = numpy.array( pim )
    return data_object


