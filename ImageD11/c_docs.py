

"""
Add documentation to compiled codes
"""


connectedpixels_doc = """
  nblobs = connectedpixels ( data=numpy.array(data, 2D, np.float32),
                             labels=Numeric.array(blob, 2D, np.int32),
                             threshold=float threshold,
                             verbose=Int verbose )
        data is the image to be labelled
        blob is an array to receive pixel -> blob assignments
        threshold is the value above which a pixel assigned to a blob
        verbose flag : 0=quiet, higher = louder printing on stdout
"""

blobproperties_doc = """
     res = blobproperties ( numpy.array(data, 2D, np.float32),
                            numpy.array(blob, 2D, np.int32)  ,
                            int32 npk, 
                            float32 omega,
                            int32 verbose )
        Computes various properties of blobs in a data image
                (created by connectedpixels or other labelling)
        data  = image data
        blob  = integer peak assignments (e.g. connectedpixels)
        omega = third dimension value (which image is this)
        npk    = number of peaks to treat
        verbose flag : 0=quiet, higher = louder printing on stdout
     Returns res : array of peak properties
"""

bloboverlaps_doc = """
 success = bloboverlaps ( numpy.array(blob1, 2D, np.int32), n1  , res1,
                          numpy.array(blob2, 2D, np.int32), n2  , res2,
                          verbose=0)
      merges the results from blob1/res1 into blob2/res2
      blob1 would be the previous image from the series, with its results
        if it does not overlap it will stay untouched in res
        if it overlaps with next image it is passed into that ones res
"""


blob_moments_doc = """
       blob_moments(numpy.array(res1), npk)
           Loop over array filling out moments from sums 
"""



closest_doc = """
    *    From python: closest(cosines_allowed, cosine_measured)
    *    Take 1 peak with hkl's assigned, compute the cosines to
    *    another powder ring (angles between this peak and peaks in
    *    that ring). This extension finds the best peak in the
    *    other ring to pair up with.
    *
    *    More generally closest(array, values) finds the closest
    *    item in values to one of the values in array.  Returns the
    *    difference and the index of the closest thing it finds.
    *
    *    Both arguments should be one dimensional Numeric arrays
    *    of type Numeric.Float
"""

score_doc = """
         From python: score(ubi, gv, tol) where ubi is an orientation
    *    matrix and gv are an array of g-vectors. Returns the number
    *    of g-vectors which have integer hkl peaks within tolerance
    *    tol. Uses the conv_double_to_int_fast function in here for 
    *    factor of !EIGHT! speed increase compared to rounding in 
    *    C. In fact this gives the nearest even integer, instead
    *    of the nearest integer, but we don't care, as a peak having
    *    hkl of 0.5 is nowhere near being indexed anyway.
    *
    *    Returns and integer - number of peaks indexed
    *    UBI is a 3x3 Numeric.Float array (figure out the transposing yourself)
    *    GV is a nx3 Numeric.Float array, and you should try to make the 3 
    *       be the fast index for best performance
"""

score_and_refine_doc = """
    *    score_and_refine(ubi, gv, tol) 
    *    From python: same as score, I hope, but ubi is overwritten
    *    with refined matrix following paciorek algorithm which is 
    *    in indexing.py
"""

score_and_assign_doc = """    
    *    score_and_assign( ubi, gv, tol, drlv, labels, (int) label)
    *    as for score, but assignments are only made if drlv is lower than
    *    the previous. Fills in label with the new label
"""

refine_assigned_doc = """
    * 5. refine_assigned( ubi, gv, labels, label, weight) 
    *    Use only the peaks where label matches labels in refinement 
    *    ...following again the Paciorek algorithm 
"""

put_incr_doc = """
    * 6. put_incr( data, indices, values)
    *    pretty much as numeric.put but as increments
"""

weighted_refine_doc = """
    * 7. weighted_refine(ubi, gv, tol, weights) 
    *    From python: same as score, but with weights, and ubi is overwritten
    *    with refined matrix following paciorek algorithm which is 
    *    in indexing.py
"""
