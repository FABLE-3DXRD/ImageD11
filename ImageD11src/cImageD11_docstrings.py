blob_moments = """fills in the reduced moments in results array.
... FIXME - this would be clearer in python, fast anyway.
"""
bloboverlaps = """determines the overlaps between labels1 and labels2
for an image series. Peaks in labels2 may be merged if they were
joined by a peak on labels1. Results in results1 are accumulated
into results2 if peaks are overlapped.
"""
blobproperties = """fills the array results with properties of each labelled
object described by data (pixel values) and labels. The omega value
is the angle for this frame.
results are FIXME
"""
compress_duplicates = """removes duplicate i,j labels. On entry then
i and j are set as the labels from two images. They are sorted
and on exit i,j hold the unique pairs and oi holds the count
for the number of overlaps. oj and tmp are temporaries.
"""
connectedpixels = """Determines which pixels in data are above the
user supplied threshold and assigns them into connected objects 
which are output in labels. Connectivity is 3x3 box (8) by default
and reduces to a +(4) is con8==0
"""
localmaxlabel = """assigns a label for each pixel so they are grouped
to the local maximum. Equal values choose to assign towards the earlier
value in memory.
cpu arg (1)0=C, (1)1=SSE2, (1)2=AVX2; if > 9 prints timing
"""
sparse_blob2Dproperties = """fills the array results with properties of each labelled
object described by v and labels (pixel values and blob) and positions i,j in the image.
results are: 
  results[ipk,s2D_1]   = sum (1), number of pixels
  results[ipk,s2D_I]   = sum (I), total intensity
  results[ipk,s2D_fI]  = sum (f*I), intensity weighted fast index
  results[ipk,s2D_sI]  = sum (s*I), intensity weighted slow index
  results[ipk,s2D_ffI] = sum (f*f*I), intensity weighted fast^2 index
  results[ipk,s2D_sfI] = sum (s*f*I), intensity weighted slow*fast index
  results[ipk,s2D_ssI] = sum (s*s*I), intensity weighted slow^2 index
"""
sparse_localmaxlabel = """assigns labels to sparse array in sorted coo format
supplied in (v,(i,j)). MV and iMV are temporaries.
single threaded 
"""
sparse_overlaps = """identifies the pixels in i1,j1 which overlap i2,j2.
The list of overlaps is returned in k1/k2 such that i1[k1]==i2[k2]
and j1[k1]==j2[k2]. Probably assumes that sparse_is_sorted was true.
"""
__all__ = ["blob_moments","bloboverlaps","blobproperties","compress_duplicates","connectedpixels","localmaxlabel","sparse_blob2Dproperties","sparse_localmaxlabel","sparse_overlaps"]