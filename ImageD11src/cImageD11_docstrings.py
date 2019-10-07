blobproperties = """fills the array results with properties of each labelled
object described by data (pixel values) and labels. The omega value
is the angle for this frame.
results are FIXME
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
__all__ = ["blobproperties","connectedpixels","localmaxlabel","sparse_blob2Dproperties","sparse_localmaxlabel"]