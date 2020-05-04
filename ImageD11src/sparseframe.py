
from __future__ import print_function, division

import time, sys
import h5py, scipy.sparse, numpy as np, pylab as pl
from ImageD11 import cImageD11


# see also sandbox/harvest_pixels.py

NAMES = {
    "filename" : "original filename used to create a sparse frame",
    "intensity" : "corrected pixel values",
    "nlabel": "Number of unique labels for an image labelling",
    "threshold" : "Cut off used for thresholding",
    }


class sparse_frame( object ):
    """
    Indices / shape mapping
    """
    def __init__(self, row, col, shape, itype=np.uint16, pixels=None):
        """ row = slow direction
        col = fast direction
        shape = size of full image
        itype = the integer type to store the indices
                our c codes currently use unsigned short...
        nnz is implicit as len(row)==len(col)
        pixels = numpy arrays in a dict to name them
                 throw in a ary.attrs if you want to save some
        """
        self.check( row, col, shape, itype )
        self.shape = shape
        self.row = np.asarray(row, dtype = itype )
        self.col = np.asarray(col, dtype = itype )
        self.nnz = len(self.row)
        # Things we could have using those indices:
        #   raw pixel intensities
        #   corrected intensities
        #   smoothed pixel intensities
        #   labelling via different algorithms
        self.pixels = {}
        if pixels is not None:
            for name, val in pixels.items():
                # print(val.data.shape, self.nnz)
                assert len(val) == self.nnz
                self.pixels[name] = val

        
    def check(self, row, col, shape, itype):
        """ Ensure the index data makes sense and fits """
        lo = np.iinfo(itype).min
        hi = np.iinfo(itype).max
        assert len(shape) == 2
        assert shape[0] >= lo and shape[0] < hi
        assert shape[1] >= lo and shape[1] < hi
        assert np.min(row) >= lo and np.max(row) < hi
        assert np.min(col) >= lo and np.max(col) < hi
        assert len(row) == len(col)
        
    def is_sorted(self):
        """ Tests whether the data are sorted into slow/fast order
        rows are slow direction
        columns are fast """
        # TODO: non uint16 cases
        assert self.row.dtype == np.uint16 and \
            cImageD11.sparse_is_sorted( self.row, self.col ) == 0

    def to_dense(self, data=None, out=None):
        """ returns the full 2D image 
        data = name in self.pixels or 1D array matching self.nnz
        Does not handle repeated indices
        e.g.  obj.to_dense( obj.pixels['raw_intensity'] )
        """
        if data in self.pixels:
            data = self.pixels[data] # give back this array
        else:
            ks = list( self.pixels.keys() )
            if len(ks)==1:
                data = self.pixels[ks[0]] # default for only one
            else:
                data = np.ones( self.nnz, np.bool ) # give a mask
        if out is None:
            out = np.zeros( self.shape, data.dtype )
        else:
            assert out.shape == self.shape
        assert len(data) == self.nnz
        adr = self.row.astype(np.intp) * self.shape[1] + self.col
        out.flat[adr] = data    
        return out

    def mask( self, msk ):
        """ returns a subset of itself """
        newpx = {}
        for name, px in self.pixels.items():
            newpx[name] = px[msk]
            if hasattr( px, "attrs" ):
                newpx[name].attrs = dict( px.attrs )
        return sparse_frame( self.row[msk],
                             self.col[msk],
                             self.shape, self.row.dtype,
                             pixels = newpx )

    def set_pixels( self, name, values ):
        """ Named arrays sharing these labels """
        assert len(values) == self.nnz
        self.pixels[name] = values

    def sort_by( self, name ):
        """ Not sure when you would do this. For sorting 
        by a peak labelling to get pixels per peak """
        assert name in self.pixels
        order = np.argsort( self.pixels[name] )
        self.reorder( self, order )

    def sort( self ):
        """ Puts you into slow / fast looping order """
        order = np.lexsort( ( self.col, self.row ) )
        self.reorder( self, order )

    def reorder( self, order ):
        """ Put the pixels into a different order (in place) """
        assert len(order) == self.nnz
        self.row[:] = self.row[order]
        self.col[:] = self.col[order]
        for name, px in self.pixels.items():
            px[:] = px[order]
        
    def threshold(self, threshold, name='intensity'):
        """
        returns a new sparse frame with pixels > threshold
        """
        return self.mask( self.pixels[name] > threshold )

    def to_hdf_group( frame, group ):
        """ Save a 2D sparse frame to a hdf group 
        Makes 1 single frame per group
        """
        itype = np.dtype( frame.row.dtype )
        meta = { "itype"  : itype.name,
                 "shape0" : frame.shape[0],
                 "shape1" : frame.shape[1] }
        for name, value in meta.items():
            group.attrs[name] = value
        opts = { "compression": "lzf",
                "shuffle" : True,
                }
        #opts = {}
        group.require_dataset( "row", shape=(frame.nnz,),
                               dtype=itype, **opts )
        group.require_dataset( "col", shape=(frame.nnz,),
                               dtype=itype, **opts )
        group['row'][:] = frame.row
        group['col'][:] = frame.col
        for pxname, px in frame.pixels.items():
            group.require_dataset( pxname, shape=(frame.nnz,),
                                   dtype=px.dtype,
                                   **opts ) 
            group[pxname][:] = px
            if hasattr( px, "attrs" ):
                group[pxname].attrs = dict( px.attrs )


    
def from_data_mask( mask, data, header ):
    """
    Create a sparse from a dense array
    """
    assert mask.shape == data.shape
    # using uint16 here - perhaps make this general in the future
    # ... but not for now
    assert data.shape[0] < pow(2,16)-1
    assert data.shape[1] < pow(2,16)-1
    nnz = (mask>0).sum()
    tmp = np.empty( data.shape[0],'i') # tmp hold px per row cumsums
    row = np.empty( nnz, np.uint16 )
    col = np.empty( nnz, np.uint16 )
    cImageD11.mask_to_coo( mask, row, col, tmp )
    intensity = data[ mask > 0 ]
    intensity.attrs = dict(header)
    return sparse_frame( row, col, data.shape, itype=np.uint16,
                  pixels = { "intensity" : intensity } )

def from_hdf_group( group ):
    itype = np.dtype( group.attrs['itype'] )
    shape = group.attrs['shape0'], group.attrs['shape1']
    row = group['row'][:] # read it
    col = group['col'][:]
    pixels = {}
    for pxname in list(group):
        if pxname in ["row", "col"]:
            continue
        data = group[pxname][:]
        header = dict( group[pxname].attrs )
        pixels[pxname] = data
        pixels[pxname].attrs = dict( header )
    return sparse_frame( row, col, shape, itype=itype, pixels=pixels )

def sparse_moments( frame, intensity_name, labels_name ):
    """ We rely on a labelling array carrying nlabel metadata (==labels.data.max())"""
    nl = frame.pixels[ labels_name ].attrs[ "nlabel" ]
    return cImageD11.sparse_blob2Dproperties(
        frame.pixels[intensity_name],
        frame.row,
        frame.col,
        frame.pixels[labels_name],
        nl )


def overlaps(frame1, labels1, frame2, labels2):
    """
    figures out which label of self matches which label of other
    Returns sparse array of:
    label in self (row)
    label in other (col)
    number of shared pixels (data)
    """
    ki = np.empty( frame1.nnz, 'i' )
    kj = np.empty( frame2.nnz, 'i' )
    npx = cImageD11.sparse_overlaps( frame1.row, frame1.col, ki,
                                     frame2.row, frame2.col, kj)
    # self.data and other.data filled during init
    row = frame1.pixels[labels1][ ki[:npx] ]  # my labels
    col = frame2.pixels[labels2][ kj[:npx] ]  # your labels
    ect = np.empty( npx, 'i')    # ect = counts of overlaps
    tj  = np.empty( npx, 'i')    # tj = temporary  for sorting
    n1  = frame1.pixels[labels1].attrs[ "nlabel" ]
    n2  = frame2.pixels[labels2].attrs[ "nlabel" ]
    tmp = np.empty( max(n1, n2)+1, 'i') # for histogram
    nedge = cImageD11.compress_duplicates( row, col, ect, tj, tmp )
    # overwrites row/col in place
    crow = row[:nedge]
    ccol = col[:nedge]
    cdata = ect[:nedge]
    cedges = scipy.sparse.coo_matrix( ( cdata, (crow, ccol)) )
    return cedges
 
    
def sparse_connected_pixels( frame,
                             label_name="connectedpixels",
                             data_name="intensity",
                             threshold=None ):
    """
    frame = a sparse frame
    label_name = the array to save labels to in that frame
    data_name = an array in that frame
    threshold = float value or take data.threshold
    """
    labels = np.zeros( frame.nnz, "i" )
    if threshold is None:
        threshold = frame.pixels[data_name].attrs["threshold"]
    nlabel = cImageD11.sparse_connectedpixels(
        frame.pixels[data_name], frame.row, frame.col,
        threshold,  labels )
    labels.attrs = { "data_name": data_name,
                     "threshold": threshold,
                     "nlabel"   : nlabel } 
    frame.set_pixels( label_name, labels )


def sparse_localmax( frame,
                     label_name="localmax",
                     data_name = "intensity" ):
    labels = np.zeros( frame.nnz, "i" )
    vmx = np.zeros( frame.nnz, np.float32 )
    imx = np.zeros( frame.nnz, 'i')
    nlabel = cImageD11.sparse_localmaxlabel(
        frame.pixels[data_name], frame.row, frame.col,
        vmx, imx, labels )
    labels.attrs = { "data_name" : data_name,
                     "nlabel" : nlabel } 
    frame.set_pixels( label_name, labels )

    
                              
        



        

