
from __future__ import print_function, division

import time, sys
import h5py, scipy.sparse, numpy as np, pylab as pl
from ImageD11 import cImageD11

# see also sandbox/collect_peak_pixels.py

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
        pixels = Various properties labelled here
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
            for k,v in pixels.items():
                print(v.shape, self.nnz)
                assert len(v) == self.nnz
                self.pixels[k] = v
        
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

    def to_dense(self, data, out=None):
        """ returns the full 2D image 
        Does not handle repeated indices
        e.g.  obj.to_dense( obj.pixels['raw_intensity'] )
        """
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
        return sparse_frame( self.row[msk],
                             self.col[msk],
                             self.shape, self.row.dtype,
                             pixels= { name: pixval[msk] for (name, pixval)
                                       in self.pixels.items()} )

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
        """ Put the pixels into a different order """
        assert len(order) == self.nnz
        self.row = self.row[order]
        self.col = self.col[order]
        self.pixels = { (name, pixval[order]) for (name, pixval)
                        in self.pixels.items() }
        
    def threshold(self, threshold, name='data'):
        """
        returns a new sparse frame with pixels > threshold
        """
        return self.mask( self.pixels[name] > threshold )


class sparse_labelling( object ):
    """ Augments a sparse frame with some labels """
    
    def __init__(self, frame, name='data'):
        self.frame = frame
        self.name  = name
        self.labels = np.zeros( self.frame.nnz, 'i' )
        
    def moments(self):
        return cImageD11.sparse_blob2Dproperties(
            self.frame.pixels[self.intensity_name],
            self.frame.row,
            self.frame.col,
            self.labels,
            self.nlabel )

    def overlaps(self, other):
        """
        figures out which label of self matches which label of other
        Returns sparse array of:
           label in self (row)
           label in other (col)
           number of shared pixels (data)
        """
        ki = np.empty(  self.frame.nnz, 'i' )
        kj = np.empty( other.frame.nnz, 'i' )
        npx = cImageD11.sparse_overlaps( self.frame.row, self.frame.col, ki,
                                         other.frame.row, other.frame.col, kj)
        # self.data and other.data filled during init
        row = self.labels[ ki[:npx] ]  # my labels
        col = other.labels[ kj[:npx] ] # your labels
        ect = np.empty( npx, 'i')    # ect = counts of overlaps
        tj  = np.empty( npx, 'i')    # tj = temporary  for sorting
        tmp = np.empty( max(self.nlabel, other.nlabel)+1, 'i') # for histogram
        
        nedge = cImageD11.compress_duplicates( row, col, ect, tj, tmp )
        # overwrites row/col in place
        crow = row[:nedge]
        ccol = col[:nedge]
        cdata = ect[:nedge]
        cedges = scipy.sparse.coo_matrix( ( cdata, (crow, ccol)) )
        return cedges
 
    
class sparse_connected_pixels( sparse_labelling ):
    
    def __init__(self, frame, threshold, name='data'):
        """
        guesses threshold as min val-1 or given threshold
        adds an array named [cplabels_tx.xxx]
        self.nclabel number of connected objects
        """
        sparse_labelling.__init__(self, frame, name=name)
        if threshold is None:
            self.threshold = self.frame.pixels[name].min()-1
        else:
            self.threshold = threshold
        self.nlabel = cImageD11.sparse_connectedpixels(
            self.frame.pixels[self.name], self.frame.row, self.frame.col,
            self.threshold,   self.labels )


class sparse_localmax( sparse_labelling ):

    def __init__(self, frame, name = 'data'):
        """
        adds an array names mxlabels
        """
        sparse_labelling.__init__(self, frame, name=name)
        vmx = np.zeros( self.nnz, np.float32 )
        imx = np.zeros( self.nnz, 'i')
        self.nlabel = cImageD11.sparse_localmaxlabel(
            self.frame.pixels[self.name], self.frame.row, self.frame.col,
            vmx, imx, self.labels )
        



        

