
from __future__ import print_function, division

import time, sys
import h5py, scipy.sparse, numpy as np #, pylab as pl
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
    
    This was developed for a single 2D frame
       See SparseScan below for something aiming towards many frames
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
        self.meta = {}
        if pixels is not None:
            for name, val in pixels.items():
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
        scipy.sparse.coo_matrix((data, (self.row, self.col)), shape=(self.shape)).todense(out=out)
        # does not handle duplicate indices if they were present:
        #        adr = self.row.astype(np.intp) * self.shape[1] + self.col
        #        out.flat[adr] = data    
        return out

    def mask( self, msk ):
        """ returns a subset of itself """
        spf = sparse_frame( self.row[msk],
                            self.col[msk],
                            self.shape, self.row.dtype )
        for name, px in self.pixels.items():
            if name in self.meta:
                m = self.meta[name].copy()
            else:
                m = None
            spf.set_pixels( name, px[msk], meta = m )
        return spf

    def set_pixels( self, name, values, meta=None ):
        """ Named arrays sharing these labels """
        assert len(values) == self.nnz
        self.pixels[name] = values
        if meta is not None:
            self.meta[name] = meta
            

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
            if pxname in frame.meta:
                group[pxname].attrs = dict( frame.meta[pxname] )


class SparseScan( object ):
    
    omeganames = ['measurement/rot_center', 'measurement/rot',
                  'measurement/diffrz_center', 'measurement/diffrz']
    dtynames   = ['measurement/dty_center', 'measurement/dty',
                  'measurement/diffty_center', 'measurement/diffty']
    
    def __init__( self, hname, scan ):
        """
        hname : file coming from a sparse segmentation
        scan : a scan within that file
        motors : which motor channels to (try) to read
        
        assumes the scan fits into memory (could be problematic)
        """
        with h5py.File(hname,"r") as hin:
            grp = hin[scan]
            self.shape = ( int(v) for v in ( grp.attrs['nframes'], 
                                            grp.attrs['shape0'], 
                                            grp.attrs['shape1'] ) )
            self.motors = {}
            for name, motors in [ ('omega',self.omeganames),
                                  ('dty',self.dtynames) ]:
                for motor in motors:
                    if motor in grp:
                        self.motors[ name ] = grp[motor][:]
                        break
                
            self.nnz = grp['nnz'][:]
            self.ipt = np.concatenate( ( (0,) , np.cumsum(self.nnz, dtype=int) ) )
            self.frame  = grp['frame'][:]
            self.row = grp['row'][:]
            self.col = grp['col'][:]
            self.intensity = grp['intensity'][:]
            
    def cplabel(self, threshold = 0 ):
        """ Label pixels using the connectedpixels assigment code
        Fills in:
           self.nlabels = number of peaks per frame
           self.labels  = peak labels (should be unique)
           self.total_labels = total number of peaks
        """
        self.nlabels = np.zeros( len(self.nnz), np.int32 )
        self.labels = np.zeros( len(self.row), "i")
        nl = 0
        # TODO: run this in parallel with threads?
        for i, npx in enumerate( self.nnz ):
            s = self.ipt[i]
            e = self.ipt[i+1]
            if npx > 0:
                self.nlabels[i] = cImageD11.sparse_connectedpixels(
                    self.intensity[ s : e ],
                    self.row[ s : e ],
                    self.col[ s : e ],
                    threshold,
                    self.labels[ s : e ] )
                # zero label is the background!
                self.labels[ s : e ] = np.where( self.labels[ s : e ] > 0, 
                                                 self.labels[ s : e ] + nl, 0 )
            else:
                self.nlabels[i] = 0
            nl += self.nlabels[i]
        self.total_labels = nl

          
    def lmlabel(self, threshold = 0 ):
        """ Label pixels using the localmax assigment code
        Fills in:
           self.nlabels = number of peaks per frame
           self.labels  = peak labels (should be unique)
           self.total_labels = total number of peaks
        """
        self.nlabels = np.zeros( len(self.nnz), np.int32 )
        self.labels = np.zeros( len(self.row), "i")
        # temporary workspaces
        npxmax = self.nnz.max()
        vmx = np.zeros( npxmax, np.float32 )
        imx = np.zeros( npxmax, 'i' )
        nl = 0
        # TODO: run this in parallel with threads?
        for i, npx in enumerate( self.nnz ):
            s = self.ipt[i]
            e = self.ipt[i+1]
            if npx > 0:
                self.nlabels[i] = cImageD11.sparse_localmaxlabel(
                    self.intensity[ s : e ],
                    self.row[ s : e ],
                    self.col[ s : e ],
                    vmx[:npx],
                    imx[:npx],
                    self.labels[s : e] )
                assert (self.labels[s:e] > 0).all()
                self.labels[ s : e ] += nl
            else:
                self.nlabels[i] = 0
            nl += self.nlabels[i]
        self.total_labels = nl
            
    def moments(self):
        """ Computes the center of mass in s/f/omega
        returns a columnfile
        """
        pks = {}
        i32 = self.intensity.astype(np.float32)
        pks['Number_of_pixels'] = np.bincount(self.labels, 
                                              weights=None,
                                              minlength = self.total_labels+1 )[1:]
        pks['sum_intensity'] = np.bincount(self.labels, 
                                           weights=i32,
                                           minlength = self.total_labels+1 )[1:]
        pks['s_raw'] = np.bincount(self.labels,
                                   weights=i32*self.row,
                                   minlength = self.total_labels+1 )[1:]
        pks['s_raw'] /= pks['sum_intensity']
        pks['f_raw'] = np.bincount(self.labels,
                                   weights=i32*self.col,
                                   minlength = self.total_labels+1 )[1:]
        pks['f_raw'] /= pks['sum_intensity']
        for name in 'omega','dty':
            if name in self.motors:
                pks[name] = np.bincount(self.labels,
                           weights=i32*self.motors[name][self.frame],
                           minlength = self.total_labels+1 )[1:]
                pks[name] /= pks['sum_intensity']
        return pks
                
    
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
    #    intensity.attrs = dict(header) # FIXME USE xarray ?
    spf = sparse_frame( row, col, data.shape, itype=np.uint16 )
    spf.set_pixels( "intensity" , intensity, dict( header ) )
    return spf

def from_hdf_group( group ):
    itype = np.dtype( group.attrs['itype'] )
    shape = group.attrs['shape0'], group.attrs['shape1']
    row = group['row'][:] # read it
    col = group['col'][:]
    spf = sparse_frame( row, col, shape, itype=itype )
    for pxname in list(group):
        if pxname in ["row", "col"]:
            continue
        data = group[pxname][:]
        header = dict( group[pxname].attrs )
        spf.set_pixels( pxname, data, header )
    return spf

def sparse_moments( frame, intensity_name, labels_name ):
    """ We rely on a labelling array carrying nlabel metadata (==labels.data.max())"""
    nl = frame.meta[ labels_name ][ "nlabel" ]
    return cImageD11.sparse_blob2Dproperties(
        frame.pixels[intensity_name],
        frame.row,
        frame.col,
        frame.pixels[labels_name],
        nl )


def overlaps(frame1, labels1, frame2, labels2):
    """
    figures out which label of self matches which label of other
    Assumes the zero label does not exist (background)
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
    n1  = frame1.meta[labels1][ "nlabel" ]
    n2  = frame2.meta[labels2][ "nlabel" ]
    tmp = np.empty( max(n1, n2)+1, 'i') # for histogram
    nedge = cImageD11.compress_duplicates( row, col, ect, tj, tmp )
    # overwrites row/col in place : ignore the zero label (hope it is not there)
    crow = row[:nedge]-1
    ccol = col[:nedge]-1
    cdata = ect[:nedge]
    cedges = scipy.sparse.coo_matrix( ( cdata, (crow, ccol)), shape=(n1, n2) )
    # really?
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
        threshold = frame.meta[data_name]["threshold"]
    nlabel = cImageD11.sparse_connectedpixels(
        frame.pixels[data_name], frame.row, frame.col,
        threshold,  labels )
    frame.set_pixels( label_name, labels, { 'nlabel' : nlabel } )
    return nlabel


def sparse_localmax( frame,
                     label_name="localmax",
                     data_name = "intensity" ):
    labels = np.zeros( frame.nnz, "i" )
    vmx = np.zeros( frame.nnz, np.float32 )
    imx = np.zeros( frame.nnz, 'i')
    nlabel = cImageD11.sparse_localmaxlabel(
        frame.pixels[data_name], frame.row, frame.col,
        vmx, imx, labels )                   
    frame.set_pixels( label_name, labels, { "nlabel" : nlabel }  )
    return nlabel


    
                              
        



        

