

import sys, os
import hdf5plugin, h5py, timeit
import fabio
import numpy as np
import scipy.sparse
from ImageD11 import sparseframe, cImageD11, columnfile, blobcorrector

timer= timeit.default_timer
T0 = START = timer()
def bench(msg):
    #flake8: global T0
    global START
    end = timer()
    print(msg, "%.4f /s %.4f /s"%( end - START, end - T0 ) )
    START = end


class sparsedatafile(object):
    
    def __init__(self, hname, folder):
        self.hname = hname
        self.h = h5py.File( hname, "r" )
        gin = self.h[folder]
        self.shape = gin.attrs['shape0'], gin.attrs['shape1']
        self.row = gin['row']
        self.col = gin['col']
        self.sig = gin['intensity']
        self.omega = gin['omega']
        self.ipt = np.cumsum( np.concatenate( ( (0,), gin['nnz'][:] ) ) )
        self.nframes = gin['nnz'].shape[0]
        
    def getframe(self, i):
        assert i>=0 and i<self.nframes
        ilo = self.ipt[i]
        ihi = self.ipt[i+1]
        f = sparseframe.sparse_frame( self.row[ilo:ihi],
                                      self.col[ilo:ihi],
                                      shape=self.shape )
        f.set_pixels("intensity", self.sig[ilo:ihi] )
        # print("%d"%(i),end=" ")
        # sys.stdout.flush()
        return f

    def peaksearch(self, frm):
        # -> f.pixels['localmax']
        sparseframe.sparse_localmax(frm)
        # -> f.pixels['connectedpixels']
        sparseframe.sparse_connected_pixels(frm, threshold = np.finfo( np.float32 ).min )
        # which localmax peaks overlap each other
        adj = sparseframe.overlaps(frm, 'localmax', frm, 'connectedpixels' )
        r, c = self.adj2overlap( adj )
        return frm, r, c

    def adj2overlap(self, adj):
        """ Locate the overlaps on the same frame """
        cts = np.bincount( adj.col, minlength=adj.shape[1] )
        pairs = np.arange( adj.shape[1], dtype=int)[ cts > 1 ]
        rows = []
        cols = []
        for k in pairs:
            overlaps = adj.row[ adj.col == k ]
            for i,p1 in enumerate(overlaps):
                for p2 in overlaps[i+1:]:
                    rows.append( p1 )
                    cols.append( p2 )
        return np.array(rows), np.array(cols)

    def peaksearch_all(self):
        rows = []
        cols = []
        
        frm, r, c = self.peaksearch(self.getframe(0))
        frames = [ frm, ]

        npks = [ frames[0].meta['localmax']['nlabel'], ]
        cpks = [0, npks[-1]] # cumulative peak count

        for i in range(1,self.nframes):
            frm, r, c  = self.peaksearch( self.getframe(i)) # labelled
            n = frm.meta['localmax']['nlabel']
            frames.append( frm )
            rows.append(r+cpks[i])
            cols.append(c+cpks[i])
            npks.append( n )
            cpks.append( cpks[-1] + n ) # cumsum
            m = sparseframe.overlaps(frames[i-1], 'localmax',
                                     frames[i  ], 'localmax')
            rows.append( m.row + cpks[i-1] )
            cols.append( m.col + cpks[i] )
        self.rows = np.concatenate( rows )
        self.cols = np.concatenate( cols )
        self.frames = frames # needed?
        self.cpks = cpks

    def overlaps(self):
        ntot = self.cpks[-1]
        g = scipy.sparse.coo_matrix( (np.ones(len(self.rows),np.uint8),
                                      (self.rows, self.cols)),
                                     shape=(ntot,ntot))
        self.ncomp, self.labels = scipy.sparse.csgraph.connected_components(
            g, return_labels=True, directed=False )
        
    def moments(self, spline=None):
        titles = "1 I fI sI".split() # for now
        npx = cImageD11.s2D_1
        I  = cImageD11.s2D_I
        sI = cImageD11.s2D_sI
        fI = cImageD11.s2D_fI
        # results
        c = { '1' : np.zeros( self.ncomp, 'f' ),
              'I' : np.zeros( self.ncomp, 'f' ),
              'sI': np.zeros( self.ncomp, 'f' ),
              'fI': np.zeros( self.ncomp, 'f' ),
              'oI': np.zeros( self.ncomp, 'f' ),
        }
        
        for i,f in enumerate(self.frames):
            m = sparseframe.sparse_moments( f, 'intensity', 'localmax' )
            #                   data  , ind, vals : data[ind] += vals
            inds = self.labels[ self.cpks[i]:self.cpks[i+1] ]
            cImageD11.put_incr( c[ '1'], inds, m[:, npx] )
            cImageD11.put_incr( c[ 'I'], inds, m[:, I]   )
            cImageD11.put_incr( c['sI'], inds, m[:, sI]  )
            cImageD11.put_incr( c['fI'], inds, m[:, fI]  )
            cImageD11.put_incr( c['oI'], inds, m[:, I]*self.omega[i] )
        s_raw = c['sI']/c['I']
        f_raw = c['fI']/c['I']
        omega = c['oI']/c['I']
        if spline is not None:
            sc, fc = spatial_correct( s_raw, f_raw, spline )
        else:
            sc, fc = s_raw, f_raw
        results = {
            'sc': sc,
            'fc': fc,
            'omega': omega,
            'Number_of_pixels' : c['1'],
            'avg_intensity' : c['I']/c['1'],
            's_raw' : s_raw,
            'f_raw' : f_raw,
            'sum_intensity': c['I'],
            }
        return results
                                

def spatial_correct( s_raw, f_raw, spline):
    si = np.round( s_raw ).clip(0,2047).astype(int)
    fi = np.round( f_raw ).clip(0,2047).astype(int)
    try:
        dx, dy = spline
        fc = f_raw + dx[si, fi]
        sc = s_raw + dy[si, fi]
    except:
        corr = blobcorrector.correctorclass(spline)
        dx, dy = corr.make_pixel_lut( (2048, 2048) )
        fc = f_raw - fi + dy[si, fi]
        sc = s_raw - si + dx[si, fi]
    return sc, fc
    

def dict_to_colfile( c ):
    titles = list(c.keys())
    nrows = len(c[titles[0]])
    for t in titles:
        assert len(c[t]) == nrows, t
    colf = columnfile.newcolumnfile( titles=titles )
    colf.nrows = nrows
    colf.set_bigarray( [ c[t] for t in titles ] )
    return colf

# /data/id11/jon/inhouse/jon/mar10/segment2d
hdffile = "fe3o4_peaks.hdf"
folder = "entry"
parfile = "avg.par"
splinefile = "frelon4m.spline"
splinefile = (fabio.open("F4M_EO_dx.edf").data,fabio.open("F4M_EO_dy.edf").data)



bench( "startup" )
obj = sparsedatafile( hdffile, folder )            ;  bench("openhdf")
obj.peaksearch_all()                               ;  bench("peaksearch")
obj.overlaps()                                     ;  bench("overlaps")
res = obj.moments(splinefile)                      ;  bench("moments")
c = dict_to_colfile( res )                         ;  bench("convert")
c.parameters.loadparameters( parfile )
c.updateGeometry()                                 ;  bench("geometry")
if os.path.exists('peaks.hdf'):
    os.remove('peaks.hdf')
    columnfile.colfileobj_to_hdf( c, "peaks.hdf" ) ; bench("save hdf")

if 1:
    import pylab as pl
    pl.plot( c.tth, c.eta, ",")
    pl.title("%d peaks"%(c.nrows))
    pl.show()

