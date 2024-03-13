
from __future__ import print_function

# Make a project file in a hdf.

# example : 
#
from __future__ import division
import os
import numpy as np
import h5py
import fabio

# how the scan was done
fmt = "/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y%03d_/interlaced_1_%d/Frelon/interlaced_1_%d_Frelon%04d.edf"
ysteps = np.arange( -60, 60.1, 1)
npts = 360 * 2
ostep = 0.5
ostart = 0.25
# We recorded a series of scans
omegas = [ np.arange( npts, dtype = np.float32 ) * ostep + ostart, ] * len( ysteps )
dtys   = [ np.full( npts, ystep, dtype = np.float32 ) for ystep in ysteps ]

def interlacedscan( ystep, fmt, npts ):
    scan = []
    for i in range( npts//2 ):
        scan.append( fmt % ( ystep, 1, 1, i ) )         # forward interlaced
        scan.append( fmt % ( ystep, 2, 2, npts//2 - 1 - i) )
    return scan

filenames = [ interlacedscan( ystep, fmt, npts ) for ystep in ysteps ] 

# Information about the edf files
def binary_info( filename ):
    """ info to read the edf files - assume flat regular files """
    im = fabio.open( filename )
    return im._frames[0].start, im._frames[0].size, im.data.shape, im.data.dtype

def writesino(h5name, omegas, dtys, filenames):
    offset, size, shape, dtype = binary_info( filenames[0][0] )
    print(offset,size,shape,dtype)
    nframes = len( omegas[0] ) * len( omegas )
    print(nframes, len(omegas), sum(len(o) for o in omegas))
    # Now create a hdf5 file:
    with h5py.File(h5name, "w", libver='latest' ) as h:
        # now create a VDS linking within the same file
        layout = h5py.VirtualLayout( shape = (nframes, shape[0], shape[1] ),
                                     dtype = dtype )
        j = 0
        graw = h.require_group('scans')
        for i, scan in enumerate(filenames):
            g = graw.require_group('scan%04d'%(i))
            g.create_dataset( "data",
                              shape = (len(scan), shape[0], shape[1]),
                              dtype = dtype,
                              external = [(fname, offset, size) for fname in scan] )
            g.create_dataset( "omega" , data = omegas[i] )
            g.create_dataset( "dty" , data = dtys[i] )
            vsource = h5py.VirtualSource( h.filename, # ok - circular?
                                          'scans/scan%04d/data'%(i),
                                          shape = (len(scan), shape[0], shape[1]) )
            layout[ j:j+len(scan), :, :] = vsource
            j += len(scan)
        g = h.require_group('sinogram')
        g.create_dataset('omega', data = np.concatenate(omegas) )
        g.create_dataset('dty', data = np.concatenate(dtys) )
        g.create_virtual_dataset( 'data', layout )
        
h5name = 'demo1.h5'
writesino( h5name, omegas, dtys, filenames )

# Now some tests ...:

for iy in range(0,len(dtys), 31):
    print('iy',iy)
    for jo in (0,1,123):
        fname = filenames[iy][jo]
        print('  ', fname, iy*npts + jo, end= ' ' )
        fab = fabio.open( fname ).data
        h5v = h5py.File( h5name, 'r' )[ 'sinogram/data' ][ iy * npts + jo ]
        # external data
        # h5e = h5py.File( h5name, 'r' )[ 'scans/scan%04d/data'%( iy ) ] [jo]
        if (h5v == fab).all():
            print('ok')
        else:
            print('debug please')
            print(fab.shape,fab)
            print(h5v.shape,h5v)
