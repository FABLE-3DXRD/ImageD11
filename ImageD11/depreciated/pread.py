
"""WARNING: work in progress - look at mpi and/or uncompressed data instead?"""
from __future__ import print_function

import os
import h5py
import sys
import numpy as np
import ImageD11.sinograms.properties
import ImageD11.cImageD11
import multiprocessing
import concurrent.futures

def read_segment( args ):
    hname, dset, shmd, start, end = args
    tbl = ImageD11.sinograms.properties.pks_table.fromSHM( shmd )
    with h5py.File( open(hname, "rb"), 'r') as hin:
        sys.stdout.flush()
        s = start
        b = 1024*1024*10
        assert end <= tbl.glabel.shape[0]
        assert end <= tbl.pk_props.shape[1]
        while s < end:
            e = min( end, s+b )
            hin[dset]['pk_props'].read_direct( tbl.pk_props, 
                                              np.s_[:, s:e], np.s_[:, s:e] )
            hin[dset]['glabel'].read_direct( tbl.glabel, 
                                            np.s_[s:e], np.s_[s:e] )
            s = e
    total_bytes = tbl.pk_props[:,start:end].nbytes + tbl.glabel[start:end].nbytes
    del tbl
    return start, end

def create( hname, dset ):
    with h5py.File(hname, 'r') as hin:
        tbl = ImageD11.sinograms.properties.pks_table()
        tbl.npk = hin['pks2d/npk'][:]
        tbl.nlabel = hin['pks2d'].attrs['nlabel']
        tbl.ipk = hin['pks2d/ipk'][:]                         
        for name in 'pk_props', 'glabel':
            dset = hin['pks2d'][name]
            setattr( tbl, name, tbl.share( name, 
                                          shape = dset.shape, 
                                          dtype = dset.dtype ) )
        del dset

    return tbl

def pread( tbl, hname, dset='pks2d', nproc=0  ):
    if nproc == 0:
        nproc = cImageD11.cores_available() 
        print("Using", nproc, 'processes to read')
    if tbl is None:
        tbl = create( hname, dset )
    shmd = tbl.export()
    n = shmd[ 'glabel' ][ 'shape' ][0]
    se = list(range(0,n,n//nproc))
    se[-1] = n
    args = [(hname, dset, shmd, start, end) for start, end in zip(se[:-1],se[1:])]
    nbytes = 0
    with multiprocessing.Pool( nproc ) as pool:
        for rb in pool.map( read_segment, args ):
            pass
    return tbl

def read2pipe( dest, hname, dset, selection ):
    with h5py.File( open( hname, 'rb' ), 'r' ) as hin:
        dest.send_bytes( hin[ dset ][ selection ][:].tobytes() )
    dest.close()
    
def readpartial( args ):
    ary, hname, dset, selection = args
    write_end, read_end = multiprocessing.Pipe()
    p = multiprocessing.Process( target = read2pipe, args = ( write_end, hname, dset, selection ) )
    p.start()
    output = ary[ selection ] 
    output[:] = np.frombuffer( read_end.recv_bytes(), dtype = output.dtype ).reshape(output.shape)
    p.join()
    read_end.close()
    return selection
    

def getse( shape, chunks, nproc ):
    dim = len(shape)
    n = [ (s+c-1) // c for s,c in zip(shape, chunks)]
    print(n, 'chunks in the dataset', shape, chunks )
    ax = np.argmax(n)
    nperproc = n[ax]//nproc
    start = 0
    block = chunks[ax]
    slices = []
    for p in range(nproc):
        slc = []
        for axis in range(len(shape)):
            if axis == ax:
                end = start + block * nperproc
                slc.append( slice( start, min(end, shape[ax] ) ) )
                start = end
            else:
                slc.append( slice(None, None, None) )
        slices.append( tuple(slc) )
    return slices
    

def readN( hname, dset, NPROC=40 ):
    with h5py.File(hname, 'r') as hin:
        ds = hin[dset]
        output = np.empty( ds.shape, ds.dtype )
        chunks = ds.chunks
    se = getse(output.shape, chunks, NPROC)
    print(se)
    with concurrent.futures.ThreadPoolExecutor( max_workers=NPROC ) as pool:
        jobs = {}
        for s in se:
            jobs[ pool.submit( readpartial, (output, hname, dset, s ) ) ] = s
        for future in concurrent.futures.as_completed(jobs):
            print(future.result())
        

if __name__=="__main__":

    hname = '/data/projects/3dipolyplast/NS_Analysis_jw/et12_z95_2/20230405/pks_ET12_7ns_slice_pz15_2.h5'
    import timeit
    
    start = timeit.default_timer()
    tbl = pread( None, hname )
    end = timeit.default_timer()
    print("Reading took {end-start:.2f}/s".format(**locals()))
    
    # Single threaded
    start = timeit.default_timer()
    chk = ImageD11.sinograms.properties.pks_table.load( hname )
    end = timeit.default_timer()
    print("One core reading took {end-start:.2f}/s".format(**locals()))

    start = timeit.default_timer()
    assert (tbl.pk_props == chk.pk_props).all()
    assert (tbl.glabel == chk.glabel).all()
    end = timeit.default_timer()
    print("Matches {end-start:.2f}/s".format(**locals()))

