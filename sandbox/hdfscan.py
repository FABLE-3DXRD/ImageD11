
"""
Wrapper to map hdf onto the old peaksearch code

Try to decide if this belongs in peaksearch or if we make something new
"""


from fabio.fabioimage import fabioimage
import hdf5plugin
import h5py

filename = "frelon.hdf"
dataset  = '/entry/measurement/data'

def fso( filename, dataset ):
    h  = h5py.File( filename, "r" )
    omega = h['/entry/measurement/omega'][:]
    nframes = len(omega)
    data = h[dataset]
    assert data.shape[0]==nframes
    header = {'Omega' : omega[0] }
    filename = "%s::%s"%(filename, dataset) 
    def frm(i):
        header['Omega'] = omega[i]
        f = fabioimage( data[i], header )
        f.filename = filename
        f.currentframe = i
        return f
    
    #
    yield( frm(0) ) # first
    for i in range(nframes):
        yield frm(i)
        
file_series_object = fso( filename, dataset )
first_image = next(file_series_object)


if __name__=="__main__":
    from timeit import default_timer as timer
    t0 = timer()
    for i in file_series_object:
        pass
    t1 = timer()
    print( "Read time", t1-t0 )



