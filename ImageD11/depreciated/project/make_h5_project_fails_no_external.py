

# Make a project file in a hdf.

# example : 
#
import os
import h5py, fabio
fmt = "/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y%03d_/interlaced_1_%d/Frelon/interlaced_1_%d_Frelon%04d.edf"
args = "ystep", "ipass", "ipass", "fnum"
ysteps = range(-60,61)
# interlaced
frames = []
for ystep in ysteps:
    for fnum in range(360):
        # forward
        frames.append( fmt % ( ystep, 1, 1, fnum ) )
        # and back
        frames.append( fmt % ( ystep, 2, 2, 359-fnum ) )
      
print("total frames:",len(frames))
for f in frames[:10]:
    print(f)
for f in frames[-10:]:
    print(f)
    
im = fabio.open( frames[0] )
# hacky
bsize = im._frames[0].size
boffset = im._frames[0].start
bshape = im.data.shape

with h5py.File("/tmp/demo.hdf" , "w" ) as h:
    g = h.require_group('sinogram')
    g.create_dataset( "data",
                      shape = (len(frames), bshape[0], bshape[1]),
                      dtype = im.data.dtype,
                      external = [(fname, boffset, bsize) for fname in frames] )

# outputs:
"""
total frames: 87120
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y-60_/interlaced_1_1/Frelon/interlaced_1_1_Frelon0000.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y-60_/interlaced_1_2/Frelon/interlaced_1_2_Frelon0359.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y-60_/interlaced_1_1/Frelon/interlaced_1_1_Frelon0001.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y-60_/interlaced_1_2/Frelon/interlaced_1_2_Frelon0358.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y-60_/interlaced_1_1/Frelon/interlaced_1_1_Frelon0002.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y-60_/interlaced_1_2/Frelon/interlaced_1_2_Frelon0357.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y-60_/interlaced_1_1/Frelon/interlaced_1_1_Frelon0003.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y-60_/interlaced_1_2/Frelon/interlaced_1_2_Frelon0356.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y-60_/interlaced_1_1/Frelon/interlaced_1_1_Frelon0004.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y-60_/interlaced_1_2/Frelon/interlaced_1_2_Frelon0355.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y060_/interlaced_1_1/Frelon/interlaced_1_1_Frelon0355.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y060_/interlaced_1_2/Frelon/interlaced_1_2_Frelon0004.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y060_/interlaced_1_1/Frelon/interlaced_1_1_Frelon0356.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y060_/interlaced_1_2/Frelon/interlaced_1_2_Frelon0003.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y060_/interlaced_1_1/Frelon/interlaced_1_1_Frelon0357.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y060_/interlaced_1_2/Frelon/interlaced_1_2_Frelon0002.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y060_/interlaced_1_1/Frelon/interlaced_1_1_Frelon0358.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y060_/interlaced_1_2/Frelon/interlaced_1_2_Frelon0001.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y060_/interlaced_1_1/Frelon/interlaced_1_1_Frelon0359.edf
/data/id11/nanoscope/Commissioning/sonja_fib_Al_z500__nobackup/difftomo_Al_y060_/interlaced_1_2/Frelon/interlaced_1_2_Frelon0000.edf

Traceback (most recent call last):
  File "make_h5_project_fails_no_external.py", line 40, in <module>
    g.create_dataset( "data",
  File "/usr/lib/python3/dist-packages/h5py/_hl/group.py", line 136, in create_dataset
    dsid = dataset.make_new_dset(self, shape, dtype, data, **kwds)
  File "/usr/lib/python3/dist-packages/h5py/_hl/dataset.py", line 167, in make_new_dset
    dset_id = h5d.create(parent.id, None, tid, sid, dcpl=dcpl)
  File "h5py/_objects.pyx", line 54, in h5py._objects.with_phil.wrapper
  File "h5py/_objects.pyx", line 55, in h5py._objects.with_phil.wrapper
  File "h5py/h5d.pyx", line 80, in h5py.h5d.create
ValueError: Unable to create dataset (object header message is too large)
"""
