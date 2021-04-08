

# Make a project file in a hdf.

# "Scans" ... dimensions ?
#     omega / dty / pz / load / temperature ... etc

# Images ...

# example : 
#
import os
dataRoot = "/data/id11/nanoscope/Commissioning/2018Apr/difftomo_al_quick"
os.chdir( dataRoot )
#dataRoot = "./"
# interlaced = 1
# roundtrip = 1
# nimg = 360
# extensions per folder
imgs = []
for i in range(360):
    imgs.append("%d_%04d.edf"%(0,i))
    imgs.append("%d_%04d.edf"%(1,359-i))

stems = ["Al_y%03d_"%(i) for i in range(-60, 61)]

import h5py, fabio, os

h = h5py.File("/tmp/demo.hdf" , "w" )
for stem in stems:
    g = h.require_group(stem)
    fname = os.path.join( stem, stem+imgs[0] )
    im = fabio.open( fname )
    # hacky
    bsize = im._frames[0].size
    bstart = im._frames[0].start
    bshape = im.data.shape
    # here : verify headers
    extlist = [ (os.path.join( stem, stem+i ), bstart, bsize)
                for i in imgs ]
    g.create_dataset( "images",
                      shape = (len(imgs), bshape[0], bshape[1]),
                      dtype = im.data.dtype,
                      external = extlist )
    h.close()
                      
    
    



