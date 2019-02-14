

# ImageD11_webgui : Feb 2019

We would like to make a graphical interface to access the pieces of ImageD11 in a way that helps new users to get started quickly while also being useful for experts and doing things like managing batch jobs. 

This is also a moment to review what is in the program in terms of functionality and what to do to improve it.

## What should it do ?

### background and peaksearching from images

- open a single image (via fabio) and show it on screen
- apply corrections (dark, flat) to see they are done right
- display 1D line cuts, numbers inside pixels, image statistics
- estimate the background (single image), e.g. via filters on radial transform
- create a file series (sequence of images) for rotation scan or difftomo or serial crystallography
- estimate the background from a sequence, e.g. via filters on image number
- single threshold blob identification (was peaksearch.py)
- multi-threshold merging (was merge_flt.py)
- localmax labeling algorithm (not yet user visible)
- saving output : ...

### Peak search output

Per frame. Should be much smaller than the frame itself. Can back reference to the frame data.

For each spot: 
-   Pixel indices?
-   slow / fast / sum_i / sum_1 / bounding_box
-   implicit "assigned to background" 
-   1D integration ? 

Typical case. Difftomo at 10 fps for 48 hours, 1.728x10^6 images, 13 Tb for raw unsigned short.

Some data storage notes for pixel indices. Per peak requirement appears to be of the order 100 bytes assuming compression:
    >>> import numpy, zlib
    >>> im = numpy.zeros((2048,2048),numpy.uint8)
    >>> im[40:80,400:450]=1
    >>> indices = numpy.arange( im.shape[0]*im.shape[1] )[ im.ravel() > 0 ]
    >>> dI = indices.copy()
    >>> dI[1:]= indices[1:]-indices[:-1]
    >>> print( "Full image",len( im.ravel().tostring() ) )
    Full image 4194304
    >>> print( "Full image",len(zlib.compress( im.ravel() ) ) )
    Full image 4177
    >>> print( "Indices as float32",len(zlib.compress( indices.astype(numpy.uint32) ) ) )
    Indices as float32 2294
    >>> print( "Diff indices",len(zlib.compress( dI.astype(numpy.uint32) ) ) )
    Diff indices 81

This is 10 bytes for some hypothetical peak. If we use the localmaxlabel code to label all pixels as belonging to a local maximum then we run into a problem (144K peaks for some typical image and labels compresses to 680 Kb). There must be a threshold to decide whether a peak should be stored or not and so there must be a way to decide what is "signal" and what is "noise". 







Per scan dimension (over omega, diffractometer_y, sample_z, etc)
- Links to say peaks are overlapping, thought to the the same



### calibration and scattering vectors

- powder style : beam centre, tilts, distance, etc. ImageD11_gui/transform or fable/transformation
- Friedel pairs : to be added to fit wedge / tilt_x
- more generalised geometry descriptions. Use imgCIF/CBF, see, e.g.:

https://sites.google.com/site/nexuscbf/mapping-draft/functional-mapping

loop_
 _axis.id 
 _axis.depends_on 
 _axis.equipment
 _axis.offset[1] 
 _axis.offset[2] 
 _axis.offset[3] 
 _axis.type
 _axis.vector[1] 
 _axis.vector[2] 
 _axis.vector[3]
y_pixel_offset GONIOMETER_OMEGA detector -38.6 39.95 0 translation 0 -1 0
x_pixel_offset y_pixel_offset detector 0 0 0 translation 1 0 0
DETECTOR_Z '.' detector 0 0 0 translation 0 0 -1
GONIOMETER_OMEGA DETECTOR_Z detector 0 0 0 rotation 0 0 -1
GONIOMETER_PHI GONIOMETER_KAPPA goniometer 0 0 0 rotation 1 0 0
GONIOMETER_KAPPA '.' goniometer 0 0 0 rotation 0 -1 0 

/detector:/NXdetector
 /transformations:NXtransformations
  DATASET "GONIOMETER_OMEGA" = [0]
  @depends_on="/entry/instrument/detector/transformations/DETECTOR_Z”
  @equipment="detector
  @long_name="GONIOMETER_OMEGA
  @transformation_type="rotation"
  @units="degrees
  @vector=[0, 0, 1] 

### 

## Historical background 

This seems to be a hard thing to make. Historically we had/have:

- Tkinter for a gui with: matplotlib for 2D plotting, pyopengl/togl for 3D plotting, and a home-made macro recorder. Still runs but now uses home-made pyopengltk for 3D which still has opengl problems.

- Fable/RCP java eclipse based gui that also hosted the rest of the fable suite. Initially use xmlrpc to call python and then changed to a jepp embedded python interpreter. Nice but no java progammers currently active on the project.

- This "webgui" folder was a proof-of-concept for a webgl widget, but was blocked web security preventing files being opened (no simple file selector).

- The silx_guis folder has a few mini utilities using PyQT based widgets. This is the toolbox being developed at ESRF.

## Options

### PyQT

This seems to be the good choice for cross platform desktop applications. Jon would need to learn QT. Easy to get expert help here at ESRF. A concern is that C++ is first class and there is more than one python binding (PyQT3, PyQT4, PyQT5, PySide, PySide2).

There are plenty of recipes for making a program that someone goes and runs. Probably it will have a web backend in the future, if not already.

There are alot of developments around silx at ESRF that seem to use pyqt for gui. Probably you can use silx without using the gui.

### Jupyter

This is well funded outside and seems to see wide adoption. Perhaps silx will target it.

There are several options:
- jupyterlab "next generation web-based user interface"
- jupyternotebook. Presumably the thing making .ipynb files
- jupyterconsole. Looks like it talks to backends
- qtconsole. Python console with inline plots

It has many nice widgets:
- pythreejs widget offering 3D
- bqplot / matplotlib / altair / bokeh  for 2D plotting
- csv viewer for large files (can scroll million rows)

For batch jobs and clusters there is:
- ipyparallel (uses zeromq)
- dask.distributed and/or dask-jobqueue can send jobs to a cluster.

It was less obvious how someone would "run a program" in the more traditional sense, but this seems to be a thing too:
 https://github.com/jupyterlab/jupyterlab/blob/master/examples/filebrowser
 https://github.com/jupyterlab/jupyterlab/tree/master/examples/app

Jupyterlab looks a bit like eclipse/rcp running in a web browser. Written in typescript, it seems to be a new version of ipython for the web.

Drawbacks / concerns:
- Installation and distribution fears : do you just need conda + recent web browser ?
- Fast moving technologies
- Can you see it running on OAR from a debian desktop ?
- Can you make a standalone exe ?
- ipynb has some drawbacks as you can edit cells that were previously run

### Bokeh

Youtube video : https://www.youtube.com/watch?v=HmI1foA0MZc

Developing Dashboard Applications Using Bokeh - Bryan Van de Ven

### Home-made web application

This got stuck and abandoned in the past, it was just too difficult for a one person effort.

### Tkinter

This still runs and is used. It seems to be very difficult to add things to it and maintain it.






