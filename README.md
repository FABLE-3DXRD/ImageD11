

ImageD11 is a python code for identifying individual grains in spotty area detector X-ray diffraction images.

Version 1.9.8, Jon Wright, wright@esrf.fr

This is the source code for ImageD11. Probably you wanted a compiled version.

If your pip is up-to-date, you can try to install it like this (numpy is needed
to compile):
```
 python -m pip install --upgrade pip
 python -m pip install numpy
 python -m pip install ImageD11
```
To get all the possible dependencies too, you can try:
 `python -m pip install ImageD11[full]`

Some (dated) documentation is here: https://imaged11.readthedocs.io/

If you are at ESRF on an old linux computer you can try "module load fable". 

To use from git, try this:

 - Download and install python 3.7+, perhaps from www.python.org but probably from conda.
 - Preload binary packages from conda (or your system package manager): 
    numpy, scipy, matplotlib, h5py, pillow, pycifrw, xfab, pyqt, pillow, silx[full] etc
 - `pip install git+https://github.com/FABLE-3DXRD/ImageD11.git`
 
If you want to work with the sources then you can try like this:
 ```
 $ python -m pip install --upgrade pip
 $ git clone https://github.com/FABLE-3DXRD/ImageD11.git && cd ImageD11
 $ python -m pip install --editable .
 ```

After it is installed, you should find a script ImageD11_gui.py, somewhere in your path.

Until 2017 this code was mostly developed on sourceforge at http://sourceforge.net/projects/fable/ 

It is now developed at http://github.com/FABLE-3DXRD/ImageD11 

Bug reports are always welcome!

Good luck!

## CI Status

Windows: [![Build status](https://ci.appveyor.com/api/projects/status/4pdlvsj2grtk0hel?svg=true)](https://ci.appveyor.com/project/jonwright/imaged11)

Linux: [![CircleCI](https://circleci.com/gh/jonwright/ImageD11.svg?style=svg)](https://circleci.com/gh/jonwright/ImageD11)

Macos + Linux [![Build Status](https://travis-ci.com/jonwright/ImageD11.svg?branch=master)](https://travis-ci.com/jonwright/ImageD11)
