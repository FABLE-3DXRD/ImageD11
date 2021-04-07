

ImageD11 is a python code for identifying individual grains in spotty area detector X-ray diffraction images.

Version 1.9.8, Jon Wright, wright@esrf.fr

This is the source code for ImageD11. Probably you wanted a compiled version.

You can try to get it via:

 `pip install ImageD11`

Some (dated) documentation is here: https://imaged11.readthedocs.io/

If you are at ESRF on an old linux computer you can try "module load fable". 

To use from git, try this:

 - Download and install python 3.7+, perhaps from www.python.org but probably from conda.
 - Preload packages from conda (or your system package manager): numpy, scipy, matplotlib, h5py, pillow, pycifrw, xfab
 - `pip install git+https://github.com/FABLE-3DXRD/ImageD11.git`
 
If you want the sources then checkout like this:
 ```
 $ git clone https://github.com/FABLE-3DXRD/ImageD11.git && cd ImageD11
 $ python setup.py build bdist_wheel
 Followed by installation in your virtual or conda environment:
 $ pip install dist/ImageD11_[version you built].whl
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
