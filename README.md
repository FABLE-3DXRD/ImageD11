

ImageD11 is a python code for identifying individual grains in spotty area detector X-ray diffraction images.

Version 2.1.2, Jon Wright, wright@esrf.fr

This is the source code for ImageD11. Probably you wanted a compiled version.

If your pip is up-to-date, you can try to install it like this:
```
 python -m pip install --upgrade pip setuptools
 python -m pip install ImageD11
```
If you want to use an existing numpy installation add a `--no-build-isolation` flag.

To get all the possible dependencies too, you can try:
 `python -m pip install ImageD11[full]`

Some (dated) documentation is here: https://imaged11.readthedocs.io/

To use from git, try this:

 - Download and install python 3.8+, perhaps from www.python.org but probably from conda.
 - Preload binary packages from conda (or your system package manager): 
    numpy, scipy, matplotlib, h5py, pillow, pycifrw, xfab, pyqt, pillow, silx[full] etc
 - `pip install git+https://github.com/FABLE-3DXRD/ImageD11.git`
 
If you want to work with the sources then you can try like this:
 ```
 $ python -m pip install --upgrade pip
 $ git clone https://github.com/FABLE-3DXRD/ImageD11.git && cd ImageD11
 $ python -m pip install --editable .
 ```

If you want multiple binaries in your home (on recent pythons) you can do and get the compiled code
for each platform in .so files that are labelled by platform. This is potentially useful for a
heterogeneous cluster (like at ESRF): 
```
  # on ppc64le:
  python3 -m pip install dist/ImageD11-1.9.8-cp38-cp38-linux_ppc64le.whl --user --ignore-installed
  # on x86_64:
  python3 -m pip install dist/ImageD11-1.9.8-cp38-cp38-linux_x86_64.whl --user --ignore-installed
  # etc
  # ~/.local/lib/python3.8/site-packages/ImageD11 % ls *.so
  _cImageD11.cpython-38-powerpc64le-linux-gnu.so  _cImageD11.cpython-38-x86_64-linux-gnu.so
```

After it is installed, you should find a script ImageD11_gui.py, somewhere in your path.

Until 2017 this code was mostly developed on sourceforge at http://sourceforge.net/projects/fable/ 

It is now developed at http://github.com/FABLE-3DXRD/ImageD11 

Bug reports are always welcome!

Good luck!

## CI Status

[![Flake, Build and PyTest](https://github.com/FABLE-3DXRD/ImageD11/actions/workflows/build_flake_pytest_ubuntu2004.yml/badge.svg)](https://github.com/FABLE-3DXRD/ImageD11/actions/workflows/build_flake_pytest_ubuntu2004.yml)

<!--

Windows: [![Build status](https://ci.appveyor.com/api/projects/status/4pdlvsj2grtk0hel?svg=true)](https://ci.appveyor.com/project/jonwright/imaged11)

Linux: [![CircleCI](https://circleci.com/gh/jonwright/ImageD11.svg?style=svg)](https://circleci.com/gh/jonwright/ImageD11)

Macos + Linux [![Build Status](https://travis-ci.com/jonwright/ImageD11.svg?branch=master)](https://travis-ci.com/jonwright/ImageD11)

-->
