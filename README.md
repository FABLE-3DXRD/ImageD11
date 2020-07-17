Windows: [![Build status](https://ci.appveyor.com/api/projects/status/4pdlvsj2grtk0hel?svg=true)](https://ci.appveyor.com/project/jonwright/imaged11)

Linux: [![CircleCI](https://circleci.com/gh/jonwright/ImageD11.svg?style=svg)](https://circleci.com/gh/jonwright/ImageD11)

Macos + Linux [![Build Status](https://travis-ci.com/jonwright/ImageD11.svg?branch=master)](https://travis-ci.com/jonwright/ImageD11)

ImageD11
Version 1.9.7
Jon Wright
wright@esrf.fr

This is the source code for ImageD11. Probably you wanted a compiled version.

Perhaps you can try to get it via:

 pip install ImageD11

Some (dated) documentation is here: https://imaged11.readthedocs.io/

If you are at ESRF on a linux computer you could try "module load fable"

To use from git, try this:

 Download and install python 3.6+, perhaps from www.python.org
 but probably from conda.
 Add the packages numpy, scipy, matplotlib, h5py, pillow, pycifrw, xfab

 Then try:
 $ git clone https://github.com/myusername/ImageD11.git && cd ImageD11
 $ python setup.py build bdist_wheel
 Followed installation in your virtual or conda environment:
 $ pip install dist/ImageD11_..[version you built]..whl


After it is installed, you should find a script ImageD11_gui.py, somewhere in your path.

Until 2017 this code was mostly developed on sourceforge at http://sourceforge.net/projects/fable/ 

It is now developed at http://github.com/FABLE-3DXRD/ImageD11 

Bug reports are always welcome!

Good luck!





