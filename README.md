Windows: [![Build status](https://ci.appveyor.com/api/projects/status/4pdlvsj2grtk0hel?svg=true)](https://ci.appveyor.com/project/jonwright/imaged11)

Linux: [![CircleCI](https://circleci.com/gh/jonwright/ImageD11.svg?style=svg)](https://circleci.com/gh/jonwright/ImageD11)

ImageD11
Version 1.9.6
Jon Wright
wright@esrf.fr

This is the source code for ImageD11. Probably you wanted a compiled version.

Perhaps you can try to get it via:

 pip install ImageD11

Some (dated) documentation is here: https://pythonhosted.org/ImageD11/

If you are at ESRF on a linux computer you could try "module load fable"

To use it, try this:

 Download and install python 2.7, perhaps from www.python.org 
 Add the packages numpy, matplotlib and pyopengltk.

 Then try: 
 $  python setup.py build
 or
 >  python setup.py build --compiler=mingw32
 Followed by:
 pip install .

After it is installed, you should find a script ImageD11_gui.py, somewhere in your path.

Until 2017 this code was mostly developed on sourceforge at http://sourceforge.net/projects/fable/ 

It is now developed at http://github.com/FABLE-3DXRD/ImageD11 




Good luck!





