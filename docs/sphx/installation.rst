===================
Installing ImageD11
===================

::

  pip install ImageD11

If you go for --user or systemwide the we suggest to skip the 
dependencies and fix them one at a time later as needed (pip can 
get into a fight with the system package manager or conda).

::

  pip install ImageD11 --no-deps -U --user

If that did not work, then read on.

ImageD11 is a python module which depends on the prior installation
python itself and several other packages. It is usually better to
get these from a package manager or in a batteries-included python
installation. Perhaps Anaconda or Python(x,y) etc.

Several dependency packages that are required are listed in the setup.py 
file of ImageD11 as `install_requires`. 

Specific packages from fable that ImageD11 needs are:

- fabio (from fable, now silx) for imageio
- xfab (from fable) for some crystallographic calculations

Recommended packages to pick up when you are installing a system:

- numpy for arrays of numbers
- matplotlib for plotting
- scipy for many different scientific computations
- h5py for access to hdf files
- pyminuit optimiser (needed for fitallb)
- pmw python metawidgets for fabian gui
- pillow (PIL) imaging library
- pyFAI for radial integration 
- pyopengl and pyopengltk for 3D plots
- six for the python2 vs python3 breakup

For development you would also need:

- C compiler for your system
- sphinx docutils to build this documentation
- latex if you want to make a pdf manual

Finally, you need to download ImageD11 itself. The file release area was 
historically at 
`sourceforge <http://sourceforge.net/projects/fable/files/ImageD11>`_.
Nowadays the sources are on  `github <http://github.com/jonwright/ImageD11>`_ 
and the package itself is on pypi. Install it in a virtualenv using:



Getting python 
--------------

Historically we had a lot of notes here about installing python.

These are removed for now. It is a bigger problam that just for us. You can
either go with the package manager for your system, or the version from
python.org, or a batteries included system like miniconda or anaconda. Or you
have a mac and all bets are off.

Whatever you do, please learn about virtualenv or conda env before you install
something. Python environments are relatively fragile and as soon as you install
one broken package the whole thing will be ruined. 

Python2 versus Python3
----------------------

A wide range of python versions are routinely tested by the CI services, including
python2.7. From about python 3.5 and upwards we hope it will be working, but there
may be difficulties with the earliest python 3 iterations.

Compiling from source
---------------------

There are some continuous integration scripts in the source which can give some
inspiration. See .circleci, .travis.yml and appveyor.yml

The C wrapper codes are interfaced using f2py in numpy. 
If you modify the C sources then the wrappers are re-built using a
script that adds some docstrings. Critical files are currently (May 2020) ::

  src/makepyf.py : builds interface declarations in 
  src/_cImageD11.pyf : interface description
  src/_cImageD11module.c : interface made by f2py
  ImageD11src/cImageD11.py : wrapper that adds docstrings
  ImageD11src/cImageD11_docstrings.py
  setup.py : drives the build

To optimise for a specific machine then setup.py should look at your CFLAGS
environment variable for anything you want to add::

  export CFLAGS=-march=native
  pip install .

Please test if you do this and remember not to copy it to another machine. 

Testing
-------

It should be possible to test using pytest::

  python -m pytest

Otherwise in the test folder there is a run_tests.py script.


Windows
-------

To develop ImageD11 on windows 64 bit you need to get hold of the right
compiler. This was freely available from microsoft, but a little non-obvious.
You need to get VC++ 2008 (msvcrt9, the free 2008 express edition is OK) and
then add on the .net sdk 3.5 to get the 64bit compiler. See the notes on the
`cython wiki <http://wiki.cython.org/64BitCythonExtensionsOnWindows>`_. Nowadays
ths mostly just works as microsoft started releasing an old compiler that was
needed for 2.7. For newer python versions you need a newer compiler.

Mac
---

ImageD11 has been installed by a number of mac owners but the author was never
sure how they did it. Recently a CI has been set up on Travis that tests For
macos as well, so hopefully this will improve in the future. The main missing 
point is the 3D opengl plots which need pyopngltk to be ported to mac, or perhaps
easier to find some other package to do then job.

Ubuntu
------
To install on Ubuntu ::

     sudo apt-get install build-essential
     sudo apt-get install git
     sudo apt-get install python-numpy 
     git clone http://github.com/jonwright/ImageD11
     cd ImageD11/
     python setup.py build bdist_wheel
     pip install dist/ImageD11-1.7.0-cp27-cp27mu-linux_x86_64.whl 
     cd ImageD11/
     cd test/
     python run_tests.py 

Installing at ESRF
------------------

Historically various ways. Currently (may 2020)::
   
  debian : module load fable
  ubuntu : make your own conda or virtual environment for now
