Getting ImageD11
================

ImageD11 is a python module which depends on the prior installation 
python itself and several other packages. It is usually better to 
get these from a package manager or in a batteries-included python 
installation. Perhaps Anaconda or Python(x,y) etc.

General purpose packages that are required:

- `Python itself <http://python.org/download/>`_ (2.5, 2.6 or 2.7)
- `numpy <http://www.scipy.org/Download>`_ for numerical arrays
- `matplotlib <http://matplotlib.org/downloads.html>`_ for plotting
- `PIL <http://www.pythonware.com/products/pil/>`_, the python imaging library
- `PyOpenGl <http://pyopengl.sourceforge.net/>`_ (with the Togl widget)

Specific packages from fable:

- fabio (from fable, now silx) for imageio
- xfab (from fable) for some crystallographic calculations

Recommended packages to pick up when you are installing a system:

- scipy for many different scientific computations
- h5py for access to hdf files
- pyminuit optimiser (needed for fitallb)
- pmw python metawidgets for fabian gui
- pyFAI for radial integration 

For development you also need:

- C compiler for your system
- opencl system 
- sphinx docutils to build this documentation
- latex for the pdf manual

Finally, you need to download ImageD11 itself. The file release area was historically at 
`sourceforge <http://sourceforge.net/projects/fable/files/ImageD11>`_.
Nowadays the sources are on  `github <http://github.com/jonwright/ImageD11>`_ and the
package itself is on pypi. Install it using:

.. 
  pip install ImageD11



Linux platforms where you have root access
------------------------------------------

Use the package manager to install pre-built versions of the main python 
packages. Python itself usually comes with the system.
There are some hints at the `scipy site
<http://scipy.github.com/download.html>`_ about how the packages are named for
different vendors. 

On debian based systems (also ubuntu etc):

..
  sudo apt-get install python-numpy python-numpy-dev python-image python-scipy python-matplotlib python-opengl build-essential

To install the bleeding edge version do an svn checkout:

.. 
  svn co http://fable.svn.sourceforge.net/projects/ImageD11/trunk fabio
  svn co http://fable.svn.sourceforge.net/projects/ImageD11/trunk xfab
  svn co http://fable.svn.sourceforge.net/projects/ImageD11/trunk ImageD11

Then you can install system wide via:

..
  cd fabio
  python setup.py build install
  cd ../xfab
  python setup.py build install
  cd ../ImageD11
  python setup.py build install
 
Or in some local area by adding that location as --prefix

..
  export PATH=$PATH:/somewhere/bin export PYTHONPATH=/somewhere/lib/python2.6/site-packages
  python setup.py build install --prefix=/somewhere 

Historically the pyopengl builds are broken for the Togl widget. 
Get it from Togl.sourceforge.net and do something like the following:

..
 tar -zxf Togl2.0-8.4-Linux.tar.gz 
 cd Togl2.0-8.4-Linux/lib 
 cp -r Togl2.0 /usr/lib/pymodules/python2.6/OpenGL/Tk/linux2-tk8.5

  



Linux platforms where you don't have root access
------------------------------------------------

This is possible, but surprisingly difficult. The best solution is to go talk
to your system administrator and ask them to help. It is really much easier to
install via the system package manager. A few keystrokes and you have everything
installed for free.

For linux there are some prebuilt versions with batteries included:

- Enthought python distribution http://www.enthought.com/products/epd.php.
- Anaconda https://store.continuum.io/cshop/anaconda

Otherwise build from source code. You download the source distribution for each
dependency and then compile it and install it. 
Usually you will need to compile further dependencies for each of the 
packages you are installing. It can be done, and has some value for 
education and character building, but it is not worth the effort.
If you can succeed to do that you'll probably now agree you were better off 
to talk to the sysadmin or buy the Enthought version.

It should be possible to install using "easy_install" or "pip" via
the pypi package index. FIXME: add the command for doing that.

Windows
-------

The easiest solution is to install a "batteries included" distribution of
python. For example:

- The fable gui http://sourceforge.net/projects/fable
- Enthought python distribution http://www.enthought.com/products/epd.php
- pythonxy www.pythonxy.com
- Anaconda https://store.continuum.io/cshop/anaconda

The author usually installs the python binaries from www.python.org and then 
add on the extension packages from each individual project if they are available. 
Otherwise Christoph Gohlke maintains a 
`handy collection <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_ of 
compiled modules, especially for 64bit windows systems, in case you 
cannot find a build from the upstream project. 

To develop ImageD11 on windows 64 bit you need to get hold of the right
compiler. This was freely available from microsoft, but a little non-obvious.
You need to get VC++ 2008 (msvcrt9, the free 2008 express edition is OK) and
then add on the .net sdk 3.5 to get the 64bit compiler. See the notes on the
`cython wiki <http://wiki.cython.org/64BitCythonExtensionsOnWindows>`_.

Mac
---

ImageD11 has been installed by a number of mac owners but the author is not 
sure how they did it. FIXME...


