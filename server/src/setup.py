


# ImageD11_v0.4 Software for beamline ID11
# Copyright (C) 2005  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

"""
Setup script for python distutils package

Can be used to generate source and binary distributions

On windows (with mingw tools installed, from www.mingw.org:
    python setup.py build --compiler=mingw32 sdist bdist bdist_wininst

In Linux
    python setup.py build sdist

Once assumed you have the "g77" compiler available for fortran and that
the usual (gnu linux/win32 x86) name mangling are done for the splines
extension to work ... now that code went through f2c, so it should just
work ??

"""




from distutils.core import setup, Extension


# Compiled extensions:

# closest is for indexing grains
cl = Extension("closest", sources=['Imaged11Functions/closest.c'])

# connectedpixels is for peaksearching images
#cp = Extension("connectedpixels", sources = ['Imaged11Functions/connectedpixels.c'])

# histogramming thing
#ch = Extension("_hist", sources = ['Imaged11Functions/hist.c'])

# Fortran hack - just need to have libsplines.a available for linking
#import os
#os.system("g77 -c src/bispev.f -o src/bispev.o")
#os.system("ar -cr src/libsplines.a src/bispev.o")

# _splines is for correcting peak positions for spatial distortion
#bl = Extension("_splines", sources = ['Imaged11Functions/splines.c', 'Imaged11Functions/bispev.c'] )
#,
#                     libraries = ["splines"], library_dirs = ["src"] )


# See the distutils docs...
setup(name='ImageD11Gui',
      version='0.0.1',
      author='Jon Wright, Jaroslaw Butanowicz',
      author_email='wright@esrf.fr jaroslaw.butanowicz@esrf.fr',
      description='ImageD11Gui',
      license = "GPL",
      ext_package = "ImageD11Gui",   # Puts extensions in the ImageD11 directory
      ext_modules = [cl],
      packages = ["example.xmlrpc.sum.pyserver"],
      scripts = ["sumserver.py"])
