


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





from distutils.core import setup, Extension

cl = Extension("closest", sources=['src/closest.c'])
cp = Extension("connectedpixels", sources = ['src/connectedpixels.c'])
# Fortran hack
import os
os.system("g77 -c src/bispev.f -o src/bispev.o")
os.system("ar -cr src/libsplines.a src/bispev.o")
bl = Extension("_splines", sources = ['src/splines.c'], 
                     libraries = ["splines"], library_dirs = ["src"] )


setup(name='ImageD11',
      version=0.4,
      author='Jon Wright',
      author_email='wright@esrf.fr',
      description='ImageD11',
      license = "GPL",
      ext_package = "ImageD11",
      ext_modules = [cl,cp,bl],
      packages = ["ImageD11"],
      scripts = ["scripts/peaksearch.py", "scripts/ImageD11_gui.py"])

      
      
"""py_modules = ["bisplev","blobcorrector","data","opendata",
                    "peakmerge","peaksearch","simplex",
                    "transform","unitcell","indexing",
                    "guiindexer","guimaker","guipeaksearch",
                    "guitransformer","listdialog","plot3d","twodplot"]
      )"""

