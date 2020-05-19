
from __future__ import print_function

# ImageD11_v1.0 Software for beamline ID11
# Copyright (C) 2005-2018  Jon Wright
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
Setup script

Do:
  python setup.py build --force
   -> it will build static libraries in the src folder
   -> compile wrappers them here
"""

# For pip / bdist_wheel etc
import setuptools
import os, sys, platform
# Need to use numpy.distutils so it works for f2py
from numpy.distutils.core import setup, Extension
from numpy import get_include
from numpy.distutils.misc_util import mingw32
from distutils.command.sdist import sdist


def get_version():
    with open("ImageD11src/__init__.py","r") as f:
        for line in f.readlines():
            if line.find("__version__")>-1:
                return eval(line.split("=")[1].strip())

print("Building version |%s|"%get_version())

# Compiled extensions:

def eca():        
    """
    If you install on a personal machine then setting
    CFLAGS=-march=native -mtune=native
    might help?
    """
    if platform.system() == "Windows" and not mingw32():
        arg = [ "/O2", "/openmp", ]
    else:
        arg=["-O2", "-fopenmp", "-fPIC", "-std=c99" ]
    if "CFLAGS" in os.environ:
        arg += os.environ.get("CFLAGS").split(" ")
    return arg

def ela():
    """ Link args - just repeats compile for now """
    return eca()


# Only one built version. Too hard to make this cross plaform. See if CFLAGS
# can be used to add extras ?

cnames =  "blobs.c cdiffraction.c cimaged11utils.c closest.c " + \
 "connectedpixels.c darkflat.c localmaxlabel.c sparse_image.c splat.c"
csources = [os.path.join('src',c) for c in cnames.split()]

ekwds = { 'sources' : [ "./src/_cImageD11.pyf" ] + csources,
          'include_dirs' : [get_include(), 'src' ],
          'extra_compile_args' : eca(),
          'extra_link_args' : ela(),
        }

# print(ekwds)

# Name using platform.machine
extensions = [ Extension( "_cImageD11",  **ekwds ), ]

# Things we depend on. This generally borks things if pip
# tries to do source installs of all of these.

needed =[
    "six",
    "numpy",
    "scipy",
    "pillow",
    "h5py",
    "matplotlib",
    "xfab>=0.0.4",
    "fabio",
    "pyopengl",
    "pyopengltk",
    ]
    # , "FitAllB", ... ]

# See the distutils docs...
setup(name='ImageD11',
      version=get_version(),
      author='Jon Wright',
      author_email='wright@esrf.fr',
      cmdclass={'sdist': sdist},
      description='ImageD11',
      license = "GPL",
      ext_package = "ImageD11",   # Puts extensions in the ImageD11 directory
      ext_modules = extensions,
      install_requires = needed,
      packages = ["ImageD11", "ImageD11.tkGui", "ImageD11.silxGui"],
      package_dir = {"ImageD11":"ImageD11src"},
      url = "http://github.com/jonwright/ImageD11",
      package_data = {"ImageD11" : ["doc/*.html", "data/*" ]},
      scripts = ["ImageD11src/rsv_mapper.py",
                 "ImageD11src/tkGui/plot3d.py",
                 "scripts/peaksearch.py",
                 "scripts/fitgrain.py",
                 "scripts/tomapper.py",
                 "scripts/ubi2cellpars.py",
                 "scripts/filtergrain.py",
                 "scripts/pars_2_sweeper.py",
                 "scripts/ImageD11_2_shelx.py",
                 "scripts/fit2dcake.py",
                 "scripts/fix_spline.py",
                 "scripts/edfheader.py",
                 "scripts/huber2bruker.py",
                 "scripts/id11_summarize.py",
                 "scripts/ImageD11_gui.py",
                 "scripts/bgmaker.py",
                 "scripts/merge_flt.py",
                 "scripts/makemap.py",
                 "scripts/plotlayer.py",
                 "scripts/plotedf.py",
                 "scripts/plotgrainhist.py",
                 "scripts/rubber.py",
                 "scripts/plotImageD11map.py",
                 "scripts/cutgrains.py",
                 "scripts/index_unknown.py",
                 "scripts/spatialfix.py",
       	         "scripts/refine_em.py",
                 "scripts/avg_par.py",
                 "scripts/powderimagetopeaks.py"])


