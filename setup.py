
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

To re-build the wrappers do:
 cd src && python make_pyf.py
"""

from io import open # this misery may never end.

# For pip / bdist_wheel etc
import setuptools
import os, sys, platform, os.path
from distutils.core import setup, Extension
from distutils.command import build_ext
from numpy import get_include

############################################################################
def get_version():
    with open("ImageD11src/__init__.py","r") as f:
        for line in f.readlines():
            if line.find("__version__")>-1:
                return eval(line.split("=")[1].strip())

print("Building version |%s|"%get_version(), "on system:", platform.system())

#############################################################################
# Set the openmp flag if needed. Also CFLAGS and LDSHARED from sys.argv ?
#
# JW https://stackoverflow.com/questions/724664/python-distutils-how-to-get-a-compiler-that-is-going-to-be-used


copt =  {
    'msvc': ['/openmp', '/O2'] , 
    'unix': ['-fopenmp', '-O2'], #, '-DF2PY_REPORT_ON_ARRAY_COPY=100'] , 
    'mingw32': ['-fopenmp', '-O2'] , 
 }

lopt =  { k : [a for a in l] for k,l in copt.items() }
lopt['msvc'] = []
if platform.system() == "Darwin":
    copt['unix'].remove("-fopenmp")
    lopt['unix'].remove("-fopenmp")


# might try:
# set CFLAGS=/arch:AVX2 for msvc
# CFLAGS=-march=native -mtune=native
# LDFLAGS=-march=native -mtune=native


class build_ext_subclass( build_ext.build_ext ):
    def build_extensions(self):
        c = self.compiler.compiler_type
        CF = [] ; LF=[]
        if "CFLAGS" in os.environ:
            CF = os.environ.get("CFLAGS").split(" ")
        if "LDFLAGS" in os.environ:
            LF = os.environ.get("LDFLAGS").split(" ")
        for e in self.extensions:
            if c in copt:
               e.extra_compile_args = copt[ c ] + CF
               e.extra_link_args = lopt[ c ] + LF
        print("Customised compiler",c,e.extra_compile_args,
                    e.extra_link_args)
        build_ext.build_ext.build_extensions(self)

cnames =  "_cImageD11module.c blobs.c cdiffraction.c cimaged11utils.c"+\
" closest.c connectedpixels.c darkflat.c localmaxlabel.c sparse_image.c "+\
" splat.c fortranobject.c"
csources = [os.path.join('src',c) for c in cnames.split()]

extension = Extension( "_cImageD11",  csources,
                include_dirs = [get_include(), 'src' ])

################################################################################


# Try to further reduce this long list
scripts = ["ImageD11src/rsv_mapper.py",
           "ImageD11src/tkGui/plot3d.py",
                 "scripts/peaksearch.py",
                 "scripts/fitgrain.py",
                 "scripts/ubi2cellpars.py",
                 "scripts/filtergrain.py",
                 "scripts/ImageD11_2_shelx.py",
                 "scripts/fix_spline.py",
                 "scripts/edfheader.py",
                 "scripts/ImageD11_gui.py",
                 "scripts/bgmaker.py",
                 "scripts/merge_flt.py",
                 "scripts/makemap.py",
                 "scripts/plotlayer.py",
                 "scripts/plotgrainhist.py",
                 "scripts/plotImageD11map.py",
                 "scripts/cutgrains.py",
                 "scripts/index_unknown.py",
                 "scripts/spatialfix.py",
       	         "scripts/refine_em.py",
                 "scripts/avg_par.py",
                 "scripts/powderimagetopeaks.py"]
                 
################################################################################
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
    "pycifrw",
    # breaks travis for macos ?? "silx",
    "pyopengl",
    "pyopengltk",
    ]
    # , "FitAllB", ... ]

 # read the contents of your README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    readme = f.read()

    
# See the distutils docs...
setup(name='ImageD11',
      version=get_version(),
      author='Jon Wright',
      author_email='wright@esrf.fr',
      cmdclass={'build_ext': build_ext_subclass},
      description='ImageD11',
      license = "GPL",
      ext_package = "ImageD11",   # Puts extensions in the ImageD11 directory
      ext_modules = [extension,],
      setup_requires = ['numpy'], # to compile
      install_requires = needed,
      packages = ["ImageD11",
                  "ImageD11.tkGui",
                  "ImageD11.silxGui",
                  "ImageD11.nbGui"],
      package_dir = {"ImageD11":"ImageD11src"},
      url = "http://github.com/jonwright/ImageD11",
      package_data = {"ImageD11" : ["doc/*.html", "data/*", "sandbox/*" ]},
      scripts = scripts,
      long_description = readme,
      long_description_content_type='text/markdown',
)


