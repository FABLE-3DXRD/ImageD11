
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

You should run src/make_pyf.py to update the pyf wrapper and docs
"""
import sys
from io import open # this misery may never end.
# For pip / bdist_wheel etc
import os, platform, os.path
from setuptools import setup, Extension, find_packages
from setuptools.command import build_ext
#
import numpy, numpy.f2py   # force wrapper re-generation

if not hasattr( numpy.f2py, 'get_include'):
    numpy.f2py.get_include = lambda : os.path.join(
        os.path.dirname(os.path.abspath(numpy.f2py.__file__)),
        'src')

############################################################################
def get_version():
    with open("ImageD11/__init__.py","r") as f:
        for line in f.readlines():
            if line.find("__version__")>-1:
                return eval(line.split("=")[1].strip())

print("Building version |%s|"%get_version(), "on system:", platform.system())
############################################################################

#############################################################################
# Set the openmp flag if needed. Also CFLAGS and LDSHARED from sys.argv ?
#
# JW https://stackoverflow.com/questions/724664/python-distutils-how-to-get-a-compiler-that-is-going-to-be-used


copt =  {
    'msvc': ['/openmp', '/O2'] ,
    'unix': ['-fopenmp', '-O2', ], # '-DF2PY_REPORT_ON_ARRAY_COPY=100', '-DNPY_DISABLE_OPTIMIZATION=1' ] ,
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
    def build_extension(self, ext):
        print('Building _cImageD11 module')
        c = self.compiler.compiler_type
        CF = [] ; LF=[]
        if "CFLAGS" in os.environ:
            CF = os.environ.get("CFLAGS").split()
        if "LDFLAGS" in os.environ:
            LF = os.environ.get("LDFLAGS").split()
        if c in copt:
            ext.extra_compile_args = copt[ c ] + CF
            ext.extra_link_args = lopt[ c ] + LF
        print("Customised compiler",c,ext.extra_compile_args,
                    ext.extra_link_args)
        if ext.sources[0].endswith('.pyf'):
            name = ext.sources[0]
            # generate wrappers
            print('Creating f2py wrapper for', name)
            numpy.f2py.run_main( [
                #'--quiet',
                name,])
            ext.sources[0] = os.path.split(name)[-1].replace('.pyf', 'module.c')
            ext.sources.append( os.path.join(numpy.f2py.get_include(), 'fortranobject.c' ) )
        build_ext.build_ext.build_extension(self, ext)

# note that the pyf must come first
cnames =  "_cImageD11.pyf blobs.c cdiffraction.c cimaged11utils.c"+\
" closest.c connectedpixels.c darkflat.c localmaxlabel.c sparse_image.c "+\
" splat.c"

csources = [os.path.join('src',c) for c in cnames.split()]

extension = Extension( "ImageD11._cImageD11",  csources,
                include_dirs = [ 'src',numpy.get_include(), numpy.f2py.get_include() ])

################################################################################

# Try to further reduce this long list
scripts = ["ImageD11/rsv_mapper.py",
           "ImageD11/tkGui/plot3d.py",
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

# sfood -I ImageD11/depreciated -I build/ -v -u >sfood.out 2>sfood.err

minimal = [  # can't compile without this
    'numpy',
    "setuptools",
    ]

useful = [   # stuff you probably want, and should be able to get easily
    'fabio==0.2.2 ; python_version < "3" and sys_platform == "win32" ',
    'fabio ; python_version >= "3" or sys_platform != "win32" ',
    "xfab>=0.0.4", #
       # comes from xfab : "PyCifRW",
    "matplotlib",  # tkGui
    'pyopengl==3.1.7 ; python_version < "3"',
    'pyopengl ; python_version >= "3"',
    "pyopengltk",  # plot3d in tkGui
    "scipy",       #
    'hdf5plugin==1.4.1 ; python_version < "3"',   #  and sys_platform == "win32" ',
    'hdf5plugin ; python_version >= "3" ',
    'h5py',
    'pyyaml',
    "pytest",       # for the CI
    'numba==0.46.0 ; python_version < "3" ',       # for some test cases
    'llvmlite==0.30.0 ; python_version < "3" ',    # for some test cases
    'numba ; python_version >= "3" ',               # for some test cases
    "bslz4_to_sparse",
    "fast_histogram",
    "scikit-image",
    "tqdm",
    'threadpoolctl ; python_version >= "3" ',
    'orix ; python_version >= "3"',  # for orix interface
    'diffpy.structure ; python_version >= "3"'  # for orix interface
]


more = [
    # Used in sandbox / test / not completely essential, but should work for CI
    "papermill",   # in test for notebook testing
    "pillow",      # in sandbox
    "lmfit",       # in sandbox
    "sympy",       # for maths
    'ipywidgets ; python_version >= "3"',  # for notebook nbGui
    'pyopencl; python_version >= "3"',    # (was?) in sandbox
    'numexpr <= 2.7.0; python_version < "3" ',
    'pyFAI ; python_version >= "3" ',   # pypi problematic
    # 'pyFAI <= 0.18.0 ; python_version  < "3" ',
    'silx[full] ; python_version >= "3" ',  # for silxGui,
    'pandas'
]

rare = [           #
    "FitAllB",     # not for python3
    "minuit",      # for fitallb
    "PyTango",     # sandbox
    "PyMca5",      # in sandbox
    ]

 # read the contents of your README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    readme = f.read()


# See the distutils docs...
setup(name='ImageD11',
      version=get_version(),
      author='Jon Wright',
      author_email='wright@esrf.fr',
      cmdclass={ 'build_ext' : build_ext_subclass },
      description='ImageD11',
      license = "GPL",
#      python_requires='<3.12',  # Numba still not working for 3.12
      ext_modules = [extension,],
      setup_requires = minimal,   # to compile
      install_requires = minimal + useful,
      extras_require = { 'full' : more, 'rare' : rare },
      packages = find_packages( include=['ImageD11'
                                        ] ),
      package_dir = {"ImageD11": "ImageD11"},
      url = "http://github.com/jonwright/ImageD11",
      include_package_data=True,
      package_data = {"ImageD11.data" : [ "data/*"],
                      "ImageD11.nbGui" : ["ImageD11.nbGui/*.ipynb", ],
                      "ImageD11.nbGui.TDXRD" : ["ImageD11.nbGui.TDXRD/*.ipynb", ],
                      "ImageD11.nbGui.S3DXRD" : ["ImageD11.nbGui.S3DXRD/*.ipynb", ],
                      "ImageD11.nbGui.calibration" : ["ImageD11.nbGui.calibration/*.ipynb", ],
                     },
      scripts = scripts,
      long_description = readme,
      long_description_content_type='text/markdown',
)


