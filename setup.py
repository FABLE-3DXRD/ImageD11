
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
"""

import setuptools
import sys, sysconfig, platform
from numpy.distutils.core import setup, Extension
from numpy import get_include
import struct

# For 32 or 64 bits
nbyte = struct.calcsize("P") # 4 or 8

# We use size_t somewhere so it can address >2Gb on 64 bit
# systems. This doesn't work on 32 bit, so we decide what
# to do according to the python version building the extension
def fix_f2py_pointer(nbyte):
    # decide if we have 32 bits or 64 bits compilation
    with open("src/cImageD11_template.pyf","r") as f:
        pyf = f.read()
    pyf = pyf.replace(  "(kind=size_t)" ,"*%d"%(nbyte))
    # we also create versions for newer and older cpu
    pyfsse = pyf.replace(  "python module cImageD11" , "python module cImageD11_sse2" )
    with open("src/cImageD11_sse2.pyf", "w") as f:
        f.write( pyfsse )
    pyfavx = pyf.replace(  "python module cImageD11" , "python module cImageD11_avx2" )
    with open("src/cImageD11_avx2.pyf", "w") as f:
        f.write( pyfavx )
# And do this:
fix_f2py_pointer(nbyte)

# We have some files with sse2 or avx2 instructions
# These can only run on machines that actually have sse2 or avx2
#
# For distribution we will want to compile with or without and choose
# at runtime. Because of the way inline and linkers work it will be 
# safer to compiler *everything* with the same compiler options and
# have a whole extension module per option. Then decide the one to
# use at runtime.
# https://randomascii.wordpress.com/2016/12/05/vc-archavx-option-unsafe-at-any-speed/
# e.g: mixing options can make a mess in the long run
# 
# We therefore build cImageD11_XXX.[so|pyd]
# in __init__.py decide which one to load (requires a decision about sse2/avx)
#
# Potentially 3 versions: nothing, SSE2, AVX2
#  ... but we will only build for the second two for now, sse2 is really old

extn_kwds = {
    "include_dirs" : [get_include(), "src"],
    "extra_compile_args" :  ["-DF2PY_REPORT_ON_ARRAY_COPY", ],
}
sources = [ "src/connectedpixels.c",
            "src/closest.c",
            "src/cdiffraction.c",
            "src/localmaxlabel.c",
            "src/sparse_image.c",
            "src/blobs.c" ]
    

# Get base args from system (mostly linux):
if sysconfig.get_config_var("CFLAGS") is not None:
    extn_kwds["extra_compile_args"] += sysconfig.get_config_var("CFLAGS").split()

# MSVC compilers
if (platform.system() == "Windows") and ("--compiler=mingw32" not in sys.argv):
    extn_kwds["extra_compile_args"] += [ "/openmp", ]
    if nbyte == 4:
        extn_kwds["extra_compile_args"] += [ "/arch:SSE2", ]
    # else sse2 is always available on 64 bit windows    
    avx2_kwds = extn_kwds.copy()
    #
    #    /arch:AVX2 option was introduced in Visual Studio 2013 Update 2, version 12.0.34567.1.
    # Visual C++    CPython
    # 14.0          3.5, 3.6
    # 10.0          3.3, 3.4
    # 9.0           2.6, 2.7, 3.0, 3.1, 3.2
    if sys.version_info[:2] > ( 3, 4 ):
        avx2_kwds["extra_compile_args"] += [ "/arch:AVX2", ]
        assert (len(avx2_kwds["extra_compile_args"])==len(extn_kwds["extra_compile_args"])+1)
    else:
        print("Warning: your compiler does not have AVX2, try mingw32 instead")
# gcc compilers
elif (platform.system() == "Linux") or ("--compiler=mingw32" in sys.argv):
    extn_kwds["extra_compile_args"] += ["-fopenmp", "-O2", "-std=c99", "-msse2"]
    extn_kwds["extra_link_args"] = extn_kwds["extra_compile_args"]
    extn_kwds["libraries"] = ("gomp","pthread")
    avx2_kwds = extn_kwds.copy()
    avx2_kwds["extra_compile_args"] = extn_kwds["extra_compile_args"] + ["-mavx2",]
    assert (len(avx2_kwds["extra_compile_args"])==len(extn_kwds["extra_compile_args"])+1)
    avx2_kwds["extra_link_args"] = avx2_kwds["extra_compile_args"]
else:
    raise Exception("Sorry, your platform/compiler is not supported")

extn_kwds['sources'] = ["src/cImageD11_sse2.pyf",] + sources
avx2_kwds['sources'] = ["src/cImageD11_avx2.pyf",] + sources

# Compiled extensions:
extensions = [ Extension( "cImageD11_sse2", **extn_kwds), 
               Extension( "cImageD11_avx2", **avx2_kwds) ]


# Removed list of dependencies from setup file
#  borked something, not sure what or why
needed =[]
# ["xfab", "fabio", "pyopengl",  "matplotlib", "numpy", "scipy", "six", "h5py",
#  "pyopengltk", "FitAllB", ... ]

# See the distutils docs...
setup(name='ImageD11',
      version='1.9.1',
      author='Jon Wright',
      author_email='wright@esrf.fr',
      description='ImageD11',
      license = "GPL",
      ext_package = "ImageD11",   # Puts extensions in the ImageD11 directory
      ext_modules = extensions,
      install_requires = needed,
      packages = ["ImageD11"],
      package_dir = {"ImageD11":"ImageD11"},
      url = "http://github.com/jonwright/ImageD11",
      package_data = {"ImageD11" : ["doc/*.html", "data/*" ]},
      scripts = ["ImageD11/rsv_mapper.py",
                 "scripts/peaksearch.py",
                 "scripts/fitgrain.py",
                 "scripts/tomapper.py",
                 "scripts/ubi2cellpars.py",
                 "scripts/filtergrain.py",
                 "scripts/filterout.py",
                 "ImageD11/plot3d.py",
                 "scripts/pars_2_sweeper.py",
                 "scripts/ImageD11_2_shelx.py",
                 "scripts/fit2dcake.py",
                 "scripts/fix_spline.py",
                 "scripts/edfheader.py",
                 "ImageD11/plot3d.py",
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


