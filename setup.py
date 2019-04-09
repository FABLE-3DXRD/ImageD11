
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
import os, sys
# Need to use numpy.distutils so it works for f2py
from numpy.distutils.core import setup, Extension
from numpy.distutils.command.build_ext import build_ext
from numpy import get_include
# Get the path for the static libraries        
import src.bldlib

if "build" in sys.argv:
    """ Ugly - requires setup.py build 
    Should learn to use build_clib at some point
    """
    os.chdir("src")
    os.system("python write_check.py")
    ok = os.system("python bldlib.py")
    os.chdir("..")
    if ok != 0:
        print("Returns was",ok)
        sys.exit(ok)

class custom_build_ext( build_ext ):    
    """ Set compiler switches per extension module (sse2/avx2 etc) 
    """
    def build_extension(self, ext):
        print("IMAGED11:using compiler",self.compiler.compiler_type)
        if ext.name.find("sse2")>0:
            ext.extra_compile_args = src.bldlib.sse2arg
            ext.extra_link_args = src.bldlib.lsse2arg
        if ext.name.find("avx2")>0:
            ext.extra_compile_args = src.bldlib.avx2arg
            ext.extra_link_args = src.bldlib.lavx2arg
        print("IMAGED11:COMPILE ARGS:", ext.extra_compile_args)
        print("IMAGED11:LINK ARGS:", ext.extra_link_args)
        build_ext.build_extension(self, ext)


ekwds = { 'include_dirs' : [get_include(), 'src' ],
          'library_dirs' : ['./src'],
        }

# Compiled extensions:
extensions = [ Extension( "cImageD11_sse2", 
                          sources = ["src/cImageD11_sse2.pyf",],
                          libraries = [src.bldlib.sse2libname],
                          **ekwds ),
               Extension( "cImageD11_avx2", 
                          sources = ["src/cImageD11_avx2.pyf",],
                          libraries = [src.bldlib.avx2libname],
                          **ekwds) ]


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
      cmdclass = {'build_ext': custom_build_ext },
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


