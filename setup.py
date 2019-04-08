
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
  python setup.py build
   -> it will build the libraries in the src folder
   -> generate sse2 and avx2 wrappers
   -> compile them here

"""

# For pip / bdist_wheel etc
import setuptools
import os, sys, sysconfig, platform, struct
# For using f2py
from numpy.distutils.core import setup, Extension
from numpy.distutils.command.build_ext import build_ext
from numpy import get_include

if "build" in sys.argv:
    os.chdir("src")
    os.system("python write_check.py")
    os.system("python bldlib.py")
    os.chdir("..")
    
import src.bldlib

# For 32 or 64 bits
nbyte = struct.calcsize("P") # 4 or 8
           
class custom_build_ext( build_ext ):    
    def build_extension(self, ext):
        print("IMAGED11:using compiler",self.compiler.compiler_type)
        if ext.name.find("sse2")>0:
            ext.extra_compiler_args = src.bldlib.sse2arg
        if ext.name.find("avx2")>0:
            ext.extra_compiler_args = src.bldlib.avx2arg    
        ext.extra_link_args = list(ext.extra_compiler_args)
        build_ext.build_extension(self, ext)

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


extn_kwds = {'sources' : ["src/cImageD11_sse2.pyf",],
             'libraries' : [src.bldlib.sse2libname],
             'include_dirs' : [get_include(), 'src' ],
             'library_dirs' : ['src'],
             }
avx2_kwds = {'sources' : ["src/cImageD11_avx2.pyf",] ,
             'libraries' : [src.bldlib.avx2libname],
             'include_dirs':[get_include(), 'src' ],
             'library_dirs' : ['src'],
             }

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


