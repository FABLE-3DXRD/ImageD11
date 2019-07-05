
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
from numpy import get_include
# Get the path for the static libraries        
import src.bldlib



if "build" in sys.argv:
    """ Ugly - requires setup.py build 
    Should learn to use build_clib at some point
    """
    os.chdir("src")
    print("Call write_check")
    os.system(sys.executable+" write_check.py")
    print("Call bldlib")
    # transmit compiler from command line
    ok = os.system(sys.executable+" bldlib.py "+ " ".join(sys.argv))
    os.chdir("..")
    if ok != 0:
        print("Return was",ok)
        sys.exit(ok)
    else:
        print("Seems to build OK")


ekwds = { 'include_dirs' : [get_include(), 'src' ],
          'library_dirs' : ['./src'],
          'extra_compile_args' : src.bldlib.arg + \
               ["-DF2PY_REPORT_ON_ARRAY_COPY=1.",],
          'extra_link_args' : src.bldlib.arg + \
               ["-DF2PY_REPORT_ON_ARRAY_COPY=1.",],
        }

# Compiled extensions:
extensions = [ Extension( "cImageD11_sse2", 
                          sources = ["src/cImageD11_sse2.pyf",],
                          libraries = [src.bldlib.sse2libname],
                          **ekwds ),
               Extension( "cImageD11_avx",
                          sources = ["src/cImageD11_avx.pyf",],
                          libraries = [src.bldlib.avx2libname],
                          **ekwds) ]


# Removed list of dependencies from setup file
#  borked something, not sure what or why
needed =[]
# ["xfab", "fabio", "pyopengl",  "matplotlib", "numpy", "scipy", "six", "h5py",
#  "pyopengltk", "FitAllB", ... ]

# See the distutils docs...
setup(name='ImageD11',
      version='1.9.2',
      author='Jon Wright',
      author_email='wright@esrf.fr',
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


