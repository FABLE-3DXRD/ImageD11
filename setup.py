
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
import sys
from numpy.distutils.core import setup, Extension
from numpy import get_include
import struct

if sys.platform == "win32" and "--compiler=mingw32" not in sys.argv:
    ecomparg = ["/openmp","-DF2PY_REPORT_ON_ARRAY_COPY", "/arch:sse2"]
    elinkarg = []
    elibs = None
else:
    ecomparg = ["-fopenmp","-O2", "-msse2", "-std=c99",
                "-flto",
                "-Wall", "-Wextra",
                "-DF2PY_REPORT_ON_ARRAY_COPY"]
    elinkarg = [a for a in ecomparg]
    elibs = ["gomp","pthread"]


nid = [get_include(),]

def fix_f2py_pointer():
    # decide if we have 32 bits or 64 bits compilation
    nbyte = struct.calcsize("P") # 4 or 8
    with open("src/cImageD11_template.pyf","r") as f:
        pyf = f.read()
    print(pyf.find( "(kind=size_t)"))
    pyf = pyf.replace(  "(kind=size_t)" ,"*%d"%(nbyte))
    with open("src/cImageD11.pyf", "w") as f:
        f.write( pyf )

fix_f2py_pointer()
# Compiled extension:
cImageD11extension = Extension( "cImageD11",
                                sources = [ "src/cImageD11.pyf",
                                            "src/connectedpixels.c",
                                            "src/closest.c",
                                            "src/cdiffraction.c",
                                            "src/localmaxlabel.c",
                                            "src/blobs.c"],
                               include_dirs = nid + ["src",],
                               extra_compile_args=ecomparg,
                               extra_link_args=elinkarg,
                               libraries = elibs
                               )
            

# Removed list of dependencies from setup file
# Do a miniconda (or something) instead...
#if sys.platform == 'win32':
#    needed = [
#        'six',
#        'numpy>=1.0.0',
#        'scipy', 
#        'xfab>=0.0.2',
#           'pycifrw'
#        'fabio>=0.0.5',
#        'matplotlib>=0.90.0',
#        ... 
#        ]

needed =[]#
# ["xfab",
#          "fabio",
#          "pyopengl",
#          "matplotlib",
#          "numpy",
#          "scipy",
#          "six",
#          "h5py",
#          ]

# See the distutils docs...
setup(name='ImageD11',
      version='1.9.0',
      author='Jon Wright',
      author_email='wright@esrf.fr',
      description='ImageD11',
      license = "GPL",
      ext_package = "ImageD11",   # Puts extensions in the ImageD11 directory
      ext_modules = [cImageD11extension,],
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


