
# ImageD11_v1.0 Software for beamline ID11
# Copyright (C) 2005-2007  Jon Wright
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




# from distutils.core import setup, Extension
# from setuptools import setup, Extension
import setuptools
from numpy.distutils.core import setup, Extension

from numpy import get_include

nid = [get_include()]


# Compiled extensions:

# closest is for indexing grains
cl = Extension("closest",
               sources=['src/closest.c'],
               include_dirs = nid )

# connectedpixels is for peaksearching images



cp = Extension("connectedpixels",
               sources = ['src/connectedpixels.c','src/blobs.c'],
               include_dirs = nid)
# No header files for distutils as sources 'src/dset.h'])


# histogramming thing
ch = Extension("_hist",
               sources = ['src/hist.c'],
               include_dirs = nid )


# _splines is for correcting peak positions for spatial distortion
bl = Extension("_splines",
               sources = ['src/splines.c', 'src/bispev.c'],
               include_dirs = nid)

import sys
                  
# New fortran code - you might regret this...
fi = Extension("fImageD11",
               sources = ['fsrc/fImageD11.f90' ],
#               extra_f90_compile_args=["-fopenmp -O2"],
               libraries = ['gomp','pthread'])
# OK, I am beginning to regret it now. Patch for older numpys
sys.argv.extend ( ['config_fc', '--fcompiler=gnu95',
                   '--f90flags="-fopenmp -O2 -shared"'])

if sys.platform == 'win32':
    needed = [
        'xfab>=0.0.2',
        'fabio>=0.0.5',
        'numpy>=1.0.0',
        'matplotlib>=0.90.0',
        ]
else: # Take care of yourself if you are on linux
    # Your package manager is inevitably f*cked
    needed = []
#        'xfab>=0.0.1',
#        'fabio>=0.0.4']

# See the distutils docs...
setup(name='ImageD11',
      version='1.7.0',
      author='Jon Wright',
      author_email='wright@esrf.fr',
      description='ImageD11',
      license = "GPL",
      ext_package = "ImageD11",   # Puts extensions in the ImageD11 directory
      ext_modules = [cl,cp,bl,ch,fi],
      install_requires = needed,
      packages = ["ImageD11"],
      package_dir = {"ImageD11":"ImageD11"},
      url = "http://fable.wiki.sourceforge.net/ImageD11",
#      download_url = ["http://sourceforge.net/project/showfiles.php?group_id=82044&package_id=147869"],
      package_data = {"ImageD11" : ["doc/*.html"]},
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

print "For windows you would need:"
print 'set LDFLAGS="-static-libgfortran -static-libgcc -static -lgomp -shared"'
print 'also gfortran/gcc installed (--compiler=mingw32)'
print 'also to patch f2py to let it run'

