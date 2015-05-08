# python f-setup.py build_ext --inplace
#   cython f.pyx -> f.cpp
#   g++ -c f.cpp -> f.o
#   g++ -c fc.cpp -> fc.o
#   link f.o fc.o -> f.so

# distutils uses the Makefile distutils.sysconfig.get_makefile_filename()
# for compiling and linking: a sea of options.

# http://docs.python.org/distutils/introduction.html
# http://docs.python.org/distutils/apiref.html  20 pages ...
# http://stackoverflow.com/questions/tagged/distutils+python

import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# from Cython.Build import cythonize

ext_modules = [Extension(
    name="cImageD11_wrap",
    sources=["cImageD11_wrap.pyx", "cImageD11.c"],
    include_dirs = [numpy.get_include()],  
    language="c",
    # libraries=
     extra_compile_args = "-fopenmp".split(),
     extra_link_args = "-fopenmp".split()
    )]

setup(
    name = 'cImageD11',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    # ext_modules = cythonize(ext_modules)  ? not in 0.14.1
    # version=
    # description=
    # author=
    # author_email=
    )

# test: import f
