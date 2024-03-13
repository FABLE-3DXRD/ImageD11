Filtering
=========

There is a new gui option, try doing:
 python -m ImageD11.silxGui.silx_colfile.py mypeaks.flt


How to clean up a peaksearch output to only get the nice peaks, or 
whatever criteria you like::

 #!/usr/bin/env fable.python
 from ImageD11.columnfile import *
 import sys
 f = columnfile(sys.argv[1])
 m = (f.fc < 1020.0) | (f.fc > 1027.0)
 f.filter(m)
 m = (f.fc < 2044.0) | (f.fc > 2051.0)
 f.filter(m)
 m = (f.sc < 1020.0) | (f.sc > 1027.0)
 f.filter(m)
 m = (f.sc < 2044.0) | (f.sc > 2051.0)
 f.filter(m)
 f.writefile(sys.argv[2])

Since version 1.0 of ImageD11 a "columnfile" facility has been added. This python object allows for easy filtering of your data by generating attributes representing the names of columns, and also a filter function. Here is an example::

 D:\wright\Grain_Stuff\sim_test\simdata_oPPA_5grains>python
 Python 2.4.4 (#71, Oct 18 2006, 08:34:43) [MSC v.1310 32 bit (Intel)] on win32
 Type "help", "copyright", "credits" or "license" for more information.
 >>> from ImageD11.columnfile import columnfile
 >>> **obj = columnfile("peaks_t11.flt")**
 >>> dir(obj)
 ['IMax_f', 'IMax_int', 'IMax_o', 'IMax_s', 'Max_f', 'Max_o', 'Max_s', 'Min_f', 'Min_o', 'Min_s', 'Number_of_pixels',
 '__doc__', '__init__', '__module__', 'avg_intensity', 'bigarray', 'covfo', 'covsf', 'covso', 'dety', 'detz', 'f_raw', 'fc',
 'filename', 'filter', 'ncols', 'nrows', 'omega', 'onfirst', 'onlast', 'parameters', 'readfile', 's_raw', 'sc',
 'set_attributes', 'sigo', 'sigs', 'spot3d_id', 'sum_intensity', 'sum_intensity^2', 'titles', 'writefile']
 >>> from matplotlib.pylab import *
 >>> plot(obj.Number_of_pixels, obj.sum_intensity)
 >>> show()
  [close plot window]
 >>>
 >>> obj.nrows
 13985
 >>> **obj.filter( obj.Number_of_pixels > 5 )**
 >>> obj.nrows
 13577
 >>> # Make a more complicated mask
 >>> mask = (obj.detz > 20) & (obj.detz < 130)
 >>> # see how many peaks will survive
 >>> mask.cumsum()
 array([  0,   0,   0, ..., 482, 482, 482])
 >>> **obj.filter( mask )**
 >>> obj.nrows
 482
 >>> **obj.writefile( "my_filtered_peaks.flt" )**
 
Some possible filtering operations clipped out of the numpy array object 
help::
 
 >>> help(mask)
 Help on ndarray object:
 class ndarray(__builtin__.object)
 |  An array object represents a multidimensional, homogeneous array
 |  of fixed-size items. [...]
 |
 |  Methods defined here:
 |
 |  __abs__(...)
 |      x.__abs__() <==> abs(x)
 |
 |  __add__(...)
 |      x.__add__(y) <==> x+y
 |
 |  __and__(...)
 |      x.__and__(y) <==> x&y
 |
 |  __contains__(...)
 |      x.__contains__(y) <==> y in x
 |
 |  __div__(...)
 |      x.__div__(y) <==> x/y
 |
 |  __divmod__(...)
 |      x.__divmod__(y) <==> divmod(x, y)
 |
 |  __eq__(...)
 |      x.__eq__(y) <==> x==y
 |
 |  __float__(...)
 |      x.__float__() <==> float(x)
 |
 |  __floordiv__(...)
 |      x.__floordiv__(y) <==> x//y
 |
 |  __ge__(...)
 |      x.__ge__(y) <==> x>=y
 |
 |  __gt__(...)
 |      x.__gt__(y) <==> x>y
 |
 |  __iadd__(...)
 |      x.__iadd__(y) <==> x+y
 |
 |  __iand__(...)
 |      x.__iand__(y) <==> x&y
 |
 |  __idiv__(...)
 |      x.__idiv__(y) <==> x/y
 |
 |  __ifloordiv__(...)
 |      x.__ifloordiv__(y) <==> x//y
 |
 |  __ilshift__(...)
 |      x.__ilshift__(y) <==> x<<y
 |
 |  __imod__(...)
 |      x.__imod__(y) <==> x%y
 |
 |  __imul__(...)
 |      x.__imul__(y) <==> x*y
 |
 |  __int__(...)
 |      x.__int__() <==> int(x)
 |
 |  __invert__(...)
 |      x.__invert__() <==> ~x
 |
 |  __ior__(...)
 |      x.__ior__(y) <==> x|y
 |
 |  __ipow__(...)
 |      x.__ipow__(y) <==> x**y
 |
 |  __irshift__(...)
 |      x.__irshift__(y) <==> x>>y
 |
 |  __isub__(...)
 |      x.__isub__(y) <==> x-y
 |
 |  __itruediv__(...)
 |      x.__itruediv__(y) <==> x/y
 |
 |  __ixor__(...)
 |      x.__ixor__(y) <==> x^y
 |
 |  __le__(...)
 |      x.__le__(y) <==> x<=y
 |
 |  __lshift__(...)
 |      x.__lshift__(y) <==> x<<y
 |
 |  __lt__(...)
 |      x.__lt__(y) <==> x<y
 |
 |  __mod__(...)
 |      x.__mod__(y) <==> x%y
 |
 |  __mul__(...)
 |      x.__mul__(y) <==> x*y
 |
 |  __ne__(...)
 |      x.__ne__(y) <==> x!=y
 |
 |  __neg__(...)
 |      x.__neg__() <==> -x
 |
 |  __nonzero__(...)
 |      x.__nonzero__() <==> x != 0
 |
 |  __or__(...)
 |      x.__or__(y) <==> x|y
 |
 |  __pow__(...)
 |      x.__pow__(y[, z]) <==> pow(x, y[, z])
 |
 |  __radd__(...)
 |      x.__radd__(y) <==> y+x
 |
 |  __rand__(...)
 |      x.__rand__(y) <==> y&x
 |
 |  __rdiv__(...)
 |      x.__rdiv__(y) <==> y/x
 |
 |  __rdivmod__(...)
 |      x.__rdivmod__(y) <==> divmod(y, x)
 |
 |  __rfloordiv__(...)
 |      x.__rfloordiv__(y) <==> y//x
 |
 |  __rlshift__(...)
 |      x.__rlshift__(y) <==> y<<x
 |
 |  __rmod__(...)
 |      x.__rmod__(y) <==> y%x
 |
 |  __rmul__(...)
 |      x.__rmul__(y) <==> y*x
 |
 |  __ror__(...)
 |      x.__ror__(y) <==> y|x
 |
 |  __rpow__(...)
 |      y.__rpow__(x[, z]) <==> pow(x, y[, z])
 |
 |  __rrshift__(...)
 |      x.__rrshift__(y) <==> y>>x
 |
 |  __rshift__(...)
 |      x.__rshift__(y) <==> x>>y
 |
 |  __rsub__(...)
 |      x.__rsub__(y) <==> y-x
 |
 |  __rtruediv__(...)
 |      x.__rtruediv__(y) <==> y/x
 |
 |  __rxor__(...)
 |      x.__rxor__(y) <==> y^x
 |
 |  __sub__(...)
 |      x.__sub__(y) <==> x-y
 |
 |  __truediv__(...)
 |      x.__truediv__(y) <==> x/y
 |
 |  __xor__(...)
 |      x.__xor__(y) <==> x^y
 |
 |  all(...)
 |      a.all(axis=None)
 |
 |  any(...)
 |      a.any(axis=None, out=None)
 |
 |  argmax(...)
 |      a.argmax(axis=None, out=None)
 |
 |  argmin(...)
 |      a.argmin(axis=None, out=None)
 |
 |  choose(...)
 |      a.choose(b0, b1, ..., bn, out=None, mode='raise')
 |
 |      Return an array that merges the b_i arrays together using 'a' as
 |      the index The b_i arrays and 'a' must all be broadcastable to the
 |      same shape.  The output at a particular position is the input
 |      array b_i at that position depending on the value of 'a' at that
 |      position.  Therefore, 'a' must be an integer array with entries
 |      from 0 to n+1.;
 |
 |  clip(...)
 |      a.clip(min=, max=, out=None)
 |
 |  nonzero(...)
 |      a.nonzero() returns a tuple of arrays
 |
 |      Returns a tuple of arrays, one for each dimension of a,
 |      containing the indices of the non-zero elements in that
 |      dimension.  The corresponding non-zero values can be obtained
 |      with
 |          a[a.nonzero()].
 |
 |      To group the indices by element, rather than dimension, use
 |          transpose(a.nonzero())
 |      instead. The result of this is always a 2d array, with a row for
 |      each non-zero element.;
 
We think it is Turing complete!
