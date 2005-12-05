ImageD11
Version 0.4
Jon Wright
wright@esrf.fr

This is the source distribution.... I assume you know what you are doing!

Try this:
Download and install python from www.python.org
Download and install the Numeric extensions from numeric.scipy.org (version >23, be sure to get Numeric, not numarray)
Download and install the matplotlib plotting package from matplotlib.sourceforge.net
Download and install the opengl extensions pyopengl.sourceforge.net

Finally, you must already have gotten ImageD11 from somewhere to be reading this, so now....

You probably want to try something like:
python setup.py build
or
python setup.py build --compiler=mingw32

Followed by:
python setup.py install



After it is installed, you should find a script ImageD11_gui.py, somewhere



Good luck!





(Thanks to Tine Knudsen for updating the links in this file!)