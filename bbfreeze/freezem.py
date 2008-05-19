

from bbfreeze import Freezer
import os

import ImageD11

# Choose a version number which climbs - should actually be the svn revision #

target = "fable_python_exe" + ImageD11.__version__

f = Freezer(target ,
         includes = 
                ("ctypes.util", 
                 "Tkinter", 
                 "OpenGL",
                 "OpenGL.Tk",
                 "Pmw",
                 "matplotlib",
                 "pytz",
                 "PIL",
                 "pyreadline",
                ))

scripts =    [
        "peaksearch.py",
        "fitgrain.py",
        "ubi2cellpars.py",
        "filtergrain.py",
        "filterout.py",
        "pars_2_sweeper.py",
        "ImageD11_2_shelx.py",
        "fit2dcake.py",
        "edfheader.py",
        "recoveromega.py",
        "id11_summarize.py",
        "bgmaker.py",
        "makemap.py",
        "plotedf.py",
        "plotgrainhist.py",
        "rubber.py",
        "edf2bruker.py",
        "index_unknown.py",
        "ImageD11Server.py",
        "powderimagetopeaks.py",
        "ImageD11_gui.py",
        "plot3d.py",
        # these are from fabian
        "collapse.py",
        "median.py",
        "fabian.py",
        # these are from fabric
        "plothst.py",
        "integrate.py",
        "fabric.py",
        # 
        # This one I just like
        "ipython.py"
        
        ]

root = """c:\python25\scripts"""

for s in scripts:
    f.addScript( os.path.join(root , s) )

f()


import os
