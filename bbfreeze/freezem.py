

from bbfreeze import Freezer
import os, sys,ImageD11

# Choose a version number which climbs - should actually be the svn revision #

target = "fable_python_" + sys.platform +"_" + ImageD11.__version__

includes =  ["ctypes.util", 
             "Tkinter", 
             "OpenGL",
             "OpenGL.Tk",
             "Pmw",
             "matplotlib",
             "matplotlib.numerix",
             "matplotlib.numerix.fft",
             "matplotlib.numerix.linear_algebra",
             "matplotlib.numerix.ma",
             "matplotlib.numerix.mlab",
             "matplotlib.numerix.random_array",
             "pytz",
             "PIL",
             "numpy",

                ]

if sys.platform == "win32":
    # linux doesnt seem to have it.
    includes.append("pytz.zoneinfo.UTC")
    # linux doesnt seem to need it.
    includes.append("pyreadline")


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
        # This is from PolyXSim
        "PolyXSim.py",
        # 
        # This one I want, but recognise might not work...
        "ipython"
        ]


if sys.version_info[0:2] == (2,5):
    if sys.platform == "win32":
        root = """c:\python25\scripts"""
        f = Freezer(target , includes = tuple(includes) )        
        for s in scripts:
            f.addScript( os.path.join(root , s) )
        f()

        import shutil
        shutil.copy( os.path.join("win32","Togl20.dll"),     
                     os.path.join(".", target) )
        shutil.copy( os.path.join("win32","pkgIndex.tcl"), 
                     os.path.join(".", target) )
        import site
        shutil.copy( site.__file__,
                     os.path.join(".", target) )
        shutil.copy( site.__file__.replace("pyc","py"),
                     os.path.join(".", target) )
        sys.exit()
    if sys.platform == "linux2" and sys.version.find("GCC 4.2.3 (Ubuntu")>0:
        # eg - an ubuntu machine...
        root = "/usr/bin"
        target = target + "_glibc_2.4"

        f = Freezer(target , includes = tuple(includes) )
        for s in scripts:
            f.addScript( os.path.join( root, s) )
        f()
        # 
        sys.exit()



print "I dont know how to do your platform / version yet, sorry"

