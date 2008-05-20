import os, glob, shutil
fl = glob.glob(os.path.join("..","Pmw_dynamic_imports","Pmw_1_2","lib","*.py"))
for f in fl:
   print f
   t = os.path.split(f)[-1].replace("Pmw","")
   print t
   shutil.copy(f, t)

# File to allow this directory to be treated as a python package.
f = open("__init__.py","w")
f.write("# Fabian specific __init__ for Pmw
#
# Apologies to anyone wanting a real Pmw - this is broken
# The dynamic imports are evil and the make package script they
# provide doesn't work as I dont have regsub.
#
#
from Pmw.Base import *
from NoteBook import *
from Color import *
from ScrolledFrame import *
from Group import *
")
