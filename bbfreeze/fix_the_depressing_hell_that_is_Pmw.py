import os, glob, shutil
fl = glob.glob(os.path.join("..","Pmw_dynamic_imports","Pmw_1_2","lib","*.py"))
for f in fl:
   print f
   t = os.path.split(f)[-1].replace("Pmw","")
   print t
   shutil.copy(f, t)

