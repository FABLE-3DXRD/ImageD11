
from __future__ import print_function, division

import os , time, sys

# generate the inp file
def measure(cmd):
    start = time.time()
    ret = os.system(cmd)
    end = time.time()
    if os.path.exists("timing.log"):
        open("timing.log","a").write("%f %s\n"%(end-start, cmd))
    else:
        open("timing.log","w+").write("%f %s\n"%(end-start, cmd))
    if ret!=0:
        print("Problem, stopping")
        sys.exit()

if os.path.exists("allgrid.map"):
    os.remove("allgrid.map")
if os.path.exists("timing.log"):
    os.remove("timing.log")    
measure(sys.executable + " make_u_t.py")
measure("PolyXSim.py -i Al1000.inp")
# measure("python grid_index.py  Al1000/Al1000.flt Al1000/Al1000.par grid 9")
measure(sys.executable + " make_u_t.py")
measure("makemap.py -u allgrid.map -U allgridfit.map -f Al1000/Al1000.flt -p Al1000/Al1000.par --omega_no_float -t 0.05")
measure("makemap.py -u allgridfit.map -U allgridfit.map -f Al1000/Al1000.flt -p Al1000/Al1000.par --omega_no_float -t 0.02")
for i in range(3):
    measure("makemap.py -u allgridfit.map -U allgridfit.map -f Al1000/Al1000.flt -p Al1000/Al1000.par --omega_slop=0.125 -t 0.0075")

os.system(sys.executable + " res2map.py Al1000/Al1000.ubi Al1000/Al1000.res ideal.map")
os.system(sys.executable + " shake.py ideal.map shaken.map")
