

import os , time

# generate the inp file
def measure(cmd):
    start = time.time()
    ret = os.system(cmd)
    end = time.time()
    if os.path.exists("timing.txt"):
        open("timing.log","a").write("%f %s\n"%(end-start, cmd))
    else:
        open("timing.log","w+").write("%f %s\n"%(end-start, cmd))

if os.path.exists("allgrid.map"):
    os.remove("allgrid.map")
if os.path.exists("timing.log"):
    os.remove("timing.log")    
measure("python make_u_t.py")
measure("PolyXSim.py -i Al1000.inp")
measure("python grid_index.py  Al1000/Al1000.flt Al1000/Al1000.par grid")
measure("makemap.py -u allgrid.map -U allgridfit.map -f Al1000/Al1000.flt -p Al1000/Al1000.par --omega_slop=0.125 -t 0.0075")
measure("makemap.py -u allgridfit.map -U allgridfit.map -f Al1000/Al1000.flt -p Al1000/Al1000.par --omega_slop=0.125 -t 0.0075")

