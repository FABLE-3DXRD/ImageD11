
import os , time

# generate the inp file
def measure(cmd):
    start = time.time()
    ret = os.system(cmd)
    end = time.time()
    open("timing.log","a").write("%f %s\n"%(end-start, cmd))

os.system("python make_u_t.py")
os.system("PolyXSim.py -i Al1000.inp")
os.system("python grid_index.py  Al1000/Al1000.flt Al1000/Al1000.par grid")
os.system("makemap.py -u allgrid.map -U allgridfit.map -f Al1000/Al1000.flt -p Al1000/Al1000.par --omega_slop=0.125 -t 0.0075")
os.system("makemap.py -u allgridfit.map -U allgridfit.map -f Al1000/Al1000.flt -p Al1000/Al1000.par --omega_slop=0.125 -t 0.0075")

