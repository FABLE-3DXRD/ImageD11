
import sys, os

args = {
    'flt' : os.path.join("Al1000","Al1000.flt"),

    'par' : os.path.join("Al1000","Al1000.par"),
    'ubi' : "shaken.map"
    }

# 1: Run makemap.py
if 0:
    cmd1 = "makemap.py -f %(flt)s -p %(par)s -s cubic --omega_no_float --no_sort "%args
    cmd2 = "-U allgrid_makemap.map -t 0.015 -u %(ubi)s"%args
    os.system(cmd1+cmd2)
    cmd2 = "-U allgrid_makemap.map -t 0.0075 -u allgrid_makemap.map"
    os.system(cmd1+cmd2)
    cmd2 = "-U allgrid_makemap.map -t 0.0075 -u allgrid_makemap.map"
    os.system("time -o makemap.log "+cmd1+cmd2)
    


args['flt']= os.path.join("ideal.flt")
if 0:
    # 2: Run sandbox/fittrans using assigned grains
    cmd = "python ../../sandbox/fittrans.py %(flt)s.new %(par)s %(ubi)s allgrid_fittrans.map "%args
    os.system("time -o fittrans.log "+ cmd)

    # 3: Run ImageD11/indexer using assigned grains
    cmd = "python ../../ImageD11/indexer.py %(par)s %(flt)s.new fit %(ubi)s allgrid_indexer.map"%args
    os.system("time -o indexer.log " + cmd)

# 4: Run sandbox/teo.py  using assigned grains
cmd = "python ../../sandbox/teo.py %(flt)s.new %(par)s %(ubi)s allgrid_teo.map"%(args)
os.system("time -o teo.log "+ cmd)
