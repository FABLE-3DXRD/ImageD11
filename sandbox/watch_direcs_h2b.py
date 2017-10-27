
from __future__ import print_function


import glob, os, sys, time, fabio

HOME = "/data/visitor/hc1181/id11"

FLOOD= "/data/visitor/hc1181/id11/F4M_irt_flood.edf"

CMD = "huber2bruker.py -n %s -f %d -l %d -F %s -d %s -w 0.158154 -D 8 --omega_zero=-1.5935 --chi_zero=-0.3705 -j 4 --maskfilename=/data/visitor/hc1181/id11/fit2d.msk"

darks = {
    # Key is time in 0.1 s units
    2 : "/data/visitor/hc1181/id11/dark02s.edf",
    5 : "/data/visitor/hc1181/id11/dark05s.edf",
    10 : "/data/visitor/hc1181/id11/dark1s.edf",
    20 : "/data/visitor/hc1181/id11/dark1s.edf",
    30 : "/data/visitor/hc1181/id11/dark3s.edf",
    50 : "/data/visitor/hc1181/id11/dark5s.edf"
    }

def convert(stem,first,last,int_time):
    nearest_int = int(10*float(int_time))
    assert nearest_int in darks
    dark = darks[nearest_int]
    mycmd = CMD%(stem, first, last, FLOOD, dark)
    print(mycmd)
    os.system("%s"%(mycmd))


def bname(s):
    o = fabio.deconstruct_filename(s)
    return "%s_0.%04d"%(o.stem,o.num)

def process(stem):
    print("Checking",stem)
    sys.stdout.flush()
    os.chdir(os.path.join(HOME, stem))
    fl = glob.glob("%s????.edf"%(stem))
    br = glob.glob("%s_0.????"%(stem))
    # edfs = filter(lambda s:s.endswith("edf"), fl)
    # print len(edfs)
    todo = []
    for f in fl:
        newname = bname(f)
        if newname not in br:
            todo.append((newname,f))
    if len(todo) > 0:
        nums = [fabio.getnum(t[0]) for t in todo]
        nums.sort()
    else:
        return
    first=nums[0]
    last = nums[-1]
    int_time = fabio.openimage.openheader(todo[0][1]).header['Integration']
    convert(stem, first, last, int_time )


checked = {0:0}
skip = [30]
skipdir="grid_edna alignKB_edna check_det_edna dark05s_edna 01_SLSCO_006_C1_9__edna".split()

while 1:
    os.chdir(HOME)
    for dirname in glob.glob("*edna"):
        if dirname in skipdir:
            continue
        stem = dirname.replace("_edna","")
        process(stem)
    print("sleeping")
    time.sleep(10)
#    break
