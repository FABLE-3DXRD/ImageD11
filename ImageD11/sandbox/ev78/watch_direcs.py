
from __future__ import print_function

import glob, os, sys, time

HOME = "/data/visitor/ev78/id11"

checked = {0:0}
skip=[]
skipdir="""01_plate8_3_2x2_alignKB_edna
01_plate8_3_align_edna
01_plate8_3_alignKB_edna
02_plate_8_3__edna
03_plate_8_3__edna
aligning_edna
ceo2calib_edna
ceo2calibKB_edna
ceo2calibKB__edna
shutterless_2x2_dark02s_edna
shutterless_dark02s_edna""".split()


def process(stem):
    print("Checking",stem, end=' ')
    sys.stdout.flush()
    os.chdir(os.path.join(HOME, stem))
    fl = glob.glob("*")
    edfs = [s for s in fl if s.endswith("edf")]
    dats = [s for s in fl if s.endswith("dat")]
    print(len(edfs))
    if len(edfs)!=len(dats):
        os.system("../process2x2.sh %s"%(stem))        

lastone=[0,]
while 1:
    os.chdir(HOME)
    for dirname in glob.glob("*edna"):
        if dirname in skipdir:
            continue
        num = dirname.split("_")[0]
        try:
            num = int(num)
        except:
            print("dirname does not start with a number",dirname)
            continue
        stem = dirname.replace("_edna","")
        process(stem)
    lastone = list(checked.keys())
    lastone.sort()
    print("sleeping")
    time.sleep(10)
    
