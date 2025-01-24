
from __future__ import print_function

import os, time, glob, sys

firstrun = True

log=open("log","w")

while 1:
    edfs = glob.glob("*.edf")
    chis = glob.glob("*.chi")
    if firstrun:
        chis = []
        firstrun = False
    todo = []
    for e in edfs:
        if e.replace("edf","chi") not in chis:
            num = int(e.split("tion")[1].split(".")[0])
            todo.append( num )
    if len(edfs) != len(chis):
        os.system("fit2dcake.py -g TbAsO4_with_transition -p newfar.par")
        os.system("rm tmp*")
        first = min(todo)
        last = max(todo)
        log.write("Fitting %d %d\n"%(first,last))
        os.system("rm fits")
        os.system("rm results")
        os.system("fitem_jon_chi TbAsO4_with_transition %d %d 5.9 6.27"%(
            first, last))
        os.system("cat results >> pkfits")
    print()
    print("Control-c to exit me\r", end=' ')
    sys.stdout.flush()
    time.sleep(30)
