

import glob, shutil 


def got_future(lines):
    for line in lines:
        if line.find("from __future__ import print_function")==0:
            return True
    return False

def got_print(lines):
    for line in lines:
        if line.find("print") > -1:
            return True
    return False

def got_input(lines):
    for line in lines:
        if line.find("input") > -1:
            return True
    return False

fl = glob.glob("*.py")

def fix(fname, lines):
    shutil.move(fname, fname+".bak2")
    g = open(fname,"w")
    g.write("\nfrom __future__ import print_function\n")
    if got_input(lines):
        g.write("from six.moves import input\n")
    for line in lines:
        g.write(line)
    g.close()


for fname in fl:
    lines = open(fname).readlines()
    p = got_print(lines) 
    f = got_future(lines)
    if p and f:
        print "OK",fname
        continue
    if p:
        print "Need to fix",fname
        fix(fname, lines)


