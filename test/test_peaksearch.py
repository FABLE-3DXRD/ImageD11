

# testing for the peaksearch command line script

# These are not really unit tests

def w(name,image):
    import Numeric
    f=open(name,"wb")
    f.write("{%1021s}\n"%("""Omega = 0.0 ;\n Dim_1 = 256 ;\n Dim_2 = 256 ; 
    DataType = UnsignedShort ;
    ByteOrder = LowByteFirst ;
    
    
    
    \n\n"""))
    f.write(Numeric.ravel(image).astype(Numeric.UInt16).tostring())
    f.close()

def makedata():
    import Numeric
    image=Numeric.ones((256,256),Numeric.UInt16)*20
    image[20:53,56:99]=20000
    for i in range(4):
        w("data%04d.edf"%(i),image)
    w("dark_good.edf",image*0+1)
    w("dark_bad.edf",image*0+65534)

def losedata():
    for i in range(4):
        name = "data%04d.edf"%(i)
        print "removing",name
        os.remove(name)
    os.remove("dark_bad.edf")
    os.remove("dark_good.edf")

cmds = ["""..\scripts\peaksearch.py -n data -f 0 """]

cmds = [cmd + "-l %d "%(last) for last in [0,1,3] for cmd in cmds]

thresholds = ["0","5500"]

cmds = [cmd + "-t %s "%(t) for t in thresholds for cmd in cmds]

cmds = cmds + [ cmd + "-t 1500 -t 25000 " for cmd in cmds]

# cmds = cmds + [cmd + "-F edf " for cmd in cmds]

cmds = cmds + [cmd + "-d %s "%(d) for d in ["dark_good.edf","dark_bad.edf"] for cmd in cmds]

# floods as dark
cmds = cmds + [cmd + "-O %s "%(d) for d in ["dark_good.edf","dark_bad.edf"] for cmd in cmds]

# one with spline , one perfect

clist = [ cmd + "-p Y " for cmd in cmds] + [cmd + "-s spatial2k.spline" for cmd in cmds]

import os
makedata()
print len(clist)
for c in clist:
    print c
    os.system(c)

losedata()
