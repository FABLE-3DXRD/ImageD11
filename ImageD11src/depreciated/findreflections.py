
"""
Goes looking for peaks...
"""


import sys

def loadfiltered(filename):
    """
    Read in file containing filtered and merged peak positions
    """
    file_obj=open(filename,"r")
    #                   0123456789012
    thisline = file_obj.readline()
    if thisline[0:12] != "# xc yc omega"[0:12]:
        print thisline
        raise Exception("Sorry That does not seem to be a filter peaks file,"+\
                        " output from the peaksearching menu option")
    # titles = line.replace("#","").split()
    bigarray=[]
    for thisline in file_obj.readlines():
        values = [float(z) for z in thisline.split()]
        bigarray.append(values)
    file_obj.close()
    return bigarray

import math

def findclosest(xin, yin, omega_in ,data):
    """ Finds the closest values in data to x,y,omega"""
    dmin = 20.
    vmin = None
    for values in data:
        xp,yp,op = values[0:3]
        dist = math.sqrt((xin-xp)*(xin-xp)+(yin-yp)*(yin-yp))
        if dist < dmin and abs(omega_in - op) < 2.:
            vmin = values
            dmin = dist
    return vmin


if __name__=="__main__":
    scans=[]
    loads={}

    for f in open(sys.argv[1]).readlines():
        name,load = f.split()
        scans.append(name)
        loads[name]=float(load)


    scandata = {}

    for scan in scans:
        scandata[scan] = loadfiltered(scan)

    results = open(sys.argv[3],"w")

    for line in open(sys.argv[2],"r"):
        items = line.split()
        h,k,l = [float(v) for v in items[2:5]]
        results.write("# hkl %f %f %f\n"%(h,k,l))
        x,y,omega =  [float(v) for v in items[7:10]]
        print h,k,l,x,y,omega
        for scan in scans:
            line = findclosest(x,y,omega,scandata[scan])
            results.write("%f %s "%(loads[scan],scan))
            if line is None:
                results.write(" Nothing found within tolerance ")
            else:
                for item in line:
                    results.write("%6g "%(item))
            results.write("\n")
        results.write("\n\n\n")
