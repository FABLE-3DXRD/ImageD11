

from ImageD11.grain import read_grain_file
import sys, os

gf = read_grain_file(sys.argv[1])
mapfile=open(sys.argv[2],"w")

def dodot(xyz,k):
    mapfile.write("%f %f %f %d\n"%(xyz[0],xyz[1],xyz[2],k))

def getmedian(s):
    items=s.split()
    j = -1
    for i in range(len(items)):
        if items[i] == "median":
            j = i
    if j == -1:
        return 0
    return float(items[j+2])
            

sf = 1.
k = 0
for g in gf:
    k+=1
    #print g.translation, g.ubi
    mapfile.write("\n\n")
    o = g.translation
    sf = pow(getmedian(g.intensity_info),0.3333)
    for u in g.ubi:
        dodot(o,k)
        dodot(o+u*sf,k)
    for u in g.ubi:
        dodot(o,k)
        dodot(o-u*sf,k)
    dodot(o,k)
    dodot(o+sf*(-g.ubi[0]-g.ubi[1]),k)
    dodot(o,k)
    dodot(o+sf*(g.ubi[0]+g.ubi[1]),k)

mapfile.close()
open("gnuplot.in","w").write("""
set ticslevel 0
set palette rgbformulae 3,2,2 model HSV
splot "%s" u 1:2:3:4 w l lw 2 lc pal z
"""%(sys.argv[2]))
os.system("gnuplot -background white gnuplot.in -")

    

