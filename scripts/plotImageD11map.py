#!/usr/bin/env python

from __future__ import print_function

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
    return abs(float(items[j+2]))
            
try:
    outersf = float(sys.argv[3])
except:
    outersf = 1.0

print("Scale factor is",outersf)
for g in gf:
    #print g.translation, g.ubi
    mapfile.write("\n\n")
    o = g.translation
    try:
        sf = pow(getmedian(g.intensity_info),0.3333)*outersf
    except:
        sf = outersf
    try:
        k = int(g.npks)
    except:
        k = 1
    for u in g.ubi:
        dodot(o,k)
        dodot(o+u*sf,int(g.npks))
    for u in g.ubi:
        dodot(o,k)
        dodot(o-u*sf,int(g.npks))
#    dodot(o,k)
#    dodot(o+sf*(-g.ubi[0]-g.ubi[1]),k)
#    dodot(o,k)
#    dodot(o+sf*(g.ubi[0]+g.ubi[1]),k)

mapfile.close()
term = " "
if "linux" in sys.platform:
    term = "set term x11"
if "win32" in sys.platform:
    term = "set term windows"
    
open("gnuplot.in","w").write("""
%s
set ticslevel 0
set title "Color proportional to number of peaks"
set palette model RGB
set palette defined ( 0 "violet", 1 "blue", 2 "green", 3 "yellow", 4 "orange", 5 "red" )
set view equal xyz
set view 75,0,1,1
#set terminal gif animate delay 10 loop 1 optimize size 1024,768
set nokey
set hidden3d
#set output "ImageD11map.gif"
splot "%s" u 1:2:3:4 w l lw 2 lc pal z
"""%(term, sys.argv[2])
# "".join(["set view 75,%d\n replot\n"%(i) for i in range(1,360,1)])
                             )


    
os.system("gnuplot -background white gnuplot.in -")

    

