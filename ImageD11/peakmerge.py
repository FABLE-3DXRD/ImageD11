



# ImageD11_v0.4 Software for beamline ID11
# Copyright (C) 2005  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA        

import sys


TOLERANCE = 1.0

class peak:
    def __init__(self,xc,yc,omega,avg,np):
        self.xc=xc
        self.yc=yc
        self.avg=avg
        self.np=np
        self.omega=omega

        self.sorter=omega
        self.forget=False
        

    def combine(self,other):
        np = self.np+other.np
        s = (self.avg*self.np + other.avg*other.np)
        avg = s/np
        omega = (self.omega*self.np*self.avg + other.omega*other.np*other.avg)/s
        xc = (self.xc*self.np*self.avg + other.xc*other.np*other.avg)/s
        yc = (self.yc*self.np*self.avg + other.yc*other.np*other.avg)/s
        newpeak = (xc,yc,omega,avg,np)
        #        print "Returning:",newpeak
        self.forget=True
        other.forget=True
        return peak(xc,yc,omega,avg,np)
        
    def __eq__(self,other):
        try:
           if abs(self.xc - other.xc) < TOLERANCE and abs(self.yc - other.yc) < TOLERANCE :
              return True
           else:
              return False
        except:
            print self,other
            raise

    def __str__(self):
        return "xc= %f yc= %f Omega= %f Avg= %f np= %d"%(self.xc,self.yc,self.omega,self.avg,self.np)
    def __repr__(self):
        return "xc= %f yc= %f Omega= %f Avg= %f np= %d"%(self.xc,self.yc,self.omega,self.avg,self.np)
    def __cmp__(self,other):
        o=cmp(self.omega,other.omega)
        if o == 0:
            x=cmp(self.xc,other.xc)
            if x ==0 :
                return cmp(self.yc,other.yc)
            else:
                return x
        else:
            return o
        

peakdict = {}

try: 
   infile = sys.argv[1]
   outfile = sys.argv[2]
except:
   print "%s infile outfile minimum_threshold"%(sys.argv[0])
try:
   th=float(sys.argv[3])
except:
   th=0.
   print "Using zero minimum threshold"

Omega = -9999.

lines = open(infile,"r").xreadlines()
i=0
t=0
print "# Opened file %s and read it"%(infile)
while i < len(lines):
    if lines[i].find(" Omega")>0:
        Omega = lines[i].split()[-1]
    if lines[i].find("Threshold level")>0:
#        print lines[i],lines[i].split()
        t = float(lines[i].split()[-1])
    if th > t:
       i=i+1
       continue
    if lines[i][0] != "#" and len(lines[i])>10:
        vals = [ float(x) for x in lines[i].split() ]
        p=peak(vals[4],vals[5],float(Omega),vals[1],vals[0])
        try:
            found=0
            for index in range(len(peakdict[Omega])):
                if p == peakdict[Omega][index]:
                    peakdict[Omega][index]=p.combine(peakdict[Omega][index])
                    found=1
            if found==0:
                peakdict[Omega].append(p)
        except KeyError:
            peakdict[Omega]=[p]
    i=i+1 # while loop

print "# Made a dictionary of peaks"
keys = peakdict.keys()
keys.sort(lambda x,y : cmp(float(x),float(y)))
print "# Sorted peaks "
peaklist=[]


#print peakdict["81"]

i=0
while i < len(keys)-1:
    # Check if this peak is on the next frame
    for k in range(len(peakdict[keys[i]])):  # peaks on this frame
        for j in range(len(peakdict[keys[i+1]])):
            if peakdict[keys[i]][k]==peakdict[keys[i+1]][j]:   # peaks on next frame
                peakdict[keys[i+1]][j]=peakdict[keys[i]][k].combine(peakdict[keys[i+1]][j])
                peakdict[keys[i]][k].forget=True
    i=i+1     

print "# Merged peak on adjacent frames"

peaklist=[]
i=0
while i < len(keys)-1:
    # Check if this peak is on the next frame
    for mypeak in peakdict[keys[i]]:  # peaks on this frame
        if not mypeak.forget:
            peaklist.append(mypeak)
#            print mypeak
    i=i+1         

print "# Made final peak list"
print "# Number of merged peaks",len(peaklist)
peaklist.sort()
outfile=open(sys.argv[2],"w")
for p in peaklist:
    if p.avg>1000:
#        print p
        outfile.write("%s\n"%(p))

print "# All done"
