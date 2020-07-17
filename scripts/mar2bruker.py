#!/usr/bin/env python

from __future__ import print_function

from Numeric import *

import sys,time
starttime=time.time()

# Convert mar to edf

def mar2edf_fit2d(infile,outfile,rebinfactor):
    macro="""INPUT DATA
MAR-PCK
%s
RE-BIN
%d
%d
1
1
EXCHANGE
OUTPUT
"KLORA"
%s
EXIT
YES
"""
    open("converter.mac","w").write(macro%(infile,rebinfactor,rebinfactor,outfile))
    import os
    os.system("fit2d -nogr -dim4096x4096 -macconverter.mac")
    
   
# Now read in a fit2d float edf file

def readedffloat(filename):
    f = open(filename,"rb")
    h = [ s.split("=") for s in f.read(1024).split(";") ]
    for n,v in h[:-1]:
        v=v.strip()
        if n.find("Dim_1")>-1: xdim = int(v)
        if n.find("Dim_2")>-1: ydim = int(v)
        if n.find("DataType")>-1:
            if v != "FLOAT":
                raise Exception("File format problem for "+filename+v)
    # seek back from end of file (independent of screwed up headers)
    f.seek(-xdim*ydim*4,2)
    return  reshape(fromstring(f.read(xdim*ydim*4),Float32),(xdim,ydim))


def fortranindex(i,size):
    index1 = i/size # rounds down
    index2 = i%size
    return index2*size+index1

def getomega(num):
    return 0.

# Get a bruker header

# from ImageD11 ...
import opendata

try:
   brukerfiletemplate=sys.argv[1]
   stem = sys.argv[2]
   first=int(sys.argv[3])
   last=int(sys.argv[4])
except:
    print("Usage: %s brukertemplate stem first last [extension=mar3450] [rebinfactor=3] [outsize=2048]"%(sys.argv[0]))
    sys.exit()

try:
   inputextn=sys.argv[5]
except:
   inputextn="mar3450"
try:
   rebinfactor=int(sys.argv[6])
except:
   rebinfactor=3
try:
   outsize=int(sys.argv[7])
except:
   outsize=2048
    
hs=opendata.readbrukerheader(brukerfiletemplate)["headerstring"]

#
hs=hs.replace("RANGE  :     0.4000000","RANGE  :     0.1000000")
hs=hs.replace("INCREME:    -0.4000000","INCREME:     0.1000000")
a=hs.find("ANGLES")
e=hs.find("ENDING")
s=hs.find("START")
x=hs.find("CENTER")
o=hs.find("NOVERFL")
r=hs.find("NROWS")
c=hs.find("NCOLS")

extn=".edf"

try:
    outsize = int(sys.argv[7])
    if outsize not in [512,1024,2048]:
        print("You should probably be aiming to make 512,1024 or 2048 bruker frames? arg 7")
except:
    outsize = 2048
omegadone =[]
skip=0
for num in range(first,last+1):
    # make edf file
    infile = stem+"%03d."%(num)+inputextn
    outfile = stem+"%03d.edf"%(num)
    mar2edf_fit2d(infile,outfile,rebinfactor)
    d = readedffloat(outfile)
    
    if d.shape != (outsize,outsize) and d.shape[0]>outsize:
        # take middle
        off = (d.shape[0]-outsize)/2
        d=d[off:off+outsize, off:off+outsize]

    if d.shape != (outsize,outsize) and d.shape[0]<outsize:
        # centre it
        off = (outsize - d.shape[0])/2
        new=zeros((outsize,outsize),Float32)
        new[off:off+d.shape[0] , off:off+d.shape[1] ] += d
        d=new

    #ANGLES :   -19.9995995  -180.0000000     0.0000000    54.9219017
    #         1234567890123012345678901230123456789012301234567890123


    omega = getomega(num)

    
    hl=hs.replace(hs[a:a+80],
                  "ANGLES :%14f%14f%14f%14f"%(0,0,omega,90)+" "*(80-14*4-8)).replace(
        hs[e:e+80],"ENDING :%14f%14f%14f%14f"%(0,0,omega+0.5,90)+" "*(80-14*4-8)).replace(
        hs[s:s+80],"START  :%14f"%(omega)+" "*(80-14-8)).replace(
        hs[x:x+80],"CENTER :%14f%14f"%(468,567)+" "*(80-14*2-8)).replace(
        hs[r:r+80],"NROWS  :%10d"%(outsize)+" "*(80-10-8)).replace(
        hs[c:c+80],"NCOLS  :%10d"%(outsize)+" "*(80-10-8))
 
    # Handle overflows
    r = ravel(d.copy())
    indices = compress(r > pow(2,16)-1, arange(r.shape[0]))
    noverfl=indices.shape[0]
    #                           123456789
    hl = hl.replace(hs[o:o+80],"NOVERFL:%10d"%(noverfl)+" "*(80-10-8) )
    #print hl[a-80:a+160]
    #print hl[e-80:e+160]
    #print "ANGLES :   -19.9995995  -180.0000000     0.0000000    54.9219017"
    r = where(r<65535,r,65535)
    out=open(opendata.makename(stem+".",num-skip,""),"wb")
    out.write(hl)
    out.write(r.astype(UInt16).tostring())
    # write overflows
    length=0
    for i in indices:
        s="%9d%7d"%(ravel(d)[i],i) 
        length+=16
        out.write(s)
    # pad to 512 length, seems optional
    out.write(" "*(512-length%512))
    out.close()
    print("Image",num,"time",(time.time()-starttime)/(num-first+1),opendata.makename(stem+".",num,""),"overfl",noverfl)
