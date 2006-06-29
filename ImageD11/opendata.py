


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

"""
ImageD11 generic file opener

Returns a data object

Try to determine file type from extension or otherwise first few
bytes???
"""

from data import data
import gzip, bz2
import Numeric 

def opendata(filename):
    """
    Guesses the file type and return a data object

    Supports edf, edf.gz, edf.bz2, xye, epf, inp, chi
     ... otherwise tries bruker
    """
    f = filename.rstrip()
    if f[-3:]in ["edf", "EDF"] :
        return openedf(filename)
    if f[-7:] in [".edf.gz"]:
        return openedf(filename)
    if f[-8:] in [".edf.bz2"]:
        return openedf(filename)
    if f[-3:] in ["xye","epf","inp"]:
        return openxye(filename)
    if f[-3:] in ["chi"]:
        return openchi(filename)
    if f[-5]=="." and f[-4:].isdigit():
        return openbruker(filename)
    s = "\n\nUnrecognised filename: %s\n"%(f)
    raise Exception(s)



def makename(stem,num,extn,width=4):
    """
    Filename creation mydata1234.edf etc
    eg stem + num + extn   ... with num being width wide,
    """
    if width>0:
        fs="%s%0"+"%d"%(width)+"d%s"
    else:
        fs="%s%d%s"
    return fs%(stem,num,extn)


def openchi(filename):
    """
    opens a chi file and returns a data object
    2 columns in resulting array (x,y)
    """
    h={}
    xy=[]
    h['columns']=2
    for line in open(filename,"r"):
        try:
            v=[float(x) for x in line.split()]
            if len(v)==2:
                xy.append(v)
        except:
            pass # titles

    a=Numeric.array(xy,Numeric.Float)
    h['rows']=a.shape[1]
    return data(a,h)

def openxye(filename):
    """
    opens an xye file and returns a data object
    3 columns in resulting array (x,y,e)
    """
    h={}
    xye=[]
    h['columns']=3
    for line in open(filename,"r"):
        try:
            xye.append(map(float,line.split()))
        except:
            pass
    a=Numeric.array(xye)
    h['rows']=a.shape[1]
    h['xye']=True
    return data(a,h)

def readbytestream(file,offset,x,y,nbytespp,datatype='int',signed='n',
                   swap='n',typeout=Numeric.UInt16):
    """
    Reads in a bytestream from a file (which may be a string indicating
    a filename, or an already opened file (should be "rb"))
    offset is the position (in bytes) where the pixel data start
    nbytespp = number of bytes per pixel
    type can be int or float (4 bytes pp) or double (8 bytes pp)
    signed: normally signed data 'y', but 'n' to try to get back the right numbers
      when unsigned data are converted to signed (python has no unsigned numeric types.)
    swap, normally do not bother, but 'y' to swap bytes
    typeout is the Numeric type to output, normally UInt16, but more if overflows occurred
    x and y are the pixel dimensions

    TODO : Read in regions of interest

    PLEASE LEAVE THE STRANGE INTERFACE ALONE - IT IS USEFUL FOR THE BRUKER FORMAT
    """
    tin="dunno"
    len=nbytespp*x*y # bytes per pixel times number of pixels
    if datatype=='int' and signed=='n':
        if nbytespp==1 : tin=Numeric.UInt8
        if nbytespp==2 : tin=Numeric.UInt16
        if nbytespp==4 : tin=Numeric.UInt32
    if datatype=='int' and signed=='y':
        if nbytespp==1 : tin=Numeric.Int8
        if nbytespp==2 : tin=Numeric.Int16
        if nbytespp==4 : tin=Numeric.Int32
    if datatype=='float':
        tin=Numeric.Float32
    if datatype=='double' :
        tin=Numeric.Float64
    if tin=="dunno" :
        raise SyntaxError, "Did not understand what type to try to read"
    opened=0
    if type(tin) == type(file):  # Did we get a string or a file pointer?
        f=open(file,'rb')
        opened=1
    else:
        f=file
    f.seek(offset)
    if swap=='y':
        ar=Numeric.array(Numeric.reshape(
           Numeric.byteswapped(Numeric.fromstring(f.read(len),tin)),(x,y)),typeout)
    else:
        ar=Numeric.array(Numeric.reshape(
                               Numeric.fromstring(f.read(len),tin) ,(x,y)),typeout)
    if(opened):f.close()
    return(ar)


def readbrukerheader(file):
    """
    Reads a Bruker file header into a Python dictionary
    file=filename or file pointer
    """
    s="string"                 # s is a string
    if type(s) == type(file):  # if arg is a string, open file, else treat as file object
        f=open(file,"rb")
        opened=1
    else:
        f=file
        opened=0                # opened var flags to close again if we open the file
    i=80
    hs=f.read(512)             # always start with a 512 byte header
    block=hs
    Header={}                  # dict to take results
    while i < 512 :            # wander along the 512 bytes
        key,val=block[i-80:i].split(":",1)   # uses 80 char lines in key : value format
        key=key.strip()         # remove the whitespace (why?)
        val=val.strip()
        if Header.has_key(key):             # append lines if key already there
            Header[key]=Header[key]+'\n'+val
        else:
            Header[key]=val
        i=i+80                  # next 80 characters
    nhdrblks=int(Header['HDRBLKS'])    # we must have read this in the first 512 bytes.
    # print "got first chunk, headerblocks is",nhdrblks
    # Now read in the rest of the header blocks, appending to what we have
    rest=f.read(512*(nhdrblks-1))
    block = block[i-80:512] + rest
    hs=hs+rest
    j=512*nhdrblks
    while i < j :
        # print i,"*",block[i-80:i].strip(),"*"
        if block[i-80:i].find(":") > 0:          # as for first 512 bytes of header
            key,val=block[i-80:i].split(":",1)
            key=key.strip()
            val=val.strip()
            if Header.has_key(key):
                Header[key]=Header[key]+'\n'+val
            else:
                Header[key]=val
        i=i+80
    Header['datastart']=f.tell()                # make a header item called "datastart"
    if(opened):f.close()
#   print hs
    Header['headerstring']=hs
#   print Header['datastart'],len(hs)
    return Header     # s



def readbruker(file):
    """
    Reads in a Bruker file, returning the data and header
    file may be a string or file object
    TODO we should later modify to take ROI ranges somehow (xmin,xmax,ymin,ymax)
    """
    s="string"                 # s is a string
    if type(s) == type(file):  # if arg is a string, open file, else treat as file object
        f=open(file,"rb")
        opened=1
    else:
        f=file
        opened=0
    Header=readbrukerheader(f) # pass in the file pointer, it stays open
    npixelb=int(Header['NPIXELB'])   # you had to read the Bruker docs to know this!
    rows   =int(Header['NROWS'])
    cols   =int(Header['NCOLS'])
    # We are now at the start of the image - assuming readbrukerheader worked
    size=rows*cols*npixelb
    data=readbytestream(f,f.tell(),rows,cols,npixelb,datatype="int",signed='n',swap='n')
    no=int(Header['NOVERFL'])        # now process the overflows
    if no>0:   # Read in the overflows
        # need at least Int32 sized data I guess - can reach 2^21
        data=data.astype(Numeric.UInt32)
        # 16 character overflows, 9 characters of intensity, 7 character position
        for i in range(no):
            ov=f.read(16)
            intensity=int(ov[0:9])
            position=int(ov[9:16])
            r=position%rows           # relies on python style modulo being always +
            c=position/rows           # relies on truncation down
            #print "Overflow ",r,c,intensity,position,data[r,c],data[c,r]
            data[c,r]=intensity
    f.close()
    Header["rows"]=rows
    Header["columns"]=cols
    return Header,data


def openbruker(filename):
    """
    Reads a bruker file into a data object
    """
    h,d=readbruker(filename)
    return data(d,h)

def edfheader(file):
    """
    Reads the header of edf files into a dictionary
    file can be fileobject or filename
    """
    if type(file)==type("string"):
        f=open(filename,"rb")
        opened=1
    else:
        f=file
        opened=0
    # Header comes in 1024 byte blocks, terminated by "}"
    fh=f.read(1024)
    i=1023
    while fh.find("}\n")<0:
        fh+=f.read(1024)
    # Interpret header
    headeritems=fh[1:-1].split(";")
    hd={}
    hd["headerstring"]=fh
    for item in headeritems:
        if item.find("=")>0:
            hd[item.split("=")[0].lstrip().rstrip()]=item.split("=")[1].lstrip().rstrip()
    hd["rows"]=int(hd["Dim_1"])
    hd["columns"]=int(hd["Dim_2"])
    #   print hd['description']
    #   line=hd['description'].split(" ")
    #   print line
    #   while 1:
    #     try:
    #      value=line.pop(-1)
    #      name=line.pop(-1)
    #      print "n,v",name,value
    #      hd[name]=value
    #     except:
    #      break
    #print hd
    # seek back from the end of the file
    #f.seek( -int(hd["Size"]) , 2  )
    #datastring=f.read( int(hd["Size"]) )
    #f.close()
    if "Size" not in hd.keys():
        if hd["DataType"]=="UnsignedShort":
            hd["Size"]=hd["rows"]*hd["columns"]*2
        if hd["DataType"]=="FLOAT":   # fit2d
            hd["Size"]=hd["rows"]*hd["columns"]*4
    if(opened): f.close()
    return hd

def openedf(filename):
    """
    Opens edf files returning data objects
    Can be gzipped or bzipped

    So far only UnsignedShort (eg frelon) or FLOAT (eg fit2d)
    """
    try:
        if filename[-3:] == ".gz":
            f=gzip.GzipFile(filename,"rb")
        elif filename[-4:] == ".bz2":
            f=bz2.BZ2File(filename,"rb")
        else:
            f=open(filename,"rb")
    except: # Attempt to find a .gz file
        try:
            f=gzip.GzipFile(filename+".gz","rb")
        except:
            try:
                f=bz2.BZ2File(filename+".bz2","rb")
            except:
                print "Cannot manage to open %s"%(filename)
                raise
    hd=edfheader(f)
    # TODO : USE READBYTESTREAM
    try:
        # seek back from the end of the file - fails on gzipped so read all
        datastring=f.read()[-int(hd["Size"]):] # all of the data
    except:
        print hd
        raise
    f.close()
    # Convert datastring to Numeric array
    numerictype = "dunno"
    if hd["DataType"]=="UnsignedShort": numerictype=Numeric.UInt16
    if hd["DataType"]=="FLOAT": numerictype=Numeric.Float32
    if numerictype == "dunno":
        raise TypeError("Unimplemented edf filetype")
    if hd["ByteOrder"]=="LowByteFirst":
        ar=Numeric.reshape(
              Numeric.fromstring(datastring,numerictype),
              (hd["columns"],hd["rows"]) )
    else:
        ar=Numeric.reshape(
              Numeric.fromstring(datastring,numerictype).byteswapped(),
              (hd["columns"],hd["rows"]) )
    return data(ar,hd)


def writebruker(filename,dataobject):
    try:
        hs = dataobject.header["headerstring"]
        o = hs.find("NOVERFL")
        b = hs.find("NPIXELB")
        assert dataobject.header["NROWS"] == dataobject.data.shape[0], "dimensions must match!"
        assert dataobject.header["NCOLS"] == dataobject.data.shape[1], "dimensions must match!"
        if dataobject.data.typecode() == Numeric.UInt16:
            # 2 bytes per pixel, no overflows needed
            noverfl = 0
            bytespp = 2
        else:
            # data must be positive
            assert Numeric.minimum.reduce(Numeric.ravel(dataobject.data)) > 0. , "data must be positive!"
        hs = hs.replace(hs[o:o+80],"NOVERFL:%10d"%(noverfl)+" "*(80-10-8))
        hs = hs.replace(hs[o:o+80],"NPIXELB:%10d"%(bytespp)+" "*(80-10-8))
        
        
    except:
        raise Exception("Sorry, no headerstring in your bruker dataobject")

def writedata(filename,dataobject):
    # identify type of dataobject
    if dataobject.has_key("FORMAT") and int(dataobject["FORMAT"]) == 86:
        writebruker(filename,dataobject)
    raise Exception("Not implemented yet, sorry")



if __name__=="__main__":
    import sys,time
    if len(sys.argv)!=2:
        print "Usage: %s filename" % (sys.argv[0])
        sys.exit()
    t1=time.time()
    testdata=opendata(sys.argv[1])
    t2=time.time()
    print "Time to read file =",t2-t2,"/s"
    print "Rows              =",testdata.header['rows']
    print "Columns           =",testdata.header['columns']
    print "Maximum           =",Numeric.maximum.reduce(Numeric.ravel(testdata.data))
    print "Minimum           =",Numeric.minimum.reduce(Numeric.ravel(testdata.data))
    t3=time.time()
    print "Time native ops   =",t3-t2,"/s"
    s =sum( Numeric.ravel(testdata.data).astype(Numeric.Float32) )  # 16 bit overflows
    sq=sum( Numeric.pow( ravel(testdata.data).astype(Numeric.Float32), 2) )
    n=testdata.header['rows']*testdata.header['columns']
    print "Sum               =",s
    print "Average           =",s/n
    print "Variance          =",Numeric.sqrt(sq/n - (s/n)*(s/n))
    t4=time.time()
    print "Time float ops    =",t4-t3,"/s"
    print "Total run time    =",t4-t1,"/s"
