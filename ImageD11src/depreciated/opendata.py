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
import numpy as np

def opendata(filename):
    """
    Guesses the file type and return a data object

    Supports edf, edf.gz, edf.bz2, xye, epf, inp, chi
     ... otherwise tries bruker
    """
    f = filename.rstrip()
    # add spd file extensions
    if f[-3:]in ["edf", "EDF","cor","COR"] :
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

def getnum(name):
    # try to figure out the file number
    # guess it starts at the back
    nl = []
    first = False
    for c in name[::-1]: # this means iterate backwards through the string
        if c.isdigit():
            first = True
            nl.append(c)
            continue
        if first: break
    num = "".join(nl[::-1])
    return int(num)

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

    a=np.array(xy,float)
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
    a=np.array(xye)
    h['rows']=a.shape[1]
    h['xye']=True
    return data(a,h)

def readbytestream(file,offset,x,y,nbytespp,datatype='int',signed='n',
                   swap='n',typeout=np.uint16):
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
        if nbytespp==1 : tin=np.uint8
        if nbytespp==2 : tin=np.uint16
        if nbytespp==4 : tin=np.uint32
    if datatype=='int' and signed=='y':
        if nbytespp==1 : tin=np.int8
        if nbytespp==2 : tin=np.int16
        if nbytespp==4 : tin=np.int32
    if datatype=='float':
        tin=np.float32
    if datatype=='double' :
        tin=np.float64
    if tin=="dunno" :
        raise SyntaxError, "Did not understand what type to try to read"
    opened=0
    if isinstance(tin, type(file)):  # Did we get a string or a file point)er?
        f=open(file,'rb')
        opened=1
    else:
        f=file
    f.seek(offset)
    ar = np.array(np.reshape(np.fromstring(f.read(len),tin),(x,y)),typeout)
    if swap=='y':
        ar.byteswap(True)
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
    try:
        npixelb=int(Header['NPIXELB'])   # you had to read the Bruker docs to know this!
    except:
        print "length",len(Header['NPIXELB'])
        for c in Header['NPIXELB']:
            print "char:",c,ord(c)
        raise
    rows   =int(Header['NROWS'])
    cols   =int(Header['NCOLS'])
    # We are now at the start of the image - assuming readbrukerheader worked
    size=rows*cols*npixelb
    data=readbytestream(f,f.tell(),rows,cols,npixelb,datatype="int",signed='n',swap='n')
    no=int(Header['NOVERFL'])        # now process the overflows
    if no>0:   # Read in the overflows
        # need at least Int32 sized data I guess - can reach 2^21
        data=data.astype(np.uint32)
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
        f=open(file,"rb")
        opened=1
    else:
        f=file
        opened=0
    # Header comes in 1024 byte blocks, terminated by "}"
    fh=f.read(1024)
    if len(fh)!=1024:
        raise Exception("File too small")
    i=1023
    j=0
    while fh.find("}\n")<0 and j<10:
        extra = f.read(1024)
        if len(extra)!=1024:
            raise Exception("File too small")
        fh+=extra
        j=j+1
        if j==9:
            raise Exception("Does not look like an edf file, header too long")
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
    if hd["DataType"]=="UnsignedShort": numerictype=np.uint16
    if hd["DataType"]=="FLOAT": numerictype=np.float32
    if hd["DataType"]=="FloatValue": numerictype=np.float32
    if numerictype == "dunno":
        raise TypeError("Unimplemented edf filetype"+hd["DataType"])
    if hd["ByteOrder"]=="LowByteFirst":
        ar=np.reshape(
              np.fromstring(datastring,numerictype),
              (hd["columns"],hd["rows"]) )
    else:
        ar=np.reshape(
              np.fromstring(datastring,numerictype).byteswap(),
              (hd["columns"],hd["rows"]) )
    pass  ## ar.savespace(1)
    return data(ar,hd)


def writebruker(filename,dataobject):
    """
    Writes 16-bit + overflow bruker images based on the headerstring
    Assumes you've got a good bruker headerstring in the data
    object
    """
    try:
        hs = dataobject.header["headerstring"]
    except:
        raise Exception("Sorry, no headerstring in your bruker dataobject")
    assert int(dataobject.header["NROWS"]) == dataobject.data.shape[0], "dimensions must match!"
    assert int(dataobject.header["NCOLS"]) == dataobject.data.shape[1], "dimensions must match!"
    minval = np.minimum.reduce(np.ravel(dataobject.data))
    assert minval >= 0 , "data must be positive! "+str(minval)
    r = np.ravel(dataobject.data)
    d = np.ravel(dataobject.data)
    bytespp = 2
    if dataobject.data.dtype.char == np.uint16:
        # 2 bytes per pixel, no overflows needed
        noverfl = 0
    else:
        # data must be positive
        twobytes = pow(2,16)-1
        indices = np.compress( r > twobytes , range(r.shape[0]) )
        r = np.where(r < twobytes, r, twobytes)
        noverfl = indices.shape[0]
    o = hs.find("NOVERFL")
    b = hs.find("NPIXELB")
    hs = hs.replace(hs[o:o+80],"%-80s"%("NOVERFL:%10d"%(noverfl)))
    hs = hs.replace(hs[b:b+80],"%-80s"%("NPIXELB:%10d"%(bytespp)))
    # for i in range(b,b+80):
    #     print "char:",hs[i],ord(hs[i]),
    # print
    out = open(filename,"wb")
    out.write(hs)
    out.write(r.astype(np.uint16).tostring())
    if noverfl > 0 :
        r = np.ravel(dataobject.data)
        length = 0
        for i in indices:
            s = "%9d%7d"%(r[i],i)
            length += 16
            out.write(s)
        out.write(" "*(512-length%512)) # pad to 512 multiple
    out.close()

def makeedfheader(dataobject):
    """
    Writes the FLOAT edf subset a la fit2d
    Feel free to improve me!
    """
    h =     """HeaderID       = %s ;
Image          = %s ;
ByteOrder      = %s ;
DataType       = %s ;
Dim_1          = %s ;
Dim_2          = %s ;
Size           = %s ;
count_time     = %s ;
point_no       = %s ;
preset         = %s ;
col_end        = %s ;
col_beg        = %s ;
row_end        = %s ;
row_beg        = %s ;
col_bin        = %s ;
row_bin        = %s ;
time           = %s ;
time_of_day    = %s ;
dir            = %s ;
suffix         = %s ;
prefix         = %s ;
run            = %s ;
title          = %s ;
DATE (scan begin)= %s ;"""
    hd = dataobject.header
    hd["DataType"]="FLOAT"
    hd["Size"]=str(np.ravel(dataobject.data).shape[0]*4)
    keys = [s.split("=")[0].lstrip().rstrip() for s in h.split(";")]
    keys.pop()
    args = tuple([str(hd[k.lstrip().rstrip()]) for k in keys])
    s = "{\n"+h%args
    keys = keys + ["headerstring","rows","columns"]
    ok = hd.keys()
    ok.sort()
    for k in ok:
        if k not in keys:
            s = "%s\n%-12s = %s ;"%(s,k,hd[k])
    return "%-4094s}\n"%(s)



def writeedf(filename,dataobject):
    """
    """
    h=makeedfheader(dataobject)
    assert(len(h)==4096) , "bad edf header size "+str(len(h))
    f=open(filename,"wb")
    f.write(h)
    f.write(dataobject.data.astype(np.float32).tostring())
    f.close()

def writedata(filename,dataobject):
    # identify type of dataobject
    if dataobject.header.has_key("FORMAT"):
        if int(dataobject.header["FORMAT"]) == 86:
            writebruker(filename,dataobject)
            return
    if dataobject.header.has_key("HeaderID"):
        if dataobject.header["HeaderID"] == "EH:000001:000000:000000":
            writeedf(filename,dataobject)
            return
    raise Exception("Not implemented yet, I can only write bruker files, sorry")





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
    print "Maximum           =",np.maximum.reduce(np.ravel(testdata.data))
    print "Minimum           =",np.minimum.reduce(np.ravel(testdata.data))
    t3=time.time()
    print "Time native ops   =",t3-t2,"/s"
    s =sum( np.ravel(testdata.data).astype(np.float32) )  # 16 bit overflows
    sq=sum( np.pow( ravel(testdata.data).astype(np.float32), 2) )
    n=testdata.header['rows']*testdata.header['columns']
    print "Sum               =",s
    print "Average           =",s/n
    print "Variance          =",np.sqrt(sq/n - (s/n)*(s/n))
    t4=time.time()
    print "Time float ops    =",t4-t3,"/s"
    print "Total run time    =",t4-t1,"/s"
