

"""
columnfile represents an ascii file with titles begining "#"
and multiple lines of data

An equals sign "=" on a "#" line implies a parameter = value pair
"""


# ImageD11_v1.0 Software for beamline ID11
# Copyright (C) 2005-2007  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  0211-1307  USA


from ImageD11 import parameters

import numpy

FLOATS = [
    "fc",
    "sc", 
    "omega",
    "f_raw",
    "s_raw",
    "sigf",
    "sigs",
    "covsf",
    "sigo",
    "covso",
    "covfo",
    "sum_intensity",
    "sum_intensity^2",
    "IMax_int",
    "IMax_o",
    "avg_intensity",
    "Min_o",
    "Max_o",
    "dety",
    "detz",
    "gx", "gy", "gz",
    "hr", "kr", "zr",
    "xl", "yl", "zl",
    "drlv2",
    "tth",
    "eta",
    "tth_hist_prob"
    ]

INTS = [
    "Number_of_pixels",
    "IMax_f",
    "IMax_s",
    "Min_f",
    "Max_f",
    "Min_s",
    "Max_s",
    "spot3d_id",
    "h", "k", "l",
    "onfirst", "onlast", "labels",
    "labels",
    "Grain",
    "grainno",
    "grain_id",
    "IKEY", 
    ]

# 9 elements
ij = ["%d%d"%(i,j) for i in range(1,4) for j in range(1,4)]
# Uij, UBIij
LONGFLOATS  = [s+v for v in ij for s in ["U","UBI"]] 





# symmetric 6 elements
ijs = [11,22,33,23,13,12]
#               'eps'+ij+'_s'
EXPONENTIALS = [ h+str(v)+t for v in ijs for h,t in [ ('eps',''),
                                                      ('eps','_s'),
                                                      ('sig',''),
                                                      ('sig','_s') ] ]

#                              'e' ij1  'e' ij2   '_s'  for triangle ij
EXPONENTIALS +=  ["%s%d%s%d%s"%(h,ijs[i],h,ijs[j],t) 
                  for i in range(6) 
                  for j in range(i,6) 
                  for h,t in [ ('e',''),('e','_s'),('s',''),('s','_s')] ]

# testing for line compression
#from ImageD11.columnfile import LONGFLOATS as oldf
#from ImageD11.columnfile import EXPONENTIALS as olde
#assert set(oldf) == set(LONGFLOATS)
#assert set(olde) == set(EXPONENTIALS)
#print "These seem to match"

   
FORMATS = {}


# Make a dictionary for formatstrings when writing files
for f in FLOATS: 
    FORMATS[f] = "%.4f" 
for f in LONGFLOATS: 
    FORMATS[f] = "%.12f" 
for f in INTS:
    FORMATS[f] = "%.0f"
for f in EXPONENTIALS:
    FORMATS[f] = "%.4e"

def clean(str_lst): 
    """ trim whitespace from titles """
    return [s.lstrip().rstrip() for s in str_lst] 

class columnfile:
    """
    Class to represent an ascii file containing multiple named columns
    """
    
    def __init__(self, filename = None, new = False):
        self.filename = filename
        self.bigarray = None
        self.titles = []
        if filename is not None:
            self.parameters = parameters.parameters(filename=filename)
        else:
            self.parameters = parameters.parameters()
        self.ncols = 0
        self.nrows = 0
        if not new:
            self.readfile(filename)

    def removerows( self, column_name, values, tol = 0 ):
        """
        removes rows where self.column_name == values
        values is a list of values to remove
        column name should be in self.titles
        tol is for floating point (fuzzy) comparisons versus integer
        """
        col = self.getcolumn( column_name )
        if tol <= 0: # integer comparisons
            col = col.astype( numpy.int )
            mskfun = lambda x, val, t: x == val
        else:        # floating point
            mskfun = lambda x, val, t: numpy.abs( x - val ) < t 
        mask  = mskfun( col, values[0], tol )
        for val in values[1:]:
            numpy.logical_or( mskfun( col, val, tol ), mask, mask)
        self.filter( ~mask )


    def writefile(self, filename):
        """
        write an ascii columned file
        """
        self.chkarray()
        fout = open(filename,"w+") # appending
        # Write as "# name = value\n"
        parnames = self.parameters.get_parameters().keys()
        parnames.sort()
        for p in parnames:
            fout.write("# %s = %s\n"%(p, str(self.parameters.get(p) ) ) )
        # self.parameters.saveparameters(filename)
        # Now titles line
        fout.write("#")
        format_str = ""
        for title in self.titles:
            fout.write("  %s"%(title))
            try:
                format_str += "  %s" % (FORMATS[title])
            except KeyError:
                format_str += "  %f"
        fout.write("\n")
        format_str += "\n"
        for i in range(self.nrows):
            fout.write(format_str % tuple( self.bigarray[:, i]) )
        fout.close()

    def readfile(self, filename):
        """
        Reads in an ascii columned file
        """
        self.titles = []
        self.bigarray = None
        self.parameters = parameters.parameters(filename=filename)
        self.ncols = 0
        self.nrows = 0
        i = 0
        try:
            # Check if this is a hdf file: magic number
            magic = open(filename,"rb").read(4)
        except:
            raise Exception("Cannot open %s for reading"%(filename))
        #              1 2 3 4 bytes
        if magic == '\x89HDF':
            print "Reading your columnfile in hdf format"
            colfile_from_hdf( filename, obj = self )
            return
        try:
            raw = open(filename,"r").readlines()
        except:
            raise Exception("Cannot open %s for reading"%(filename))
        header = True
        while header and i < len(raw):
             if len(raw[i].lstrip())==0:
                 # skip blank lines
                 i += 1    
                 continue
             if raw[i][0] == "#":
                 # title line
                 if raw[i].find("=") > -1:
                     # key = value line
                     name, value = clean(raw[i][1:].split("="))
                     self.parameters.addpar(
                         parameters.par( name, value ) )
                 else:
                     self.titles = raw[i][1:].split()
                 i += 1    
             else:
                 header = False
        try:
            cc = [numpy.fromstring(v,sep=' ') for v in raw[i:]]
            self.bigarray = numpy.array(cc).transpose()
        except:
            raise Exception("Non numeric data on all lines\n")
        (self.ncols, self.nrows) = self.bigarray.shape

#         data = []
#         try:
#             fileobj = open(filename,"r").readlines()
#         except:
#             raise Exception("Cannot open %s for reading"%(filename))
#         for line in fileobj:
#             i += 1
#             if len(line.lstrip())==0:
#                 # skip blank lines
#                 continue
#             if line[0] == "#":
#                 # title line
#                 if line.find("=") > -1:
#                     # key = value line
#                     name, value = clean(line[1:].split("="))
#                     self.parameters.addpar(
#                         parameters.par( name, value ) )
#                 else:
#                     self.titles = clean(line[1:].split())
#                     self.ncols = len(self.titles)
#                 continue
#             # Assume a data line
#             try:
#                 vals = [ float(v) for v in line.split() ]
#             except:
#                 raise Exception("Non numeric data on line\n"+line)
#             if len(vals) != self.ncols:
#                 raise Exception("Badly formatted column file\n"\
#                                     "expecting %d columns, got %d\n"\
#                                     " line %d in file %s"%
#                                 (self.ncols, len(vals),
#                                  i, self.filename))
#             self.nrows += 1
#             data.append(vals)
#         self.bigarray = numpy.transpose(numpy.array(data))
#         if self.nrows > 0:
#             assert self.bigarray.shape == (self.ncols, self.nrows)
        self.set_attributes()

    def set_attributes(self):
        """
        Set object vars to point into the big array
        """
        if self.nrows == 0:
            # use empty arrays for now... not sure why this was avoided in the past?
            pass
            #return
        for title, i in zip(self.titles, range(len(self.titles))):
            setattr(self, title, self.bigarray[i])
            assert getattr(self, title).shape == (self.nrows,)

    def filter(self, mask):
        """
        mask is an nrows long array of true/false
        """
        self.chkarray()
        if len(mask) != self.nrows:
            raise Exception("Mask is the wrong size")
        self.nrows = int(numpy.sum(
            numpy.compress(mask, numpy.ones(len(mask)))))
        self.bigarray = numpy.compress(mask, self.bigarray, axis = 1)
        assert self.bigarray.shape == (self.ncols, self.nrows)
        self.set_attributes()
 
    def copy(self):
        """
        Returns a (deep) copy of the columnfile
        """
        cnw = columnfile(self.filename, new = True)
        self.chkarray()
        cnw.bigarray = self.bigarray.copy()
        cnw.titles = [t for t in self.titles ]
        cnw.parameters = self.parameters
        (cnw.ncols, cnw.nrows) = cnw.bigarray.shape
        cnw.set_attributes()
        return cnw

    def chkarray(self):
        """
        Ensure self.bigarray holds our attributes
        """
        for i, title in enumerate(self.titles):
            self.bigarray[i] = getattr(self, title)

    def addcolumn(self, col, name):
        """
        Add a new column col to the object with name "name"
        Overwrites in the case that the column already exists
        """
        if len(col) != self.nrows:
            raise Exception("Wrong length column")
        if name in self.titles:
            # Make this overwrite instead of throwing an exception
            self.bigarray[self.titles.index(name)] = col
            # raise Exception("Already got a column called "+name)
            setattr(self, name,             
                    self.bigarray[self.titles.index(name)] )
        else:
            self.titles.append(name)
            self.ncols += 1
            lis = list(self.bigarray)
            lis.append(col)
            self.bigarray = numpy.array( lis )
            assert self.bigarray.shape == (self.ncols, self.nrows)
            self.set_attributes()

    # Not obvious, but might be a useful alias
    setcolumn = addcolumn

    def getcolumn(self, name):
        """
        Gets data, if column exists
        """
        if name in self.titles:
            return self.bigarray[self.titles.index(name)]
        raise Exception("Name "+name+" not in file")
    
class newcolumnfile(columnfile):
    """ Just like a columnfile, but for creating new
    files """
    def __init__(self, titles):
        columnfile.__init__(self, filename=None, new=True)
        self.titles = titles
        self.ncols = len(titles)


try:
    import h5py, os
    def colfile_to_hdf( colfile, hdffile, name=None, compression=None,
                        compression_opts=None):
        """
        Copy a columnfile into hdf file
        FIXME TODO - add the parameters somewhere (attributes??)
        """
        if isinstance(colfile, columnfile):
            c = colfile
        else:
            c = columnfile( colfile )
        if isinstance(hdffile, h5py.File):
            h = hdffile
            opened = False
        else:
            h = h5py.File( hdffile , 'a') # Appending if exists
            opened = True
        if name is None:
            # Take the file name
            name = os.path.split(c.filename)[-1]
        if name in h.keys():
            g = h[name]
        else:
            g = h.create_group( name )
        g.attrs['ImageD11_type'] = 'peaks'
        for t in c.titles:
            if t in INTS:
                ty = numpy.int32
            else:
                ty = numpy.float32
            # print "adding",t,ty
            dat = getattr(c, t).astype( ty )
            if t in g.keys():
                g[t][:] = dat
            else:
                g.create_dataset( t, data = dat,
                                  compression=compression,
                                  compression_opts=compression_opts )
        if opened:
            h.close()

    def colfileobj_to_hdf( cf, hdffile, name=None):
        """
        Save a columnfile into hdf file format
        FIXME TODO - add the parameters somewhere (attributes??)
        """
        h = h5py.File(hdffile, 'a' )
        if name is None:
            name = str(cf.filename)
        try:
            g = h.create_group( name )
        except:
            print name, h
            raise
        g.attrs['ImageD11_type'] = 'peaks'
        for t in c.titles:
            if t in INTS:
                ty = numpy.int32
            else:
                ty = numpy.float32
            g.create_dataset( t, data = getattr(c, t).astype( ty ) )
        h.close()

    def colfile_from_hdf( hdffile , name=None, obj=None ):
        """
        Read a columnfile from a hdf file
        FIXME TODO - add the parameters somewhere (attributes??)
        """
        import time
        h = h5py.File( hdffile, 'r' )
        if hasattr(h, 'listnames'):
            groups = h.listnames()
        else: # API changed
            groups = h.keys()
        if name is not None:
            if name in groups:
                g = h[name]
            else:
                print groups
                raise Exception("Did not find your "+str(name)+" in "+hdffile)
        else:
            assert len(groups) == 1, "Your hdf file has many groups. Which one??"
            g = h[groups[0]]
            name = groups[0]
        if hasattr(g, 'listnames'):
            titles = g.listnames()
        else: # API changed
            titles = g.keys()
        otitles = [t for t in titles]
        otitles.sort()
        newtitles = []
        # Put them back in the order folks might have hard wired in their
        # programs
        for t in ['sc', 'fc', 'omega' , 'Number_of_pixels',  'avg_intensity',
        's_raw',  'f_raw',  'sigs',  'sigf',  'covsf' , 'sigo',  'covso',
        'covfo',  'sum_intensity',  'sum_intensity^2',  'IMax_int',  'IMax_s',
        'IMax_f',  'IMax_o',  'Min_s',  'Max_s',  'Min_f',  'Max_f',  'Min_o',
        'Max_o',  'dety',  'detz',  'onfirst',  'onlast',  'spot3d_id',  'xl',
        'yl',  'zl',  'tth',  'eta',  'gx',  'gy',  'gz']:
            if t in otitles:
                newtitles.append(t)
                otitles.remove(t)
        # Anything else goes in alphabetically
        [newtitles.append(t) for t in otitles]
        assert len(newtitles) == len( titles )
        if obj is None:
            col = columnfile( filename=name, new=True )
        else:
            col = obj
        col.titles = newtitles
        dat = g[newtitles[0]][:]
        col.bigarray = numpy.zeros( (len(newtitles), len(dat) ), numpy.float32)
        col.bigarray[0] = dat
        col.ncols = len(newtitles)
        col.nrows = len(dat)
        i = 1
        for t in newtitles[1:]:
            col.bigarray[i] = g[t][:]
            i += 1
        col.set_attributes()        
        return col
    
except ImportError:
    def hdferr():
        raise Exception("You do not have h5py installed!")

    def colfile_to_hdf( a,b,name=None):
        hdferr()
    
    def colfile_from_hdf( hdffile , name=None ):
        hdferr()

    def colfileobj_to_hdf( cf, hdffile, name=None):
        hdferr()

def bench():
    """
    Compares the timing for reading with columfile
    versus numpy.loadtxt
    """
    import sys, time
    start = time.time()
    colf = columnfile(sys.argv[1])
    print colf.bigarray.shape
    print "ImageD11", time.time() - start
    start = time.time()
    nolf = numpy.loadtxt(sys.argv[1])
    print nolf.shape
    print "numpy", time.time() - start
    
    # os.system("time -p ./a.out")




import sqlite3 as database_module
# Perhaps other modules follow the same api.
# Doubtless one does not get away with using a filename?


def colfile2db( colfilename, dbname ):
    """
    Read the columnfile into a database
    Ignores parameter metadata (not used yet)
    """
    colf = columnfile( colfilename )
    dbo = database_module.connect( dbname )
    curs = dbo.cursor()
    # Build up columnames and types to make table
    tablecols = []
    # Not allowed for sql to have ^ in string
    colf.titles = [t.replace("^","_pow_") for t in colf.titles]
    for name in colf.titles:
        if name in INTS:
            tablecols.append(name + " INTEGER")
            continue
        if name in FLOATS:
            tablecols.append(name + " REAL")
            continue
        tablecols.append(name + " REAL")
    curs.execute("create table peaks \n( " + \
                 " , ".join(tablecols)     + " ) ; \n" )
    # Make a format string for inserting data
    ins = "insert into peaks values ("  + \
          ",".join(["?"]*colf.ncols)    +") ;"
    # insert the data
    for i in range(colf.nrows):
        curs.execute( ins , tuple(colf.bigarray[:, i]) )
    curs.close()
    dbo.commit()
    dbo.close()

                                       
if __name__ == "__main__":
    bench()
 
