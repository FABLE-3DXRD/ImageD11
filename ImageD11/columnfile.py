

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

import numpy.oldnumeric as n

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
    "sum_intensity^2"
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
    "tth_hist_prob", 
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
    "Grain"
    ]
FORMATS = {}


# Make a dictionary for formatstrings when writing files
for f in FLOATS: 
    FORMATS[f] = "%.4f" 
for f in INTS:
    FORMATS[f] = "%.0f"

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
            col = col.astype( n.Int )
            mskfun = lambda x, v, t: x == v
        else:        # floating point
            mskfun = lambda x, v, t: n.abs( x - v ) < t 
        mask  = mskfun( col, values[0], tol )
        for v in values[1:]:
            n.logical_or( mskfun( col, v, tol ), mask, mask)
        self.filter( ~mask )


    def writefile(self, filename):
        """
        write an ascii columned file
        """
        self.chkarray()
        self.parameters.saveparameters(filename)
        fout = open(filename,"w+") # appending
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
        data = []
        i = 0
        try:
            fileobj = open(filename,"r")
        except:
            raise Exception("Cannot open %s for reading"%(filename))
        for line in fileobj:
            i += 1
            if len(line.lstrip())==0:
                # skip blank lines
                continue
            if line[0] == "#":
                # title line
                if line.find("=") > -1:
                    # key = value line
                    name, value = clean(line[1:].split("="))
                    self.parameters.addpar(
                        parameters.par( name, value ) )
                else:
                    self.titles = clean(line[1:].split())
                    self.ncols = len(self.titles)
                continue
            # Assume a data line
            try:
                vals = [ float(v) for v in line.split() ]
            except:
                raise Exception("Non numeric data on line\n"+line)
            if len(vals) != self.ncols:
                raise Exception("Badly formatted column file\n"\
                                    "expecting %d columns, got %d\n"\
                                    " line %d in file %s"%
                                (self.ncols, len(vals),
                                 i, self.filename))
            self.nrows += 1
            data.append(vals)
        self.bigarray = n.transpose(n.array(data))
        assert self.bigarray.shape == (self.ncols, self.nrows)
        self.set_attributes()

    def set_attributes(self):
        """
        Set object vars to point into the big array
        """
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
        self.nrows = n.sum(n.compress(mask, n.ones(len(mask))))
        self.bigarray = n.compress(mask, self.bigarray, axis = 1)
        assert self.bigarray.shape == (self.ncols, self.nrows)
        self.set_attributes()
 
    def copy(self):
        """
        Returns a (deep) copy of the columnfile
        """
        cn = columnfile(self.filename, new = True)
        self.chkarray()
        cn.bigarray = self.bigarray.copy()
        cn.titles = [t for t in self.titles ]
        cn.parameters = self.parameters
        (cn.ncols, cn.nrows) = cn.bigarray.shape
        cn.set_attributes()
        return cn

    def chkarray(self):
        """
        Ensure self.bigarray holds our attributes
        """
        for i,t in enumerate(self.titles):
            self.bigarray[i] = getattr(self, t)

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
            self.bigarray = n.array( lis )
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
        
if __name__=="__main__":
    import numpy, sys, time, os
    from ImageD11.columnfile import columnfile
    start = time.time()
    c = columnfile(sys.argv[1])
    print "ImageD11",time.time()-start
    start = time.time()
    c = numpy.loadtxt(sys.argv[1])
    print c.shape
    print "numpy",time.time()-start
    os.system("time -p ./a.out")
