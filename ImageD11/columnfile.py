

from __future__ import print_function

import io
from functools import partial

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

import warnings

from ImageD11 import parameters, transform
import numpy as np

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
    "spot4d_id",
    "h", "k", "l",
    "onfirst", "onlast", "labels",
    "labels",
    "Grain",
    "grainno",
    "grain_id",
    "IKEY",
    "npk2d",
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


def fillcols(lines, cols):
    for i,line in enumerate(lines): # could be parallel
        for j,item in enumerate(line.split()):
            cols[j][i] = float(item)

class columnfile(object):
    """
    Class to represent an ascii file containing multiple named columns
    """

    def __init__(self, filename = None, new = False):
        self.titles = []
        self.filename = filename
        self.__data = []
        if filename is not None:
            self.parameters = parameters.parameters(filename=filename)
        else:
            self.parameters = parameters.parameters()
        self.ncols = 0
        self.nrows = 0
        if not new:
            self.readfile(filename)

    def get_bigarray(self):
        # if someone uses this we have to go back to the old
        # representation
        if not hasattr(self,"__bigarray") or len(self.__data) != len(self.__bigarray):
            self.__bigarray = np.asarray( self.__data )
        self.__data = self.__bigarray
        return self.__bigarray

    def set_bigarray(self, ar):
#        print("setting bigarray",len(ar),len(ar[0]))
#        warnings.filter("once")
#        warnings.warn("Setting bigarray on colfile", stacklevel=2)
        assert len(ar) == len(self.titles), \
            "Wrong length %d to set bigarray"%(len(ar))+\
            " ".join(self.titles)
        nrows = len(ar[0])
        for col in ar:
            assert len(col) == nrows, "ar is not rectangular"
        self.nrows = nrows
        # use a list of arrays
        self.__bigarray = ar
        self.__data = self.__bigarray
        self.set_attributes()

    bigarray = property(fget=get_bigarray, fset=set_bigarray)

    def set_attributes(self):
        """
        Set object vars to point into the big array
        """
        if self.nrows == 0:
            # use empty arrays for now...
            # not sure why this was avoided in the past?
            pass
            #return
        for i, name in enumerate(self.titles):
            setattr(self, name, self.__data[i])
            a  = getattr(self, name)
            assert len(a) == self.nrows, "%s %d %d"%(name,len(a),self.nrows)

    def __getitem__(self, key):
        if key in self.titles:
            return self.getcolumn( key )
        else:
            raise KeyError

    def __setitem__(self, key, value):
        if key in self.titles:
            self.getcolumn(key)[:] = value
        else:
            self.addcolumn(value, key)

    def __setattr__(self, key, value):
        if key == 'titles':
            super(columnfile, self).__setattr__(key, value)
            return
        if key in self.titles:
            if np.isscalar(value): # broadcast
                self.__data[self.titles.index(key)][:] = value
            else:
                assert len(value) == self.nrows
                self.__data[self.titles.index(key)] = value
        super(columnfile, self).__setattr__(key, value)

    def keys(self):
        return self.titles

    def removerows( self, column_name, values, tol = 0 ):
        """
        removes rows where self.column_name == values
        values is a list of values to remove
        column name should be in self.titles
        tol is for floating point (fuzzy) comparisons versus integer
        """
        col = self.getcolumn( column_name )
        if tol <= 0: # integer comparisons
            col = col.astype( int )
            mskfun = lambda x, val, t: x == val
        else:        # floating point
            mskfun = lambda x, val, t: np.abs( x - val ) < t
        mask  = mskfun( col, values[0], tol )
        for val in values[1:]:
            np.logical_or( mskfun( col, val, tol ), mask, mask)
        self.filter( ~mask )

    def sortby( self, name ):
        """
        Sort arrays according to column named "name"
        """
        col = self.getcolumn( name )
        order = np.argsort( col )
        self.reorder( order )

    def reorder( self, indices ):
        """
        Put array into the order given by indices
        ... normally indices would come from np.argsort of something
        """
        for col in self.__data:
            col[:] = col[indices]
        self.set_attributes()

    def writefile(self, filename):
        """
        write an ascii columned file
        """
        self.chkarray()
        fout = open(filename,"w") # appending
        # Write as "# name = value\n"
        parnames = list(self.parameters.get_parameters().keys())
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
            fout.write(format_str % tuple( [col[i] for col in self.__data] ) )
        fout.close()

    def readfile(self, filename):
        """
        Reads in an ascii columned file
        """
        self.titles = []
        self.parameters = parameters.parameters(filename=filename)
        self.ncols = 0
        self.nrows = 0
        i = 0
        # Check if this is a hdf file: magic number
        with open(filename,"rb") as f:
            magic = f.read(4)
        #              1 2 3 4 bytes
        if magic == b'\x89HDF':
            print("Reading your columnfile in hdf format")
            colfile_from_hdf( filename, obj = self )
            return
        with open(filename,"r") as f:
            raw = f.readlines()
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
                     name, value = clean(raw[i][1:].split("=",1))
                     self.parameters.addpar(
                         parameters.par( name, value ) )
                 else:
                     self.titles = raw[i][1:].split()
                 i += 1
             else:
                 header = False
        try:
            row0 = [ float( v ) for v in raw[i].split() ]
            lastrow = [ float( v ) for v in raw[-1].split() ]
            if len(row0) == len(lastrow ):
                nrows = len(raw)-i
                last = len(raw)
            else:
                nrows = len(raw)-i-1 # skip the last row
                last = len(raw)-1
            cols = [ np.empty( nrows , float ) for _ in range(len(row0))]
            fillcols( raw[i:last], cols )
            self.__data=cols
        except:
            raise # Exception("Problem interpreting your colfile")
        self.ncols, self.nrows = len(row0), nrows
        self.parameters.dumbtypecheck()
        self.set_attributes()


    def filter(self, mask):
        """
        mask is an nrows long array of true/false
        """
        self.chkarray()
        if len(mask) != self.nrows:
            raise Exception("Mask is the wrong size")
        msk = np.array( mask, dtype=bool )
        # back to list here
        # previously makes a 2X memory footprint here:
        # self.__data = [col[msk] for col in self.__data]
        #
        if not isinstance( self.__data, list ):
            self.__data = list( self.__data )
        for i, col in enumerate( self.__data ):   # could be parallel
            self.__data[i] = col[msk]
        self.nrows = len(self.__data[0])
        self.set_attributes()

    def copy(self):
        """
        Returns a (deep) copy of the columnfile
        """
        cnw = columnfile(self.filename, new = True)
        self.chkarray()
        cnw.titles = [t for t in self.titles ]
        cnw.parameters = parameters.parameters( **self.parameters.parameters.copy() )
        cnw.set_bigarray( [col.copy() for col in self.__data] )
        cnw.ncols = self.ncols
        cnw.set_attributes()
        return cnw

    def copyrows(self, rows):
        """
        Returns a copy of select rows of the columnfile
        """
        self.chkarray()
        cnw = columnfile(self.filename, new = True)
        cnw.titles = [t for t in self.titles ]
        cnw.parameters = parameters.parameters( **self.parameters.parameters.copy() )
        cnw.bigarray = [col[rows] for col in self.__data]
        #cnw.ncols, cnw.nrows = cnw.bigarray.shape
        #cnw.set_attributes()
        return cnw

    def chkarray(self):
        """
        Ensure self.bigarray holds our attributes
        """
        for i, name in enumerate(self.titles):
            a = getattr(self, name)
            self.__data[i] = a

    def addcolumn(self, col, name):
        """
        Add a new column col to the object with name "name"
        Overwrites in the case that the column already exists
        """
        if len(col) != self.nrows:
            raise Exception("Wrong length column")
        if name in self.titles:
            idx = self.titles.index(name)
            # Make this overwrite instead of throwing an exception
            self.__data[idx] = col
            # raise Exception("Already got a column called "+name)
        else:
            # assert self.bigarray.shape == (self.ncols, self.nrows)
            data = np.asanyarray( col )
            assert  data.shape[0] == self.nrows
            self.titles.append(name)
            idx = len(self.titles)-1
            self.ncols += 1
            self.__data.append( data )
        setattr(self, name, self.__data[idx] )

    def setcolumn(self, col, name):
        assert name in self.titles
        self.addcolumn( col, name )

    def getcolumn(self, name):
        """
        Gets data, if column exists
        """
        if name in self.titles:
            return self.__data[self.titles.index(name)]
        raise KeyError("Name "+name+" not in file")

    def setparameters( self, pars ):
        """
        update the parameters
        """
        self.parameters = pars
        self.parameters.dumbtypecheck()

    def updateGV(self, pars=None, translation=None, fast=True ):
        """ This only computes gvectors """
        if pars is not None:
            self.setparameters( pars )
        pars = self.parameters
        assert "omega" in self.titles
        if "sc" in self.titles and "fc" in self.titles:
            sc, fc = self.sc, self.fc
        elif "xc" in self.titles and "yc" in self.titles:
            sc, fc = self.xc, self.yc
        else:
            raise Exception("columnfile file misses xc/yc or sc/fc")
        if translation is not None:
            tx, ty, tz = translation
        else:
            tx = pars.get('t_x')
            ty = pars.get('t_y')
            tz = pars.get('t_z')
        if fast:
            computer = transform.Ctransform( pars.parameters )
            gve = computer.sf2gv( sc, fc, self.omega, tx, ty, tz )
            for i, name in enumerate(('gx','gy','gz')):
                self.addcolumn( gve[:,i], name )
        else:
            self.updateGeometry( pars = pars, translation=translation, fast=fast )

    def updateGeometry(self, pars=None, translation=None, fast=True ):
        """
        changing or not the parameters it (re)-computes:
           xl,yl,zl = ImageD11.transform.compute_xyz_lab
           tth, eta = ImageD11.transform.compute_tth_eta
           gx,gy,gz = ImageD11.transform.compute_g_vectors
           ds = length of scattering vectors

        The option "fast=True" means using the newer C code from
        ImageD11.transform.Ctransform rather than the older python
        ImageD11.transform code. We keep it for testing.
        """
        if pars is not None:
            self.setparameters( pars )
        pars = self.parameters
        assert "omega" in self.titles
        if "sc" in self.titles and "fc" in self.titles:
            sc, fc = self.sc, self.fc
        elif "xc" in self.titles and "yc" in self.titles:
            sc, fc = self.xc, self.yc
        else:
            raise Exception("columnfile file misses xc/yc or sc/fc")
        if fast:
            if translation is not None:
                tx, ty, tz = translation
            else:
                tx = pars.get('t_x')
                ty = pars.get('t_y')
                tz = pars.get('t_z')
            computer = transform.Ctransform( pars.parameters )
            xyz = computer.sf2xyz( sc, fc )
            for i, name in enumerate(('xl','yl','zl')):
                self.addcolumn( xyz[:,i], name )
            out = computer.xyz2geometry( xyz, self.omega, tx, ty, tz )
            for i,name in enumerate(("tth", "eta", "ds", "gx", "gy", "gz")):
                self.addcolumn( out[:,i], name )
        else:
            pars = { name : pars.parameters[name] for name in pars.parameters }
            if translation is not None:
                pars['t_x'] = translation[0]
                pars['t_y'] = translation[1]
                pars['t_z'] = translation[2]
            xl,yl,zl = transform.compute_xyz_lab( (sc, fc),
                                                  **pars)
            self.addcolumn(xl,"xl")
            self.addcolumn(yl,"yl")
            self.addcolumn(zl,"zl")
            peaks_xyz = np.array((xl,yl,zl))
            assert "omega" in self.titles,"No omega column"
            om = self.omega *  float( pars.get("omegasign") )
            tth, eta = transform.compute_tth_eta_from_xyz(
                peaks_xyz, om,
                **pars)
            self.addcolumn(tth, "tth", )
            self.addcolumn(eta, "eta", )
            gx, gy, gz = transform.compute_g_vectors(
                tth, eta, om,
                wvln  = pars.get("wavelength"),
                wedge = pars.get("wedge"),
                chi   = pars.get("chi") )
            self.addcolumn(gx, "gx")
            self.addcolumn(gy, "gy")
            self.addcolumn(gz, "gz")
            modg =  np.sqrt( gx * gx + gy * gy + gz * gz )
            self.addcolumn(modg, "ds") # dstar



class newcolumnfile(columnfile):
    """ Just like a columnfile, but for creating new
    files """
    def __init__(self, titles):
        columnfile.__init__(self, filename=None, new=True)
        self.titles = titles
        self.ncols = len(titles)
        
        
def colfile_from_dict( c ):
    """ convert from a dictonary of numpy arrays """
    titles = list(c.keys())
    nrows = len(c[titles[0]])
    for t in titles:
        assert len(c[t]) == nrows, t
    colf = newcolumnfile( titles=titles )
    colf.nrows = nrows
    colf.set_bigarray( [ c[t] for t in titles ] )
    return colf


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
            try:
                name = os.path.split(c.filename)[-1]
            except:
                name = 'peaks'
        if name in list(h.keys()):
            g = h[name]
        else:
            g = h.create_group( name )
        g.attrs['ImageD11_type'] = 'peaks'
        for t in c.titles:
            if t in INTS:
                ty = np.int64
            else:
                ty = np.float64
            # print "adding",t,ty
            dat = getattr(c, t).astype( ty )
            if t in list(g.keys()):
                if g[t].shape != dat.shape:
                    try:
                        g[t].resize( dat.shape )
                    except TypeError:
                        raise TypeError("Columnfile already exists on disk with a different length, cannot write HDF5!")
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
            print(name, h)
            raise
        g.attrs['ImageD11_type'] = 'peaks'
        for t in cf.titles:
            if t in INTS:
                ty = np.int64
            else:
                ty = np.float64
            g.create_dataset( t, data = getattr(cf, t).astype( ty ) )
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
            groups = list(h.keys())
        if name is not None:
            if name in groups:
                g = h[name]
            else:
                print(groups)
                raise Exception("Did not find your "+str(name)+" in "+hdffile)
        else:
            groups = [g for g in groups 
                                if 'ImageD11_type' in h[g].attrs and 
                                   h[g].attrs['ImageD11_type'] in ('peaks', b'peaks') ]
            assert len(groups) == 1, "Your hdf file has many groups. Which one??"+str(groups)
            g = h[groups[0]]
            name = groups[0]
        if hasattr(g, 'listnames'):
            titles = g.listnames()
        else: # API changed
            titles = list(g.keys())
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
        col.nrows = len( g[newtitles[0]] )
        for name in newtitles:
            col.addcolumn( g[name][:].copy(), name )
        h.close()
        return col


    def mmap_h5colf( fname, path='peaks', mode = 'r' ):
        """
        From:
        https://gist.github.com/anonymous/19814198398da86b8f98

        The columnfile must have non-compressed data

        Memory maps the underlying array for read only data.
        You will need to use something like copyrows on the result.
        """
        if mode != 'r':
            # The purpose of this code is to read columnfiles across
            # many processes. The implications of writing are messy.
            # You can't have lots of processes writing into the hdf5
            # files as you can't resize. So just force read only.
            raise Exception("Sorry, read only for now")
        with h5py.File(fname,'r') as hin:
            if path not in hin:
                for path in list(hin['/']):
                     if ('ImageD11_type' in hin[path].attrs and 
                         hin[path].attrs['ImageD11_type'] in ('peaks', b'peaks')):
                            break
            grp = hin[path]
            names = list(grp)
            offsets = {}
            dtypes = {}
            shapes = {}
            for name in names:
                ds = grp[name]
                # We get the dataset address in the HDF5 file.
                offsets[name] = ds.id.get_offset()
                # We ensure we have a non-compressed contiguous array.
                assert ds.chunks is None
                assert ds.compression is None
                assert offsets[name] > 0
                dtypes[name] = ds.dtype
                shapes[name] = ds.shape
        ardict = { name : np.memmap(fname,
                                    mode='r',
                                    shape=shapes[name],
                                    offset=offsets[name],
                                    dtype=dtypes[name])
                   for name in names }
        return colfile_from_dict( ardict )



except ImportError:
    def hdferr():
        raise Exception("You do not have h5py installed!")

    def colfile_to_hdf( a,b,name=None):
        hdferr()

    def colfile_from_hdf( hdffile , name=None ):
        hdferr()

    def colfileobj_to_hdf( cf, hdffile, name=None):
        hdferr()


try:
    import pandas as pd
    class PandasColumnfile(columnfile):
        def __init__(self, filename=None, new=False):
            self._df = pd.DataFrame()
            self.filename = filename
            if filename is not None:
                self.parameters = parameters.parameters(filename=filename)
            else:
                self.parameters = parameters.parameters()
            if not new:
                self.readfile(filename)

        def __getattribute__(self, item):
            # called whenever getattr(item, 'attr') or item.attr is called
            # for speed, look in self._df FIRST before giving up
            if item == '_df':
                return object.__getattribute__(self, item)
            if item in self._df.columns:
                return super(PandasColumnfile, self).__getattribute__('_df')[item].to_numpy()
            return super(PandasColumnfile, self).__getattribute__(item)

        def __setattr__(self, key, value):
            if key == "_df":
                object.__setattr__(self, key, value)
                return
            if key in self.titles:
                self._df[key] = value
                return
            super(PandasColumnfile, self).__setattr__(key, value)

        @property
        def nrows(self):
            return self._df.index.size

        @nrows.setter
        def nrows(self, value):
            # we don't set the number of rows explicitly for Pandas dataframes
            pass

        @property
        def ncols(self):
            return len(self._df.columns)

        @property
        def titles(self):
            return self._df.columns.tolist()

        def removerows(self, column_name, values, tol=0):
            """
            removes rows where self.column_name == values
            values is a list of values to remove
            column name should be in self.titles
            tol is for floating point (fuzzy) comparisons versus integer
            """

            if tol <= 0:
                # exact match comparisons
                matching_row_indices = self._df[self._df[column_name].isin(values)].index
                self._df.drop(matching_row_indices, inplace=True)
            else:
                # floating point
                matching_row_indices = self._df.loc[np.abs(self._df[column_name] - values) < tol].index
                self._df.drop(matching_row_indices, inplace=True)

        def sortby(self, name):
            """
            Sort arrays according to column named "name"
            """
            self._df.sort_values(name, inplace=True)

        def reorder(self, indices):
            """
            Put array into the order given by indices
            ... normally indices would come from np.argsort of something
            """
            self._df = self._df.loc[indices]

        def writefile(self, filename):
            """
            write an ascii columned file
            """

            # TODO: Can be greatly sped up with pd.write_csv

            # self.chkarray()

            # create a text buffer
            s_buf = io.StringIO()

            # Write as "# name = value\n"
            parnames = list(self.parameters.get_parameters().keys())
            parnames.sort()
            for p in parnames:
                s_buf.write("# %s = %s\n" % (p, str(self.parameters.get(p))))

            s_buf.write('# ')

            def format_func(value, fmt):
                return fmt % value

            # get a dictionary of format strings

            format_funcs = {}
            for title in self.titles:
                try:
                    format_str = FORMATS[title]
                except KeyError:
                    format_str = "%f"
                format_funcs[title] = partial(format_func, fmt=format_str)

            self._df.to_string(s_buf, index=False, formatters=format_funcs)

            s_buf.seek(0)

            with open(filename, 'w') as file:
                file.write(s_buf.getvalue())

        def readfile(self, filename):
            """
            Reads in an ascii columned file
            """
            self.parameters = parameters.parameters(filename=filename)

            i = 0
            # Check if this is a hdf file: magic number
            with open(filename, "rb") as f:
                magic = f.read(4)
            #              1 2 3 4 bytes
            if magic == b'\x89HDF':
                print("Reading your columnfile in hdf format")
                colfile_from_hdf(filename, obj=self)
                return

            # pandas read from CSV method
            # key-value stuff

            with open(filename, "r") as f:
                raw = f.readlines()
            header = True
            while header and i < len(raw):
                if len(raw[i].lstrip()) == 0:
                    # skip blank lines
                    i += 1
                    continue
                if raw[i][0] == "#":
                    # title line
                    if raw[i].find("=") > -1:
                        # key = value line
                        name, value = clean(raw[i][1:].split("=", 1))
                        self.parameters.addpar(
                            parameters.par(name, value))
                    else:
                        pass
                    i += 1
                else:
                    first_data_row = i
                    header = False

            # print(f"Skipping {first_data_row-1} header rows")

            column_titles = raw[first_data_row - 1].replace("# ", "").lstrip(" ").split()
            # print(f"Column titles are {column_titles}")

            self._df = pd.read_csv(filename, skiprows=range(first_data_row), sep=r'\s+', header=None, names=column_titles)

            self.parameters.dumbtypecheck()
            # self.set_attributes()

        def filter(self, mask):
            """
            mask is an nrows long array of true/false
            """
            if len(mask) != self.nrows:
                raise Exception("Mask is the wrong size")
            msk = np.array(mask, dtype=bool)
            # back to list here
            self._df.drop(self._df[~msk].index, inplace=True)

        def copy(self):
            """
            Returns a (deep) copy of the columnfile
            """
            cnw = PandasColumnfile(self.filename, new=True)
            cnw.parameters = parameters.parameters( **self.parameters.parameters.copy() )
            cnw._df = self._df.copy(deep=True)
            return cnw

        def copyrows(self, rows):
            """
            Returns a copy of select rows of the columnfile
            """

            cnw = PandasColumnfile(self.filename, new=True)
            cnw.parameters = parameters.parameters( **self.parameters.parameters.copy() )
            cnw._df = self._df.iloc[rows]

            return cnw

        def addcolumn(self, col, name):
            """
            Add a new column col to the object with name "name"
            Overwrites in the case that the column already exists
            """
            self._df[name] = col

        # Not obvious, but might be a useful alias
        setcolumn = addcolumn

        def getcolumn(self, name):
            """
            Gets data, if column exists
            """
            return self._df[name].to_numpy()


    class NewPandasColumnfile(PandasColumnfile):
        """ Just like a columnfile, but for creating new
        files """

        def __init__(self, titles):
            super().__init__(self, filename=None, new=True)
            self._df.columns = titles

except ImportError:
    def pderr():
        raise Exception("You do not have pandas installed!")


    class PandasColumnfile(columnfile):
        def __init__(self, filename=None, new=False):
            pderr()

    class NewPandasColumnfile(PandasColumnfile):
        def __init__(self, titles):
            pderr()


def bench():
    """
    Compares the timing for reading with columfile
    versus np.loadtxt
    """
    import sys, time
    start = time.time()
    import cProfile, pstats
    pr = cProfile.Profile()
    pr.enable()
    colf = columnfile(sys.argv[1])
    pr.disable()
    ps = pstats.Stats(pr, stream=sys.stdout )
    ps.sort_stats('tottime')
    ps.reverse_order()
    print(colf.bigarray.shape)
    print("ImageD11", time.time() - start)
    start = time.time()
    nolf = np.loadtxt(sys.argv[1])
    print(nolf.shape)
    print("np", time.time() - start)
    ps.print_stats()
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
