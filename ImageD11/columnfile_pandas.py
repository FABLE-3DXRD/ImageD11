from __future__ import print_function

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
import io
from functools import partial

from ImageD11 import parameters, transform
from ImageD11 import columnfile
import numpy as np

import pandas as pd


class PandasColumnfile(columnfile.columnfile):
    def __init__(self, filename=None, new=False):
        self.filename = filename
        if filename is not None:
            self.parameters = parameters.parameters(filename=filename)
        else:
            self.parameters = parameters.parameters()
        self._df = pd.DataFrame()
        if not new:
            self.readfile(filename)

    def __getattribute__(self, item):
        # called whenever getattr(item, 'attr') or item.attr is called
        # for speed, look in self._df FIRST before giving up
        try:
            # we can't use self.getcolumn() here without retriggering self.__getattribute__()
            # so we have to directly get the column from the dataframe
            return super(PandasColumnfile, self).__getattribute__('_df')[item].to_numpy()
        except KeyError:
            # not in the dataframe, so just call getattribute normally (e.g for cf.nrows)
            return super(PandasColumnfile, self).__getattribute__(item)
        except AttributeError as e:
            # didn't find it
            raise e

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
                format_str = columnfile.FORMATS[title]
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
            columnfile.colfile_from_hdf(filename, obj=self)
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
                    name, value = columnfile.clean(raw[i][1:].split("=", 1))
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

        self._df = pd.read_csv(filename, skiprows=range(first_data_row), sep='\s+', header=None, names=column_titles)

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
        cnw.parameters = parameters.parameters(**self.parameters.parameters)
        cnw._df = self._df.copy(deep=True)
        return cnw

    def copyrows(self, rows):
        """
        Returns a copy of select rows of the columnfile
        """

        cnw = PandasColumnfile(self.filename, new=True)
        cnw.parameters = self.parameters
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
    colf = PandasColumnfile(sys.argv[1])
    pr.disable()
    ps = pstats.Stats(pr, stream=sys.stdout)
    ps.sort_stats('tottime')
    ps.reverse_order()
    # print(colf.bigarray.shape)
    print("ImageD11", time.time() - start)
    start = time.time()
    nolf = np.loadtxt(sys.argv[1])
    print(nolf.shape)
    print("np", time.time() - start)
    ps.print_stats()
    # os.system("time -p ./a.out")


if __name__ == "__main__":
    bench()
