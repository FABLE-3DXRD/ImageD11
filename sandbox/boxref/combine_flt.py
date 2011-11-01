

# ImageD11 Software for beamline ID11
# Copyright (C) 2011  Jon Wright
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



import numpy, columnfile, sys

class combine_flt(object):
    def __init__(self, fltfile, samtx=0, samty=0, samtz=0):
        """
        First flt file (defines columns, should be same for all)
        Translations to be supplied if not zero
        """
        print "reading", fltfile
        sys.stdout.flush()
        self.fltfile = columnfile.columnfile( fltfile )
        onearray = numpy.ones(self.fltfile.nrows, numpy.float32)
        self.fltfile.addcolumn( onearray*samtx, "samtx" )
        self.fltfile.addcolumn( onearray*samty, "samty" )
        self.fltfile.addcolumn( onearray*samtz, "samtz" )
        # Convert to list for easier appending later
        self.datalist = [ list(col) for col in self.fltfile.bigarray ]
        print len(self.datalist),len(self.datalist[0])
    def addfltfile(self, fltfile, samtx=0, samty=0, samtz=0):
        """
        Add another flt file, probably with different translations
        """
        c = columnfile.columnfile(fltfile)
        onearray = numpy.ones(c.nrows, numpy.float32)
        c.addcolumn( onearray*samtx, "samtx" )
        c.addcolumn( onearray*samty, "samty" )
        c.addcolumn( onearray*samtz, "samtz" )
        # Check columnfiles have equivalent titles
        for i,j in zip(c.titles, self.fltfile.titles):
            assert i == j
        for i in range(len(self.datalist)):
            self.datalist[i] += list(c.bigarray[i])
        print len(c.bigarray[0]), len(self.datalist[0])
    def write_hdf(self, name, datagroup):
        """ Save the output, the hdf code is in columnfile.py """
        self.datalist = numpy.array(self.datalist)
        self.fltfile.bigarray = self.datalist
        self.fltfile.ncols, self.fltfile.nrows = self.fltfile.bigarray.shape
        self.fltfile.set_attributes()
        columnfile.colfile_to_hdf(self.fltfile, name, datagroup)



def combine_flt_from_file( fname, hdfname, hdfgroup ):
    """ Read one flt per line with x,y,z translations """
    o = None
    for line in open(fname).readlines():
        name = line.split()[0]
        x,y,z = [float(v) for v in line.split()[1:]]
        print name, x, y, z
        if o is None:
            o = combine_flt( name, x, y, z )
        else:
            o.addfltfile( name, x, y, z)
    o.write_hdf( hdfname, hdfgroup )


if __name__=="__main__":
    try:
        xyzfile = sys.argv[1]
        hdfile =  sys.argv[2]
        hdfgroup = sys.argv[3]
    except:
        print "Usage: %s xyzfile hdffile hdfgroup"%( sys.argv[0] )
        print """ ... where:
   xyzfile = ascii file with lines containing:
        fltfilename   samtx   samty   samtz
   hdffile = output hdffile with combined data in it
   hdfgroup = group name (directory) in the hdf file """
 
    combine_file_from_file( sys.argv[1], sys.argv[2], sys.argv[3] )
