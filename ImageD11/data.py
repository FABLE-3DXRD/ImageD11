
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
ImageD11 data object

Try not to blown away, it is just a Numeric array with a
dictionary tagged on so that pseudodata (eg motor positions etc)
can follow around together. Throw whatever it is that you need
in the dictionary
"""

class data:
    """ Generic datatype for handling 2D images
        Just a wrapper for Numeric arrays now, so that information
        from file headers can be transported around in a dict
    """
    def __init__(self,array,header=None):
        """ Minimal things - dimensions and a numeric array """
        self.rows=array.shape[0]
        self.columns=array.shape[1]
        self.data=array
        if header==None:
            self.header={}
            self.header["rows"]=array.shape[0]
            self.header["columns"]=array.shape[1]
        else:
            self.header=header
    def __add__(self,other):
        if type(other) in [type(1),type(1.0)]:
            return data(self.header,self.data+other)
        if self.header["rows"]==other.header["rows"] and \
           self.header["columns"]==other.header["columns"]:
            return data(self.header,self.data+other.data)
        else:
            raise Exception("Dimensions not compatible")
    def __sub__(self,other):
        if type(other) in [type(1),type(1.0)]:
            return data(self.header,self.data-other)
        if self.header["rows"]==other.header["rows"] and \
           self.header["columns"]==other.header["columns"]:
            return data(self.header,self.data-other.data)
        else:
            raise Exception("Dimensions not compatible")
    def __mul__(self,other):
        if type(other) in [type(1),type(1.0)]:
            return data(self.header,self.data*other)
        if self.header["rows"]==other.header["rows"] and \
           self.header["columns"]==other.header["columns"]:
            return data(self.header,self.data*other.data)
        else:
            raise Exception("Dimensions not compatible")
    def __div__(self,other):
        if type(other) in [type(1),type(1.0)]:
            return data(self.header,self.data/other)
        if self.header["rows"]==other.header["rows"] and \
           self.header["columns"]==other.header["columns"]:
            return data(self.header,self.data/other.data)
        else:
            raise Exception("Dimensions not compatible")
