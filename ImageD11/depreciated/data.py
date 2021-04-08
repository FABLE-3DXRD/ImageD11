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

TODO : rebinning - is this generally useful??
TODO : PILimage methods ??
TODO : ROI integration - need to define the slice mapping better to
     : be compatible with fabian

Added some methods from:
fabian:edfimage.py by:
Authors: Henning O. Sorensen & Erik Knudsen
         Center for Fundamental Research: Metal Structures in Four Dimensions
         Risoe National Laboratory
         Frederiksborgvej 399
         DK-4000 Roskilde
         email:erik.knudsen@risoe.dk
"""

import numpy as np
import logging

class data: #IGNORE:R0902
    """ Generic datatype for handling 2D images
        Just a wrapper for Numeric arrays now, so that information
        from file headers can be transported around in a dict
    """
    def __init__(self, array, header=None):
        """ Minimal things - dimensions and a numeric array """
        self.rows = array.shape[0]
        self.columns = array.shape[1]
        self.data = array
        self.m = self.maxval = self.stddev = self.minval = None
        if header == None:
            self.header = {}
            self.header["rows"] = array.shape[0]
            self.header["columns"] = array.shape[1]
        else:
            self.header = header
        self.header_keys = self.header.keys()

    def getheader(self):
        """ return the header"""
        return self.header

    def getmax(self):
        """ Return data maximum value (type matches data) """
        if self.maxval == None:
            self.maxval = np.maximum.reduce(np.ravel(self.data))
        return self.maxval

    def getmin(self):
        """ Return data minimum value (type matches data) """
        if self.minval == None:
            self.minval = np.minimum.reduce(np.ravel(self.data))
        return self.minval

    def getmean(self):
        """ return mean """
        if self.m == None:
            d = np.ravel(self.data.astype(float))
            self.m = np.sum(d) / len(d)
        return float(self.m)

    def getstddev(self):
        """ return standard deviation of image """
        if self.m == None:
            self.getmean()
            logging.debug("recalc mean")
        if self.stddev == None:
            #  wikipedia method:
            # "In other words, the standard deviation of a discrete
            # uniform random variable X can be calculated as follows:
            #
            #   1. For each value xi calculate the difference
            # x_i - <x> between xi and the average value <x>.
            d = np.ravel(self.data.astype(float)) - self.m
            #   2. Calculate the squares of these differences.
            S = d*d
            #   3. Find the average of the squared differences.
            # This quantity is the variance sigma2.
            N = len(S)
            S2 = np.sum(S)/N
            #  4. Take the square root of the variance.
            import math # ensure not an array here
            self.stddev = math.sqrt(S2)
        return float(self.stddev)

    def resetvals(self):
        """ resets properties in the event of data changing """
        self.m = self.stddev = self.maxval = self.minval = None


    def __add__(self, other):
        if type(other) in [type(1), type(1.0)]:
            return data(self.header, self.data + other)
        if self.header["rows"] == other.header["rows"] and \
           self.header["columns"] == other.header["columns"]:
            return data(self.header , self.data + other.data)
        else:
            raise Exception("Dimensions not compatible")
    def __sub__(self, other):
        if type(other) in [type(1), type(1.0)]:
            return data(self.header, self.data - other)
        if self.header["rows"] == other.header["rows"] and \
           self.header["columns"] == other.header["columns"]:
            return data(self.header, self.data - other.data)
        else:
            raise Exception("Dimensions not compatible")
    def __mul__(self, other):
        if type(other) in [type(1), type(1.0)]:
            return data(self.header, self.data * other)
        if self.header["rows"] == other.header["rows"] and \
           self.header["columns"] == other.header["columns"]:
            return data(self.header, self.data * other.data)
        else:
            raise Exception("Dimensions not compatible")
    def __div__(self, other):
        if type(other) in [type(1), type(1.0)]:
            return data(self.header, self.data / other)
        if self.header["rows"] == other.header["rows"] and \
           self.header["columns"] == other.header["columns"]:
            return data(self.header, self.data / other.data)
        else:
            raise Exception("Dimensions not compatible")
