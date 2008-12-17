#!/usr/bin/python


## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py

#!/bliss/users/blissadm/python/bliss_python/suse82/bin/python


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
Script for making a pseudo-dark file as the minimum of a file series

Defines one class (minimum_image) which might perhaps be reused
"""


from fabio.openimage import openimage
from fabio import file_series
import numpy

class minimum_image:
    """
    Used for forming the minimum of a series of images
    """
    def __init__(self, filename = None, image = None):
        """
        file - the initial image to use a minimum
        image - a second image to make the min of the two args (??)
        """
        self.minimum_image = None
        if image is not None:
            # if an image is supplied we take that as the initial minimum
            self.minimum_image = image
        if filename != None:
            # if a filename is supplied we take the minimum of the image in that
            # file and the one stored
            self.add_file(filename)

    def add_file(self, filename):
        """
        Include another file
        """
        try:
            data_object = openimage(filename)
        except IOError:
            print filename, "not found"
            return
        picture = data_object.data
        if self.minimum_image is None:
            self.minimum_image = picture.copy()
        else:
            # Check dimensions match
            if self.minimum_image.shape == picture.shape:
                self.minimum_image = numpy.minimum(self.minimum_image, 
                                                   picture)
            else:
                raise Exception("Incompatible image dimensions")



def get_options(parser):
    """ add the command line options to parser """
    parser.add_option("-n", "--namestem", action = "store", 
                      type = "string", dest = "stem",
     help="Name of the files up the digits part, eg mydata in mydata0000.edf" )
    parser.add_option("-f", "--first", action = "store", type = "int", 
                      dest = "first",default = 0,
                      help = "Number of first file to process, default=0")
    parser.add_option("-l", "--last", action = "store", type = "int", 
                      dest = "last",default = 0, 
                      help = "Number of last file to process")
    parser.add_option("-o", "--outfile", action = "store", type = "string", 
                          dest = "outfile", default = "bkg.edf", 
                      help = "Output filename, default=bkg.edf")
    parser.add_option("-F", "--Format", action = "store", 
                      type = "string", dest = "format", default = ".edf",
                      help = "File format [edf|bruker]")
    parser.add_option("-s", "--step", action = "store", type = "int", 
                      dest = "step",default = 1,
                      help = "step - every nth image")
    parser.add_option("--ndigits", action = "store", type = "int",
                      dest = "ndigits", default = 4,
                      help = "Number of digits in file numbering [4]")
    return parser

def check_options( options ):
    """ Validate the command line options """
    for att in [ "stem", "first", "last"]:
        if getattr( options, att) is None:
            raise Exception("You must supply an option for "+att)


def bgmaker( options ):
    """ execute the command line script """
    # Generate list of files to proces
    if options.format in ['bruker', 'BRUKER', 'Bruker', 'GE']:
        extn = ""
    else:
        extn = options.format
    file_series_object = file_series.numbered_file_series(
        options.stem,
        options.first,
        options.last,
        extn,
        options.ndigits,
        step = options.step)
    minim = minimum_image(file_series_object[0])
    for filein in file_series_object:
        print filein
        minim.add_file(filein)
    # finally write out the answer
    # model header + data
    obj = openimage(file_series_object[0])
    obj.data = minim.minimum_image
    try:
        obj.write(options.outfile, force_type = obj.data.dtype)
    except TypeError:
        obj.write(options.outfile)
        
if __name__ == "__main__":

    # If we are running from a command line:
    import time
    # For benchmarking
    START = time.time()

    try:
        from optparse import OptionParser
        MYPARSER = OptionParser()
        MYPARSER = get_options( MYPARSER )
        OPTS , DUMMY = MYPARSER.parse_args()
        check_options( OPTS )
        bgmaker( OPTS )

    except:
        MYPARSER.print_help()
        raise

    END = time.time()

    print "Total time = %f /s" % (END - START)
