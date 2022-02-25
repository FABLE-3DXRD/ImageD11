#!/usr/bin/env python
from __future__ import print_function



## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py



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

import fabio
import fabio.edfimage
from fabio.openimage import openimage
import numpy
import random # to do images in random order
import logging
from ImageD11 import ImageD11options

class minimum_image(object):
    """
    Used for forming the minimum of a series of images
    """
    def __init__(self, filename = None, image = None):
        """
        file - the initial image to use a minimum
        image - a second image to make the min of the two args (??)
        """
        self.bkg = None
        if image is not None:
            # if an image is supplied we take that as the initial minimum
            self.bkg = image
        if filename is not None:
            # if a filename is supplied we take the minimum of the
            # image in that
            # file and the one stored
            self.add_file(filename)

    def add_image(self, picture):
        """
        """
        if self.bkg is None:
            self.bkg = picture.copy()
        else:
            # Check dimensions match
            if self.bkg.shape == picture.shape:
                numpy.minimum(self.bkg, picture, self.bkg)
            else:
                raise Exception("Incompatible image dimensions")

    def add_file(self, filename):
        """
        Include another file
        """
        try:
            data_object = openimage(filename)
        except IOError:
            print(filename, "not found")
            return
        self.add_image(data_object.data)

class kbg(object):
    """
    Kalman style background filtering
    """
    def __init__(self, x0, p0, Q=0.01):
        """
        x0 = initial background estimate
        p  = current noise estimate
        p0 = initial noise estimate
        Q  = drift rate [fixed here]
        """
        self.bkg = x0.astype(numpy.float32)
        self.p = p0
        self.p0= numpy.float32(p0)
        self.Q = Q

    def add_image(self, z):
        """
        Update the current estimates of background and noise
        """
        assert z.shape == self.bkg.shape, "Incompatible dimensions"
        self.p = self.p + self.Q
        # R = noise estimate for this image
        #  ...either read noise or obs-calc for downweighting peaks
        # Misfit of the current image
        err = z - self.bkg
        R = numpy.where( err > self.p0, err, self.p0)
        # Usual Kalman updates
        K = self.p/(self.p + R)
        self.bkg = self.bkg + K * err
        self.p = ( 1 - K ) * self.p

    def add_file(self, filename):
        """
        Include another file
        """
        try:
            data_object = openimage(filename)
        except IOError:
            print(filename, "not found")
            return
        self.add_image(data_object.data)


def get_options(parser):
    """ add the command line options to parser """
    parser.add_argument("-n", "--namestem", action = "store", 
                      type = str, dest = "stem",
     help="Name of the files up the digits part, eg mydata in mydata0000.edf" )
    parser.add_argument("-f", "--first", action = "store", type = int, 
                      dest = "first",default = 0,
                      help = "Number of first file to process, default=0")
    parser.add_argument("-l", "--last", action = "store", type = int, 
                      dest = "last",default = 0, 
                      help = "Number of last file to process")
    parser.add_argument("-o", "--outfile", action = "store", 
                        type = ImageD11options.ImageFileType(mode='w'), 
                          dest = "outfile", default = "bkg.edf", 
                      help = "Output filename, default=bkg.edf")
    parser.add_argument("-F", "--Format", action = "store", 
                      type = str, dest = "format", default = ".edf",
                      help = "File format [edf|bruker]")
    parser.add_argument("-s", "--step", action = "store", type = int,
                      dest = "step",default = 1,
                      help = "step - every nth image")
    parser.add_argument("--ndigits", action = "store", type = int,
                      dest = "ndigits", default = 4,
                      help = "Number of digits in file numbering [4]")
    parser.add_argument("-k", "--kalman-error", action="store", type = float,
            dest = "kalman_error", default = 0,
            help = "Error value to use Kalman style filter (read noise)" )

    return parser

def check_options( options ):
    """ Validate the command line options """
    for att in [ "stem", "first", "last"]:
        if getattr( options, att) is None:
            logging.error("You must supply an option for %s",att)
            raise Exception("Missing "+att)


def bgmaker( options ):
    """ execute the command line script """
    # Generate list of files to proces

    if options.format in ['bruker', 'BRUKER', 'Bruker', 'GE']:
        extn = ""
    else:
        extn = options.format

    first_image_name = fabio.filename_object(
        options.stem,  
        num = options.first,
        extension = extn,
        digits = options.ndigits)


    first_image = openimage( first_image_name )
    print(first_image.filename)

    allimagenumbers = list(range(options.first,
                                 options.last + 1 - options.step,
                                 options.step))

    if options.kalman_error <= 0:
        print("Using minimum image algorithm")
        bko = minimum_image( image = first_image.data )
    else:
        print("Using Kalman algorithm with error =",options.kalman_error)
        bko = kbg( first_image.data, options.kalman_error*options.kalman_error )
        print("Taking images in random order")
        random.seed(42) # reproducible
        random.shuffle( allimagenumbers )

    for current_num in allimagenumbers:
        try:
            im = first_image.getframe( current_num )
            print(im.filename)
            bko.add_image( im.data )
        except KeyboardInterrupt:
            print("Got a keyboard interrupt")
            break
        except:
            import traceback
            traceback.print_exc()
            print("Failed for",current_num)
        
    # finally write out the answer
    # model header + data


    # write as edf - we should actually have a way to flag
    # which fabioimage formats know how to write themselves
    if options.outfile[-3:] == "edf":
        print("writing",options.outfile,"in edf format")
        im = fabio.edfimage.edfimage( data = bko.bkg )
    else:
        im = first_image
        
    im.data = bko.bkg
    try:
        im.write(options.outfile, force_type = im.data.dtype)
    except TypeError: # WTF?
        im.write(options.outfile)
    except:
        print("problem writing")
        print("trying to write",options.outfile,"in edf format")
        im = fabio.edfimage.edfimage( data = bko.minimum_image )
        try:
            im.write(options.outfile, force_type = im.data.dtype)
        except TypeError:
            im.write(options.outfile)
        
if __name__ == "__main__":

    # If we are running from a command line:
    import time
    # For benchmarking
    START = time.time()

    try:
        from argparse import ArgumentParser
        MYPARSER = ArgumentParser()
        MYPARSER = get_options( MYPARSER )
        OPTS = MYPARSER.parse_args()
        check_options( OPTS )
        bgmaker( OPTS )

    except:
        MYPARSER.print_help()
        raise

    END = time.time()

    print("Total time = %f /s" % (END - START))
