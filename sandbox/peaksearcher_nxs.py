
from __future__ import print_function



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



"""

Adaptation of ImageD11src/peaksearcher.py by James Ball, Mar 2020


Script for peaksearching images from a supplied nexus file

Uses the connectedpixels extension for finding blobs above a threshold
and the blobcorrector(+splines) for correcting them for spatial distortion

Defines one function (peaksearch) which might be reused
"""

# For benchmarking
import time
reallystart = time.time()

from six.moves import queue
# import threading
import sys , os.path
import numpy

# Generic file format opener from fabio
from fabio.openimage import openimage

from ImageD11 import blobcorrector, ImageD11options
from ImageD11.correct import correct
from ImageD11.labelimage import labelimage
from ImageD11 import ImageD11_thread
ImageD11_thread.stop_now = False


class timer:
    def __init__(self):
        self.start = time.time()
        self.now = self.start
        self.msgs = []
    def msg(self ,msg):
        self.msgs.append(msg)
    def tick(self ,msg=""):
        now = time.time()
        self.msgs.append("%s %.2f/s " %(msg ,now-self.now))
        self.now = now
    def tock(self ,msg=""):
        self.tick(msg)
        print(" ".join(self.msgs) ,"%.2f/s "% (self.now-self.start))
        sys.stdout.flush()


def peaksearch( filename ,
                data_object ,
                corrector ,
                thresholds ,
                labims ):
    """
    filename  : The name of the image file for progress info
    data_object : Fabio object containing data and header
    corrector : spatial and dark/flood , linearity etc

    thresholds : [ float[1], float[2] etc]

    labims : label image objects, one for each threshold

    """
    t = timer()
    picture = data_object.data.astype(numpy.float32)

    assert "Omega" in data_object.header, "Bug in peaksearch headers"

    for lio in list(labims.values()):
        f = lio.sptfile
        f.write("\n\n# File %s\n" % (filename))
        f.write("# Frame %d\n" % (data_object.currentframe) )
        f.write("# Processed on %s\n" % (time.asctime()))
        try:
            f.write("# Spatial correction from %s\n" % (corrector.splinefile))
            f.write("# SPLINE X-PIXEL-SIZE %s\n" % (str(corrector.xsize)))
            f.write("# SPLINE Y-PIXEL-SIZE %s\n" % (str(corrector.ysize)))
        except:
            pass
        for item in list(data_object.header.keys()):
            if item == "headerstring": # skip
                continue
            try:
                f.write("# %s = %s\n" % (item,
                                         str(data_object.header[item]).replace("\n" ," ")))
            except KeyError:
                pass

    # Get the rotation angle for this image
    ome = float(data_object.header["Omega"])
    # print "Reading from header"
    #
    # Now peaksearch at each threshold level
    t.tick(filename)
    for threshold in thresholds:
        labelim = labims[threshold]
        f = labelim.sptfile
        if labelim.shape != picture.shape:
            raise "Incompatible blobimage buffer for file %s" %(filename)
        #
        #
        # Do the peaksearch
        f.write("# Omega = %f\n" % (ome))
        labelim.peaksearch(picture, threshold, ome)
        f.write("# Threshold = %f\n" % (threshold))
        f.write("# npks = %d\n" % (labelim.npk))
        #
        if labelim.npk > 0:
            labelim.output2dpeaks(f)
        labelim.mergelast()
        t.msg("T=%-5d n=%-5d;" % (int(threshold), labelim.npk))
        # Close the output file
    # Finish progress indicator for this file
    t.tock()
    sys.stdout.flush()
    return None


def peaksearch_driver(options, args):
    """
    To be called with options from command line
    """
    ################## debugging still
    for a in args:
        print("arg: " + str(a) + "," + str(type(a)))
    for o in list(options.__dict__.keys()):  # FIXME
        if getattr(options, o) in ["False", "FALSE", "false"]:
            setattr(options, o, False)
        if getattr(options, o) in ["True", "TRUE", "true"]:
            setattr(options, o, True)
        print("option:", str(o), str(getattr(options, o)), ",", \
              str(type(getattr(options, o))))
    ###################
    print("This peaksearcher is from", __file__)

    if options.killfile is not None and os.path.exists(options.killfile):
        print("The purpose of the killfile option is to create that file")
        print("only when you want peaksearcher to stop")
        print("If the file already exists when you run peaksearcher it is")
        print("never going to get started")
        raise ValueError("Your killfile " + options.killfile + " already exists")

    if options.thresholds is None:
        raise ValueError("No thresholds supplied [-t 1234]")

    if len(args) == 0 and options.nexusfile is None:
        raise ValueError("No files to process")

        # What to do about spatial

    if options.perfect == "N" and os.path.exists(options.spline):
        print("Using spatial from", options.spline)
        corrfunc = blobcorrector.correctorclass(options.spline)
    else:
        print("Avoiding spatial correction")
        corrfunc = blobcorrector.perfect()

    # Get list of filenames to process
    # if len(args) > 0 :
    #     # We no longer assume unlabelled arguments are filenames
    #     file_series_object = file_series.file_series(args)

    # This is always the case now
    corrfunc.orientation = "edf"

    import h5py

    # Read list of files and list of motor positions from Nexus file:
    nexus_path = options.nexusfile
    nexus_file = h5py.File(nexus_path, "r")
    group = nexus_file[options.group_path]
    omega_dset = group.get(options.omega_dset)
    image_dset = group.get(options.image_dset)
    omega_list = [x for x in omega_dset[..., :]]
    image_list = [x.decode("utf-8") for x in image_dset[..., :]]

    import fabio

    # Output files:

    import fabio.file_series
    # Use traceback = True for debugging
    first_image = openimage(image_list[0])

    file_series_object = fabio.file_series.new_file_series(
        first_image,
        nimages=len(image_list),
        traceback=True)

    if options.outfile[-4:] != ".spt":
        options.outfile = options.outfile + ".spt"
        print("Your output file must end with .spt, changing to ", options.outfile)

    # Make a blobimage the same size as the first image to process

    # List comprehension - convert remaining args to floats
    # must be unique list so go via a set
    thresholds_list = list(set([float(t) for t in options.thresholds]))
    thresholds_list.sort()

    li_objs = {}  # label image objects, dict of

    s = first_image.data.shape  # data array shape

    # Create label images
    for t in thresholds_list:
        # the last 4 chars are guaranteed to be .spt above
        mergefile = "%s_t%d.flt" % (options.outfile[:-4], t)
        spotfile = "%s_t%d.spt" % (options.outfile[:-4], t)
        li_objs[t] = labelimage(shape=s,
                                fileout=mergefile,
                                spatial=corrfunc,
                                sptfile=spotfile)
        print("make labelimage", mergefile, spotfile)
    if options.dark is not None:
        print("Using dark (background)", options.dark)
        darkimage = openimage(options.dark).data.astype(numpy.float32)
    else:
        darkimage = None
    if options.darkoffset != 0:
        print("Adding darkoffset", options.darkoffset)
        if darkimage is None:
            darkimage = options.darkoffset
        else:
            darkimage += options.darkoffset
    if options.flood is not None:
        floodimage = openimage(options.flood).data
        cen0 = int(floodimage.shape[0] / 6)
        cen1 = int(floodimage.shape[0] / 6)
        middle = floodimage[cen0:-cen0, cen1:-cen1]
        nmid = middle.shape[0] * middle.shape[1]
        floodavg = numpy.mean(middle)
        print("Using flood", options.flood, "average value", floodavg)
        if floodavg < 0.7 or floodavg > 1.3:
            print("Your flood image does not seem to be normalised!!!")

    else:
        floodimage = None
    start = time.time()
    print("File being treated in -> out, elapsed time")
    # If we want to do read-ahead threading we fill up a Queue object with data
    # objects
    # THERE MUST BE ONLY ONE peaksearching thread for 3D merging to work
    # there could be several read_and_correct threads, but they'll have to get the order right,
    # for now only one
    if options.oneThread:
        # Wrap in a function to allow profiling (perhaps? what about globals??)
        def go_for_it(file_series_object, darkimage, floodimage,
                      corrfunc, thresholds_list, li_objs):
            for inc, data_object in enumerate(file_series_object):
                t = timer()
                if not isinstance(data_object, fabio.fabioimage.fabioimage):
                    # Is usually an IOError
                    if isinstance(data_object[1], IOError):
                        sys.stdout.write(data_object[1].strerror + '\n')
                        # data_object[1].filename
                    else:
                        import traceback
                        traceback.print_exception(data_object[0], data_object[1], data_object[2])
                        sys.exit()
                    continue
                filein = data_object.filename

                data_object.header["Omega"] = float(omega_list[inc])

                data_object = correct(data_object, darkimage, floodimage,
                                      do_median=options.median,
                                      monitorval=options.monitorval,
                                      monitorcol=options.monitorcol,
                                      )
                t.tick(filein + " io/cor")
                peaksearch(filein, data_object, corrfunc,
                           thresholds_list, li_objs)
            for t in thresholds_list:
                li_objs[t].finalise()

        go_for_it(file_series_object, darkimage, floodimage,
                  corrfunc, thresholds_list, li_objs)


    else:
        print("Going to use threaded version!")
        try:
            # TODO move this to a module ?

            class read_only(ImageD11_thread.ImageD11_thread):
                def __init__(self, queue, file_series_obj, myname="read_only"):
                    """ Reads files in file_series_obj, writes to queue """
                    self.queue = queue
                    self.file_series_obj = file_series_obj
                    ImageD11_thread.ImageD11_thread.__init__(self,
                                                             myname=myname)
                    print("Reading thread initialised", end=' ')

                def ImageD11_run(self):
                    """ Read images and copy them to self.queue """
                    for inc, data_object in enumerate(self.file_series_obj):
                        if self.ImageD11_stop_now():
                            print("Reader thread stopping")
                            break
                        if not isinstance(data_object, fabio.fabioimage.fabioimage):
                            # Is usually an IOError
                            if isinstance(data_object[1], IOError):
                                #                                print data_object
                                #                                print data_object[1]
                                sys.stdout.write(str(data_object[1].strerror) + '\n')
                            #                                  ': ' + data_object[1].filename + '\n')
                            else:
                                import traceback
                                traceback.print_exception(data_object[0], data_object[1], data_object[2])
                                sys.exit()
                            continue
                        ti = timer()
                        filein = data_object.filename + "[%d]" % (data_object.currentframe)
                        try:
                            data_object.header["Omega"] = float(omega_list[inc])
                        except KeyboardInterrupt:
                            raise
                        except:
                            continue
                        ti.tick(filein)
                        self.queue.put((filein, data_object), block=True)
                        ti.tock(" enqueue ")
                        if self.ImageD11_stop_now():
                            print("Reader thread stopping")
                            break

                    # Flag the end of the series
                    self.queue.put((None, None), block=True)

            class correct_one_to_many(ImageD11_thread.ImageD11_thread):
                def __init__(self, queue_read, queues_out, thresholds_list,
                             dark=None, flood=None, myname="correct_one",
                             monitorcol=None, monitorval=None,
                             do_median=False):
                    """ Using a single reading queue retains a global ordering
                    corrects and copies images to output queues doing
                    correction once """
                    self.queue_read = queue_read
                    self.queues_out = queues_out
                    self.dark = dark
                    self.flood = flood
                    self.do_median = do_median
                    self.monitorcol = monitorcol
                    self.monitorval = monitorval
                    self.thresholds_list = thresholds_list
                    ImageD11_thread.ImageD11_thread.__init__(self,
                                                             myname=myname)

                def ImageD11_run(self):
                    while not self.ImageD11_stop_now():
                        ti = timer()
                        filein, data_object = self.queue_read.get(block=True)
                        if filein is None:
                            for t in self.thresholds_list:
                                self.queues_out[t].put((None, None),
                                                       block=True)
                            # exit the while 1
                            break
                        data_object = correct(data_object, self.dark,
                                              self.flood,
                                              do_median=self.do_median,
                                              monitorval=self.monitorval,
                                              monitorcol=self.monitorcol,
                                              )
                        ti.tick(filein + " correct ")
                        for t in self.thresholds_list:
                            # Hope that data object is read only
                            self.queues_out[t].put((filein, data_object),
                                                   block=True)
                        ti.tock(" enqueue ")
                    print("Corrector thread stopping")

            class peaksearch_one(ImageD11_thread.ImageD11_thread):
                def __init__(self, q, corrfunc, threshold, li_obj,
                             myname="peaksearch_one"):
                    """ This will handle a single threshold and labelimage
                    object """
                    self.q = q
                    self.corrfunc = corrfunc
                    self.threshold = threshold
                    self.li_obj = li_obj
                    ImageD11_thread.ImageD11_thread.__init__(
                        self,
                        myname=myname + "_" + str(threshold))

                def run(self):
                    while not self.ImageD11_stop_now():
                        filein, data_object = self.q.get(block=True)
                        if not isinstance(data_object, fabio.fabioimage.fabioimage):
                            break
                        peaksearch(filein, data_object, self.corrfunc,
                                   [self.threshold],
                                   {self.threshold: self.li_obj})
                    self.li_obj.finalise()

            # 8 MB images - max 40 MB in this queue
            read_queue = queue.Queue(5)
            reader = read_only(read_queue, file_series_object)
            reader.start()
            queues = {}
            searchers = {}
            for t in thresholds_list:
                print("make queue and peaksearch for threshold", t)
                queues[t] = queue.Queue(3)
                searchers[t] = peaksearch_one(queues[t], corrfunc,
                                              t, li_objs[t])
            corrector = correct_one_to_many(read_queue,
                                            queues,
                                            thresholds_list,
                                            dark=darkimage,
                                            flood=floodimage,
                                            do_median=options.median,
                                            monitorcol=options.monitorcol,
                                            monitorval=options.monitorval)
            corrector.start()
            my_threads = [reader, corrector]
            for t in thresholds_list[::-1]:
                searchers[t].start()
                my_threads.append(searchers[t])
            nalive = len(my_threads)

            def empty_queue(q):
                while 1:
                    try:
                        q.get(block=False, timeout=1)
                    except:
                        break
                q.put((None, None), block=False)

            while nalive > 0:
                try:
                    nalive = 0
                    for thr in my_threads:
                        if thr.isAlive():
                            nalive += 1
                    if options.killfile is not None and \
                            os.path.exists(options.killfile):
                        raise KeyboardInterrupt()
                    time.sleep(1)
                except KeyboardInterrupt:
                    print("Got keyboard interrupt in waiting loop")
                    ImageD11_thread.stop_now = True
                    try:
                        time.sleep(1)
                    except:
                        pass
                    empty_queue(read_queue)
                    for t in thresholds_list:
                        q = queues[t]
                        empty_queue(q)
                    for thr in my_threads:
                        if thr.isAlive():
                            thr.join(timeout=1)
                    print("finishing from waiting loop")
                except:
                    print("Caught exception in waiting loop")
                    ImageD11_thread.stop_now = True
                    time.sleep(1)
                    empty_queue(read_queue)
                    for t in thresholds_list:
                        q = queues[t]
                        empty_queue(q)
                    for thr in my_threads:
                        if thr.isAlive():
                            thr.join(timeout=1)
                    raise


        except ImportError:
            print("Probably no threading module present")
            raise


def get_options(parser):
    """ Add our options to a parser object """
    parser.add_argument("-n", "--nexusfile", action="store",
                        dest="nexusfile", type=str, default=None,
                        help="The Nexus file path")
    parser.add_argument("-g", "--group_path", action="store",
                        dest="group_path", type=str, default=None,
                        help="Internal NXS path to datasets group")
    parser.add_argument("--image_dset", action="store",
                        dest="image_dset", type=str, default=None,
                        help="Name of image files dataset")
    parser.add_argument("--omega_dset", action="store",
                        dest="omega_dset", type=str, default=None,
                        help="Name of omegas dataset")
    parser.add_argument("-o", "--outfile", action="store",
                        dest="outfile", default="peaks.spt", type=str,
                        help="Output filename, default=peaks.spt")
    parser.add_argument("-d", "--darkfile", action="store",
                        dest="dark", default=None, type=ImageD11options.ImageFileType(mode='r'),
                        help="Dark current filename, to be subtracted, default=None")
    dod = 0
    parser.add_argument("-D", "--darkfileoffset", action="store",
                        dest="darkoffset", default=dod, type=float,
                        help=
                        "Constant to subtract from dark to avoid overflows, default=%d" % (dod))
    parser.add_argument("-s", "--splinefile", action="store",
                        dest="spline", default=None, type=ImageD11options.SplineFileType(mode='r'),
                        help="Spline file for spatial distortion, default=None")
    parser.add_argument("-p", "--perfect_images", action="store",
                        choices=["Y", "N"], default="Y", dest="perfect",
                        help="Ignore spline Y|N, default=N")
    parser.add_argument("-O", "--flood", action="store",
                        type=ImageD11options.ImageFileType(mode='r'),
                        default=None, dest="flood",
                        help="Flood file, default=None")
    parser.add_argument("-t", "--threshold", action="append", type=float,
                        dest="thresholds", default=None,
                        help="Threshold level, you can have several")
    parser.add_argument("--singleThread", action="store_true",
                        dest="oneThread", default=False,
                        help="Do single threaded processing")
    parser.add_argument("-k", "--killfile", action="store",
                        dest="killfile", default=None,
                        type=ImageD11options.FileType(),
                        help="Name of file to create stop the peaksearcher running")
    parser.add_argument("-m", "--median1D", action="store_true",
                        default=False, dest="median",
                        help="Computes the 1D median, writes it to file .bkm and" \
                             + " subtracts it from image. For liquid background" \
                             + " on radially transformed images")
    parser.add_argument("--monitorcol", action="store", type=str,
                        dest="monitorcol",
                        default=None,
                        help="Header value for incident beam intensity")
    parser.add_argument("--monitorval", action="store", type=float,
                        dest="monitorval",
                        default=None,
                        help="Incident beam intensity value to normalise to")
    parser.add_argument("--interlaced", action="store_true",
                        dest="interlaced", default=False,
                        help="Interlaced DCT scan")
    parser.add_argument("--iflip", action="store_true",
                        dest="iflip", default=False,
                        help="Reverse second half of interlaced scan")
    return parser


def get_help(usage=True):
    """ return the help string for online help """
    try:
        import StringIO as io
    except:
        # python3
        import io
    import argparse
    if usage:
        o = get_options(argparse.ArgumentParser())
    else:
        o = get_options(argparse.ArgumentParser(usage=argparse.SUPPRESS))
    f = io.StringIO()
    o.print_help(f)
    return f.getvalue()


if __name__ == "__main__":
    raise Exception("Please use the driver script peaksearch.py")


