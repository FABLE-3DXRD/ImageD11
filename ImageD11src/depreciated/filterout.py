#!/usr/bin/env python

from __future__ import print_function

import sys, logging

from ImageD11 import refinegrains, indexing


def filtergrain(options):
    """
    Filter a peaks file according to which peaks are indexed
    """
    o = refinegrains.refinegrains(options.tol)
    o.loadparameters(options.parfile)
    o.readubis(options.ubifile)
    o.loadfiltered(options.fltfile)
    o.generate_grains()
    for gn in o.grainnames:
        o.reset_labels(options.fltfile)
        o.set_translation(gn, options.fltfile)
        o.compute_gv(gn, options.fltfile)
        drlv2 = indexing.calc_drlv2(o.grains[(gn,options.fltfile)].ubi,
                                    o.gv )
        o.scandata[options.fltfile].filter(drlv2 > options.tol*options.tol)
        print(gn, o.scandata[options.fltfile].nrows)
    o.scandata[options.fltfile].writefile(options.newfltfile)



if __name__=="__main__":

    console = logging.StreamHandler(sys.stdout)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(levelname)-8s : %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    console.setLevel(logging.DEBUG)
    root = logging.getLogger('')
    root.addHandler(console)
    root.setLevel(logging.DEBUG) # should we process everything...?

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-p",  "--parfile", action="store",
                      dest="parfile", type="string",
                      help="Name of parameter file")
    parser.add_option("-u",  "--ubifile", action="store",
                      dest="ubifile", type="string",
                      help="Name of ubi file")
    parser.add_option("-U", "--newubifile", action="store",
                      dest="newubifile", type="string",
                      help="Name of new ubi file")
    parser.add_option("-f",  "--fltfile", action="store",
                      dest="fltfile", type="string",
                      help="Name of flt file")
    parser.add_option("-F",  "--newfltfile", action="store",
                      dest="newfltfile", type="string",
                      help="Name of new flt file")
    parser.add_option("-t", "--tol", action="store",
                      dest="tol", type="float",
                      default = None,
                      help="Tolerance to use in peak assignment")
    parser.description = """
todo
    """

    options, args = parser.parse_args()

    if None in [options.parfile,
                options.ubifile,
                options.newfltfile,
                options.tol,
                options.fltfile]:
        parser.print_help()
        logging.error("You have not filled in all the required options")
        import sys
        sys.exit()

    try:
        filtergrain(options)
    except:
        parser.print_help()
        import traceback
        logging.error("An error occurred, details follow")
        traceback.print_exc()



