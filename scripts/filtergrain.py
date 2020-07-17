#!/usr/bin/env python

from __future__ import print_function

"""
Utility script to pick out the peaks belonging to a certain
grain

Used to stabilise refinements
"""

import sys, logging

import numpy as np
from argparse import ArgumentParser

from ImageD11 import refinegrains, indexing, grain, ImageD11options



def filtergrain(options):
    """
    Filter a peaks file according to which peaks are indexed
    """
    o = refinegrains.refinegrains(tolerance=0.9,
                                  OmFloat=options.omega_float,
                                  OmSlop=options.omega_slop)
    o.loadparameters(options.parfile)
    o.readubis(options.ubifile)

    if options.grain is None:
        if len(o.grainnames) == 1:
            gn = o.grainnames[0]
        else:
            for i,gn in zip(list(range(len(o.grainnames))),o.grainnames):
                logging.info("Choose %d for grain %s"%(i, gn))
            gn = o.grainnames[int(input("select which grain "))]
    else:
        gn = o.grainnames[int(options.grain)]
    o.grainnames = [ gn, ]

    o.loadfiltered(options.fltfile)
    o.generate_grains()
    o.assignlabels()
    o.set_translation( gn, options.fltfile )
    o.compute_gv(o.grains[ (gn, options.fltfile) ] )
    if options.tol is None:
        for tol in [0.01,0.025,0.05,0.1,0.15,0.25,0.5]:
            o.tolerance = tol
            logging.info("tol %f"%(o.tolerance))
            o.refineubis(quiet=False, scoreonly=True)
        o.tolerance = float(input("Enter tolerance "))
        options.tol = o.tolerance
    else:
        o.tolerance = options.tol
    matrix = o.grains[(gn,options.fltfile)].ubi

    o.assignlabels()
    drlv2 = indexing.calc_drlv2( matrix, o.gv )
    logging.info("Total peaks before filtering %d"%
                     o.scandata[options.fltfile].nrows)
    gotpks = o.scandata[options.fltfile].copy()
    gotpks.filter(gotpks.labels == gn)
    gotpks.writefile(options.newfltfile)
    logging.info("Peaks which were indexed %d written to %s"%(
                gotpks.nrows, options.newfltfile))
    # don't bother to copy here as we can overwrite
    if options.notindexed is not None:
        notpks = o.scandata[options.fltfile].copy()
        notpks.addcolumn(o.tth, "tth")
        notpks.addcolumn(o.eta, "eta")
        notpks.addcolumn(o.gv[:,0] , "gx")
        notpks.addcolumn(o.gv[:,1] , "gy")
        notpks.addcolumn(o.gv[:,2] , "gz")
        notpks.filter(drlv2 > o.tolerance*o.tolerance)
        notpks.writefile(options.notindexed)
        logging.info("Peaks which were not indexed %d written to %s"%(
                notpks.nrows, options.notindexed))
    if options.newubifile is not None:
        o.scandata[options.fltfile] = gotpks
#        matrix = o.refine(matrix,quiet=True)
        grain.write_grain_file(options.newubifile,[o.grains[(gn,options.fltfile)]])
        logging.info("Refined ubi in %s "%(
                         options.newubifile))

def get_options(parser):
    parser=refinegrains.get_options(parser)
    parser.add_argument("-N",  "--notindexed", action="store",
                      dest="notindexed", 
                      type=ImageD11options.ColumnFileType(mode='w'),
                      help="Name of flt file for unindexed peaks")
    parser.add_argument("-g","--grain", action="store",
                      dest = "grain", type=int, default=None,
                      help = "Which grain to choose")

    parser.description = """
Filtergrain should choose the peaks from a filtered
peaks output file according to those which are closest
to a particular grain
    """
    return parser

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


    parser = get_options( ArgumentParser() )

    options = parser.parse_args()

    if None in [options.parfile,
                options.ubifile,
                options.newfltfile,
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



