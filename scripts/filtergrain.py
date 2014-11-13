#!/usr/bin/env python

"""
Utility script to pick out the peaks belonging to a certain
grain

Used to stabilise refinements
"""

import sys, logging

import numpy as np

from ImageD11 import refinegrains, indexing, grain


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
            for i,gn in zip(range(len(o.grainnames)),o.grainnames):
                logging.info("Choose %d for grain %s"%(i, gn))
            gn = o.grainnames[int(raw_input("select which grain "))]
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
        o.tolerance = float(raw_input("Enter tolerance "))
        options.tol = o.tolerance
    else:
        o.tolerance = options.tol
    matrix = o.grains[(gn,options.fltfile)].ubi

    o.assignlabels()
    drlv2 = indexing.calc_drlv2( matrix, o.gv )
    logging.info("Total peaks before filtering %d"%
                     o.scandata[options.fltfile].nrows)
    gotpks = o.scandata[options.fltfile].copy()
    gotpks.filter(gotpks.labels == 0)
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
    parser.add_option("-N",  "--notindexed", action="store",
                      dest="notindexed", type="string",
                      help="Name of flt file for unindexed peaks")
    parser.add_option("-t", "--tol", action="store",
                      dest="tol", type="float",
                      default = None,
                      help="Tolerance to use in peak assignment")
    parser.add_option("-g","--grain", action="store",
                      dest = "grain", type="string", default=None,
                      help = "Which grain to choose")
    parser.add_option( "--omega_no_float", action="store_false",
                      dest = "omega_float",
                      default = True,
                      help= "Use exact observed omega values")

    parser.add_option( "--omega_slop", action="store", type="float",
                      dest = "omega_slop",
                      default = 0.5,
                      help= "Omega slop (step) size")

    parser.description = """
Filtergrain should choose the peaks from a filtered
peaks output file according to those which are closest
to a particular grain
    """

    options, args = parser.parse_args()

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



