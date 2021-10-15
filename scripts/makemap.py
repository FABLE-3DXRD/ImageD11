#!/usr/bin/env python

from __future__ import print_function


from ImageD11.indexing import readubis, write_ubi_file
from ImageD11.refinegrains import refinegrains
import ImageD11.refinegrains
from ImageD11 import ImageD11options
import sys, os, argparse


def makemap(options):
    try:
        if options.tthrange is None:
            tthr = (0.,180.)
        else:
            tthr = options.tthrange
            if len(tthr) == 1:
                tthr = ( 0, tthr[0] )
            print("Using tthrange",tthr)
        func = getattr(ImageD11.refinegrains, options.latticesymmetry )
        o = refinegrains(intensity_tth_range = tthr,
                         latticesymmetry = func, 
                         OmFloat=options.omega_float,
                         OmSlop=options.omega_slop)
    except:
        raise
        o = refinegrains()
    o.loadparameters(options.parfile)
    print("got pars")
    o.loadfiltered(options.fltfile)
    print("got filtered")
    o.readubis(options.ubifile)
    if options.symmetry != "triclinic":
        # Grainspotter will have already done this
        print("transform to uniq")
        o.makeuniq(options.symmetry)
    print("got ubis")
    o.tolerance = float(options.tol)
    print("generating")
    o.generate_grains()
    print("Refining posi too")
    # o.refineubis(quiet = False , scoreonly = True)
    print("Refining positions too")
    o.refinepositions()
    print("Done refining positions too")    
    # o.refineubis(quiet = False , scoreonly = True)
    o.savegrains(options.newubifile, sort_npks = options.sort_npks)
    col = o.scandata[options.fltfile].writefile( options.fltfile+".new")
    if hasattr(options, "newfltfile") and options.newfltfile is not None:
        print("re-assignlabels")
        o.assignlabels()
        col = o.scandata[options.fltfile].copy()
        print("Before filtering",col.nrows)
        col.filter(col.labels < -0.5)
        # print col.labels[:10]
        print("After filtering",col.nrows)
        col.writefile(options.newfltfile)
        

def get_options(parser):
    parser = ImageD11.refinegrains.get_options(parser)
    parser.add_argument("--no_sort", action="store_false",
                      dest="sort_npks", default = True,
                      help="Sort grains by number of peaks indexed")
    parser.add_argument( "--tthrange", action="append",
                      dest = "tthrange", type=float,
                      default = None,
                      help= "Two theta range for getting median intensity")
    return parser


if __name__ == "__main__":
    import logging, sys
    
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


    parser = get_options( argparse.ArgumentParser() )
    
    options = parser.parse_args()

    for name in ["parfile" , 
                 "ubifile", 
                 "newubifile",
                 "fltfile"]:
        if getattr(options, name) is None:
            parser.print_help()
            logging.error("Missing option "+name)
            import sys
            sys.exit()

        
    try:
        makemap(options)
    except:
        parser.print_help()
        import traceback
        logging.error("An error occurred, details follow")
        traceback.print_exc()

