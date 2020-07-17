#!/usr/bin/env python

from __future__ import print_function


"""
Wrapper script to refine a single grain using 
all peaks in a dataset (so that hkl assignments
are not a problem).
"""
from argparse import ArgumentParser
from ImageD11 import refinegrains, indexing, ImageD11options
import logging, sys

def fitgrain(options):
    """
    Fits grains to a dataset using all peaks within tol
    """
    func = getattr(refinegrains, options.latticesymmetry )
    o = refinegrains.refinegrains(tolerance = options.tol,
                                  latticesymmetry = func, 
                                  OmFloat=options.omega_float,
                                  OmSlop=options.omega_slop)
    o.loadparameters(options.parfile)
    o.readubis(options.ubifile)
    o.loadfiltered(options.fltfile)
    o.generate_grains()
    #o.refineubis(quiet = False)
    #print o.grains.keys()
    #print "#pks",o.grains[(0,"0.flt")].x.shape
    o.parameterobj.varylist = options.varylist
    for p in options.fixlist:
        try:
            o.parameterobj.varylist.remove(p)
        except:
            pass
    logging.info("Varying " + str(o.parameterobj.varylist))
    #print "#pks",o.grains[(0,"0.flt")].x.shape
    o.fit(maxiters = options.steps)
    #print "#pks",o.grains[(0,"0.flt")].x.shape
    o.refineubis(quiet = False)
    o.saveparameters(options.newparfile)
    # ul = [g.ubi for g in o.grains.values()]
    # indexing.write_ubi_file(options.newubifile, ul)
    # Keep the original ordering and add translation information
    o.savegrains(options.newubifile, sort_npks=False)

def get_options(parser):
    parser = refinegrains.get_options( parser )
    parser.add_argument("-P", "--newparfile", action="store",
                      dest="newparfile", 
                      type=ImageD11options.ParameterFileType(mode='w'),
                      help="Name of new parameter file")
    parser.add_argument("-v", "--vary", action="append",
                      dest="varylist", type=str,
                      default =    [ "y_center","z_center",
                                     "tilt_y","tilt_x","tilt_z","wedge",
                                     "t_x","t_y","distance"],
                      help="Parameters to vary"  )
    parser.add_argument("-x", "--fiX", action="append",
                      dest="fixlist", type=str, default = [],
                      help="Parameters to fix (overrides vary)")
    parser.add_argument("-S", "--steps", action="store",
                      dest="steps", type=int,
                      default =   1000,
                      help="Number of simplex iterations")
    
    parser.description = """
Fitgrain should attempt to fit one or more grains to a dataset
using the parameters specified on the command line.
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

    for name in ["parfile" , "newparfile",
                 "ubifile", "newubifile",
                    "fltfile"]:
        if getattr(options, name) is None:
            parser.print_help()
            logging.error("Missing option "+name)
            import sys
            sys.exit()
        
    try:
        fitgrain(options)
    except:
        parser.print_help()
        import traceback
        logging.error("An error occurred, details follow")
        traceback.print_exc()
    
