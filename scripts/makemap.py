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
    if options.symmetry is not "triclinic":
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
    parser.add_argument("-p",  "--parfile", action="store",
                      dest="parfile", 
                      type=ImageD11options.ParameterFileType(mode='r'),
                      help="Name of parameter file")
    parser.add_argument("-u",  "--ubifile", action="store",
                      dest="ubifile", 
                      type=ImageD11options.UbiFileType(mode='r'),
                      help="Name of ubi file")
    parser.add_argument("-U",  "--newubifile", action="store",
                      dest="newubifile", 
                      type=ImageD11options.UbiFileType(mode='w'),
                      help="Name of new ubi file to output")
    parser.add_argument("-f",  "--fltfile", action="store",
                      dest="fltfile", 
                      type=ImageD11options.ColumnFileType(mode='r'),
                      help="Name of flt file")
    parser.add_argument("-F",  "--newfltfile", action="store",
                      dest="newfltfile", 
                      type=ImageD11options.ColumnFileType(mode='w'),
                      help="Name of flt file containing unindexed peaks")
    parser.add_argument("-t", "--tol", action="store",
                      dest="tol", type=float,
                      default =   0.5,
                      help="Tolerance to use in peak assignment")
    parser.add_argument("--no_sort", action="store_false",
                      dest="sort_npks", default = True,
                      help="Sort grains by number of peaks indexed")
    lattices = ["cubic", "hexagonal", "trigonal","rhombohedralP",
                "tetragonal", "orthorhombic", "monoclinic_a",
                "monoclinic_b","monoclinic_c","triclinic"]
    parser.add_argument("-s", "--sym", action="store",
                      dest="symmetry", # type="choice",
                      default = "triclinic",
                      choices = lattices,
                      help="Lattice symmetry for choosing orientation")

    parser.add_argument("-l", "--lattice", action="store",
                      dest="latticesymmetry", #type="choice",
                      default = "triclinic",
                      choices = lattices,
                      help="Lattice symmetry for constraints "+
                      "|".join(lattices))

    parser.add_argument( "--tthrange", action="append",
                      dest = "tthrange", type=float,
                      default = None,
                      help= "Two theta range for getting median intensity")

    parser.add_argument( "--omega_no_float", action="store_false",
                      dest = "omega_float",
                      default = True,
                      help= "Use exact observed omega values")

    parser.add_argument( "--omega_slop", action="store", type=float,
                      dest = "omega_slop",
                      default = 0.5,
                      help= "Omega slop (step) size")
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

