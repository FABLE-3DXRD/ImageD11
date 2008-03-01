#!/usr/bin/python
## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py

#!/bliss/users/blissadm/python/bliss_python/suse82/bin/python

"""
Utility script to pick out the peaks belonging to a certain
grain

Used to stabilise refinements
"""

from ImageD11 import refinegrains, indexing
import numpy.oldnumeric as Numeric
import sys, logging


def filtergrain(options):
    """
    Filter a peaks file according to which peaks are indexed
    """   
    o = refinegrains.refinegrains(tolerance=0.1)
    o.loadparameters(options.parfile)
    o.readubis(options.ubifile)
    o.loadfiltered(options.fltfile)
    o.generate_grains()
    if options.grain is None:
        if len(o.grainnames) == 1:
            gn = o.grainnames[0]
        else:
            for i,gn in zip(range(len(o.grainnames)),o.grainnames):
                logging.info("Choose %d for grain %s"%(i, gn))
            gn = o.grainnames[int(raw_input("select which grain "))]
    else:
        gn = o.grainnames[int(options.grain)]  
    o.compute_gv(gn, options.fltfile)
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
    drlv2 = indexing.calc_drlv2( matrix, o.gv )    
    
    logging.info("Total peaks before filtering %d"%
                     o.scandata[options.fltfile].nrows)    
    gotpks = o.scandata[options.fltfile].copy()
    # Output some derived information (h,k,l etc)
    gotpks.addcolumn(o.tth, "tth")
    gotpks.addcolumn(o.eta, "eta")
    gotpks.addcolumn(o.gv[:,0] , "gx")
    gotpks.addcolumn(o.gv[:,1] , "gy")
    gotpks.addcolumn(o.gv[:,2] , "gz")
    # Should these calcs go to refinegrains...?
    hkl = Numeric.dot(matrix, Numeric.transpose(o.gv))
    gotpks.addcolumn(hkl[0,:] , "hr")
    gotpks.addcolumn(hkl[1,:] , "kr")
    gotpks.addcolumn(hkl[2,:] , "lr")
    hkli = Numeric.floor(hkl+0.5).astype(Numeric.Int)
    gotpks.addcolumn(hkli[0,:] , "h")
    gotpks.addcolumn(hkli[1,:] , "k")
    gotpks.addcolumn(hkli[2,:] , "l")
    gotpks.addcolumn(drlv2, "drlv2")
    gotpks.filter(drlv2 < o.tolerance*o.tolerance)
    gotpks.writefile(options.newfltfile)
    logging.info("Peaks which were indexed %d written to %s"%(
                gotpks.nrows, options.newfltfile))
    # don't bother to copy here as we can overwrite
    if options.notindexed is not None:
        notpks = o.scandata[options.fltfile]
        notpks.filter(drlv2 > o.tolerance*o.tolerance)
        notpks.writefile(options.notindexed)
        logging.info("Peaks which were not indexed %d written to %s"%(
                notpks.nrows, options.notindexed))
    if options.newubifile is not None:
        o.scandata[options.fltfile] = gotpks
        matrix = o.refine(matrix,quiet=True)
        indexing.write_ubi_file(options.newubifile,[matrix])
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
    


