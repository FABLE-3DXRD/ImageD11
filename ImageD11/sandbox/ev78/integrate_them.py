#!/usr/bin/python 
from __future__ import print_function
import sys
#sys.path.append('/users/wright/software/lib/python')

import os, time, fabio, numpy
import pyFAI
print(pyFAI.__file__)
SOLID_ANGLE = True
#print "PATH:", sys.path

from pyFAI.azimuthalIntegrator import AzimuthalIntegrator

class darkflood(object):
    """ 
    For doing dark/flood subtractions
    Could instead be done as a wrapper or box around a file series
    to supply a series of clean data
    """
    def __init__(self, darkfile=None, floodfile=None):
        if darkfile is None:
            print("Warning: no dark supplied")
            self.darkobj = None
        else: 
            self.darkobj = fabio.open(darkfile)
        if floodfile is None:
            print("Warning: no flood supplied")
            self.floodmult = None
        else: 
            self.floodobj = fabio.open(floodfile)
            self.floodmult = 1.0 / fabio.open(floodfile).data
    def correct( self, datao ):
        if self.darkobj is not None:
            datao.data = datao.data.astype(numpy.float32)
            numpy.subtract( datao.data, self.darkobj.data, datao.data)
        if self.floodmult is not None:
            numpy.multiply( datao.data, self.floodmult, datao.data)
        return datao

def edffilenameseries(stem, first=None, last=None, glob=False, extn=".edf"):
    """
    Should be in fabio?
    """
    if glob:
        import glob, os
        files = glob.glob("%s*%s"%(stem,extn))
        filesdone = glob.glob("%s*%s"%(stem,".dat"))
        filesdone = [f.replace(".dat",extn) for f in filesdone]
        for f in files:
            if f not in filesdone:
                yield f
    else:
        if first is not None and last is not None:
            num = first
            while num <= last:
                fname = "%s%04d%s"%( stem, num, extn )
                yield fname
                num += 1
        if first is not None and last is None:
            print("Starting with image",first,"then incrementing")
            num = first
            while 1:
                fname = "%s%04d%s"%( stem, num, extn )
                yield fname
                num += 1
        if first is None and last is None:
            raise Exception("File series not specified")



def calcfrom1d( integrator, tth, I, shape):
    """
    Computes a 2D image from a 1D integrated profile
    """
    ttha = integrator.twoThetaArray(shape)
    print(ttha.ravel().shape)
    print(tth.shape)
    print(I.shape)
    calcimage = numpy.interp( ttha.ravel(),
                              tth*numpy.pi/180,
                              I )
    calcimage.shape = shape
    # Solid angle correction
    # flake8: global SOLID_ANGLE
    if SOLID_ANGLE:
        numpy.multiply( calcimage, integrator._dssa, calcimage )
    return calcimage


def determineparfile( fname ):
    """
    Guess if we have a fit2dcake.py parameter file or a pyFAI 
    parameter file
    """
    bytes = open(fname).read()
    if bytes.find("Poni1") != -1:
        return "pyfai"
    if bytes.find("X-BEAM CENTRE") != -1:
        return "fit2d"
    print(bytes)
    raise Exception("Parameters not recognised in file: %s"%(fname))

class fit2dcakepars:
    """
    For loading fit2dcake parameters
    """
    option_names = [
        "DARK CURRENT",
        "DC FILE",
        "FLAT-FIELD",
        "FF FILE",
        "FF SCALE",
        "FF MULTIPLIER",
        "SPATIAL DIS.",
        "SD FILE",
        "X-PIXEL SIZE",
        "Y-PIXEL SIZE",
        "DISTANCE",
        "WAVELENGTH",
        "X-BEAM CENTRE",
        "Y-BEAM CENTRE",
        "TILT ROTATION",
        "ANGLE OF TILT",
        "START AZIMUTH",
        "END AZIMUTH",
        "INNER RADIUS",
        "OUTER RADIUS",
        "SCAN TYPE",
        "1 DEGREE AZ",
        "AZIMUTH BINS",
        "RADIAL BINS",
        "CONSERVE INT.",
        "POLARISATION",
        "GEOMETRY COR.",
        "USE MASK",
        "MASK FILE",
        "DIM1_DATA",
        "DIM2_DATA",
        "input_extn", "saving_format", "output_extn"
        ]
    options = {}
    def __init__(self, fname ):
        """Read parfile"""
        for o in self.option_names: 
            self.options[o]=None
        pf = open(fname, "r").readlines()
        for line in pf:
            if line[0] == "#" or \
               len(line) < 3 or \
               line.find("=") < 0 :
                continue
            try:
                key, value = line.split(" = ")
                value = value.rstrip() # get rid of trailing \n
            except:
                print("$$$%s$$$"% (line))
                raise
            if value == "None":
                value = None
            self.setparameter(key, value)
        self.parsread = True
    def setparameter( self, key, value ):
        if key in self.option_names:
            self.options[key] = value
        else:
            print("Failure to set {%s}=%s"%( key, value))
    def __getitem__(self, key):
        return self.options[key]
    def __setitem__(self, key, value):
        self.options[key] = value



def display(tth, I, img):
    from matplotlib.pylab import imshow, show, subplot, colorbar,figure,plot
    figure()
    plot(tth,I,"-")
    figure()
    imshow(img,vmin=-.21,vmax=.21)
    colorbar()
    show()
    input("continue?")

def integrate_them( o ):
    """ 
    Process a series of files
    o = options object 
    o.parfile gives name of parameter file (fit2d or poni format)
    o.dark overrides to supply dark filename
    o.flood overrides to supply flood filename
    o.mask  overrides to supply mask filename
    o.backcalc asks for the back computation of the image
    o.npts gives number of output points
    """
    # pyFAI.load( ponifile )
    integrator = AzimuthalIntegrator()
    #integrator.tth = integrator.newtth
    # integrator.setChiDiscAtZero()
    ptype = determineparfile( o.parfile )
    if ptype == "pyfai":
        integrator.load( o.parfile )
        if o.dark is not None:
            print("Using dark from command line",o.dark)
        if o.flood is not None:
            print("Using dark from command line",o.flood)
    elif ptype == "fit2d":
        f2d = fit2dcakepars( o.parfile )
        if f2d["SPATIAL DIS."][0] not in ["Y","y"]:
            # Set to None. Spatial is from parfile
            f2d["SD FILE"] = None
        integrator.setFit2D(
            float(f2d["DISTANCE"]),
            float(f2d["X-BEAM CENTRE"]),
            float(f2d["Y-BEAM CENTRE"]),
            tilt=float(f2d["ANGLE OF TILT"]),
            tiltPlanRotation=float(f2d["TILT ROTATION"]),
            pixelX=float(f2d["X-PIXEL SIZE"]),
            pixelY=float(f2d["Y-BEAM CENTRE"]),
            splineFile=f2d["SD FILE"]
            )
        integrator.rot3=0
        integrator.reset()
        print(integrator.param, integrator.detector.pixel1)
        # First choice is command line. Then from pars if supplied
        if o.dark is None:
            if f2d["DARK CURRENT"][0] in ["Y","y"]:
                o.dark = f2d["DC FILE"]
                print("Using dark from fit2d parameter file",o.dark)
        else:
            print("Using dark from command line",o.dark)
        if o.flood is None:
            if f2d["FLAT-FIELD"][0] in ["Y","y"]:
                o.flood = f2d["FF FILE"]
                print("Using flood from fit2d parameter file",o.flood)
        else:
            print("Using flood from command line",o.flood)
    # Should be in fabio utilities
    df = darkflood( o.dark, o.flood )
    # Should be in fabio
    fs = edffilenameseries( o.stem, o.first, o.last, o.glob, o.extn )
#    integrator.polarization( factor = 1, shape=(2048,2048) )
    # Command line is first priority for make
    if o.mask is not None:
        mask = fabio.open( o.mask ).data
        print("Using mask", o.mask)
        # assume poni file deals with this independently?
    elif ptype == "fit2d": # try in fit2d parfile
        if f2d["USE MASK"][0] in ['y','Y']:
            mask = fabio.open( f2d["MASK FILE"] ).data
            print("Using mask",f2d["MASK FILE"])
        else:
            mask = None
    if mask is not None:
        print("mask mean:",mask.mean())
    integrator.write(os.path.splitext(o.parfile)[0]+".poni")
    for f in fs:
        print("Processing",f, end=' ')
        try:
            fo = df.correct( fabio.open(f) )
        except:
            
            continue
        if ptype == "fit2d":
            outFile = f.replace(f2d["input_extn"], f2d["output_extn"])
        else:
            outFile = f.replace(o.extn,".dat")
        #flake8: global SOLID_ANGLE
        if 0:
            from matplotlib.pylab import imshow, figure, show, log, plot
            #imshow(log(fo.data))
            #figure()
            if mask is not None:
                imshow(log(fo.data * (1-mask)),vmin=log(10),vmax=log(30000))
            else:
                imshow(log(fo.data),vmin=log(100),vmax=log(3000))
            #        show()
        if o.npts is None:
            npts = min(fo.data.shape)
        else:
            npts = int(o.npts)
        tth, I = integrator.integrate1d(fo.data,
                                 nbPt=npts,
                                 filename=outFile,
                                 correctSolidAngle=SOLID_ANGLE,
                                 mask=mask, # 1 for valid
                                 unit="q_A^-1",
                                 #dummy=dummy, # mask pixels == dummy
                                 #delta_dummy=delta_dummy # precision of dummy
                                 )

        print("wrote",outFile)
        if o.backcalc:
            calcimage = calcfrom1d( integrator, tth, I, fo.data.shape ) * integrator._polarization
            err = (calcimage - fo.data) * (1-mask)/(calcimage+mask)
            e = fabio.edfimage.edfimage( data = err.astype(numpy.float32) )
            e.write( outFile+".edf")
            fitcen( fo.data, calcimage, (1-mask) )
#            from matplotlib.pylab import imshow, show
#            imshow( integrator._polarization )       
#            show()
    
    if o.display:
        if mask is not None:
            display( tth, I, (calcimage - fo.data) * (1-mask)/(calcimage+1))
        else:
            display( tth, I, (calcimage - fo.data)/(calcimage+1))


def fitcen( data, cim , mask=None ):
    """ Do a fit to the beam center in pixel coordinates """
    dcdx = numpy.zeros(data.shape, numpy.float32)
    dcdx[1:-1] = cim[2:]-cim[:-2]
    dcdy = numpy.zeros(data.shape, numpy.float32)
    dcdy[:,1:-1] = cim[:,2:]-cim[:,:-2]
    diff = (cim - data)
    if mask is not None:
        numpy.multiply( dcdx, mask, dcdx)
        numpy.multiply( dcdy, mask, dcdy)
        numpy.multiply( diff, mask, diff)
    myy = (dcdy*dcdy).sum()
    mxy = (dcdx*dcdy).sum()
    mxx = (dcdx*dcdx).sum()
    ry = (dcdy*diff).sum()
    rx = (dcdx*diff).sum()
    try:
        imat = numpy.linalg.inv( [ [ mxx, mxy], [mxy, myy] ] )
        shft = numpy.dot( imat, [rx, ry] )
        print("Image offset estimate (pixels)", shft)
        print("Intensity mean, min, max errors", (diff*diff).mean(), diff.min(), diff.max())
    except:
        print("Cannot estimate beam centre error")
    

def get_options(parser):
    parser.add_option("-n", "--namestem", action="store",
                      dest="stem", type="string", default=None,
                      help="Name of the files up the digits part  "+
                      "eg mydata in mydata0000.edf" )
    parser.add_option("-F", "--extn", action="store",
                      dest="extn", type="string", default=".edf",
                      help="Filename extension  "+
                      "eg .edf in mydata0000.edf" )
    parser.add_option("-f", "--first", action="store",
                      dest="first", default=None, type="int",
                      help="Number of first file to process, default=0")
    parser.add_option("-l", "--last", action="store",
                      dest="last", type="int",default=None,
                      help="Number of last file to process")
    parser.add_option("-m", "--mask", action="store",
                      default=None, dest="mask",
                      help="Mask file (from fit2d)")
    parser.add_option("-g", "--glob", action="store_true",
                      default=False, dest="glob",
                      help="Keep watching for new edf files with stem and process them")
    parser.add_option("-p", "--parfile", action="store",
                      default=None, dest="parfile",
                      help="Parameter file from fit2dcake.py or pyFAI")
    parser.add_option("-d", "--darkfile", action="store",
                      dest="dark", default=None,  type="string",
                      help="Dark current filename, to be subtracted, default=None")
    parser.add_option("-O", "--flood", action="store", type="string",
                      default=None, dest="flood",
                      help="Flood file, default=None")
    parser.add_option("-N", "--Npoints", action="store", type="int",
                      default=None, dest="npts",
                      help="Number of points in output profile")
    parser.add_option("-c","--backcalc",action="store_true", default=False,
                      dest = "backcalc", 
                      help="Compute ideal image and fit parameters")
    parser.add_option("--display",action="store_true", default=False,
                      dest = "display", 
                      help="Display difference image after computing ideal")

    return parser

def check_pars( o ):
    assert o.stem is not None, "Stem name for images is not supplied"
    assert o.parfile is not None and os.path.exists( o.parfile ), \
        "Need a parameter file, -p option"






if __name__ == "__main__":
    from optparse import OptionParser
    o = get_options( OptionParser() )
    options , args = o.parse_args()
    try:
        check_pars( options )
    except:
        if o is not None:
            o.print_help()
        print("\nHere is the problem:\n")
        raise
    integrate_them( options )

    # calcimage = calcfrom1d( integrator, tth, I, data.shape, mask)

    print(pyFAI.__file__)


    

#    show()
#    raw_input()
