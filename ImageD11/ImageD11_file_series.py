
from __future__ import print_function

"""
To be moved to fabio sometime
"""

import fabio.file_series
import fabio.fabioimage
import fabio.openimage
import numpy, h5py
import gzip, bz2
from ImageD11 import ImageD11options

def get_options(parser):

    parser.add_argument("-5","--hdf5",action="store", type=str,
                      dest = "hdf5", default = None,
                      help = "hdf file containing input image series")
    # or, eventually:
    # stem, first, last, format, (omegas better be in the headers)
    parser.add_argument("-n","--stem",action="store", type=str,
                      dest = "stem", default = None,
                      help = "stem name for input image series")
    parser.add_argument("-f","--first",action="store", type=int,
                      dest = "first", default = None,
                      help = "first number for input image series")
    parser.add_argument("-l","--last",action="store", type=int,
                      dest = "last", default = None,
                      help = "last number for input image series")
    parser.add_argument("--ndigits", action="store", type=int,
                      dest = "ndigits", default = 4,
                      help = "Number of digits in file numbering [4]")
    parser.add_argument("-P", "--padding", action="store",
                      choices=["Y","N"], 
                      default="Y", dest="padding",
                      help="Is the image number to padded Y|N, e.g. "\
                          "should 1 be 0001 or just 1 in image name, default=Y")
    parser.add_argument("-F","--format",action="store", type=str,
                      dest = "format", default = ".edf",
                      help = "format [.edf] for input image series")

    parser.add_argument("-O", "--flood", action="store", 
                        type=ImageD11options.ImageFileType(mode='r'),
                      dest = "flood", default = None,
                      help = "Flood")
    
    parser.add_argument("-d", "--dark", action="store", 
                      dest = "dark", default = None,
                      type=ImageD11options.ImageFileType(mode='r'),
                      help = "Dark image")
    parser.add_argument("-S", "--step", action="store", type=float,
                      dest = "OMEGASTEP", default = None,
                      help = "omega step size")
    parser.add_argument("-T", "--start", action="store", type=float,
                      dest = "OMEGA", default = None,
                      help = "start omega")
    parser.add_argument("--omega_motor", action="store", type=str,
                      dest = "omegamotor", default = "Omega",
       help = "Header value to use for rotation motor position [Omega]")
    parser.add_argument("--omega_motor_step", action="store", type=str,
                      dest = "omegamotorstep", default = "OmegaStep",
       help = "Header value to use for rotation width [OmegaStep]")

    return parser


def get_series_from_hdf( hdf_file, dark = None, flood = None ):
    groups = hdf_file.listnames()
    for group in groups:
        imagenames = hdf_file[group].listnames()
        for image in imagenames:
            im = hdf_file[group][image]
            om = float(im.attrs['Omega'])
            data = im[:,:]
            if (dark, flood) != (None, None):
                data = data.astype(numpy.float32)
            if dark is not None:
                numpy.subtract( data, dark, data )
            if flood is not None:
                numpy.divide( data, flood, data )
            yield fabio.fabioimage.fabioimage( data = data,
                                               header = {
                    'Omega': om } )

def series_from_fabioseries( fabioseries, dark, flood, options ):
    for filename in fabioseries:
        try:
            fim = fabio.openimage.openimage(filename)
        except:
            print("Missing image",filename)
            continue
        if (dark, flood) != (None, None):
            fim.data = fim.data.astype(numpy.float32)
        if dark is not None:
            numpy.subtract( fim.data, dark, fim.data )
        if flood is not None:
            numpy.divide( fim.data, flood, fim.data )
        if options.omegamotor in fim.header:
            fim.header['Omega'] = float(fim.header[options.omegamotor])
            try:
                fim.header['OmegaStep'] = float(fim.header[options.omegamotorstep])
            except:
                fim.header['OmegaStep'] = float(options.OMEGASTEP)            
        else:
            fim.header['Omega'] = float(options.OMEGA)
            fim.header['OmegaStep'] = float(options.OMEGASTEP)            
            options.OMEGA = float(options.OMEGA) + float(options.OMEGASTEP)
        yield fim
        


def get_series_from_stemnum( options, args, dark = None, flood = None ):
    """
    Returns a file series thing - not a fabio one
    """
    if options.format in ['bruker', 'BRUKER', 'Bruker']:
        extn = ""
    elif options.format == 'GE':
        extn = ""
    else:
        extn = options.format
        
    fso = fabio.file_series.numbered_file_series(
        options.stem,
        options.first,
        options.last,
        extn,
        digits = options.ndigits,
        padding = options.padding )
    return series_from_fabioseries( fso , dark, flood, options )
    

def get_series_from_options( options, args ):
    """
    Returns a file series thing - not a fabio one

    This gives back a fabioimage object with dark and flood
    corrected data
    """

    try:
        if options.dark is not None:
            dark = fabio.openimage.openimage( options.dark ).data
        else:
            dark = None
    except:
        print("Problem with your dark",options.dark)
        raise
    
    try:
        if options.flood is not None:
            flood = fabio.openimage.openimage( options.flood ).data
        else:
            flood = None
    except:
        print("Problem with your flood",options.flood)
        raise
        

    if len(args) > 0 :
        # We assume unlabelled arguments are filenames 
        fso = fabio.file_series.file_series(args)
        return series_from_fabioseries( fso, dark, flood, options )

    if options.hdf5 is not None:
        hf = h5py.File(options.hdf5)
        # print "Getting images from",options.hdf5
        return get_series_from_hdf( hf, dark, flood )
    
    return get_series_from_stemnum( options, args,
                                     dark, flood) 
    


def getedfheader(filename):
    """
    Reads a header from an edf file in 1024 byte chunks.
    Assumes enclosing { }
    Returns string enclosed
    Adds a filename key at the top
    """
    h = "filename = "
    if filename[-3:]==".gz":
        fp=gzip.GzipFile(filename,"rb")
    elif filename [-4:]==".bz2":
        fp=bz2.BZ2File(filename,"rb")
    else:
        try:
            fp=open(filename,"rb")
        except IOError:
            return ""
    h=h+filename+";\n"
    s=fp.read(1024)
    if s.find("{")==-1:
        raise Exception("Not an edf file")
    while 1:
        if s.find("}")>=0:
            h=h+s[0:s.find("}")+2]
            break
        else:
            h=h+s
        s=fp.read(1024)
    return h

def motor_mne(hd):
    """
    expands the _mne and _pos header items of edf headers
    """
    h = {}
    order = []
    for line in hd.split(";"):
        try:
            key,vals = line.split("=")
        except ValueError:
            continue
        key = key.lstrip().rstrip()
        h[key] = vals.split(";")[0]
        order.append( key )
    for k in order:
        if k.endswith("_mne"):
            stem = k.split("_")[0]
            p = k.replace("_mne","_pos")
            newkeys = h[k].split()
            newvals = h[p].split()
            for ik, iv in zip(newkeys, newvals):
                kk = stem+":"+ik
                h[kk]=iv
                order.append( kk )
    return h, order
