

"""
Make an image for masking a powder diagram to convert it to a set of
isolated peaks
"""


from ImageD11 import transform, parameters, blobcorrector
import numpy as np
from fabio.openimage import openimage

def make_powder_mask( parfile,
                      ndeg = 1,
                      splinefile=None,
                      dims=(2048, 2048) ):
    """
    Compute a two theta and azimuth image
    """
    pars = parameters.parameters()
    pars.loadparameters( parfile )
    if splinefile is None:
        spatial = blobcorrector.perfect()
    else:
        spatial = blobcorrector.correctorclass( splinefile )
    xim, yim = spatial.make_pixel_lut ( dims )
    peaks = [ np.ravel( xim ) , np.ravel( yim ) ]
    tth , eta = transform.compute_tth_eta( peaks , **pars.get_parameters() )
    tth.shape = dims
    eta.shape = dims
    # Assume a circle geometry for now
    # tth * eta ~ length on detector
    # lim = tth * eta
    # need some idea how to cut it up...
    #  degree bins
    m =  (eta.astype(int) % 2)==0
    return m

if __name__=="__main__":
    import sys
    parfile = sys.argv[1]
    m = make_powder_mask( parfile )
    obj = openimage(sys.argv[2])
    np.multiply( obj.data , m , obj.data )
    obj.write( sys.argv[3])

