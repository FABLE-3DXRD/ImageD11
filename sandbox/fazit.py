
"""
fast radial regrouping method proposed by Peter Boesecke
"""
import numpy as np, pylab as pl

# 2018 code:
#   compute radius / azimuth
#        integer bin for radius
#        sort on radial bin, then azimuth (or vice versa)
#
#   Re-order an input array to write output
#    no splits or sums n->n mapping
#

def compute_r_chi( dims, pars ):
    """
    Replace this with pyFAI / ImageD11 / whatever
      dims = array dimension
      pars = center in slow/fast
    returns r/chi (float)
    """
    i, j = np.mgrid[0:dims[0], 0:dims[1]]
    i = i - pars[0]
    j = j - pars[1]
    r = np.sqrt( i*i + j*j )
    c = np.arctan2( i, j )
    return r, c

dims = ( 2048, 2048 )
pars = 1022.1, 1018.5, 1.25
rad, chi = compute_r_chi( dims, pars )

irad = np.round( rad / pars[2] ).astype( int )
import time
start=time.time()
lut = np.lexsort( (chi.ravel(), irad.ravel() ) )
print(time.time()-start)
start=time.time()
radout = rad.ravel()[lut]
iradout = irad.ravel()[lut]
chiout = chi.ravel()[lut]
# this is where the pixels read from
adrout = np.arange(dims[0]*dims[1])[lut]
# required to save
#   length of radial rows, r, IR and chi value for each
# TODO: optimise the order to apply lut
#      assert len(set(lut)) == len(lut)


print(time.time()-start)



pl.subplot(221)
pl.imshow(rad)
pl.title("radius")
pl.colorbar()
pl.subplot(222)
pl.imshow(chi)
pl.title("azimuth")
pl.subplot(223)
pl.imshow(adrout.reshape(dims),aspect='auto',
          cmap='prism', interpolation='nearest')

pl.title("Source address to come from")
pl.subplot(224)
pl.imshow(np.argsort(lut).reshape(dims),aspect='auto',
          cmap='prism', interpolation='nearest')

pl.title("Output address to go to")
pl.show()
