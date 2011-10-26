# ImageD11 Software for beamline ID11
# Copyright (C) 2011  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


import numpy, sys

from ImageD11 import transform, indexing, parameters
from ImageD11 import grain, closest
import columnfile
print "You are using box-refine grains"


def update_det_pars( peakdata, pars ):
    """
    Update the xl, yl, zl columns
    """
    
    print "Update xl, yl, zl for current parameters"
    peakdata.xl,peakdata.yl,peakdata.zl = transform.compute_xyz_lab(
        [peakdata.sc, peakdata.fc],
        **pars.get_parameters() )
    
    # Imagine sample translates along +x. Detector is at x=+distance
    # To get the difference vector right we add samtx to xl for omega=0
    # Same for samty, yl and samtz, zl
    # when omega rotates to + 90 degrees (right handed)
    #   samtx is along +yl
    #   samty is along -xl
    rad = float(pars.parameters['omegasign'])*numpy.pi/180.0
    peakdata.sinomega = numpy.sin(peakdata.omega * rad) # == 0 at omega = 0
    peakdata.cosomega = numpy.cos(peakdata.omega * rad) # == 1 at omega = 0
    
    peakdata.xl += peakdata.samtx * peakdata.cosomega - peakdata.samty * peakdata.sinomega
    peakdata.yl += peakdata.samtx * peakdata.sinomega + peakdata.samty * peakdata.cosomega
    peakdata.zl += peakdata.samtz
    print "lab x shape",peakdata.xl.shape
    return peakdata
    
    


def potentialpeaks( grainlist, peakdata, pars, tol ):
    """
    For each grain, see which peaks are potentially interesting
    """
    i = 0
    assign = {}
    print "potential peak finding"
    #print grainlist
    indices = numpy.arange( len( peakdata.omega ) )
    for gc in grainlist:
        #print i
        tx, ty, tz = gc.translation
        t = transform.compute_grain_origins(-peakdata.omega,
                                             t_x=tx, t_y=ty, t_z=tz )
        t[0,:] = peakdata.xl - t[0,:]
        t[1,:] = peakdata.yl - t[1,:]
        t[2,:] = peakdata.zl - t[2,:]
        #print t
        # unit vectors
        fac = numpy.sqrt( (t*t).sum(axis=0) )
        numpy.divide( t, fac, t )
        #print t
        # t is now the k vector
        beam_direction = (1,0,0)
        for j in range(3):
            numpy.subtract( t[j,:], beam_direction[j], t[j,:] )
        #print "k",t
        g = numpy.zeros(t.shape, numpy.float32)
        g[0,:]=  peakdata.cosomega * t[0,:] + peakdata.sinomega*t[1,:]
        g[1,:]= -peakdata.sinomega * t[0,:] + peakdata.cosomega*t[1,:]
        g[2,:]=  t[2,:]
        #print "g",g
        hklr = numpy.dot( gc.ubi, g ) / float(pars.parameters['wavelength'])
        #print hklr
        hkli = numpy.floor( hklr + 0.5 )
        drlv = (hklr - hkli)
        drlv = numpy.sqrt( (drlv*drlv).sum(axis=0))
        #print i,drlv.max(), drlv.min(),drlv.mean(),drlv.shape,g.shape,hklr.shape
        inds = numpy.compress( drlv < tol, indices )
        # Keep only hkl indexing which is fixed, all other things can change
        assign[i] = inds, numpy.take(hkli, inds, axis=1)
        print i, len(inds)#, tol, (drlv<tol).sum(), len(indices)
        i+=1
    return assign
            
            
            



if __name__=="__main__":

    # Data object
    peakdata = columnfile.colfile_from_hdf( sys.argv[1] )

    # Grain collection object
    grainlist = grain.read_grain_file( sys.argv[2] )

    # Detector parameters
    pars = parameters.parameters()
    pars.loadparameters( sys.argv[3] )

    peakdata = update_det_pars( peakdata, pars )


    tol = 0.1
    assignments = potentialpeaks( grainlist, peakdata, pars, tol)

    
    
    
