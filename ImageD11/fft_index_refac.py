

"""
Classes for finding lattice vectors

We assume a very large unit cell
Compute the Patterson Function using this cell
Origin peaks show up for our crystals
These will be real space vectors
"""

import numpy.oldnumeric as n
import logging, time
from numpy.oldnumeric.fft import fftnd , inverse_fftnd

from ImageD11 import labelimage

class grid:
    def __init__(self, np=128, mr=1.0 ,nsig =5):
        """
        Set up the grid to use (the large unit cell)
        np - number of points in the grid
        mr - maximum resolution limit to consider
        """
        self.np = np
        self.nsig = nsig
        self.grid = n.zeros((np,np,np),n.Float32)
        self.old_grid = n.zeros((np,np,np),n.Float32)
        self.cell_size = np*mr/2.
        print "Using an FFT unit cell of ",self.cell_size

    

    def gv_to_grid_new(self, gv):
        """
        Put gvectors into our grid
        gv - ImageD11.indexing gvectors
        """
        # Compute hkl indices in the fft unit cell
        print "Gridding data"
        hrkrlr = self.cell_size * gv
        hkl = n.floor(hrkrlr + 0.5).astype(n.Int)
        # Filter to have the peaks in asym unit
        # ... do we need to do this? What about wrapping?
        hmx = n.where( hkl.max(axis=1) <  self.np/2 - 2.,
                       1, 0 )
        hmn = n.where( hkl.min(axis=1) > -self.np/2 + 2.,
                       1, 0 )        
        my_g = n.compress( hmx & hmn, gv, axis=0 )
        
        # Compute hkl indices in the fft unit cell
        hrkrlr = self.cell_size * gv
        hkl = n.floor(hrkrlr).astype(n.Int)

        # Fractional part of indices
        remain = hrkrlr - hkl
        grid = self.grid
        import time
        start = time.time()
        # print hkl
        for cor in [ (0,0,0),   #1
                     (1,0,0),   #2
                     (0,1,0),   #3
                     (0,0,1),   #4
                     (1,1,0),   #5
                     (1,0,1),   #6
                     (0,1,1),   #7
                     (1,1,1) ]: #8
            
            thkl = hkl + cor
            fac = 1 - n.absolute(remain - cor)
            vol = n.absolute(fac[:,0]*fac[:,1]*fac[:,2])
            thkl = n.where(thkl < 0, self.np + thkl, thkl)
            ind = thkl[:,0]*grid.shape[1]*grid.shape[2] + \
                  thkl[:,1]*grid.shape[1] + \
                  thkl[:,2]
            # print thkl,vol
            vals = n.take(grid.flat, ind)
            # print "vals",vals
            vals = vals + vol
            # print "vol",vol
            # print "ind",ind
            # This is wrong - if ind repeats several times we only take the last value
            n.put(grid, ind, vals)
            # print n.argmax(n.argmax(n.argmax(grid)))
            # print grid[0:3,-13:-9,-4:]
        logging.warn("Grid filling loop takes "+str(time.time()-start)+" /s")

        
    def gv_to_grid_old(self, gv):
        g = self.old_grid
        n_use = 0
        from math import floor
        for p in gv:
            # compute hkl indices in g grid
            # units start for a 1x1x1 unit cell
            hrkrlr= self.cell_size * p
            hkl = n.floor(self.cell_size * p+0.5).astype(n.Int)
            assert n.maximum.reduce(abs(hrkrlr - hkl)) <= 0.5, str(hrkrlr)+" "+str(hkl)
    
            if n.maximum.reduce(hkl)<self.np/2-2 and n.minimum.reduce(hkl)>-self.np/2+2:
                # use the peak
                n_use += 1
    
                sm = 0.
                # print "oldhrkrlr",hrkrlr
                
                for hs in [0,1]:
                    h=int(floor(hrkrlr[0]))+hs
                    hfac = 1-abs(hrkrlr[0] - h)
                    for ks in [0,1]:
                        k=int(floor(hrkrlr[1]))+ks    
                        kfac = 1-abs(hrkrlr[1] - k)
                        for ls in [0,1]:
                            l=int(floor(hrkrlr[2]))+ls
                            lfac = 1-abs(hrkrlr[2] - l)

                            v = hfac*kfac*lfac
                            #                    print v,h,k,l,hfac,kfac,lfac
                            sm += v
                            # print h,k,l,v
                            g[h,k,l] = g[h,k,l] +  v
                            # print h,k,l,v
                            # g[-h,-k,-l] = g[-h,-k,-l] + v
                            #                assert np == 8, "np"
                assert abs(sm-1)<1e-5, "total peak, hkl "+str(hrkrlr)+" "+str(sm)
            else:
                pass # n_ignore +=1


    def fft(self):
        """ Compute the Patterson """
        start = time.time()
        self.patty = abs(fftnd(self.grid))
        logging.warn("Time for fft "+str(time.time()-start))
        self.origin = self.patty[0,0,0]
        logging.warn("Patterson origin height is :"+str(self.origin))

    def props(self):
        """ Print some properties of the Patterson """
        print "Patterson info"
        print self.patty.shape, type(self.patty)
        p = n.ravel(self.patty)
        m = n.mean(p)
        print "Average:",m
        p2 = p*p
        self.mean = m
        v = n.sqrt( (n.sum(p2) - m*m*len(p) ) /(len(p)-1) )
        print "Sigma",v
        self.sigma = v

    def peaksearch(self, peaksfile):
        """ Peaksearch in the Patterson """
        lio = labelimage.labelimage( (self.np, self.np),
                                     peaksfile )
        print "Peaksearching at",self.nsig,"sigma"
        thresh = self.mean + self.nsig  * self.sigma
        for i in range(self.np):
            lio.peaksearch(self.patty[i],
                           thresh,
                           i)
            lio.mergelast()
            if i%2 == 0:
                print ".",
                sys.stdout.flush()
        print
        lio.finalise()
                                     


#####
#        To Do
#        Convert the labelimage/columnfile to write to memory
#        process / score / group the peaks after peaksearching
#
#        --- Separate program needed now?
#            foreach peak - is it a new peak, or sum of previous ones
#                         - how many gvecs does it make integers of
#                         - which gvecs does it make integers of w.r.t others        
#             ngvecs * npeaks == 175 * 200,000 K - 10^7 = 10Mpixel image ?
#

def test(options):
    gvfile = options.gvfile
    mr = options.max_res
    np = options.np
    nsig = options.nsig
    from ImageD11 import indexing
    go = grid(np = np , mr = mr , nsig = nsig)
    io = indexing.indexer()
    io.readgvfile(gvfile)
    go.gv_to_grid_new(io.gv )
    go.fft()
    go.props()
    go.peaksearch(open(gvfile+".patterson_pks","w"))
    sys.exit()

    
    go.gv_to_grid_old(io.gv)
    
    diff = n.ravel(go.grid - go.old_grid)
    print "All OK if this is zero",n.add.reduce(diff)
    print "error at",n.argmax(diff),diff.max(),n.argmin(diff),diff.min()


if __name__=="__main__":
    
    import sys
    class options:
        max_res = 3.0
        np = 128
        nsig = 10
    
    options.gvfile = sys.argv[1]
    
    test(options)
