

"""
Classes for finding lattice vectors

We assume a very large unit cell
Compute the Patterson Function using this cell
Origin peaks show up for our crystals
These will be real space vectors
"""

import numpy.oldnumeric as n
import logging, time, sys
from numpy.oldnumeric.fft import fftnd , inverse_fftnd

from ImageD11 import labelimage, closest

def get_options(parser):
    parser.add_option( '-n', '--ngrid', 
                       action = 'store',
                       dest = 'np',
                       type = 'int',
                       help = 'number of points in the fft grid [128]',
                       default = 128 )
    parser.add_option( '-r', '--max_res',
                       action = 'store',
                       dest = 'mr',
                       type = 'float',
             help = 'Maximum resolution limit for fft (d-spacing) [1.0]',
                       default = 1.0)
    parser.add_option( '-s', '--nsig',
                       action = 'store',
                       dest = 'nsig',
                       type = 'float',
             help = 'Number of sigma for patterson peaksearch threshold [5]',
                       default = 5)
    return parser

class grid:
    def __init__(self, np = 128, mr = 1.0 ,nsig = 5):
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
        logging.info("Using an FFT unit cell of %s"%(str(self.cell_size)))

    

    def gv_to_grid_new(self, gv):
        """
        Put gvectors into our grid
        gv - ImageD11.indexing gvectors
        """
        # Compute hkl indices in the fft unit cell
        logging.info("Gridding data")
        self.gv = gv
        hrkrlr = self.cell_size * gv
        hkl = n.floor(hrkrlr + 0.5).astype(n.Int)
        # Filter to have the peaks in asym unit
        # ... do we need to do this? What about wrapping?
        hmx = n.where( hkl.max(axis=1) <  self.np/2 - 2.,
                       1, 0 )
        hmn = n.where( hkl.min(axis=1) > -self.np/2 + 2.,
                       1, 0 )        
        my_g = n.compress( hmx & hmn, gv, axis=0 )
        
        # Compute hkl indices in the fft unit cell using filtered peaks
        hrkrlr = self.cell_size * my_g

        # Integer part of hkl (eg 1.0 from 1.9)
        hkl = n.floor(hrkrlr).astype(n.Int)

        # Fractional part of indices ( eg 0.9 from 1.9)
        remain = hrkrlr - hkl
        grid = self.grid
        import time
        start = time.time()
        # Loop over corners with respect to floor corner
        ng = grid.shape[0]*grid.shape[1]*grid.shape[2]
        flatgrid = grid.reshape( ng )
        for cor in [ (0,0,0),   #1
                     (1,0,0),   #2
                     (0,1,0),   #3
                     (0,0,1),   #4
                     (1,1,0),   #5
                     (1,0,1),   #6
                     (0,1,1),   #7
                     (1,1,1) ]: #8
            # The corner
            thkl = hkl + cor
            fac = 1 - n.absolute(remain - cor)
            vol = n.absolute(fac[:,0]*fac[:,1]*fac[:,2]).astype(n.float32)
            thkl = n.where(thkl < 0, self.np + thkl, thkl)
            ind = thkl[:,0]*grid.shape[1]*grid.shape[2] + \
                  thkl[:,1]*grid.shape[1] + \
                  thkl[:,2]

            closest.put_incr( flatgrid , ind.astype(n.int32), vol )
            # print thkl,vol

            # vals = n.take(grid.flat, ind)
            # print "vals",vals
            # vals = vals + vol
            # print "vol",vol
            # print "ind",ind
            # FIXME
            # This is wrong - if ind repeats several times we 
            # only take the last value
            # n.put(grid, ind, vals)
            # print n.argmax(n.argmax(n.argmax(grid)))
            # print grid[0:3,-13:-9,-4:]
        logging.warn("Grid filling loop takes "+str(time.time()-start)+" /s")



    def fft(self):
        """ Compute the Patterson """
        start = time.time()
        self.patty = abs(fftnd(self.grid))
        logging.info("Time for fft "+str(time.time()-start))
        self.origin = self.patty[0,0,0]
        logging.info("Patterson origin height is :"+str(self.origin))

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
        logging.info("Peaksearching at %f sigma"%(self.nsig))
        thresh = self.mean + self.nsig  * self.sigma
        for i in range(self.np):
            lio.peaksearch(self.patty[i],
                           thresh,
                           i)
            lio.mergelast()
            #if i%2 == 0:
            #    print ".",
            #    sys.stdout.flush()
        #print
        lio.finalise()
                                     
    def pv(self, v): return ("%8.4f "*3)%tuple(v)

    def reduce(self, vecs):
        raise Exception("You want lattice_reduction instead")
        from ImageD11 import sym_u
        g =  sym_u.trans_group(tol = self.rlgrid*2) 
        # print self.rlgrid, len(vecs)
        for i in range(len(vecs)):
            print self.pv(vecs[i]),
            assert len(vecs[i]) == 3
            t = g.reduce(vecs[i])
            if n.dot(t, t) < self.minlen * self.minlen:
                print "Short ",self.pv(vecs[i]), i, self.pv(t)
                continue
            if not g.isMember(t):
                g.additem(t)
                print "Adding",self.pv(vecs[i]), i, self.pv(t)
                print len(g.group), g.group
            else:
                print "Got   ",self.pv(vecs[i]), i , self.pv(t)
        print g.group
        return g.group[1:]

    def read_peaks(self, peaksfile):
        """ Read in the peaks from a peaksearch """
        from ImageD11.columnfile import columnfile
        import time
        start = time.time()
        self.colfile = columnfile(peaksfile)
        print "reading file",time.time()-start
        # hmm - is this the right way around?
        self.rlgrid = 1.0*self.cell_size/self.np
        self.px = self.colfile.omega
        self.px = n.where(self.px > self.np/2 ,
                          self.px - self.np  ,
                          self.px)*self.rlgrid
        self.py = self.colfile.sc
        self.py = n.where(self.py > self.np/2 ,
                          self.py - self.np  ,
                          self.py)*self.rlgrid
        self.pz = self.colfile.fc
        self.pz = n.where(self.pz > self.np/2 ,
                          self.pz - self.np   ,
                          self.pz)*self.rlgrid
        self.UBIALL = n.transpose(n.array( [self.px, self.py, self.pz] ))
        # self.UBIALL = self.reduce(self.UBIALL)
        print "Number of peaks found",self.px.shape[0],time.time()-start
        #        print self.UBIALL.shape, self.gv.shape

    def slow_score(self):
        import time
        start = time.time()
        scores = n.dot( self.UBIALL, n.transpose( self.gv ) )
        scores_int = n.floor( scores + 0.5).astype(n.Int)
        diff = scores - scores_int
        nv = len(self.UBIALL)
        print "scoring",nv,time.time()-start
        self.tol = 0.1
        scores = n.sqrt(n.average(diff*diff, axis = 1))
        n_ind = n.where(n.absolute(diff)< self.tol, 1 , 0)
        nind = n.sum(n_ind, axis=1)
        order = n.argsort(nind)[::-1]
        mag_v = n.sqrt( self.px * self.px +
                        self.py * self.py +
                        self.pz * self.pz )
        f = open("fft.pks","w")
        for i in range(nv):
            j = order[i]
            f.write("%d %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %7d "%(
                i,
                self.px[j],
                self.py[j],
                self.pz[j],
                mag_v[j],
                scores[j],
                self.colfile.sum_intensity[j],
                nind[j]
                ))
            for k in range(i+1,nv):
                l = order[k]
                nij = n.sum( n_ind[j] * n_ind[l] )
                f.write("%4d : %-7d "%(k,nij))
            f.write("\n")
                
        f.close()
        print diff.shape
        return diff
###################################################
## sum_sq_x = 0
## sum_sq_y = 0
## sum_coproduct = 0
## mean_x = x[1]
## mean_y = y[1]
## for i in 2 to N:
##     sweep = (i - 1.0) / i
##     delta_x = x[i] - mean_x
##     delta_y = y[i] - mean_y
##     sum_sq_x += delta_x * delta_x * sweep
##     sum_sq_y += delta_y * delta_y * sweep
##     sum_coproduct += delta_x * delta_y * sweep
##     mean_x += delta_x / i
##     mean_y += delta_y / i 
## pop_sd_x = sqrt( sum_sq_x / N )
## pop_sd_y = sqrt( sum_sq_y / N )
## cov_x_y = sum_coproduct / N
## correlation = cov_x_y / (pop_sd_x * pop_sd_y)
###################################################





        
        sm = n.zeros( (len(diff), len(diff)), n.Float)
        for k in range(len(diff)):
            i = order[k]
            sm[i,i] = n.dot(diff[i], diff[i])
        for k in range(len(diff)-1):
            i = order[k]
            for l in range(i+1, len(diff)):
                j = order[l]
                sm[i,j] = n.dot(diff[i],diff[j])/sm[i,i]/sm[j,j]
                sm[j,i] = sm[i,j]
        for i in range(len(diff)):
            sm[i,i] = 1.
        print sm[:5,:5]
        print "Scoring takes",time.time()-start
        return sm

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
    print np, mr, nsig
    go = grid(np = np , mr = mr , nsig = nsig)
    io = indexing.indexer()
    io.readgvfile(gvfile)
    go.gv_to_grid_new(io.gv )
    go.fft()
    go.props()
    go.peaksearch(open(gvfile+".patterson_pks","w"))
    im = go.read_peaks(gvfile+".patterson_pks")
    #print "Before reduction", len(go.UBIALL)
    #sys.stdout.flush()
    #ubinew = go.reduce(go.UBIALL)
    #print "After reduction", len(ubinew)
    go.UBIALL = ubinew
    go.slow_score()
    #from matplotlib.pylab import imshow, show
    #imshow(im)
    #show()
    #return im, go
    sys.exit()

    
    go.gv_to_grid_old(io.gv)
    
    diff = n.ravel(go.grid - go.old_grid)
    print "All OK if this is zero",n.add.reduce(diff)
    print "error at",n.argmax(diff),diff.max(),n.argmin(diff),diff.min()





if __name__=="__main__":
    
    import sys
    class options:
        max_res = 0.7
        np = 128
        nsig = 10
    
    options.gvfile = sys.argv[1]
    
    test(options)
