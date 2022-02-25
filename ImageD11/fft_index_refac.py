

from __future__ import print_function


"""
Classes for finding lattice vectors

We assume a very large unit cell
Compute the Patterson Function using this cell
Origin peaks show up for our crystals
These will be real space vectors
"""

import logging, time, sys
import numpy as np

from ImageD11 import labelimage, cImageD11, columnfile

def get_options(parser):
    parser.add_argument( '-n', '--ngrid',
                       action = 'store',
                       dest = 'npx',
                       type = int,
                       help = 'number of points in the fft grid [128]',
                       default = 128 )
    parser.add_argument( '-r', '--max_res',
                       action = 'store',
                       dest = 'mr',
                       type = float,
             help = 'Maximum resolution limit for fft (d-spacing) [1.0]',
                       default = 1.0)
    parser.add_argument( '-s', '--nsig',
                       action = 'store',
                       dest = 'nsig',
                       type = float,
             help = 'Number of sigma for patterson peaksearch threshold [5]',
                       default = 5)
    return parser


def refine_vector( v, gv, tol=0.25, ncycles=25, precision=1e-6 ):
    """
    Refine a (single) real space lattice vector
    Input
      v  :  real space lattice vector
      gv : scattering vectors (no errors)
      tol : tolerance for initial assignment to peak
      ncycles : how much to refine
      precision : when to stop refining in v += s when |s|/|v| < precision
    Returns
      vref : refined real space lattice vector

    Uses:
      hr = v[0]*gv[i,0] + v[0]*gv[i,1] + v[1]*gv[i,2]  1 per peak
      hi = round(hr) = calc                            1 per peak
      diff = calc - hr                                 1 per peak
      gradient = d(diff)/dv[i] = gv[j,i]               3 per peak
      matrix = gradient[i].gradient[j]                 3x3 peaks^2
      rhs = gradient . diff                            3

    Solve for shifts and iterates reducing tolerance as 1/(ncycles+1)
    Peaks that change index during refinement are removed
    Stops on : 
       no peaks reassigned
       ncycles runs out
       |shift| / |vec| < precision
    """
    assert gv.shape[1] == 3
    vref = v.copy()
    mg = (gv*gv).sum(axis=1)
    mold = None
    wt = 1./(mg + mg.max()*0.01)
    gvt = gv.T.copy()
    for i in range(ncycles):
        hr = np.dot( vref, gvt )   # npks
        hi = np.round( hr )         # npks
        diff = hi - hr              # npks
        m = abs(diff) < tol/(i+1)   #  keep or not ?
        diff = (diff*wt)[m]         # select peaks used
        grad = (gvt*wt).T[m]
        gvt = gvt[:,m]
        wt = wt[m]
        # lsq problem:
        rhs = np.dot( grad.T, diff ) 
        mat = np.dot( grad.T, grad )
        # avoid problems if singular:
        U,s,V = np.linalg.svd( mat )
        one_over_s = np.where( s/s.max() < 1e-6, 1, 1./s )
        mati = np.dot( U, np.dot(np.diag( one_over_s ), V ) )
        vshft = np.dot( mati, rhs )
        vref = vref + vshft
        if m.sum()==0:
            break
        slen = np.sqrt((vshft * vshft).sum())
        vlen = np.sqrt((vref*vref).sum())
        if vlen < precision or abs(slen/vlen) < precision:
            break
    return vref
    

class grid:
    def __init__(self, npx = 128, mr = 1.0 ,nsig = 5):
        """
        Set up the grid to use (the large unit cell)
        npx - number of points in the grid
        mr - maximum resolution limit to consider
        """
        self.npx = npx
        self.nsig = nsig
        self.grid = np.zeros((npx,npx,npx),np.float32)
        self.old_grid = np.zeros((npx,npx,npx),np.float32)
        self.cell_size = npx * mr / 2.
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
        hkl = np.round( hrkrlr ).astype(int)
        # Filter to have the peaks in asym unit
        # ... do we need to do this? What about wrapping?
        hmx = hkl.max(axis=1) <  (self.npx/2. - 2.)
        hmn = hkl.min(axis=1) >  (2. - self.npx/2.)
        my_g = np.compress( hmx & hmn, gv, axis=0 )
        # Compute hkl indices in the fft unit cell using filtered peaks
        hrkrlr = self.cell_size * my_g
        # Integer part of hkl (eg 1.0 from 1.9)
        hkl = np.floor(hrkrlr).astype(int)
        # Fractional part of indices ( eg 0.9 from 1.9)
        remain = hrkrlr - hkl
        grid = self.grid
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
            fac = 1 - abs(remain - cor)
            vol = abs(fac[:,0]*fac[:,1]*fac[:,2]).astype(np.float32)
            thkl = np.where(thkl < 0, self.npx + thkl, thkl)
            ind = thkl[:,0]*grid.shape[1]*grid.shape[2] + \
                  thkl[:,1]*grid.shape[1] + \
                  thkl[:,2]

            cImageD11.put_incr( flatgrid , ind.astype(np.intp), vol )
        logging.info("Grid filling loop takes "+str(time.time()-start)+" /s")



    def fft(self):
        """ Compute the Patterson """
        start = time.time()
        self.patty = np.ascontiguousarray(abs(np.fft.fftn(self.grid)),
                                          np.float32)
        logging.info("Time for fft "+str(time.time()-start))
        self.origin = self.patty[0,0,0]
        logging.info("Patterson origin height is :"+str(self.origin))

    def props(self):
        """ Print some properties of the Patterson """
        logging.info("Patterson info "+str(self.patty.shape)
                     +str(type(self.patty)))
        p = np.ravel(self.patty)
        m = np.mean(p)
        logging.info("Average: %f"%(m))
        p2 = p*p
        self.mean = m
        v = np.sqrt( (np.sum(p2) - m*m*len(p) ) /(len(p)-1) )
        logging.info("Sigma: %f"%(v))
        self.sigma = v

    def peaksearch(self, peaksfile):
        """ Peaksearch in the Patterson """
        lio = labelimage.labelimage( (self.npx, self.npx),
                                     peaksfile )
        logging.info("Peaksearching at %f sigma"%(self.nsig))
        thresh = self.mean + self.nsig  * self.sigma
        for i in range(self.npx):
            lio.peaksearch(self.patty[i],
                           thresh,
                           i)
            lio.mergelast()
        lio.finalise()

    def pv(self, v):
        """ print vector """
        return ("%8.4f "*3)%tuple(v)

    def reduce(self, vecs):
        raise Exception("You want lattice_reduction instead")

    def read_peaks(self, peaksfile):
        """ Read in the peaks from a peaksearch """
        start = time.time()
        colf = columnfile.columnfile(peaksfile)
        logging.info("reading file %f/s"%(time.time()-start))
        # hmm - is this the right way around?
        self.rlgrid = 1.0*self.cell_size/self.npx
        self.px = colf.omega
        self.px = np.where(self.px > self.npx/2 ,
                           self.px - self.npx  ,
                           self.px)*self.rlgrid
        self.py = colf.sc
        self.py = np.where(self.py > self.npx/2 ,
                           self.py - self.npx  ,
                           self.py)*self.rlgrid
        self.pz = colf.fc
        self.pz = np.where(self.pz > self.npx/2 ,
                           self.pz - self.npx   ,
                           self.pz)*self.rlgrid
        self.UBIALL = np.array( [self.px, self.py, self.pz] ).T
        logging.info("Number of peaks found %d  %f/s, now fit some"%(
                     self.px.shape[0],time.time()-start))
        for i in range(colf.nrows):
            print(".",end="")
            self.UBIALL[i] = refine_vector( self.UBIALL[i], self.gv )
            
        logging.info("Fitting vectors %f /s"%(time.time()-start))
        self.colfile = colf
        
    def slow_score(self):
        logging.info("running slow_score")
        import time
        start = time.time()
        scores = np.dot( self.UBIALL, np.transpose( self.gv ) )
        scores_int = np.floor( scores + 0.5).astype(int)
        diff = scores - scores_int
        nv = len(self.UBIALL)
        print("scoring",nv,time.time()-start)
        self.tol = 0.1
        scores = np.sqrt(np.average(diff*diff, axis = 1))
        n_ind = np.where(np.absolute(diff)< self.tol, 1 , 0)
        nind = np.sum(n_ind, axis=1)
        order = np.argsort(nind)[::-1]
        mag_v = np.sqrt( self.px * self.px +
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
                nij = np.sum( n_ind[j] * n_ind[l] )
                f.write("%4d : %-7d "%(k,nij))
            f.write("\n")

        f.close()
        print(diff.shape)
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






        # sm = np.zeros( (len(diff), len(diff)), float)
        # for k in range(len(diff)):
        #     i = order[k]
        #     sm[i,i] = np.dot(diff[i], diff[i])
        # for k in range(len(diff)-1):
        #     i = order[k]
        #     for l in range(i+1, len(diff)):
        #         j = order[l]
        #         sm[i,j] = np.dot(diff[i],diff[j])/sm[i,i]/sm[j,j]
        #         sm[j,i] = sm[i,j]
        # for i in range(len(diff)):
        #     sm[i,i] = 1.
        # print(sm[:5,:5])
        # print("Scoring takes",time.time()-start)
        # return sm

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
    npx = options.npx
    nsig = options.nsig
    from ImageD11 import indexing
    print(npx, mr, nsig)
    go = grid(npx = npx , mr = mr , nsig = nsig)
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
    #go.UBIALL = ubinew
    go.slow_score()
    #from matplotlib.pylab import imshow, show
    #imshow(im)
    #show()
    #return im, go
    sys.exit()


    go.gv_to_grid_old(io.gv)

    diff = np.ravel(go.grid - go.old_grid)
    print("All OK if this is zero",np.add.reduce(diff))
    print("error at",np.argmax(diff),diff.max(),np.argmin(diff),diff.min())





if __name__=="__main__":

    import sys
    class options:
        max_res = 0.7
        npx = 128
        nsig = 10

    options.gvfile = sys.argv[1]

    test(options)
