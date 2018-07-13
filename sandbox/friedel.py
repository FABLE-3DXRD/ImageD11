
from __future__ import print_function, division
from ImageD11 import columnfile, parameters, transform

import numpy as np
import sys, time

class binned_array(object):
    def __init__(self, ar, nbins=255 ):
        """
        Digitise an array onto nbins precision
        """
        self.ar = np.asarray( ar )
        self.nbins = nbins
        lo = self.ar.min()
        hi = self.ar.max()
        self.step = (hi - lo) / (nbins - 1)
        self.low = lo - self.step/2
        self.ibin = self.bin( self.ar )
        self.bc = np.bincount( self.ibin )
        cs = np.cumsum( self.bc )
        self.bin_inds = np.concatenate( ([0,], cs) )
        
        
    def bin(self, x):
        """Assign x values into bins """
        return np.round( (x - self.low) / self.step ).astype( int )
        


def absanglediff( a1 , a2 ):
    ra1 = np.radians( a1 )
    ra2 = np.radians( a2 )
    er = np.arctan2( np.sin( ra1 - ra2 ), np.cos( ra1 - ra2 ) )
    return abs(np.degrees( er ))

def findpairs_2d( cf, dangle, omestep, fac = 1.5):
    """
    Search for Friedel pairs (and harmonics) as parallel (or nearly)
    vectors

    Steps:
      1) Compute g-vectors 
      2) Normalise them to get directions
      3) locate them onto the 2D surface of directions (sphere in 3d)
      4) "flip sign" to get pairs of (g, -g) to overlap
      5) foreach 2D neighborhood look for |g|, 2|g|, 3|g|, etc.

    cf contains gx, gy, gz already
    dangle is the twotheta and eta precision
    omestep is the omega step
    """
    # Normalise the g-vectors to direction cosines
    nve = ( np.array((cf.gx, cf.gy, cf.gz)) / cf.ds )
    # We cut on y > 0 (could equally be x... )
    nve = np.where( nve[1] > 0, nve, -nve )
    # World map geometry ...
    #     angle to z-axis   n.z = |n||z|cos( latitude )
    az = np.degrees( np.arccos( nve[2] ) )
    #     angle to xy-axis  x/y = tan( longitude )
    ay = np.degrees( np.arctan2( nve[0], nve[1]) )
    # Make bins reflect the estimated precision (not ideal for now)
    # ... guess latitude error is also eta error, make bins accordingly
    nz = int( 180/(dangle*fac) )
    zstep = 180.0/nz
    zedges = np.arange( az.min() - zstep/2, az.max() + zstep, zstep )
    # ... and longitude errors is roughly omega error
    ny = int( 180/(omestep*fac) )
    ystep = 180.0/ny
    yedges = np.arange( ay.min() - ystep/2, ay.max() + ystep, ystep )
    #
    iz = (az - zedges[0])*( len(zedges) - 1 ) / ( zedges[-1] - zedges[0]  )
    iy = (ay - yedges[0])*( len(yedges) - 1 ) / ( yedges[-1] - yedges[0]  )
    iz = np.floor(iz).astype( int )
    iy = np.floor(iy).astype( int )

    # Collect all the spots into a histogram in 2D
    start = time.time()
    H, zHedges, yHedges = np.histogram2d( az, ay, ( zedges, yedges ) )
    H = H.astype( np.int )
    NZ, NY = H.shape
    if __debug__: # checkng how binning worked
        print("H2D:",time.time()-start)
        assert (zHedges == zedges).all()
        assert (yHedges == yedges).all()
        #     v   - minbin
        # Find which bin the sqpots are in
        assert (iz == np.digitize( az, zedges )-1).all(),(iz,np.digitize( az, zedges ))
        assert (iy == np.digitize( ay, yedges )-1).all()
        # H is ordered as    H.shape[0]   H.shape[1]
        #print(H.shape   ,zedges.shape, yedges.shape )
        #print( iz.min(), iz.max(), iy.min(), iy.max() )
        #b = np.argmax(iz)
        #print(az[b],iz[b],zedges[-1])
        #b = np.argmax(iy)
        #print(ay[b],iy[b],yedges[-1])
    #
    # Find the storage locations to sort the peaks onto the 2D histogram
    pointers = np.cumsum( H.ravel(), dtype = np.int64 )
    # Not sure we have a compiled code for this ?
    # loop over peaks, make location be given by pointer
    #  ... increment pointer to make space for next
    # first bin contains 2 peaks. Place these at 0,1. Next starts at 2.
    ptemp = np.zeros(pointers.shape, np.int )
    ptemp[0] = 0
    ptemp[1:] = pointers[:-1]
    # To index these pointers
    pidx = iz*NY + iy
    # resulting sort
    order = np.zeros( cf.nrows, int )
    start = time.time()
    for i in range(cf.nrows):
        p = pidx[i]           # p = which pixel is this peak
        order[i] = ptemp[p]   # output position
        ptemp[p] += 1
    # sort the columnfile and histogram indices
    indices = np.arange( cf.nrows, dtype=np.int )
    indices.put( order, indices.copy())
    cf.reorder( indices )
    iy = iy[indices]
    iz = iz[indices]
    pidx = iz*NY + iy
    if __debug__:
        print('splatting', time.time()-start )        
        assert len(np.unique(order)) == cf.nrows
        assert order.min() == 0
        assert order.max() == cf.nrows-1
    allpointers = np.concatenate( ( [0,], pointers ) )
    pimage = np.reshape( pointers, H.shape )
    # Now the peaks are sorted according to which direction they have
    # each peak has an i,j co-ordinate showing where it is on the world map
    # we want to select peaks for each pixel
    izimage, iyimage = np.mgrid[ 0:H.shape[0], 0:H.shape[1] ]
#    print(iyimage.shape,NY,iyimage.max(),H.shape)
#    print(izimage.shape,NZ,izimage.max(),H.shape)
    mask = H > 0
    k=0
    hits = []
    for iiy, iiz in zip( iyimage[mask], izimage[mask] ) :
        # print("Pixel",iiy, iiz)
        # In the y direction we will wrap around
        iys = (iiy-1)%NY , iiy, (iiy+1)%NY
        # On the z direction we do not wrap as this is a pole (fable geometry)
        if iiz > 0:
            if iiz < H.shape[0]-2: # can do iz - 1
                izs = iiz - 1, iiz, iiz + 1
            else:
                izs = iiz - 1, iiz
        else:
                izs =            iiz, iiz + 1
        npk = 0
        newhits = []
        for this_z in izs:
            for this_y in iys:
                npk += H[this_z, this_y]
                # location of these peaks is:
                plo = allpointers[ this_z * NY + this_y ]
                phi = allpointers[ this_z * NY + this_y + 1]
                newhits += range(plo, phi )
                continue
            #                print(" ", this_z, this_y, H[this_z, this_y])
            #               for k in range(plo, phi):
            #                  print("    ",k,cf.ds[k], cf.omega[k], cf.eta[k], cf.gx[k], cf.gy[k], cf.gz[k] )
        if len(newhits)>1:
            hits.append(newhits)
#    print([len(h) for h in hits] )
    print(len(hits))
        
    
    
    
    if __debug__:
        # assert order is a permutation ...
        assert len( np.unique(order)  ) == cf.nrows, (len( np.unique(order)  ), cf.nrows)
        assert order.min() == 0
        assert order.max() == cf.nrows-1
    


    
    if 1:
        import pylab as pl
        pl.figure(1)
        pl.imshow( H, interpolation='nearest',
                   origin='low',
                   extent=( yedges[0], yedges[-1], zedges[0], zedges[-1] )
                   )
        pl.plot( ay, az, "+")
        pl.show()

    

def findpairs( cf, dangle, omestep, fac = 1.5 ):
    """
    cf = columnfile
    dangle = average angle error on detector (from sample size)
    omestep = acceptable omega error
    """
    cf.sortby( 'ds' )
    allinds = np.arange( cf.nrows, dtype=np.int )
    paired = allinds < 0
    allpairs = []
    
    dang = dangle * fac
    dome = omestep * fac

    omega = cf.omega * cf.parameters.get( "omegasign" ) 
    
    # Find eta/omega for friedel pairs, and self ?
    gve = np.array( (cf.gx,cf.gy,cf.gz))
    tth, (fe1, fe2), (fo1, fo2) = transform.uncompute_g_vectors(
        -gve,
        cf.parameters.get("wavelength"),
        wedge = cf.parameters.get("wedge"),
        chi   = cf.parameters.get("chi") )
    tth, (se1, se2), (so1, so2) = transform.uncompute_g_vectors(
        gve,
        cf.parameters.get("wavelength"),
        wedge = cf.parameters.get("wedge"),
        chi   = cf.parameters.get("chi") )
    
    cf.sortby( "tth" )
    # Outer loop on tth 
    itth = binned_array( cf.tth, 1000 )

    lo_tth_bin = np.clip( itth.bin( cf.tth - dang ) , 1, itth.nbins )
    hi_tth_bin = np.clip( itth.bin( cf.tth + dang ) + 1, 0, itth.nbins )
    for i in range(cf.nrows):
        lc = itth.bin_inds[lo_tth_bin[i]]
        hc = itth.bin_inds[hi_tth_bin[i]]
        if (hc - lc) > 0:
            # check against eta and omega
            e = cf.eta[lc:hc]
            o = omega[lc:hc]
            m1 = (absanglediff( e, fe1[i]) < dang ) & \
                 (absanglediff( o, fo1[i]) < dome )
            m2 = (absanglediff( e, fe2[i]) < dang ) & \
                 (absanglediff( o, fo2[i]) < dome )
            s1 = (absanglediff( e, se1[i]) < dang ) & \
                 (absanglediff( o, so1[i]) < dome )
            s2 = (absanglediff( e, se2[i]) < dang ) & \
                 (absanglediff( o, so2[i]) < dome )
            #print(  (absanglediff( o, so1[i]) < dome) & (absanglediff( e, se1[i]) < dang ))
            #print(  m1 | m2 | s1 | s2 )
            pairs = allinds[lc:hc][ m1 | m2 | s1 | s2 ]
            if len(pairs) < 1:
                raise Exception("Missing self!")
            if len(pairs) == 1:
                continue
                pass
            print("#",i,lc,hc,hc - lc,len(pairs),list(pairs))
            #print(cf.tth[lc],cf.tth[hc],cf.tth[i])
            allpairs.append( tuple( [i,]+list(pairs) ) )
            continue
        
            print("%-6d     % 9.4f % 7.2f % 7.2f %9d %9d"%(
                i,cf.tth[i], cf.eta[i], cf.omega[i],
                cf.sum_intensity[i], cf.IMax_int[i]),
                  cf.s_raw[i], cf.f_raw[i])

            for k in pairs:
                print("    %-6d % 9.4f % 7.2f % 7.2f %9d %9d"%(
                    k, cf.tth[k],cf.eta[k], cf.omega[k],
                    cf.sum_intensity[k], cf.IMax_int[k]),
                      cf.s_raw[k], cf.f_raw[k])
                sol, _, _ = fitpair( cf, (i, k), cf.parameters )
                origin, vec, w = fitpair( cf, (i, k), cf.parameters, sol0=sol )
                print("          ( % 6.2f % 6.2f % 6.2f ) +/- ( % 6.2f % 6.2f % 6.2f )"%(
                    origin[0], origin[1], origin[2],
                    vec[0], vec[1], vec[2]), w )
                # A pair is two peaks and has an origin and direction vector
                
            print()

def unit(vec):
    return vec/np.sqrt(np.dot(vec,vec))


def skewlines( o1, v1, o2, v2 ):
    # distance between two skew lines
    n = np.cross(v1, v2)
    mag = np.sqrt( np.dot( n, n ) )
    dv = o1 - o2
    if mag > 0:
        un = n / mag
        dv = np.dot( un, dv )
    else:
        # lines are parallel, distance is constant
        pass
    return np.sqrt( np.dot( dv, dv ) )
            

def fitpair( cf, pair, pars, sol0=None ):
    # pair is a pair of peaks with g, -g
    i, j = pair
    if sol0 is None:
        tx0 = 0.
        ty0 = 0.
        tz0 = 0.
        t0 = np.zeros(3)
    else:
        tx0 , ty0, tz0 = sol0
        t0 = sol0
    xyz = np.array(( ( cf.xl[i], cf.yl[i], cf.zl[i]),
                     ( cf.xl[j], cf.yl[j], cf.zl[j]) ) ).T
    om = np.array((cf.omega[i], cf.omega[j]))
    tth0, eta0 = transform.compute_tth_eta_from_xyz(
        xyz, om, t_x=tx0, t_y=ty0, t_z=tz0,
        wedge=pars.get('wedge'), chi=pars.get('chi'))
    g0 = transform.compute_g_vectors(
        tth0, eta0, om, pars.get("wavelength"),
        wedge=pars.get('wedge'), chi=pars.get('chi'))
    tth1, eta1 = transform.compute_tth_eta_from_xyz(
        xyz, om, t_x=tx0+1, t_y=ty0, t_z=tz0,
        wedge=pars.get('wedge'), chi=pars.get('chi'))
    g1 = transform.compute_g_vectors(
        tth1, eta1, om, pars.get("wavelength"),
        wedge=pars.get('wedge'), chi=pars.get('chi'))
    dg_dtx = g1 - g0
    tth1, eta1 = transform.compute_tth_eta_from_xyz(
        xyz, om, t_x=tx0, t_y=ty0+1, t_z=tz0,
        wedge=pars.get('wedge'), chi=pars.get('chi'))
    g1 = transform.compute_g_vectors(
        tth1, eta1, om, pars.get("wavelength"),
        wedge=pars.get('wedge'), chi=pars.get('chi'))
    dg_dty = g1 - g0
    tth1, eta1 = transform.compute_tth_eta_from_xyz(
        xyz, om, t_x=tx0, t_y=ty0, t_z=tz0+1,
        wedge=pars.get('wedge'), chi=pars.get('chi'))
    g1 = transform.compute_g_vectors(
        tth1, eta1, om, pars.get("wavelength"),
        wedge=pars.get('wedge'), chi=pars.get('chi'))
    dg_dtz = g1 - g0
    #      kx = kx          ,  ky = -ky        , kz = -kz
    gerr = g0[:,0] + g0[:,1]
    de_dtx = dg_dtx[:,0] + dg_dtx[:,1]
    de_dty = dg_dty[:,0] + dg_dty[:,1]
    de_dtz = dg_dtz[:,0] + dg_dtz[:,1]
    grad = ( de_dtx, de_dty, de_dtz )
    rhs = np.zeros(3)
    mat = np.zeros((3,3))
    if 0 and sol0 is not None:
        print(tth0)
        print(eta0)
        print(g0)
        print(gerr)
    for i in range(3):
        rhs[i] = np.dot( grad[i], gerr )
        for j in range(i,3):
            mat[i,j] = mat[j,i] = np.dot( grad[i], grad[j] )
    w , v = np.linalg.eigh( mat )
    ### mat   = v.T .  w .  v.T-1
    ### mat   = v.T .  w .  v      since v is like a U matrix (orthonormal)
    ### mat-1 = v.T . 1/w . v
    ww = np.eye(3)
    ans = np.zeros(3)
    sol = t0
    # Take the two largest eigenvalues, eigh says sorted
    assert w[0] < w[1] < w[2], "check eigh returns sorted?"
    for i in range(2,0,-1):
        ww[i,i] = 1/w[i]
        ans += np.dot(v, ww[i]*np.dot(v.T, rhs) )
        sol -= ans/w[i]
        print(i,ans,sol,w[i])
    ans = -unit(np.dot(v, ww[2]*np.dot(v.T, rhs) ))
    sol -= ans*w[2]
    print(ans,sol,w)
    imat = np.dot(np.dot(v, ww),v.T)
    #print('ima',np.dot(np.linalg.inv( mat ) , rhs ))
    return sol, ans, (w[0]/w[1],w[1]/w[2])


def guess_omega_step( colf ):
    if "Max_o" in colf.titles and "Min_o" in colf.titles:
        mask = colf.Max_o == colf.Min_o
        if mask.sum() == 0:
            return 0
        uniqomega =  np.array(list(set(colf.omega[mask])))
    else:
        uniqomega =  np.array(list(set(colf.omega)))
    uniqomega.sort()
    domegauniq = uniqomega[1:] - uniqomega[:-1]
    # mode would be better statistic ?
    med = np.median( domegauniq )
    return med
    
    
    
        

def main():
    colf = columnfile.columnfile( sys.argv[1] )
    print("Read",colf.nrows,"peaks")
    p = parameters.read_par_file( sys.argv[2] )
    # input a sample radius to use for computing ranges
    radius = int( sys.argv[3] )
    omegastep = guess_omega_step( colf )
    print("Guess omega step as", omegastep)
    colf.updateGeometry( p )

    
    for i,t in enumerate(colf.titles):
        print("\n", t, end=" ")
        for j in (0,1,2, colf.nrows//2,colf.nrows//2+1,
                  colf.nrows//2+2, colf.nrows-3, colf.nrows-2, colf.nrows-1):
            print( colf.bigarray[i,j], end=" ")
    print()

    modXL = np.sqrt(colf.xl*colf.xl + colf.yl*colf.yl + colf.zl*colf.zl)
    modXLmean = modXL.mean()
    dangle = np.degrees( np.arctan2( radius , modXLmean ) )
    print ("Angle precision is about ",dangle,"for average distance",modXLmean,
           "and sample size ",radius)

    findpairs_2d( colf, dangle, omegastep )
    1/0
    
    findpairs( colf, dangle , omestep=omegastep )


if __name__=="__main__":
        main()















