
import numpy, columnfile
from ImageD11 import grain, parameters, transform

from utils_numpy import outerprod_res, vec3mod
class axis(object):
    def __init__(self, n):
        """
        n is the axis normal, eg, direction of the rotation axis
        """
        arn = numpy.asarray( n )
        assert arn.shape == (3,)
        # check normalised
        assert abs(numpy.dot( n, n ) - 1) < 1e-6
        self.n = arn

    def rotate_one_v( self, r , cost, sint, result=None, tmp=None, back=False):
        """
        r = vector, one to many angles
        cost, sint = cosine and sin of angles
        rotates a single vector to multiple angles
        if back then invert sign (sin(-x)==-sin(x), cos(-x)==cos(x))
        http://mathworld.wolfram.com/RotationFormula.html
        r' = r cos(t) + n(n.r)(1-cos(t)) + rxn sin(t)
        """
        ar = numpy.asarray( r )
        assert ar.shape == (3,)
        if result is None:
            result = numpy.zeros( (3, len(cost) ), numpy.float32)
        if tmp is None:
            tmp =  numpy.zeros( (3, len(cost) ), numpy.float32)
        # 3x1     3x1      1         3x1   3x1
        nndotr = self.n * numpy.dot(self.n, ar)
        # 3x1             3x1     3x1  
        rxn = numpy.cross( ar, self.n )
        # inplace to skip mallocs
        outerprod_res( ar, cost, result )
        outerprod_res( nndotr, (1-cost), tmp)
        numpy.add( result,  tmp, result )
        outerprod_res( rxn, sint, tmp )
        # reverse rotation, sin-> -sin
        if back:
            numpy.subtract( result,  tmp, result )
        else:
            numpy.add( result, tmp, result )
        return result


    def rotate_vs( self, vs, cost, sint, result, tmp, back=False):
        """
        rotate many vectors to many different angles
        
        R = I3cos(t) + ww^T (1-cos(t)) - W sin(t)
        t = angle
        W = vector with hat operator = [ 0  -wz  wy 
                                         wz  0  -wx
                                        -wy  wx  0  ]
        """
        if 1:
            # Rotate the vectors into sample frame
            assert (self.n == (0,0,1)).all()
            numpy.multiply(  peakdata.cosomega , vs[0], result[0] )
            numpy.multiply(  peakdata.sinomega , vs[1], tmp )
            if back:
                numpy.subtract( result[0], tmp, result[0] )
            else:
                numpy.add( result[0], tmp, result[0] )
            numpy.multiply(  peakdata.sinomega, vs[0], tmp )
            numpy.multiply(  peakdata.cosomega, vs[1], result[1] )
            if back:
                numpy.add(  result[1] , tmp, result[1] )
            else:
                numpy.subtract(  result[1] , tmp, result[1] )
            result[2]=  vs[2]
        else:
            dx, dy, dz = self.direction
            e = n.array([self.direction]*3)
            w = n.transpose(n.array( [ [  0, -dz,  dy]  , 
                                       [ dz,   0, -dx]  , 
                                       [-dy,  dx,   0] ], n.float))
            # self.matrix = n.identity(3, n.float)*ct - st * w  + \
            #               (1 - ct)*e*n.transpose(e)
            #
            # identity term:
            result[i] = cost
            # 1-cost terms:
            omc = 1 - cost
            for i in range(3):
                for j in range(3):
                    numpy.add(w[i,j]*w[i,j]*omc, result[i], result[i] )
                for j in range(i+1,3):
                    numpy.add(w[i,j]*w[i,j]*omc, result[i], result[i] )
                
            raise Exception("FIXME")     



 
def update_det_pars( peakdata, pars ):
    """
    Update the xl, yl, zl columns
    peakdata = columnfile
    pars = ImageD11 parameters object
    """
    print "Update xl, yl, zl for current parameters"
    for label in ["xl","yl","zl"]:
        if label not in peakdata.titles:
            peakdata.titles.append(label)
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
    peakdata.xl += peakdata.samtx * peakdata.cosomega - \
                   peakdata.samty * peakdata.sinomega
    peakdata.yl += peakdata.samtx * peakdata.sinomega + \
                   peakdata.samty * peakdata.cosomega
    peakdata.zl += peakdata.samtz
    print "lab x shape",peakdata.xl.shape
    return peakdata
    




def potentialpeaks( grainlist, peakdata, pars, tol ):
    """
    For each grain, see which peaks are potentially interesting
    grainlist = python list of grain objects
    peakdata  = columnfile object with xl,yl,zl,cosomega,sinomega
    pars      = ImageD11 parameters object
    tol       = Indexing integer tolerance
    """
    assign = {} # Dict to hold output
    print "potential peak finding"
    # Pre-allocate necessary work arrays
    # Logically these could be part of peakdata, but that is bad
    # for doing something in parallel
    # ... note ... npks is of the order 1e6
    npks = len(peakdata.omega)
    indices = numpy.arange( npks, dtype = numpy.intp )
    fac = numpy.zeros( npks, numpy.float32)
    tmp = numpy.zeros( npks , numpy.float32)
    g = numpy.zeros( (3, npks), numpy.float32)
    t = numpy.zeros( (3, npks), numpy.float32)
    hkli = numpy.zeros( g.shape, numpy.int32)
    hklr = numpy.zeros( g.shape, numpy.float32)
    drlv = numpy.zeros( g.shape, numpy.float32)
    i = 0
    for gc in grainlist:
        # Find the scattering vectors in lab frame:
        tx, ty, tz = gc.translation
        omega_axis.rotate_one_v( [tx,ty,tz],
                                 peakdata.cosomega,
                                 peakdata.sinomega,
                                 result = t,
                                 tmp = g, 
                                 back=True )
        numpy.subtract( peakdata.xl , t[0], t[0] )
        numpy.subtract( peakdata.yl , t[1], t[1] )
        numpy.subtract( peakdata.zl , t[2], t[2] )
        # unit vectors
        vec3mod( t, fac, tmp)
        numpy.divide( t, fac, t )
        # t is now the k vector
        beam_direction = (1,0,0)
        for j in range(3):
            if beam_direction[j] > 0:
                numpy.subtract( t[j,:], beam_direction[j], t[j,:] )
        #
        # Rotate the vectors into sample frame
        omega_axis.rotate_vs( t,
                              peakdata.cosomega,
                              peakdata.sinomega,
                              g, tmp, back=False)
        # Index the vectors using h,k,l indices
        ubioy = gc.ubi / float(pars.parameters['wavelength'])
        hklr = numpy.dot( ubioy, g )
        # Score to match up how good they are
        numpy.floor( hklr + 0.5, hkli )
        drlv = (hklr - hkli)
        drlv = numpy.sqrt( (drlv*drlv).sum(axis=0))
        inds = numpy.compress( drlv < tol, indices )
        assign[i] = inds, numpy.take(hkli, inds, axis=1), numpy.take(drlv, 
                                                                    inds)
        print i, len(inds)
        i+=1
    peakdata.drlv = drlv
    return assign





def pick_peaks( peakdata, inds, priority="drlv"):
    """
    Given drlv and whatever you like from peakdata, choose your
    favourites.
    We will go in the order of drlv, the getting uniq omega,
    samtx, samty and samtz

    Another option would be to pick the more intense peaks
    """
    priority_array = numpy.take( getattr( peakdata, priority), inds)
    order = numpy.argsort( priority_array )
    assert len(priority_array) == len(inds)
    j = order[0]
    uniq = [ inds[j] ]
    for i in order[1:]:
        j = inds[i]
        # Decide if we have already see this peak in uniq
        attrs = ['samtx', 'samty', 'samtz', 'omega' ]
        tols  = numpy.array([ 0.1,     0.1,     0.1,     10 ])
        new = True
        for p in uniq:
            diffs = [abs(getattr(peakdata,a)[p] - getattr(peakdata,a)[j])
                     for a,t  in zip(attrs, tols)]
            if ((tols < diffs) == False).all():
                new = False
#            print diffs, new, 
        if new:
            uniq.append( inds[i] )
    return uniq




def grain_choose_single_peaks( peakdata, agrain, assign, magic = 32,
                               proirity = 'sum_intensity' ):
    """
    peakdata = all the data
    agrain = one single grain
    assign = current list of indices and hkl indexing

    Choose the best peak for each possible hkl of this grain
    Assuming the grain has got a list of potential peaks, we want to
    whittle this down so it only has one observation per calculated
    point. It is for the case of a grain that indexes into a cloud
    of data
    """
    inds, hkls, drlv = assign
    numpy.put(peakdata.drlv, inds, drlv)
    # print inds.shape, hkls.shape
    assert magic > abs(hkls.ravel()).max()
    key = (hkls[0]*magic + hkls[1])*magic + hkls[2]
    order = numpy.argsort(key)
    # Several cases:
    #   Single observation of this hkl, keep it
    #   Multiple observations of this hkl:
    #       group into "equivalence", eg: samtx, omega etc
    #       pick best in equivalent group

    uniq = []
    j = 0
    doubles = []
    
    while j < len(drlv)-1:
        i = j + 1
        pks = [ order[j] ] # local indexing for this grain
        kj = key[ order[j] ]
        while i<len(drlv) and key[order[i]] == kj:
            pks.append( order[i])
            i = i + 1
        j = i + 1
        if len(pks) == 1:
            # Single observation for this hkl
            uniq.append( pks[0] )
            continue
        else:
            allpks = inds[pks] # global selection into dataset
            #uniq += pick_peaks( peakdata , allpks, priority = 'sum_intensity')
            uniq += pick_peaks( peakdata , allpks, priority = 'drlv')
    # print
    # s = set(key)
    print len(uniq)
    return     inds[uniq], hkls[uniq], drlv[uniq]




  
def refine_sparse( peakdata, gr, indices, hkls ):
    hcalc = gr.ubi
    pass



if __name__=="__main__":
    import sys
    if 0:
     print """fuckit - you need a proper grain object that
     just copies all the data for the peak it is interested in"""
     sys.exit()

    omega_axis = axis( [ 0, 0, 1] ) # z
    chi_axis   = axis( [ 1, 0, 0] ) # x
    wedge_axis = axis( [ 0, 1, 0] ) # y
    


    # Data object
    peakdata = columnfile.colfile_from_hdf( sys.argv[1] )
    if not hasattr( peakdata, 'sum_intensity'):
        peakdata.sum_intensity = peakdata.npixels * peakdata.avg_intensity

    # Grain collection object
    grainlist = grain.read_grain_file( sys.argv[2] )

    if sys.argv[2][0:2] == "x_":
        items = sys.argv[2].split("_")
        tx = float( items[1] )*1500 - 4500
        ty = float( items[3] )*1500 - 4500
        tz = float( items[5].split(".")[0] )*500 -750
        for g in grainlist:
            #print tx,ty,tz,g.translation
            g.translation += [tx,ty,tz]
            #print tx,ty,tz,g.translation
        print "Sand translations added"
    # Detector parameters
    pars = parameters.parameters()
    pars.loadparameters( sys.argv[3] )

    peakdata = update_det_pars( peakdata, pars )


    tol = 0.05
    
    assignments = potentialpeaks( grainlist, peakdata, pars, tol )

    print "Choose uniq within potentials"
    narrowassign = {}
    for i in range(len(grainlist)):
        narrowassign[i] =  grain_choose_single_peaks( peakdata,
                                                      grainlist[i],
                                                      assignments[i])

    
    

