

from __future__ import print_function
from ImageD11 import peakmerge, indexing, transformer, cImageD11
from ImageD11 import grain, unitcell, refinegrains, sym_u
import xfab.tools
import sys, os, numpy as np, time, random
import multiprocessing, traceback
from multiprocessing import Pool
from multiprocessing import Queue as PQueue

if sys.version_info[0] < 3:
    import Queue # for exception
else:
    import queue as Queue

if "win" in sys.platform:
    nulfile = "NUL"
else:
    nulfile = "/dev/null"






def domap(  pars,
            colfile,
            grains,
            gridpars):
    """
    mapping function - does what makemap.py does, but in a function
    """
    if 'FITPOS' not in gridpars:
        gridpars['FITPOS']=True
        
    OmSlop = gridpars['OMEGAFLOAT']
    OmFloat= OmSlop > 0
    #
    ss = sys.stdout # turns off printing
    if gridpars['NUL']:
        NUL = open(nulfile,"w")
        sys.stdout = NUL
    for tol in gridpars['TOLSEQ']:
        o = refinegrains.refinegrains( OmFloat = OmFloat, OmSlop = OmSlop,
                                       tolerance = tol,
                                       intensity_tth_range = (0,180),
                                       )
        o.parameterobj = pars
        # o.loadfiltered ...
        o.scannames = ["internal"]
        o.scantitles = colfile.titles
        o.scandata["internal"] = colfile
        o.tolerance = tol
        # o.readubis( grainsfile )
        for i, g in enumerate(grains):
            name = i
            o.grainnames.append(i)
            o.ubisread[name] = g.ubi
            o.translationsread[name] = g.translation
        if gridpars['SYMMETRY'] != "triclinic":
            o.makeuniq( gridpars['SYMMETRY'] )
        o.generate_grains()
        if gridpars['FITPOS']:
            o.refinepositions()
        else:
            o.assignlabels()
            for key in o.grains.keys():
                g = o.grains[key]
                g.set_ubi( o.refine( g.ubi, quiet=False ) )
        # This fills in the uniq for each grain
        o.savegrains( nulfile, sort_npks = False)
        if 'NUNIQ' in gridpars:
            keep = lambda g: g.nuniq > gridpars['NUNIQ'] and g.npks > gridpars['NPKS']
        else:
            keep = lambda g: g.npks > gridpars['NPKS']
        gl = [ g for g in o.grains.values() if keep(g) ]
        if len(gl) == 0:
            break
        grains = gl
    if gridpars['NUL']:
        sys.stdout = ss
    return gl



def doindex( gve, x, y, z, w, gridpars):
    """
    Try to index some g-vectors
    """
    NPKS = gridpars['NPKS']
    UC =      gridpars['UC']
    TOLSEQ = gridpars['TOLSEQ']
    COSTOL = gridpars[ 'COSTOL']
    DSTOL  = gridpars[ 'DSTOL']
    if "2RFIT" in gridpars:
        DOFIT = gridpars[ '2RFIT' ]
    else:
        DOFIT = False
    ss = sys.stdout # turns off printing
    if gridpars['NUL']:
        NUL = open(nulfile,"w")
        sys.stdout = NUL
    myindexer = indexing.indexer(
        wavelength = w,
        unitcell = UC,
        gv = gve.T
        )
    # added in indexer.__init__
    #myindexer.ds = np.sqrt( (gve * gve).sum(axis=0) )
    #myindexer.ga = np.zeros(len(myindexer.ds),int)-1 # Grain assignments
    for ring1 in gridpars['RING1']:
        for ring2 in gridpars['RING2']:
            myindexer.parameterobj.set_parameters(  {
                'cosine_tol' : COSTOL,
                'ds_tol': DSTOL, 
                'minpks': NPKS, 
                'max_grains': 1000, 
                'hkl_tol': TOLSEQ[0], 
                'ring_1': ring1,
                'ring_2': ring2
                } )
            myindexer.loadpars( )
            myindexer.assigntorings( )
            try:
                myindexer.find( )
                myindexer.scorethem( fitb4 = DOFIT )
            except:
                pass
    # filter out crap
    vol = 1/np.linalg.det( UC.B )
    grains = [ grain.grain(ubi, [x,y,z]) for ubi in myindexer.ubis
               if np.linalg.det(ubi) > vol*0.5 ]
    if gridpars['NUL']:
        sys.stdout = ss
    return grains

def test_many_points( args ):
    """
    Grid index - loop over points
    Places the results in a multiprocessing Queue
    """
    colfile, parameters, translations, gridpars = args
    s = "Hello from %s %d"%(multiprocessing.current_process().name ,os.getpid())
    s += " %d to do"%(len(translations))
    s += "%s %s"%(colfile, parameters)
    print(s)
    mytransformer = transformer.transformer()
    mytransformer.loadfiltered( colfile )
    mytransformer.loadfileparameters(  parameters )
    w =  mytransformer.parameterobj.get("wavelength")
    first=True
    ni = len(translations)/100.0
    for i,(t_x, t_y, t_z) in enumerate(translations):
        mytransformer.updateparameters( )
        mytransformer.parameterobj.set_parameters(  {
            't_x':t_x, 't_y':t_y, 't_z':t_z
        } )
        mytransformer.compute_tth_eta( )
        mytransformer.computegv( )
        #    mytransformer.savegv( tmp+".gve" )
        gve = np.vstack(( mytransformer.colfile.gx,mytransformer.colfile.gy,mytransformer.colfile.gz ))
        if first:
            first=False
        grains = doindex(gve, t_x,t_y,t_z, w, gridpars)
        ng = len(grains)
        if ng > 0:
            grains = domap( mytransformer.parameterobj ,
                            mytransformer.colfile ,
                            grains,
                            gridpars)
            nk = len(grains)
            if len(grains) > 0:
                test_many_points.q.put( grains, False ) # do not wait
        else:
            nk = 0
        sys.stderr.write("         % 6.2f%% Position %d %d %d"%(i/ni,t_x, t_y, t_z)+
                         " grains found %d kept %d\n"%(ng, nk))

    
class uniq_grain_list(object):
    """
    Cope with finding the same grain over and over...
    """
    def __init__(self, symmetry, toldist, tolangle, grains=None):
        self.grp = getattr( sym_u, symmetry )()
        self.dt2 = toldist*toldist
        self.tar  = np.radians(tolangle)
        self.uniqgrains = []
        if grains is not None:
            self.add( grains )
            
    def add(self, grains):
        for i,gnew in enumerate(grains):
            newgrain = True
            for gold in self.uniqgrains:
                dt = gnew.translation - gold.translation
                dt2  =np.dot(dt,dt)
                if dt2 > self.dt2:
                    continue
                aumis = np.dot(gold.asymusT, gnew.U)
                arg = (aumis[:,0,0]+aumis[:,1,1]+aumis[:,2,2] - 1. )/2.
                angle = np.arccos(np.clip(arg, -1, 1)).min()
                if angle < self.tar:
                    # too close in angle and space
                    print( "           matched",i,np.degrees(angle),np.sqrt(dt2) )
                    gold.nfound += 1
                    newgrain = False
                    break
            if newgrain:
                self.append_uniq( gnew )
                
    def append_uniq( self, g ):
        symubis = [np.dot(o, g.ubi) for o in self.grp.group]
        g.asymusT   = np.array([xfab.tools.ubi_to_u_b(ubi)[0].T for ubi in symubis])
        g.nfound = 1
        self.uniqgrains.append( g )

        
def initgrid( fltfile, parfile, tmp, gridpars ):
    """
    Sets up a grid indexing by preparing the unitcell for indexing
    and checking the columns we want are in the colfile
    """
    mytransformer = transformer.transformer()
    mytransformer.loadfiltered( fltfile )
    mytransformer.loadfileparameters(  parfile )
    gridpars[ 'UC' ] = unitcell.unitcell_from_parameters( mytransformer.parameterobj )
    col = mytransformer.colfile
    if not "drlv2" in col.titles:
        col.addcolumn( np.ones(col.nrows, float),
                       "drlv2" )
    if not "labels" in col.titles:
        col.addcolumn( np.ones(col.nrows, float)-2,
                       "labels" )
    if not "sc" in col.titles:
        assert "xc" in col.titles
        col.addcolumn( col.xc.copy(), "sc")
    if not "fc" in col.titles:
        assert "yc" in col.titles
        col.addcolumn( col.yc.copy(), "fc")
    mytransformer.colfile.writefile( "%s.flt"%(tmp))
    if ('output_filename' not in gridpars) or (gridpars['output_filename'] is None):
        gridpars['output_filename'] = "all"+tmp+".map"
    try: # verify writable
        with open( gridpars['output_filename'], 'a' ) as fout:
            fout.write('#\n')
        print("Writing output to", gridpars['output_filename'])
    except Exception as e:
        print("Error writing to", gridpars['output_filename'])
        raise
    return gridpars


# Debugging multiprocessing - print the exceptions
def wrap_test_many_points(x):
    try:
        test_many_points( x ) 
    except Exception as e:
        print('Caught exception in worker thread')
        # This prints the type, value, and stack trace of the
        # current exception being handled.
        traceback.print_exc()
        raise e

def wrap_test_many_points_init(q):
    """
    This passes the q to the function
    Something from stackoverflow.
    It happens during the Pool initializer function
    """
    test_many_points.q = q


def grid_index_parallel( fltfile, parfile, tmp, gridpars, translations ):
    """
    fltfile containing peaks
    parfile containing instrument geometry and unit cell
    tmp - base name for scratch files and results
    gridpars : dictionary of control parameters (rings to use, etc)
    translations : list of translation positions to try

    Runs a grid index algorithm using pythons multiprocessing module
    splits workload over processes (blocks of translations to each process)
    This thread should catch results via a queue
    """
    gridpars = initgrid( fltfile, parfile, tmp, gridpars )
    print( "Done init" )
    if 'NPROC' not in gridpars or gridpars['NPROC'] is None:
        NPR = max( cImageD11.cores_available() - 1, 1 )
        cImageD11.cimaged11_omp_set_num_threads(1)
    else:
        NPR = int(gridpars['NPROC'])
    if 'NTHREAD' in gridpars:
        cImageD11.cimaged11_omp_set_num_threads(int(gridpars['NTHREAD']))
    elif NPR > 1:
        cImageD11.cimaged11_omp_set_num_threads(1)
    tsplit = [ translations[i::NPR] for i in range(NPR) ]
    args = [("%s.flt"%(tmp), parfile, t, gridpars) for i,t in enumerate(tsplit) ]
    q = PQueue()
    p = Pool(processes=NPR, initializer=wrap_test_many_points_init, initargs=[q])
    print( "Using a pool of",NPR,"processes" )
    pa = p.map_async( wrap_test_many_points, args )
    ul = uniq_grain_list( gridpars['SYMMETRY'],
                          gridpars['toldist'],
                          gridpars['tolangle'] )
    lastsave = 0

    while True:
        # If we go more than 30 seconds without something, die
        try:
            grs = q.get(True, 10)
            gb4 = len(ul.uniqgrains)
            ul.add( grs )
            gnow =  len(ul.uniqgrains)
            print( "Got % 5d new %d from %d"%(gnow, gnow-gb4, len(grs) ) )
            if len(ul.uniqgrains) > lastsave:
                lastsave = len( ul.uniqgrains )
                grain.write_grain_file( gridpars['output_filename'], ul.uniqgrains )
            if pa._number_left == 0:
                break
        except Queue.Empty:
            sys.stderr.write(" Caught queue empty exception\n")
            if pa._number_left == 0:
                break
        except KeyboardInterrupt:
            break
    # write here to be on the safe side .... 
    grain.write_grain_file( gridpars['output_filename'], ul.uniqgrains )
    p.close()
    p.join()


if __name__=="__main__":
    print("#Here is an example script")
    print("""
import sys, random
from ImageD11.grid_index_parallel import grid_index_parallel

if __name__=="__main__":
    # You need this idiom to use multiprocessing on windows (script is imported again)
    gridpars = {
        'DSTOL' : 0.004,
        'OMEGAFLOAT' : 0.13,
        'COSTOL' : 0.002,
        'NPKS' : int(  sys.argv[4] ),
        'TOLSEQ' : [ 0.02, 0.015, 0.01],
        'SYMMETRY' : "cubic",
        'RING1'  : [5,10],
        'RING2' : [5,10],
        'NUL' : True,
        'FITPOS' : True,
        'tolangle' : 0.25,
        'toldist' : 100.,
        'NPROC' : None, # guess from cpu_count
        'NTHREAD' : 2 ,
        'output_filename' : None, # defaults to "all"+tmp+".map"
    }
            
    # grid to search
    translations = [(t_x, t_y, t_z) 
        for t_x in range(-500, 501, 50)
        for t_y in range(-500, 501, 50) 
        for t_z in range(-500, 501, 50) ]
    # Cylinder: 
    # translations = [( x,y,z) for (x,y,z) in translations if (x*x+y*y)< 500*500 ]
    #
    random.seed(42) # reproducible
    random.shuffle(translations)

    fltfile = sys.argv[1]
    parfile = sys.argv[2]
    tmp     = sys.argv[3]
    grid_index_parallel( fltfile, parfile, tmp, gridpars, translations )
""")
