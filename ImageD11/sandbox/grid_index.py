
from __future__ import print_function

from ImageD11 import peakmerge, indexing, transformer
from ImageD11 import grain, unitcell, refinegrains
import sys, os, numpy as np, time, random

mytransformer = transformer.transformer()
#
# Your work starts here:
#    
mytransformer.loadfiltered( sys.argv[1] )
mytransformer.loadfileparameters(  sys.argv[2] )
tmp = sys.argv[3]
DSTOL = 0.003
OMEGAFLOAT = 0.13
NPKS = 80
TOLSEQ = [ 0.04, 0.01]
SYMMETRY = "cubic"
RING1 = [1,0]
RING2 = [1,]
# grid to search
translations = [(t_x, t_y, t_z) 
        for t_x in range(-500, 501, 200)
        for t_y in range(-500, 501, 200) 
        for t_z in range(-500, 501, 200) ]
# Cylinder: 
# filter( lambda x,y,z : (x*x+y*y)< 500*500, translations )
random.seed(42) # reproducible
random.shuffle(translations)


# printing or quiet:
if 1:
    NUL = open("/dev/null","w")
else:
    NUL = sys.stdout
    
# Get the unit cell
UC = unitcell.unitcell_from_parameters( mytransformer.parameterobj )


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


class tick:
    def __init__(self):
        self.last = time.time()
        self.start = self.last
    def tick(self, msg=None):
        if msg is not None:
            sys.stdout.write(msg)
        now = time.time()
        sys.stdout.write(" %.2f total %.2f /s\n"%(now-self.last, now-self.start))
        self.last=now

ticker = tick()


def domap(  OmFloat,  OmSlop,
            pars,
            colfile,
            grainsfile,
            tolseq = [ 0.03, 0.02, 0.01],
            symmetry = "triclinic"):
    """mapping function - does what makemap.py does"""
    #flake8: global NPKS
    ss = sys.stdout # turns off printing
    for tol in tolseq:
        sys.stdout = NUL
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
        o.readubis( grainsfile )
        if symmetry != "triclinic":
            o.makeuniq( symmetry )
        o.generate_grains()
        o.refinepositions()
        gl = [x for x in list(o.grains.values()) if x.npks > NPKS]
        sys.stdout = ss
        if len(gl) == 0:
            print("I killed all your grains!")
            break
        else:
            print("Keeping",len(gl),"from",len(list(o.grains.values())),"grains with at least",NPKS,"peaks",tol)
            grain.write_grain_file( grainsfile ,  gl )
    return len(gl)



def doindex( gve, x, y, z, w):
    #flake8: global NUL
    #flake8: global tmp
    #flake8: global NPKS
    #flake8: global UC
    #flake8: global TOLSEQ
    ss = sys.stdout # turns of printing
    sys.stdout = NUL
    myindexer = indexing.indexer(
        wavelength = w,
        unitcell = UC,
        gv = gve.T
        )
    myindexer.ds = np.sqrt( (gve * gve).sum(axis=0) )
    myindexer.ga = np.zeros(len(myindexer.ds),int)-1 # Grain assignments
    #   myindexer.readgvfile( gve )

    for ring1 in RING1:
        for ring2 in RING2:
            myindexer.parameterobj.set_parameters(  {
                'ds_tol': DSTOL, 
                'minpks': NPKS, 
                'max_grains': 1000, 
                'hkl_tol': TOLSEQ[0], 
                'ring_1': ring1,
                'ring_2': ring2
                } )
            myindexer.loadpars( )
            myindexer.assigntorings( )
            myindexer.find( )
            myindexer.scorethem( )
    grains = [grain.grain(ubi, [x,y,z]) for ubi in myindexer.ubis]
    grain.write_grain_file("%s.ubi"%(tmp),grains)
    sys.stdout = ss
    return len(grains)

    
### Computation starts here


os.system("cp %s %s.flt"%(sys.argv[1],tmp))

first = True

w =  mytransformer.parameterobj.get("wavelength")

totalgrains = 0

for t_x, t_y, t_z in translations:           
    mytransformer.updateparameters( )
    mytransformer.parameterobj.set_parameters(  {
                  't_x':t_x, 't_y':t_y, 't_z':t_z
                  } )
    mytransformer.compute_tth_eta( )
    if first:
        mytransformer.addcellpeaks( )
    mytransformer.computegv( )
    #    mytransformer.savegv( tmp+".gve" )
    gve = np.vstack(( mytransformer.colfile.gx,mytransformer.colfile.gy,mytransformer.colfile.gz ))
    if first:
        first=False
    ng = doindex(gve, t_x,t_y,t_z, w)
    print("Position",t_x, t_y, t_z,"Grains", ng, end=' ')
    ticker.tick()
    if ng > 0:
        #ret = os.system("makemap.py -u %s.ubi -U %s.map -s cubic --omega_slop=0.13 "%(tmp,tmp) +
        #                    "-f %s.flt  -t 0.02 -p %s  "%(
        #                       tmp,sys.argv[2]))
        nfit = domap( True,  OMEGAFLOAT,
                      mytransformer.parameterobj ,
                      mytransformer.colfile ,
                      "%s.ubi"%(tmp),
                      tolseq = TOLSEQ,
                      symmetry = SYMMETRY)
        if nfit > 0:
            open("all%s.map"%(tmp),"a").write(open("%s.ubi"%(tmp)).read())
            mytransformer.colfile.filter( mytransformer.colfile.labels < -0.1 )
            totalgrains += nfit
        ticker.tick()
    if  mytransformer.colfile.nrows == 0:
        print("All peaks indexed")
        break
    else:
        print("still got",mytransformer.colfile.nrows," to index, found",totalgrains,"in total")
            
        


