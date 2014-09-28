
from ImageD11 import peakmerge, indexing, transformer
from ImageD11 import grain
import sys, os, numpy as np

mytransformer = transformer.transformer()
#
# Your work starts here:
#    
mytransformer.loadfiltered( sys.argv[1] )
mytransformer.loadfileparameters(  sys.argv[2] )

tmp = sys.argv[3]

NUL = open("/dev/null","w")

def doindex( gve, x, y, z):
    global NUL
    global tmp
    ss = sys.stdout # turns of printing
    #sys.stdout = NUL
    myindexer = indexing.indexer()
    myindexer.readgvfile( gve )
    for ring1 in [1,0,2]:
        for ring2 in [1,0]:
            myindexer.parameterobj.set_parameters(  {
                'ds_tol': 0.004, 
                'minpks': 20, 
                'max_grains': 1000, 
                'hkl_tol': 0.02, 
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

    
import random
random.seed(42)

translations = [(t_x, t_y, t_z) 
        for t_x in range(-250, 251, 25)
        for t_y in range(-250, 251, 25) 
        for t_z in range(-200, 201, 25) ]

random.shuffle(translations)
import time
start = time.time()

os.system("cp %s %s.flt"%(sys.argv[1],tmp))

first = True

for t_x, t_y, t_z in translations:           
    mytransformer.updateparameters( )
    mytransformer.parameterobj.set_parameters(  {
                  't_x':t_x, 't_y':t_y, 't_z':t_z
                  } )
    mytransformer.compute_tth_eta( )
    if first:
        mytransformer.addcellpeaks( )
    mytransformer.computegv( )
    mytransformer.savegv( tmp+".gve" )
    if first:
        first=False
    ng = doindex( tmp+".gve", t_x,t_y,t_z)
    print t_x, t_y, t_z, ng,time.time()-start
    if ng > 0:
        ret = os.system("makemap.py -u %s.ubi -U %s.map -s cubic --omega_slop=0.13 "%(tmp,tmp) +
                            "-f %s.flt  -t 0.02 -p %s  "%(
                               tmp,sys.argv[2]))
        if ret !=0 :
            print "bad 1"
            raise
        ret = os.system("makemap.py -u %s.ubi -U %s.map -s cubic --omega_slop=0.13 "%(tmp,tmp) +
                            "-f %s.flt  -t 0.01 -p %s  "%(
                               tmp,sys.argv[2]))
        if ret !=0 :
            print "bad 1"
            raise
        os.system("cutgrains.py %s.map %s.map 15"%(tmp,tmp))
        ret = os.system("makemap.py -u %s.map -U %s.map -s cubic --omega_slop=0.13  "%(tmp,tmp) +
                            "-f %s.flt -F %s.flt -t 0.01 -p %s "%(
                                tmp,tmp,sys.argv[2]))
        if ret !=0 :
            print "bad 3"
            raise
        open("all%s.map"%(tmp),"a").write(open("%s.map"%(tmp)).read())
        mytransformer.loadfiltered("%s.flt"%(tmp))
        


