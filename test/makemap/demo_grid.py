#Here is an example script

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
        'RING1'  : [1, 5,],
        'RING2' : [1, 5],
        'NUL' : True,
        'FITPOS' : True,
        'tolangle' : 0.25,
        'toldist' : 100.,
        'NPROC' : None, # guess from cpu_count
        'NTHREAD' : 2 ,
    }
            
    # grid to search
    translations = [(t_x, t_y, t_z) 
        for t_x in range(-350, 351, 100)
        for t_y in range(-350, 351, 100) 
        for t_z in range(-10, 11, 10) ]
    # Cylinder: 
    # translations = [( x,y,z) for (x,y,z) in translations if (x*x+y*y)< 500*500 ]
    #
    random.seed(42) # reproducible
    random.shuffle(translations)

    fltfile = sys.argv[1]
    parfile = sys.argv[2]
    tmp     = sys.argv[3]
    
    if len(sys.argv)>5:
        gridpars['output_filename'] = sys.argv[5]
    
    grid_index_parallel( fltfile, parfile, tmp, gridpars, translations )

    
