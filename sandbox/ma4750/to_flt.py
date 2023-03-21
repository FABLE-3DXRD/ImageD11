
import sys
    

"""
Script for repairing use of incorrect spline file
"""
import h5py, tqdm, numpy as np
from ImageD11 import columnfile, blobcorrector

spline = sys.argv[1]
cor = blobcorrector.correctorclass( spline )
h5in = sys.argv[2]
pks = sys.argv[3]
outname = sys.argv[4]

with h5py.File( h5in, 'r') as hin:
    print(list(hin),pks)
    g = hin[pks]
    cd = { name:g[name][:] for name in list(g) }
    cd['sc'] = np.zeros_like( cd['s_raw'] )
    cd['fc'] = np.zeros_like( cd['f_raw'] )

inc = columnfile.colfile_from_dict( cd )
       
for i in tqdm.tqdm(range( inc.nrows )):
    inc.sc[i], inc.fc[i] = cor.correct( inc.s_raw[i], inc.f_raw[i] )
    
columnfile.colfile_to_hdf( inc, outname )
