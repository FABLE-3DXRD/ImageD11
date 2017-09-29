
from __future__ import print_function


from numpy import arange, ones, concatenate
from ImageD11.columnfile import columnfile
import glob, sys
n = 0
allc=None


for f in sys.argv[2:]:
        c = columnfile(f)
        c.filter(c.Number_of_pixels > 4)
        #print f.split("_")
        nums = int(f.split("_")[-3])
        print(f,nums)
        c.addcolumn(ones(c.nrows)*nums,"dty")
        if allc is None:
            allc = c.copy()
        else:
            allc.bigarray = concatenate( (allc.bigarray, c.bigarray), axis=1 )
            allc.nrows = allc.bigarray.shape[1]
            allc.set_attributes()
        
allc.writefile(sys.argv[1])
