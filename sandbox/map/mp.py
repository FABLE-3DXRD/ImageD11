
from __future__ import print_function

from ImageD11 import sym_u
from ImageD11.grain import read_grain_file

c = sym_u.cubic()
c.makegroup()

import glob

fl = glob.glob("x*.ubi")

for f in fl:
    gl = read_grain_file( f )
    x = int(f.split("_")[0][1:])
    y = int(f.split("_")[1].split(".")[0][1:])
    for g in gl:
        ubi = sym_u.find_uniq_u( g.ubi, c )
        g.set_ubi( ubi )
        print(x, y, end=' ')
        for i in range(3):
            for j in range(3):
                print(g.u[i][j], end=' ')
        print()
        
        
