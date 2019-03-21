

import numpy as np
from ImageD11 import cImageD11


m = np.zeros((10,10),np.int8)
mm = m.copy()
m[5,5:7]=1
m[2,2]=1
cImageD11.clean_mask(m, mm)
print(mm)
