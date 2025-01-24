

from ImageD11.columnfile import columnfile
import sys, numpy as np


cf = columnfile(sys.argv[1])
x = float(sys.argv[2])
y = float(sys.argv[3])
cor = float(sys.argv[4])
tol = float(sys.argv[5])

r = cf.omega * np.pi/ 180.0

xcalc = np.sin(r)*x + np.cos(r)*y + cor

cf.filter( abs (cf.xpos - xcalc) < tol )

cf.writefile( sys.argv[6] )
