
from __future__ import print_function
import sympy
sympy.init_printing()
x,y,z=sympy.symbols("x,y,z")
gvec = sympy.Matrix( [x,y,z] )
axis=sympy.Matrix([0,0,1])
omegadir=axis.cross(gvec)
etadir = omegadir.cross(gvec)

print("# gvec:",gvec)
print("# axis:",axis)
print("# omegadir:",omegadir)
print("# etadir:",etadir)
print("# |etadir|2:",sympy.simplify(etadir.dot(etadir)))
print("# |omegadir|2:",sympy.simplify(omegadir.dot(omegadir)))

# gvec: Matrix([[x], [y], [z]])
# axis: Matrix([[0], [0], [1]])
# omegadir: Matrix([[-y], [x], [0]])
# etadir: Matrix([[x*z], [y*z], [-x**2 - y**2]])
# |etadir|2: x**2*z**2 + y**2*z**2 + (x**2 + y**2)**2
# |omegadir|2: x**2 + y**2
