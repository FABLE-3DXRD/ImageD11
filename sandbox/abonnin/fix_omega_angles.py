from ImageD11.columnfile import columnfile
from pylab import plot
c = columnfile("mainphase.flt")
c.addcolumn(c.omega.copy(), "image_number")
c.omega[:] = (c.image_number/147).astype(int)
c.addcolumn(0.2*(c.image_number - 147*c.omega),"xpos")
plot(c.omega,c.sum_intensity,",")
c.writefile("mainphase_fixed_omega.flt")
