

from ImageD11 import refinegrains
import sys, os

parfile = sys.argv[1]
ubifile = sys.argv[2]
fltfile = sys.argv[3]
newpars = sys.argv[4]

if os.path.exists(newpars):
    ok = raw_input("OK to overwrite %s ? [y/n]"%(newpars))
    if ok not in ["y","Y","Yes","yes","YES"]:
        sys.exit()

o = refinegrains.refinegrains(tolerance=0.5)
o.loadparameters(parfile)
o.readubis(ubifile)
o.loadfiltered(fltfile)
o.generate_grains()
o.refineubis(quiet=False)
o.parameterobj.varylist = [ "y_center","z_center",
        "tilt_y","tilt_x","tilt_z","wedge",
#        "wavelength",
        "t_x","t_y","distance"]

o.fit(maxiters=1000)
o.refineubis(quiet=False)

o.saveparameters(newpars)

from matplotlib.pylab import plot, show
plot(o.drlv)
show()

