
from ImageD11 import refinegrains
import Numeric
import sys
fltfile = "peaks.out_merge_t100"

o = refinegrains.refinegrains(tolerance=0.02)

o.loadparameters("g3.pars")
o.readubis("g3.ubi")
o.loadfiltered(fltfile)
o.generate_grains()
for tol in range(1,11):
    o.tolerance = tol*0.01
    print "tol",o.tolerance,
    o.refineubis(quiet=False)
o.tolerance = float(raw_input("Enter tolerance"))
for i,gn in zip(range(len(o.grainnames)),o.grainnames):
    print i,gn
gn = o.grainnames[int(raw_input("select which grain"))]
o.compute_gv(gn,fltfile)
matrix = o.refine(o.grains[(gn,fltfile)].ubi)

h=Numeric.matrixmultiply(matrix,Numeric.transpose(o.gv))
hint=Numeric.floor(h+0.5).astype(Numeric.Int) # rounds down
diff=h-hint
drlv=Numeric.sqrt(Numeric.sum(diff*diff,0))

print o.scandata[fltfile].shape
out = open(sys.argv[2],"w")
out.write("# xc yc omega npixels avg_intensity x_raw y_raw sigx sigy covxy\n")
for i in range(o.scandata[fltfile].shape[0]):
    if drlv[i] < o.tolerance:
        for v in o.scandata[fltfile][i,:]:
            out.write("%f "%(v))
        out.write("\n")
out.close()

