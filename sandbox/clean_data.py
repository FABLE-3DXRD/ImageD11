
from __future__ import print_function
from six.moves import input
from ImageD11 import peakmerge, indexing, transformer
from ImageD11 import grain
import sys, os, numpy as np

mytransformer = transformer.transformer()
#
# Your work starts here:
#    
print("select only peaks of certain rings ")

mytransformer.loadfiltered( sys.argv[1] )
mytransformer.loadfileparameters(  sys.argv[2] )
outfile = sys.argv[3]

mytransformer.updateparameters( )
tth, eta = mytransformer.compute_tth_eta( )

from matplotlib.pylab import plot, show, figure
figure(1)
plot(tth,eta,",")

mytransformer.colfile.filter( mytransformer.colfile.Number_of_pixels > 3 )        
tth =    mytransformer.colfile.tth

mytransformer.addcellpeaks()

rh = mytransformer.unitcell.ringhkls
peaks = list(rh.keys())
peaks.sort()
m = 0
print("Ring ds  tth  (h k l) mult  npks  npks/mult")
for i,d in enumerate(peaks):
    thistth = mytransformer.theorytth[i]
    mult = len(rh[d])
    print("Ring %3d %2.4f %8.4f"%(i,d,thistth),"(%2d %2d %2d)"%rh[d][0],"%2d"%(mult), end=' ')
    msk = ((thistth+0.2) > tth)&(tth > (thistth-0.2))
    npk = msk.sum()
    print("%8d %8d"%(npk,1.0*npk/mult))
    m += msk.sum()
print(m, tth.shape)

figure(2)
plot(tth,np.log(mytransformer.colfile.sum_intensity),  ",")
show()

if 1:
    while 1:
        rings_to_use = input("Enter comma separated list of rings to use")
        try:
            rs = [int(v) for v in rings_to_use.split(',')]
        except:
            print("I don't understand, please try again")
            continue
        print("You entered"," ".join([str(r) for r in rs]))
        if input("OK?") in ['Y','y']:
            break
else:
    rs = [ 0, 2, 6, 10]

thistth = mytransformer.theorytth[rs[0]]
msk = ((thistth+0.2) > tth)&(tth > (thistth-0.2))
for r in rs[1:]:
    thistth = mytransformer.theorytth[r]
    msk |= ((thistth+0.2) > tth)&(tth > (thistth-0.2))

mytransformer.colfile.filter(msk)

mytransformer.write_colfile(sys.argv[3])
t,e=mytransformer.compute_tth_eta( )
figure(1)
plot(t,e,"o")
show()
