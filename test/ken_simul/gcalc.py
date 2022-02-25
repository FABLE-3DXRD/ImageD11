from __future__ import print_function
from ImageD11.indexing import readubis
from ImageD11.transform import uncompute_g_vectors
from ImageD11 import unitcell
import numpy as np

ul = readubis("test.ubi")
ui = [ np.linalg.inv(ubi) for ubi in ul ]

r = range(-4,5)
hkls = np.array([ (h,k,l) for h in r for k in r for l in r ])


gcalc = [ np.dot( ub, hkls.T).T for ub in ui]

# 30 keV
energy = 30
wvln = 12.3985 / energy
ng = 0
for g, ubi in zip(gcalc, ul):
    ng +=1
    print ("# Grain",ng, "\n# Energy",energy,"keV, wavelength %.5f"%(wvln))
    tth, eta, omega = uncompute_g_vectors( g.T, wvln )
    order = np.argsort(tth)
    # hkls = np.dot( ubi, g.T).T
    print ("#   h    k    l   tth     eta    omega")
    for j in range(len(order)):
        i = order[j]
        for s in (0,1):
            h,k,l = [int(v) for v in hkls[i]]
            if tth[i]>0  and tth[i] < 20 and \
                    abs(eta[s][i])< 45 and not unitcell.F(h,k,l):

                
                print ((" %4d"*3 + " %7.2f"*3) % tuple([h,k,l,
                        tth[i],
                        eta[s][i], omega[s][i]]))
    print()
