
from __future__ import print_function

from ImageD11 import unitcell, parameters, transform
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

pars = parameters.read_par_file( sys.argv[1] )

# Input unit cell from par file
u = unitcell.unitcell_from_parameters(pars)

# Generate diffraction spots for unit cell:
u.makerings(0.7)
hkls = []
for d in u.ringds:
    hkls += u.ringhkls[d]
hkls = np.array(hkls).T

# A grain orientation
U0 = transform.detector_rotation_matrix(0.1, 0.2, 0.3 )
UB0 = np.dot( U0, u.B)


# Ideal diffraction vectors:
gve = np.dot( UB0, hkls )
modg= np.sqrt((gve*gve).sum(axis=0))

# print scattering vectors in order
order = np.argsort(modg)
for i in order[10::-1]:
    print (i,hkls.T[i],modg[i],gve.T[i])


def make_uniq( UBIlist ):
    # consider two UBI as equal if dot(UBIa, inv(UBIb)) == {1|0}
    # Sort by trace before to get back the ones with biggest traces
    dsu = [ (np.trace(x)+i*1e-10, x) for i,x in enumerate(UBIlist) ]
    dsu.sort()
    UL = [x[1] for x in dsu[::-1]]
    uniqs = [( UL[0] , np.linalg.inv( UL[0] ) ), ]
    for u1 in UL[1:]:
        new = True
        for u2,u2i in uniqs:
            testmat = np.dot( u1, u2i )
            itestmat = np.round( testmat )
            scor = abs(testmat - itestmat).sum()
            if scor < 1e-8:
                # Seen already
                new = False
                break
        if new:
            uniqs.append( (u1, np.linalg.inv( u1 ) ) )
    return [x[0] for x in uniqs]

def get_ubis( h1, h2, UB0, u ):
    # Compute g-vectors for these peaks in current unit cell
#    print "Using h1",h1,"h2",h2 
    g1 = np.dot( UB0, h1 )
    g2 = np.dot( UB0, h2 )
    modg1 = np.sqrt( np.dot( g1, g1 ) )
    modg2 = np.sqrt( np.dot( g2, g2 ) )
    # find unitcell powder rings for these
    r1 = np.argmin( abs(u.ringds-modg1) )
    r2 = np.argmin( abs(u.ringds-modg2) )
    # Generate orientations for these:
    u.orient( r1, g1, r2, g2)
#    print "Found",len(u.UBIlist),"orientations" 
    return u.UBIlist


# peak pair to use for twinning:
UBIs=[]
for h1, h2 in [ ( (1,0,0), (0,1, 1) ),
                ( (1,0,0), (0,1,-1) ),
                ( (0,1,0), (1,0, 1) ),
                ( (0,1,0), (1,0,-1) ),
                ( (0,0,1), (1,1, 0) ),
                ( (0,0,1),(-1,1, 0) ) ]:
        UBIs+=[x for x in get_ubis( h1, h2, UB0, u )]
uniqs = make_uniq( UBIs)

# Now do that again for each uniq one:
UBIs=[x for x in uniqs]

for UBi1 in uniqs:
#    print UBi1
    for h1, h2 in [ ( (1,0,0), (0,1, 1) ),
                ( (1,0,0), (0,1,-1) ),
                ( (0,1,0), (1,0, 1) ),
                ( (0,1,0), (1,0,-1) ),
                ( (0,0,1), (1,1, 0) ),
                ( (0,0,1),(-1,1, 0) ) ]:
        UB1 = np.linalg.inv( UBi1 )
        UBIs+=[x for x in get_ubis( h1, h2, UB1, u )]

print( len(uniqs))
uniqs2 = make_uniq( UBIs)
print( len(uniqs2))

# Now do that again for each uniq one:
UBIs=[x for x in uniqs2]

for UBi1 in uniqs2:
#    print UBi1
    for h1, h2 in [ ( (1,0,0), (0,1, 1) ),
                ( (1,0,0), (0,1,-1) ),
                ( (0,1,0), (1,0, 1) ),
                ( (0,1,0), (1,0,-1) ),
                ( (0,0,1), (1,1, 0) ),
                ( (0,0,1),(-1,1, 0) ) ]:
        UB1 = np.linalg.inv( UBi1 )
        UBIs+=[x for x in get_ubis( h1, h2, UB1, u )]
uniqs3 = make_uniq( UBIs)
print( len(uniqs3))


# Now do that again for each uniq one:
UBIs=[x for x in uniqs2]

for UBi1 in uniqs3:
#    print UBi1
    for h1, h2 in [ ( (1,0,0), (0,1, 1) ),
                ( (1,0,0), (0,1,-1) ),
                ( (0,1,0), (1,0, 1) ),
                ( (0,1,0), (1,0,-1) ),
                ( (0,0,1), (1,1, 0) ),
                ( (0,0,1),(-1,1, 0) ) ]:
        UB1 = np.linalg.inv( UBi1 )
        UBIs+=[x for x in get_ubis( h1, h2, UB1, u )]
uniqs4 = make_uniq( UBIs)
print( len(uniqs4))


uniqs= uniqs4

# Plot out peak positions on detector for the uniqs
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(223)
ax3 = fig.add_subplot(222, projection='3d')
fig2= plt.figure(2)
ax4 = fig2.add_subplot(111)
colors = 'bgryck'

hkls = np.array([(0,0,2),(0,0,-2),(2,0,0),(-2,0,0),(0,2,0),(0,-2,0)]).T
XX=True
def doplot(ubi, col):
    g = np.dot( np.linalg.inv( ubi ), hkls )
    tth, (eta0,eta1), (omega0,omega1) = transform.uncompute_g_vectors( g, pars.get('wavelength') )
    ax1.plot(tth, eta0, col)
    ax1.plot(tth, eta1, col)
    sc, fc = transform.compute_xyz_from_tth_eta( tth, eta0, omega0, **pars.parameters )
    ax2.plot( sc, fc, col)
    sc, fc = transform.compute_xyz_from_tth_eta( tth, eta1, omega1, **pars.parameters )
    ax2.plot( sc, fc, col)
    ax3.plot( g[0], g[1], g[2], col)
    ax4.scatter( eta0, omega0, c=tth)
    cax = ax4.scatter( eta1, omega1, c=tth)
    global XX
    if XX:
        plt.colorbar( cax )
        XX=False
for i,ubi in enumerate(uniqs[1:]):
    col = "%s+"%(colors[i%len(colors)])
    print( i,col)
    print( ubi)
    doplot( ubi, col )

doplot( uniqs[0], "ko" )

plt.show()




