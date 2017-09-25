
import matplotlib.pylab as pl
import os, sys

if sys.platform == "win32":
    try:
        os.remove("connectedpixels.pyd")
    except WindowsError:
        pass
    f2py = "f2py.py"
else:
    f2py = "f2py"
os.system(f2py+" -I%s -c connectedpixels.pyf connectedpixels.c blobs.c closest.c -DF2PY_REPORT_ON_ARRAY_COPY"%(os.getcwd()))

import connectedpixels, numpy as np
import ImageD11.closest
np.random.seed(42)

cosines = np.linspace( -1, 1, 3 ).astype(np.float)
print "cosines"
artest = np.array([ 0.99, 0.2, 0.3, 0.7, 0.8],np.float)
oldy = ImageD11.closest.closest(artest, cosines )
newy = connectedpixels.closest( artest, cosines )
assert oldy == newy
print "Closest.closest seems OK"

ubitest = (np.random.random((3,3))-0.5)*10
gvtest = np.random.random( (4000, 3) )
def pyscore( u, g, tol ):
    h = np.dot( u, g.T )
    diff = (h-np.round(h))
    drlv2 = (diff*diff).sum(axis=0)
    n = (drlv2 < (tol*tol)).sum()
    return n
for tol in [0.01,0.1,0.2,0.4,0.8]:
    ninew = connectedpixels.score( ubitest, gvtest, tol )
    niold = ImageD11.closest.score( ubitest, gvtest, tol )
    npyth = pyscore( ubitest, gvtest, tol )
    assert ninew == niold, (ninew, niold)
print "Seems that score is OK"


print connectedpixels.s_1

BTYPE = np.intc


a = np.arange(3*4).reshape(3,4).astype(np.float32)
r = np.zeros( a.shape, BTYPE )
npk = connectedpixels.connectedpixels( a, r, 1.01, 1)
print a
print "Pks",npk
print r
# npk = connectedpixels.connectedpixels( a.astype(float), r, 1, 1)
npk = connectedpixels.connectedpixels( a, r, 1.01, 1, 0)
print "Pks",npk
print r

a = np.array([
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,1,0,1,0,0,1,1,1,1,0,0,0,0,0 ],
    [ 0,0,1,0,0,0,1,0,0,0,1,0,0,0,0 ],
    [ 0,0,0,0,0,0,1,0,1,0,1,0,0,0,0 ],
    [ 0,0,0,0,0,0,1,0,0,0,1,0,0,0,0 ],
    [ 0,0,0,0,0,0,1,1,1,1,0,0,0,0,1 ],
    [ 0,0,1,0,1,0,0,0,0,0,0,0,0,1,0 ],
    [ 0,0,1,0,1,0,1,0,0,0,0,0,1,0,0 ],
    [ 0,0,1,0,1,1,1,0,0,0,0,0,0,0,0 ]], np.float32)
r = np.zeros( a.shape, BTYPE)

npk = connectedpixels.connectedpixels( a, r, 0.5, 1, 0)
print a
print npk
print r

npk = connectedpixels.connectedpixels( a, r, 0.5, 1, 1)
print a
print npk
print r


import ImageD11.connectedpixels , time
N=3096
M=4096
t=8.

mytimer = time.clock

a = np.random.random((N*M)).reshape((N,M)).astype(np.float32)*10
bnew = np.zeros( a.shape, BTYPE )
start = mytimer()

npknew = connectedpixels.connectedpixels( a, bnew, t, verbose=1 )
newtime   = mytimer() - start

bold = np.zeros( a.shape, BTYPE )
start = mytimer() 
npkold = ImageD11.connectedpixels.connectedpixels( a, bold, t )
oldtime   = mytimer()  - start



if (bnew == bold).all() and npkold == npknew:
    print "Matches old code", npkold
    print "Old timer",oldtime,"New timer", newtime
    print "Speed up is a factor of %.2f"%(oldtime/newtime)
else:
    print 'old' , npkold
    print 'new' , npknew
    print bold
    print bnew
    print "here: looks bad, one is wrong"
    pl.subplot(221)
    pl.imshow(a>t)
    pl.subplot(224)
    pl.imshow(bnew)
    pl.colorbar()
    pl.subplot(223)
    pl.imshow(bold)
    pl.colorbar()
    pl.subplot(222)
    pl.imshow(bold-bnew)
    pl.colorbar()
    pl.show()
    

oldmoments = ImageD11.connectedpixels.blobproperties( a, bold, npkold, 30.0  )
newmoments = connectedpixels.blobproperties( a, bold, npkold, 30.0  )

if oldmoments.shape == newmoments.shape and (oldmoments == newmoments).all():
    print "blobproperties seems OK"

ImageD11.connectedpixels.blob_moments( oldmoments )
connectedpixels.blob_moments( newmoments )
if oldmoments.shape == newmoments.shape and (oldmoments == newmoments).all():
    print "blob moments seems OK"


N=40
M=20
b0=[] ; n0=[] ; m0=[]; b1=[]; n1=[]; m1=[]
t = 8
NI = 3

for i in range(NI):
    print i,
    aimg = np.random.random((N*M)).reshape((N,M)).astype(np.float32)*10
    # New code
    bnew = np.zeros( aimg.shape, BTYPE )
    npnew = connectedpixels.connectedpixels( aimg, bnew, t, verbose=0 )
    mnew = connectedpixels.blobproperties(   aimg, bnew, npnew, i  )
    print "b4new",i, npnew, mnew[:,0].astype(int)
    if i>1:
        npknew = connectedpixels.bloboverlaps( bnew  ,npnew ,mnew,
                                               b0[-1],n0[-1],m0[-1],
                                               verbose=0)
    print "Afternew",i, npnew, mnew[:,0].astype(int)
    b0.append(bnew)
    n0.append(npnew)
    m0.append(mnew)
    # old code
    bold = np.zeros( aimg.shape, BTYPE )
    npold = ImageD11.connectedpixels.connectedpixels( aimg, bold, t, verbose=0 )
    mold = ImageD11.connectedpixels.blobproperties(   aimg, bold, npold, i  )
    print "b4old",i, npold, mold[:,0].astype(int)
    if i>1:
        npkold = ImageD11.connectedpixels.bloboverlaps( bold, npold, mold,
                                                      b1[-1],n1[-1],m1[-1],
                                                      verbose=0)
    print "afterold",i, npold, mold[:,0].astype(int)
    b1.append(bold)
    n1.append(npold)
    m1.append(mold)
    if npknew == npkold:
        print "OK",

def check(x,y,s):
    print "Checking",s
    sys.stdout.flush()
    if (x == y).all():
        print "OK",
    else:
        pl.subplot(221)
        pl.imshow(x)
        pl.title("old")
        pl.colorbar()
        pl.subplot(222)
        pl.title("new")
        pl.imshow(y)
        pl.colorbar()
        pl.subplot(223)
        pl.imshow(x-y)
        pl.show()
        raise Exception("Stop here")

for i in range(NI):
    print i,
    check( b0[i], b1[i], "b" )
    check( np.array(n0[i]) , np.array(n1[i]), "n")
    check( m0[i] , m1[i], "m")




