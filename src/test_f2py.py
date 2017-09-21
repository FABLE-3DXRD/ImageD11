

import os, sys

if sys.platform == "win32":
    f2py = "f2py.py"
else:
    f2py = "f2py"
os.system(f2py+" -I%s -c connectedpixels.pyf connectedpixels.c blobs.c  -DF2PY_REPORT_ON_ARRAY_COPY"%(os.getcwd()))

import connectedpixels, numpy as np

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
np.random.seed(42)
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
    import matplotlib.pylab as pl
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


N=200
M=300
b=[]
n=[]
m=[]
for i in range(5):
    anew = np.random.random((N*M)).reshape((N,M)).astype(np.float32)*10
    bnew = np.zeros( anew.shape, BTYPE )
    npnew = connectedpixels.connectedpixels( anew, bnew, t, verbose=0 )
    mnew = connectedpixels.blobproperties( anew, bnew, npnew, i  )
    if i>1:
        connectedpixels.bloboverlaps( bnew, npnew, mnew,
                                      b[-1],n[-1],m[-1] )
    b.append(bnew)
    n.append(npnew)
    m.append(mnew)
    


