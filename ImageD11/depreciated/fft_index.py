

from ImageD11 import indexing
import sys
from Numeric import *
from LinearAlgebra import inverse, determinant
from FFT import fftnd , inverse_fftnd
from math import  cos,sin,radians,degrees

# test case
if 1:
    # Read in the g-vectors
    io = indexing.indexer()
    io.readgvfile(sys.argv[1])

# Decide on max resolution to use from detector
# higher resolution in reciprocal space = lower in real space
mr = float(sys.argv[2]) # raw_input("Max resolution ?")

cutoff  = float(sys.argv[3])
# Array to hold g-vectors on a regular grid
s = 128
g = zeros((s,s,s),Float)

cell_size = s*mr/2.
print "FFT cell size is 128*mr/2",cell_size

n_use = 0
n_ignore = 0
import time
start = time.time()
for p in io.gv:
    # compute hkl indices in g grid
    # units start for a 1x1x1 unit cell
    hrkrlr= cell_size * p
    hkl = floor(cell_size * p+0.5).astype(Int)
    assert maximum.reduce(abs(hrkrlr - hkl)) <= 0.5, str(hrkrlr)+" "+str(hkl)
    
    if maximum.reduce(hkl)<s/2-2 and minimum.reduce(hkl)>-s/2+2:
        # use the peak
        n_use += 1
        np =0
        sm = 0.
        #        print hrkrlr
        
        for hs in [0,1]:
            h=int(floor(hrkrlr[0]))+hs
            hfac = 1-abs(hrkrlr[0] - h)
            for ks in [0,1]:
                k=int(floor(hrkrlr[1]))+ks
                kfac = 1-abs(hrkrlr[1] - k)
                for ls in [0,1]:
                    l=int(floor(hrkrlr[2]))+ls
                    lfac = 1-abs(hrkrlr[2] - l)
                    np+=1
                    v = hfac*kfac*lfac
#                    print v,h,k,l,hfac,kfac,lfac
                    sm += v
                    g[h,k,l] = g[h,k,l] +  v
                    g[-h,-k,-l] = g[-h,-k,-l] + v
        assert np == 8, "np"
        assert abs(sm-1)<1e-5, "total peak, hkl "+str(hrkrlr)+" "+str(sm)
    else:
        n_ignore +=1
print "grid filling takes",time.time()-start
print "Used",n_use,"peaks and ignored",n_ignore,"peaks"
print "Starting FFT",g.shape, g.typecode()

#r = abs(inverse_fftnd(g))*multiply.reduce(g.shape)
r = abs(fftnd(g))

print "Done FFT", r.shape, r.typecode(),r[0,0,0]

# result is complex numbers
# peak search?

sh = r.shape
origin = r[0,0,0]
print "Origin peak",origin,"array shape", maximum.reduce(ravel(r))
l = multiply.reduce(r.shape)
indices = compress(r.flat > origin*cutoff, arange(0,l,1,Int))
#indices = array([0,1,2,67,12345,44563])
vals = take(r.flat,indices)
rlgrid = cell_size/s
print "center",r[64,64,64]
order = argsort(vals)
pks = []

class peak:
    grid_x = 128
    grid_y = 128
    grid_z = 128
    def __init__(self,x,y,z,h):
        self.x=x ; self.y = y; self.z = z
        self.h=h
    def add(self,other):
        h = self.h + other.h
        x = (self.x*self.h + other.x*other.h) / h
        y = (self.y*self.h + other.y*other.h) / h
        z = (self.z*self.h + other.z*other.h) / h
        #print self,other,peak(x,y,z,h)
        return peak(x,y,z,h)
    def dist(self,other):
        try:
            dx = self.x - other.x
            dy = self.y - other.y
            dz = self.z - other.z
            return sqrt(dx*dx  + dy*dy + dz*dz)
        except:
            print d2,self,other
            raise
                  
        dx = min( (self.x - other.x)%self.grid_x,
                  (other.x - self.x)%self.grid_x )
        dy = min( (self.y - other.z)%self.grid_y,
                  (other.y - self.y)%self.grid_y )
        dz = min( (self.z - other.z)%self.grid_z,
                  (other.z - self.z)%self.grid_z )
        return sqrt(dx*dx+dy*dy+dz*dz)
    def str(self):
        return "x=%d,y=%d,z=%d,h=%d"%(self.x,self.y,self.z,self.h)
    def __str__(self):
        return "x=%d,y=%d,z=%d,h=%d"%(self.x,self.y,self.z,self.h)
    def repr(self):
        return "x=%d,y=%d,z=%d,h=%d"%(self.x,self.y,self.z,self.h)
    def __repr__(self):
        return "x=%d,y=%d,z=%d,h=%d"%(self.x,self.y,self.z,self.h)

for o in order[::-1][:1000]:
    i,v = indices[o],vals[o]
    z = i%r.shape[2]
    y = (i/r.shape[2])%r.shape[1]
    x = (i/r.shape[2]/r.shape[1])%r.shape[0]
    new = True
    n = peak(x,y,z,r[x,y,z])
    for i in range(len(pks)):
        p=pks[i]
        d = p.dist(n)
        if d<4:
            pks[i] = p.add(n)
            new = False
            break
    if new: 
        pks.append(n)
origin = peak(0,0,0,0)
ooo = peak(s-1,s-1,s-1,0)

print len(pks),"peaks found, real space grid size",rlgrid
print "   i    j    k    height     length    x y z "
dsu = [ (p.h,p) for p in pks]
dsu.sort()
pks = [p[1] for p in dsu[::-1]]

qks = []
for p in pks[:40]:
    if p.x > s/2: p.x = p.x-s
    if p.y > s/2: p.y = p.y-s
    if p.z > s/2: p.z = p.z-s
    
    x,y,z = p.x*rlgrid , p.y*rlgrid, p.z*rlgrid
    l = sqrt(x*x+y*y+z*z)
    if l > 2:
        qks.append(p)
pks = qks
 
def dotp(p,q):
    sm = p.x*q.x + p.y*q.y + p.z*q.z
    lp = p.x*p.x + p.y*p.y + p.z*p.z
    lq = q.x*q.x + q.y*q.y + q.z*q.z
    return sm/sqrt(lp*lq)

for p in pks[:40]:
    if p.x > s/2: p.x = p.x-s
    if p.y > s/2: p.y = p.y-s
    if p.z > s/2: p.z = p.z-s
    
    x,y,z = p.x*rlgrid , p.y*rlgrid, p.z*rlgrid
    l = sqrt(x*x+y*y+z*z)
    print "%4d %4d %4d %9.2f %9.2f"%(p.x, p.y, p.z, p.h,l),
    print " %8.3f %8.3f %8.3f"%(x,y,z),
    for q in pks[:40]:        
        print "%.2f "%dotp(p,q),
    print 



sys.exit()
for g in io.gv:
    for i in g: print "%6.4f "%(i),
    h = dot(true_orientation,g)
    hint = floor(h+0.5).astype(Int)
    drlv = sqrt(dot(h-hint,h-hint))
    print "%8.4f "%(drlv),
    for i in h: print "%8.4f "%(i),
    print

    

