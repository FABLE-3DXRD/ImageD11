from __future__ import print_function, division
from ImageD11 import cImageD11, indexing
from PIL import Image, ImageTk
from timeit import default_timer as timer
import numpy as np, os, sys
import nbsplat
import cProfile
try:
    import tkinter as Tk
except:
    import Tkinter as Tk

i = indexing.indexer()
try:
    i.readgvfile(sys.argv[1])
except:
    i.readgvfile(os.path.join(os.path.split(__file__)[0],"eu3.gve"))
gve = i.gv.copy().astype(float)/3.

r = Tk.Tk(sync=1)
w,h = r.winfo_screenwidth()*2//3,r.winfo_screenheight()*2//3

rgba = np.empty((h,w,4), np.uint8)
u = np.eye(3, dtype=float).ravel()
t0 = timer()
cImageD11.splat( rgba, gve, u, 2 )
t1 = timer()

p = ImageTk.PhotoImage( Image.fromarray(rgba.copy(), "RGBA") ) 
l = Tk.Label( r, width=w, height=h , background='black', image = p )
l.photo = p
l.pack()

s = np.sin(np.radians(2))
c = np.cos(np.radians(2))
ry = np.array( [[c,0,s],[0,1,0],[-s,0,c]] )
s = np.sin(np.radians(2/np.sqrt(2)))
c = np.cos(np.radians(2/np.sqrt(2)))
rx = np.array( [[1,0,0],[0,c,s],[0,-s,c]] )
rz = np.dot( ry, rx)

colors = (np.random.random( (len(gve),4))*255).astype( np.uint8 )


def rotate():
    #flake8: global u, rz, rgba, gve, l, nframes
    global u
    t0 = timer()
    u = np.dot( rz, np.reshape(u,(3,3) ))
    nbsplat.nbsplat( rgba, gve.T, u*(w+h)/2, 2 , colors=colors )
    x = Image.fromarray( rgba, "RGBA")
    t1 = timer()
   # import pdb;pdb.set_trace()
    l.photo.paste( x )
    t2 = timer()
    r.update()
    print("splat %8.3f blit %8.3f total %8.3f ms per frame"%(
        (t1-t0)*1e3,(t2-t1)*1e3,(t2-t0)*1e3), end="\n")
    r.after(1, rotate)

r.after(60, rotate)
r.mainloop()
