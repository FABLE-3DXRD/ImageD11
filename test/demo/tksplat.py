from __future__ import print_function, division
from ImageD11 import cImageD11, indexing
from time import perf_counter_ns
import numpy as np
from PIL import ImageTk, Image
import tkinter as Tk

i = indexing.indexer()
i.readgvfile("eu3.gve")
gve = i.gv.copy().astype(np.float)

w,h = 512,512
rgba = np.zeros( (w, h, 4), np.uint8 )
u = np.eye(3, dtype=np.float).ravel()
t0 = perf_counter_ns()
cImageD11.splat( rgba, gve, u, 2 )
t1 = perf_counter_ns()
im = Image.fromarray( rgba, mode="RGBA")
r = Tk.Tk()
myimg = ImageTk.PhotoImage( im )
l = Tk.Label( master=r, image = myimg )
l.pack()

s = np.sin(np.radians(2))
c = np.cos(np.radians(2))
ry = np.array( [[c,0,s],[0,1,0],[-s,0,c]] )
s = np.sin(np.radians(2/np.sqrt(2)))
c = np.cos(np.radians(2/np.sqrt(2)))
rx = np.array( [[1,0,0],[0,c,s],[0,-s,c]] )
rz = np.dot( ry, rx)
def rotate():
    global u, rz, rgba, gve, myimg, l
    t0 = perf_counter_ns()
    u = np.dot( rz, np.reshape(u,(3,3) )).ravel()
    cImageD11.splat( rgba, gve, u, 2 )
    im = Image.fromarray( rgba, mode="RGBA") 
    myimg.paste(im)
    r.after( 60, rotate )
    t1 = perf_counter_ns()
    print("%8.3f ms per frame"%((t1-t0)/1e6),end="\r")

r.after(60, rotate)
r.mainloop()
