#!/usr/bin/env python

from __future__ import print_function

"""
from example by Tarn Weisner Burton <twburton@users.sourceforge.net> in pyopengl
"""

__author__ = 'Jon Wright <jpwright@users.sourceforge.net> from example by Tarn Weisner Burton <twburton@users.sourceforge.net>'

import numpy
import sys
import os
from pyopengltk import Opengl
import OpenGL.GL as GL
import OpenGL.GLU as GLU
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk


class myOpengl(Opengl):


    # Make a parallel projection
    # mostly copied from Tk.Opengl class with small mods
    def tkRedraw(self, *dummy):
        """Cause the opengl widget to redraw itself."""
        if not self.initialised: 
            return
        self.activate()
        #print self.distance
        GL.glPushMatrix()			# Protect our matrix
        self.update_idletasks()
        self.activate()
        w = self.winfo_width()
        h = self.winfo_height()
        GL.glViewport(0, 0, w, h)
        # Clear the background and depth buffer.
        GL.glClearColor(self.r_back, self.g_back, self.b_back, 0.)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        r = 1.*w/h
        GL.glOrtho( -self.distance*r, self.distance*r, -self.distance, self.distance, 
                     -self.distance*3, self.distance*3)
#        GLU.gluPerspective(self.fovy, float(w)/float(h), self.near, self.far)
        GL.glMatrixMode(GL.GL_MODELVIEW)
        self.redraw(self)
        GL.glFlush()				# Tidy up
        GL.glPopMatrix()			# Restore the matrix
#        self.tk.call(self._w, 'swapbuffers')
        self.tkSwapBuffers()



class plot3d(Tk.Toplevel):
    def __init__(self,parent,data=None,lines=None,
                 ubis=None,image=None,pars=None,spline=None):
        """
        Data would be your observed g-vectors. Lines will
        be a computed lattice
        """
        Tk.Toplevel.__init__(self,parent)
        self.parent=parent
        if data is not None:
            xyz=data.copy()
        else:
            xyz=numpy.array([0,0,0])
        self.ps=Tk.StringVar()
        self.ps.set('1.')
        self.pointsize=1.
        self.npeaks=xyz.shape[0]

        self.o = myOpengl(self, width = 400, height = 400)
        self.o.redraw = self.redraw
        self.o.autospin_allowed = 1
        self.o.fovy=5
        self.o.near=1e6
        self.o.far=1e-6
        import math
        self.o.distance=3.
#numpy.maximum.reduce(numpy.ravel(xyz))*4 / \
#            math.tan(self.o.fovy*math.pi/180)
        print(type(xyz),xyz.dtype.char,xyz.shape)
        self.xyz=xyz
        f=Tk.Frame(self)
        Tk.Button(f,text="Help",command=self.o.help).pack(side=Tk.LEFT)
        Tk.Button(f,text="Reset",command=self.o.reset).pack(side=Tk.LEFT)
        Tk.Button(f,text="Pointsize",command=self.setps).pack(side=Tk.LEFT)
        Tk.Entry(f,textvariable=self.ps).pack(side=Tk.LEFT)
        Tk.Button(f,text="Quit",command=self.goaway).pack(side=Tk.RIGHT)
        self.dataoff=0
        self.o.pack(side = 'top', expand = 1, fill = 'both')
        f.pack(side=Tk.BOTTOM,expand=Tk.NO,fill=Tk.X)
        Tk.Label(self,text="Red=[1,0,0] Green=[0,1,0] Blue=[0,0,1]").pack(
            side=Tk.BOTTOM,expand=Tk.NO,fill=Tk.X)
        self.ubis=ubis
        self.color=numpy.ones((xyz.shape[0],3),float)
        print(self.color.shape)
        self.tex=False
        if ubis is not None:
           self.ubis = self.readubis(ubis)
           self.scorecolor(0)
        if pars is not None:
           self.tex=True
           self.readspline(spline)
           self.readprms(pars)
           self.readimage(image)
        self.after(100, self.changedata)

    def readspline(self,spline):
        from ImageD11 import blobcorrector
        self.corrector = blobcorrector.correctorclass(spline)

    def readubis(self,ubis):
        from ImageD11 import indexing
        return indexing.readubis(ubis)

    def readprms(self,prms):
        from ImageD11 import parameters
        o = parameters.parameters()
        o.loadparameters(prms)
        self.pars=o.get_parameters()

    def readimage(self,image):
        from ImageD11 import transform
        from fabio import openimage
        self.imageobj=openimage.openimage(image)
        # map from 2048x2048 to 1024x1024
        d = self.imageobj.data.astype(numpy.float32)
        mi= d.mean() - d.std()*2
        mx= d.mean() * d.std()*2
        shape=self.imageobj.data.shape
        d=numpy.reshape(numpy.clip(self.imageobj.data,mi,mx),shape) # makes a clipped copy
        d=(255.*(d-mi)/(mx-mi)) # scale intensity
        print(d.min(),d.max(),d.mean())
        self.image=numpy.zeros((1024,1024),numpy.uint8)
        if d.shape==(2048,2048):
            # rebin 2x2
            im=(d[::2,::2]+d[::2,1::2]+d[1::2,::2]+d[1::2,1::2])/4
            self.image=(255-im).astype(numpy.uint8).tostring()
        self.imageWidth=1024
        self.imageHeight=1024
        # make a 2D array of x,y
        p=[]
        pk=[]
        step = 64
        r=[ [ 0,0 ], [0,step], [step,step], [step,0] ]
        for i in range(0,1024,step):
            for j in range(0,1024,step):
                # i,j 1024x1024 texture coords
                # x,y spatially corrected
                for v in r:
                    pk.append([i+v[0],j+v[1]])
                    x,y = self.corrector.correct((i+v[0])*2 , (j+v[1])*2) # corrected
                    p.append([x,y])
        p=numpy.array(p).T 
        pk=numpy.array(pk).T
        omega=float(self.imageobj.header['Omega'])
        self.pars['distance']=float(self.pars['distance'])*1000
        tth,eta=transform.compute_tth_eta(p,**self.pars)
        gve = transform.compute_g_vectors(tth,eta,omega*self.pars['omegasign'],self.pars['wavelength'])
        self.pts = []
        print("Setting up image mapping",p.shape,gve.shape)
        for i in range(pk.shape[1]):
            self.pts.append([pk[1,i]/1024.,pk[0,i]/1024.,gve[0,i],gve[1,i],gve[2,i]])
        #for p in self.pts:
        #    print p
        self.setupTexture()

    def setupTexture(self):
        GL.glDisable(GL.GL_TEXTURE_2D)
        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        GL.glTexImage2D(GL.GL_TEXTURE_2D,#target
                    0,#level
                    3,#internalformat
                    self.imageWidth, self.imageHeight,
                    0,#border
                    GL.GL_LUMINANCE,#format
                    GL.GL_UNSIGNED_BYTE,# type
                    self.image)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_REPEAT)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_REPEAT)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST)
        GL.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_DECAL)
        GL.glEnable(GL.GL_TEXTURE_2D)
        GL.glEnable(GL.GL_NORMALIZE)
        GL.glShadeModel(GL.GL_FLAT)



    def scorecolor(self,i=0):
        cc = [ [ 1,0,0] , [0,1,0] , [0,0,1], [1,1,0], [1,0,1], [0,1,1],
                [ 0.5,0,0] , [0,0.5,0] , [0,0,0.5], [0.5,0.5,0], [0.5,0,0.5],
                [0,0.5,0.5]]
        if self.ubis is not None:
            from ImageD11 import indexing
            for u,i in zip(self.ubis,list(range(len(self.ubis)))):
                scores=indexing.calc_drlv2(u,self.xyz)
                print(self.xyz.shape,scores.shape)
                ind = numpy.compress( numpy.less(scores,0.05*0.05) , 
                                      numpy.arange(self.xyz.shape[0]) )
                print("Grain",i,scores.shape,ind.shape)
                for j in range(3):
                    c=numpy.ones(self.color.shape[0])
                    numpy.put(c,ind,cc[i%len(cc)][j])
                    self.color[:,j]*=c

    def go(self):
        """
        Allow the toplevel to return a handle for changing data
        """
        self.o.mainloop()

    def goaway(self):
        print("Called goaway")
        self.o.destroy()
        self.destroy()
        if self.parent is None: sys.exit()
        print("Ought to be gone now...")

    def changedata(self,xyz=None):
        if xyz is not None:
            self.xyz=xyz.copy()
            self.npeaks=xyz.shape[0]
        GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
        GL.glDisableClientState(GL.GL_COLOR_ARRAY)
        GL.glVertexPointer( 3, GL.GL_FLOAT, 0, self.xyz.astype(numpy.float32).tostring() )
        GL.glColorPointer( 3,  GL.GL_FLOAT, 0, self.color.astype(numpy.float32).tostring() )
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
        GL.glEnableClientState(GL.GL_COLOR_ARRAY)
        self.o.tkRedraw()

    def setps(self):
        self.pointsize=float(self.ps.get())
        self.o.tkRedraw()



    def redraw(self,o):
        
        GL.glDisable(GL.GL_LIGHTING)
        GL.glClearColor(0., 0., 0., 0)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glColor3f(1.0, 1.0, 1.0) # white
        GL.glPointSize(self.pointsize)
        GL.glDrawArrays(GL.GL_POINTS, 0, self.npeaks )

        if self.ubis is not None and len(self.ubis)==1:
            hkl = numpy.dot(numpy.linalg.inv(self.ubis[0]), 
                            numpy.identity(3,float)).T
            # print hkl
        else:
            hkl = numpy.identity(3,float)
            # print hkl
                            
        GL.glBegin(GL.GL_LINE_LOOP)
        GL.glColor3f(1.0, 0.0, 0.0) # red
        GL.glVertex3f(0.,0.,0.)
        GL.glVertex3f(hkl[0][0],hkl[0][1],hkl[0][2])
        GL.glEnd()

        GL.glBegin(GL.GL_LINE_LOOP)
        GL.glColor3f(0.0, 1.0, 0.0) # green
        GL.glVertex3f(0.,0.,0.)
        GL.glVertex3f(hkl[1][0],hkl[1][1],hkl[1][2])
        GL.glEnd()

        GL.glBegin(GL.GL_LINE_LOOP)
        GL.glColor3f(0.0, 0.0, 1.0) # blue
        GL.glVertex3f(0.,0.,0.)
        GL.glVertex3f(hkl[2][0],hkl[2][1],hkl[2][2])
        GL.glEnd()

        if self.tex:
#            print "drawing images"
            GL.glEnable(GL.GL_TEXTURE_2D)
            GL.glColor4f(.0, 1.0, .0, 1.0) # red        
            GL.glBegin(GL.GL_QUADS)
            # generate a grid of squares to map the texture in 3D
            # opengl has better "map" methods to do this
            for i,j,g1,g2,g3 in self.pts:
#                print i,j,g1,g2,g3
                GL.glTexCoord2f(i,j)
                GL.glVertex3f(g1, g2, g3) 
            GL.glEnd()
#            GL.glDisable(GL.GL_TEXTURE_2D)
            
        GL.glFlush()
        GL.glEnable(GL.GL_LIGHTING)







if __name__=="__main__":

    try:
        lines=open(sys.argv[1],"r").readlines()
    except:
        print("Usage %s gvector_file [ubifile] [image parfile]"%(sys.argv[0]))
        raise
        # sys.exit()
   
    on=0
    xyz=[]
    for line in lines:
        if on==1:
            try:
                vals=[float(x) for x in line.split()]
                xyz.append( [ vals[0],vals[1],vals[2] ])
            except:
                pass
        if line.find("xr yr zr")>0 or line.find("gx ")>0:
            on = 1
    xyz=numpy.array(xyz)
    if len(xyz) == 0 and lines[0][0]=="#":
        from ImageD11 import columnfile
        c = columnfile.columnfile( sys.argv[1] )
        xyz = numpy.array( (c.gx, c.gy, c.gz )).T
    npeaks = len(xyz)            
    if len(sys.argv)==3:
       o=plot3d(None,data=xyz,ubis=sys.argv[2])
    elif len(sys.argv)==6:
       o=plot3d(None,data=xyz,ubis=sys.argv[2],image=sys.argv[3],pars=sys.argv[4],spline=sys.argv[5])
    else:
       o=plot3d(None,data=xyz,ubis=None)
    def runit():
        o.changedata(o.xyz)
    o.after(100, runit )
    o.mainloop()
