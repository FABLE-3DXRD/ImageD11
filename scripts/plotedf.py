#!/usr/bin/env python

from __future__ import print_function

## Automatically adapted for numpy.oldnumeric Sep 06, 2007 by alter_code1.py



# This is statement is required by the build system to query build info
if __name__ == '__build__':
    raise Exception


## example python/pyopengl script to do tiled texturemapping.
## By david konerding (dek@cgl.ucsf.edu)

import sys
#from Image import *
import OpenGL.GL as OGL
from OpenGL.GLU.projection import gluProject, gluUnProject
import OpenGL.Tk as OTk
import numpy as np
import math
import fabio

class edfFile:
    """
    Absolute minimum edf file format, UInt16, rows and columns
    """
    def __init__(self,filename):
        #f=open(filename,"rb")
        #h=f.read(1024)
        #rows=None
        #cols=None
        #for item in h.split(";"):
        #       if item.find("Dim_1")>=0:
        #               rows=int(item.split()[-1])
        #       if item.find("Dim_2")>=0:
        #               cols=int(item.split()[-1])
        #f.seek(-rows*cols*2,2)
        #self.rows=rows
        #self.cols=cols
        #self.data=array(fromstring(f.read(rows*cols*2),UInt16),savespace=1)
#        self.data=fabio.open(filename).data[:1024,:1024].copy()
        self.data=fabio.open(filename).data#[:1024,:1024]
        self.rows=self.data.shape[0]
        self.cols=self.data.shape[1]
        self.data=np.ravel(self.data)
        self.minI=np.minimum.reduce(np.ravel(self.data))
        self.maxI=np.maximum.reduce(np.ravel(self.data))
        print("Opened",filename,"max=",self.maxI,"min=",self.minI)

class myOpengl(OTk.Opengl):

    def __init__(self, master=None, cnf={}, **kw):
        OTk.Opengl.__init__(*(self, master, cnf), **kw)

    def StartRotate(self,event):
        """
        Clear old selection box
        Start new one
        """
        pass
#               Opengl.StartRotate(self,event)


    def tkRotate(self, event):
        """
        Draw selection box ??? Not working
        """
        win_height = max( 1, self.winfo_height() )

        obj_c   = ( 0., 0., 0. )
        win     = gluProject( obj_c[0], obj_c[1], obj_c[2])
        obj     = gluUnProject( win[0], win[1] + 0.5 * win_height, win[2])
        dist    = math.sqrt( (obj_c[0]-obj[0])**2 + (obj_c[1]-obj[1])**2 + (obj_c[2]-obj[2])**2 )
        scale     = abs( dist / ( 0.5 * win_height ) )
        realy = self.winfo_height() - event.y
        p1 = gluUnProject(event.x, realy, 0.) # Image is at z = 0
        p2 = gluUnProject(event.x, realy, 1.) # Image is at z = 0
        print(p1[0],p1[1],p1[2])
#               Opengl.tkRotate(self,event)

    def tkAutoSpin(self, event):
        """
        Finish drawing selection box
        """
        pass
#               Opengl.tkAutoSpin(self,event)


class checker:

    def makeImage(self):
        try:
            mi=int(self.minI.get())
            mx=int(self.maxI.get())
        except:
            mi=self.edfFile.minI
            mx=self.edfFile.maxI
        shape=(self.edfFile.rows, self.edfFile.cols)
        d=np.reshape(np.clip(self.edfFile.data,mi,mx),shape) # makes a clipped copy
        print("makeImage",mx,mi,np.maximum.reduce(np.ravel(d)),np.minimum.reduce(np.ravel(d)),d.dtype.char, end=' ')
        newshape = []
        for i in shape:
            j=4
            print(j,pow(2,j),i,i<pow(2,j))
            while i > pow(2,j):
                j+=1
            newshape.append(j)
        newshape = tuple([pow(2,v) for v in newshape])
        print("newshape",newshape)
        d=255.*(d-mi)/(mx-mi)
        self.image=np.zeros((newshape[0],newshape[1],3),np.uint8)
        print(self.image.shape,d.shape)
        self.image[:shape[0],:shape[1],0] = d
        self.image[:shape[0],:shape[1],1] = d
        self.image[:shape[0],:shape[1],2] = d
        print(self.image.shape)
#        import pylab as pl
#        pl.imshow(self.image)
#        pl.show()
        # self.image = self.image # .tostring()
        self.imageWidth = newshape[1]
        self.imageHeight = newshape[0]
        print("Returning")


    def display(self, event=None):
        OGL.glClearColor( .7, 0.8, 0.9, 0)
        OGL.glClear(OGL.GL_COLOR_BUFFER_BIT | OGL.GL_DEPTH_BUFFER_BIT)
        OGL.glBegin(OGL.GL_QUADS)

        w=self.imageWidth/2
        h=self.imageHeight/2
        
        OGL.glTexCoord2f(0.0, 0.0);         OGL.glVertex3f(-w, -h, 0.0)
        OGL.glTexCoord2f(0.0, 1.0);         OGL.glVertex3f(-w,  h, 0.0)
        OGL.glTexCoord2f(1.0, 1.0);         OGL.glVertex3f( w,  h, 0.0)
        OGL.glTexCoord2f(1.0, 0.0);         OGL.glVertex3f( w, -h, 0.0)

        OGL.glEnd()
        OGL.glFlush()

    def change(self, event=None):
        self.SetupTexture()
        self.ogl.tkRedraw()

    def SetupWindow(self):

        self.OglFrame = OTk.Frame()
        self.OglFrame.pack(side = 'top', expand=1 ,fill='both')

        self.ogl = myOpengl(master=self.OglFrame, width = 500, height = 500, double = 1)


        self.ogl.pack(side = 'top', expand = 1, fill = 'both')
        self.ogl.distance=max(self.imageWidth+10,self.imageHeight+10)*10
        self.ogl.near=max(self.imageWidth+10,self.imageHeight+10)*100.
        self.ogl.far=max(self.imageWidth+10,self.imageHeight+10)/100.
        self.ogl.fovy=10.
        self.ogl.autospin_allowed=1
        self.ogl.redraw = self.display


        # Control buttons for scaling
        self.bf=OTk.Frame()
        self.QuitButton = OTk.Button(self.bf, {'text':'Quit'})
        self.QuitButton.bind('<ButtonRelease-1>', sys.exit)
        self.QuitButton.pack(side=OTk.RIGHT)

        OTk.Label(self.bf,text="MIN:").pack(side=OTk.LEFT)
        self.minI=OTk.StringVar()
        try:
            self.minI.set(str(int(sys.argv[2])))
        except:
            self.minI.set(str(self.edfFile.minI))
        self.minIentry=OTk.Entry(self.bf, textvariable=self.minI)
        self.minIentry.bind('<KeyPress-Return>', self.change)
        self.minIentry.pack(side=OTk.LEFT)

        OTk.Label(self.bf,text="   MAX:").pack(side=OTk.LEFT)
        self.maxI=OTk.StringVar()
        try:
            top=int(sys.argv[3])
        except:
            top=self.edfFile.minI+(self.edfFile.maxI-self.edfFile.minI)/10.
        self.maxI.set(str(top))
        self.maxIentry=OTk.Entry(self.bf, textvariable=self.maxI)
        self.maxIentry.bind('<KeyPress-Return>', self.change)
        self.maxIentry.pack(side=OTk.LEFT)

        self.Update = OTk.Button(self.bf, text="Update", command=self.change).pack(side=OTk.LEFT)
        OTk.Button(self.bf, text="Reset", command=self.ogl.reset).pack(side=OTk.LEFT)
        self.bf.pack(side=OTk.TOP,expand=0,fill=OTk.X)
        helpframe=OTk.Frame()
        OTk.Label(helpframe,text="Left mouse button to translate, Right to zoom").pack(side=OTk.BOTTOM)
        helpframe.pack(side=OTk.BOTTOM)




    def SetupTexture(self):
        self.makeImage()
        OGL.glPixelStorei(OGL.GL_UNPACK_ALIGNMENT, 1)
##              glTexImage2D(GL_TEXTURE_2D, 0, 3, self.imageWidth, self.imageHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE,  self.image)
        s = self.image.tostring()
        print(len(s),self.imageWidth*self.imageHeight)
        print(self.image.min(),self.image.max())
        OGL.glTexImage2D(OGL.GL_TEXTURE_2D, 0, OGL.GL_RGB, self.imageWidth, self.imageHeight, 0, OGL.GL_RGB ,OGL.GL_UNSIGNED_BYTE,  self.image)
##              glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP)
##              glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP)
        OGL.glTexParameterf(OGL.GL_TEXTURE_2D, OGL.GL_TEXTURE_WRAP_S, OGL.GL_REPEAT)
        OGL.glTexParameterf(OGL.GL_TEXTURE_2D, OGL.GL_TEXTURE_WRAP_T, OGL.GL_REPEAT)
        OGL.glTexParameterf(OGL.GL_TEXTURE_2D, OGL.GL_TEXTURE_MAG_FILTER, OGL.GL_NEAREST)
        OGL.glTexParameterf(OGL.GL_TEXTURE_2D, OGL.GL_TEXTURE_MIN_FILTER, OGL.GL_NEAREST)
        OGL.glTexEnvf(OGL.GL_TEXTURE_ENV, OGL.GL_TEXTURE_ENV_MODE, OGL.GL_DECAL)
        OGL.glEnable(OGL.GL_TEXTURE_2D)
        OGL.glShadeModel(OGL.GL_FLAT)







    def __init__(self):
        try:
            self.edfFile = edfFile(sys.argv[1])
        except:
            sys.stderr.write("usage: %s edf_file\n"%(sys.argv[0]))
            raise
        self.imageWidth = self.edfFile.rows
        self.imageHeight = self.edfFile.cols

        self.SetupWindow()
        self.SetupTexture()
        self.ogl.mainloop()

if __name__ == '__main__':
    checker()
