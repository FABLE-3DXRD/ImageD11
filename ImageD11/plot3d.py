#!

# This is statement is required by the build system to query build info
if __name__ == '__build__':
	raise Exception


import string
__version__ = string.split('$Revision$')[1]
__date__ = string.join(string.split('$Date$')[1:3], ' ')
__author__ = 'Jon Wright <jpwright@users.sourceforge.net> from example by Tarn Weisner Burton <twburton@users.sourceforge.net>'

try:
    from Numeric import *
except:
    import sys
    print "This demo requires the Numeric extension, sorry."
    sys.exit()


from OpenGL.GL import *
from OpenGL.Tk import *


import sys


class plot3d(Toplevel):
   def __init__(self,parent,data=None):
      Toplevel.__init__(self,parent)
      if data!=None:
         xyz=data.copy()
      else:
         xyz=array([0,0,0])
      self.ps=StringVar()
      self.ps.set('1.')
      self.pointsize=1.
      self.npeaks=xyz.shape[0]

      self.o = Opengl(self, width = 400, height = 400, double = 1)
      self.o.redraw = self.redraw
      self.o.autospin_allowed = 1
      self.o.fovy=5
      self.o.near=1e6
      self.o.far=1e-6
      import math  
      self.o.distance=maximum.reduce(ravel(xyz))*4/math.tan(self.o.fovy*math.pi/180)
      glVertexPointerd(xyz)
      glEnableClientState(GL_VERTEX_ARRAY)
      f=Frame(self)
      Button(f,text="Help",command=self.o.help).pack(side=LEFT)
      Button(f,text="Reset",command=self.o.reset).pack(side=LEFT)
      Button(f,text="Pointsize",command=self.setps).pack(side=LEFT)
      Entry(f,textvariable=self.ps).pack(side=LEFT)
      Button(f,text="Quit",command=self.goaway).pack(side=RIGHT)
      self.dataoff=0
      self.o.pack(side = 'top', expand = 1, fill = 'both')
      f.pack(side=BOTTOM,expand=NO,fill=X)
      Label(self,text="Red=[1,0,0] Green=[0,1,0] Blue=[0,0,1]").pack(side=BOTTOM,expand=NO,fill=X)


   def go(self):
      """
      Allow the toplevel to return a handle for changing data
      """
      self.o.mainloop()       

   def goaway(self):
      print "Called goaway"
      self.o.destroy()
      self.destroy()
      print "Ought to be gone now..."

   def changedata(self,xyz):
      self.xyz=xyz.copy()
      self.npeaks=xyz.shape[0]
      glDisableClientState(GL_VERTEX_ARRAY)
      glVertexPointerd(self.xyz)
      glEnableClientState(GL_VERTEX_ARRAY)
      self.o.tkRedraw()

   def setps(self):
      self.pointsize=float(self.ps.get())
      self.o.tkRedraw()
      
      
   def redraw(self,o):
	glClearColor(0., 0., 0., 0)
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	glOrtho(-1,1,-1,1,-1,1)
	glDisable(GL_LIGHTING)
        glColor3f(1.0, 1.0, 1.0) # white
        glPointSize(self.pointsize)
#	glEnable(GL_POINT_SMOOTH)
        glDrawArrays(GL_POINTS, 0, self.npeaks )
        glBegin(GL_LINE_LOOP)
        glColor3f(1.0, 0.0, 0.0) # red
        glVertex3f(0.,0.,0.)
        glVertex3f(1.,0.,0.)
        glEnd()
        glBegin(GL_LINE_LOOP)
        glColor3f(0.0, 1.0, 0.0) # green
        glVertex3f(0.,0.,0.)
        glVertex3f(0.,1.,0.)
        glEnd()
        glBegin(GL_LINE_LOOP)
        glColor3f(0.0, 0.0, 1.0) # blue
        glVertex3f(0.,0.,0.)
        glVertex3f(0.,0.,1.)
        glEnd()
	glEnable(GL_LIGHTING)


      


   

if __name__=="__main__":

   try:
      lines=open(sys.argv[1],"r").readlines()
   except:
      print "Usage %s gvector_file"%(sys.argv[0])
      sys.exit()

   on=0
   xyz=[]
   for line in lines:
      if on==1:
         try:
            vals=[float(x) for x in line.split()]
            xyz.append( [ vals[0],vals[1],vals[2] ])
         except:
            pass
      if line.find("xr yr zr")>0:
         on=1

   npeaks = len(xyz)

   xyz=array(xyz)
   main(xyz)

