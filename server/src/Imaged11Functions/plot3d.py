#!

"""
from example by Tarn Weisner Burton <twburton@users.sourceforge.net> in pyopengl
"""
# This is statement is required by the build system to query build info
if __name__ == '__build__':
    raise Exception


import string
__version__ = string.split('$Revision$')[1]
__date__ = string.join(string.split('$Date$')[1:3], ' ')
__author__ = 'Jon Wright <jpwright@users.sourceforge.net> from example by Tarn Weisner Burton <twburton@users.sourceforge.net>'

try:
    import Numeric
except:
    import sys
    print "This demo requires the Numeric extension, sorry."
    sys.exit()


import OpenGL.GL as GL
import OpenGL.Tk as Tk


import sys


class plot3d(Tk.Toplevel):
    def __init__(self,parent,data=None,lines=None):
        """
        Data would be your observed g-vectors. Lines will
        be a computed lattice
        """
        Tk.Toplevel.__init__(self,parent)
        if data!=None:
            xyz=data.copy()
        else:
            xyz=Numeric.array([0,0,0])
        self.ps=Tk.StringVar()
        self.ps.set('1.')
        self.pointsize=1.
        self.npeaks=xyz.shape[0]

        self.o = Tk.Opengl(self, width = 400, height = 400, double = 1)
        self.o.redraw = self.redraw
        self.o.autospin_allowed = 1
        self.o.fovy=5
        self.o.near=1e6
        self.o.far=1e-6
        import math
        self.o.distance=Numeric.maximum.reduce(Numeric.ravel(xyz))*4/math.tan(self.o.fovy*math.pi/180)
        print type(xyz),xyz.typecode(),xyz.shape
        GL.glVertexPointerd(xyz)
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
        f=Tk.Frame(self)
        Tk.Button(f,text="Help",command=self.o.help).pack(side=Tk.LEFT)
        Tk.Button(f,text="Reset",command=self.o.reset).pack(side=Tk.LEFT)
        Tk.Button(f,text="Pointsize",command=self.setps).pack(side=Tk.LEFT)
        Tk.Entry(f,textvariable=self.ps).pack(side=Tk.LEFT)
        Tk.Button(f,text="Quit",command=self.goaway).pack(side=Tk.RIGHT)
        self.dataoff=0
        self.o.pack(side = 'top', expand = 1, fill = 'both')
        f.pack(side=Tk.BOTTOM,expand=Tk.NO,fill=Tk.X)
        Tk.Label(self,text="Red=[1,0,0] Green=[0,1,0] Blue=[0,0,1]").pack(side=Tk.BOTTOM,expand=Tk.NO,fill=Tk.X)


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
        GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
        GL.glVertexPointerd(self.xyz)
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
        self.o.tkRedraw()

    def setps(self):
        self.pointsize=float(self.ps.get())
        self.o.tkRedraw()


    def redraw(self,o):
        GL.glClearColor(0., 0., 0., 0)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glOrtho(-1,1,-1,1,-1,1)
        GL.glDisable(GL.GL_LIGHTING)
        GL.glColor3f(1.0, 1.0, 1.0) # white
        GL.glPointSize(self.pointsize)
#       GL.glEnable(GL.GL_POINT_SMOOTH)
        GL.glDrawArrays(GL.GL_POINTS, 0, self.npeaks )
        GL.glBegin(GL.GL_LINE_LOOP)
        GL.glColor3f(1.0, 0.0, 0.0) # red
        GL.glVertex3f(0.,0.,0.)
        GL.glVertex3f(1.,0.,0.)
        GL.glEnd()
        GL.glBegin(GL.GL_LINE_LOOP)
        GL.glColor3f(0.0, 1.0, 0.0) # green
        GL.glVertex3f(0.,0.,0.)
        GL.glVertex3f(0.,1.,0.)
        GL.glEnd()
        GL.glBegin(GL.GL_LINE_LOOP)
        GL.glColor3f(0.0, 0.0, 1.0) # blue
        GL.glVertex3f(0.,0.,0.)
        GL.glVertex3f(0.,0.,1.)
        GL.glEnd()
        GL.glEnable(GL.GL_LIGHTING)







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
    xyz=Numeric.array(xyz)
    o=plot3d(None,data=xyz)
    o.mainloop()
