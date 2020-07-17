
from __future__ import print_function


import OpenGL.GL as GL
import OpenGL.Tk as Tk

class plotvol(Tk.Toplevel):
    def __init__(self, im ):
        Tk.Toplevel.__init__(self)


        self.f = Tk.Frame(self)
        self.bf = Tk.Frame(self)
        self.max = Tk.StringVar() ; self.max.set('60000')
        Tk.Label(self.bf,text="Max").pack(side=Tk.LEFT)
        e = Tk.Entry(self.bf,textvariable=self.max)
        e.pack(side=Tk.LEFT)
        e.bind("<Enter>",self.redraw)
        self.min = Tk.StringVar() ; self.min.set('0')
        Tk.Label(self.bf,text="Min").pack(side=Tk.LEFT)
        e=Tk.Entry(self.bf,textvariable=self.min)
        e.pack(side=Tk.LEFT)
        e.bind("<Enter>",self.redraw)
        self.zoom = Tk.StringVar() ; self.zoom.set('1')
        Tk.Label(self.bf,text="Zoom").pack(side=Tk.LEFT)
        e=Tk.Entry(self.bf,textvariable=self.zoom)
        e.pack(side=Tk.LEFT)
        e.bind("<Enter>",self.redraw)
        self.ow = Tk.StringVar() ; self.ow.set('0')
        Tk.Label(self.bf,text="Ow").pack(side=Tk.LEFT)
        e=Tk.Entry(self.bf,textvariable=self.ow)
        e.pack(side=Tk.LEFT)
        e.bind("<Enter>",self.redraw)
        self.oh = Tk.StringVar() ; self.oh.set('0')
        Tk.Label(self.bf,text="Oh").pack(side=Tk.LEFT)
        e=Tk.Entry(self.bf,textvariable=self.oh)
        e.pack(side=Tk.LEFT)
        e.bind("<Enter>",self.redraw)
        
        self.f.bind("<Enter>",self.redraw)
        self.f.bind("<Button-1>",self.redraw)
        self.ogl = Tk.Opengl(self.f,width = 500, height = 500, double = 1)


        self.bf.pack(side=Tk.TOP)
        self.f.pack(fill=Tk.BOTH,expand=1)
        self.ogl.pack(side=Tk.TOP, fill = Tk.BOTH, expand = 1)
        self.volslice = im


        self.ogl.redraw = self.redraw
        self.ogl.tkTranslate = self.tkTranslate

    def tkTranslate(self, event):
        self.activate()
        
            
        

    def redraw(self, o):
        GL.glClearColor(0., 0., 0., 0)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        GL.glDisable(GL.GL_LIGHTING)
        
        # Set the zoom factor
        z = float(self.zoom.get())
        GL.glPixelZoom(z, z)

        h, w = self.volslice.shape

        vp = GL.glget.glGetFloatv(GL.GL_VIEWPORT)        

        #  [    0.     0.  1215.   855.]
        #  [ origin,       w , h ] 

        # the bottom left corner is the place we start drawing
        # the top right is size * zoom
              
        ow = float(self.ow.get())
        oh = float(self.oh.get())
        GL.glWindowPos2f( ow, oh )


        # Set the color mapping - user is in data space
        mx = float(self.max.get())
        mn = float(self.min.get())


        # RGBA = SCALE*DATA + BIAS = [0,1]
        # 0 = SCALE*mn + BIAS
        # 1 = SCALE*mx + BIAS
        # 1 = SCALE*(mx - mn)
        scale = 1.0/(mx - mn)
        bias = -mn*scale

        #print "scale=",scale
        #print "bias=",bias
        #print "scale*mx + bias",scale*mx+bias
        #print "scale*mn + bias",scale*mn+bias

        GL.glPixelTransferf( GL.GL_RED_SCALE, scale)
        GL.glPixelTransferf( GL.GL_BLUE_SCALE, scale)
        GL.glPixelTransferf( GL.GL_GREEN_SCALE, scale)
        GL.glPixelTransferf( GL.GL_RED_BIAS, bias)
        GL.glPixelTransferf( GL.GL_BLUE_BIAS, bias)
        GL.glPixelTransferf( GL.GL_GREEN_BIAS, bias)
        GL.glDrawPixels( self.volslice.shape[1],
                         self.volslice.shape[0],
                         GL.GL_LUMINANCE,
                         GL.GL_FLOAT,
                         self.volslice )

        GL.glEnable(GL.GL_LIGHTING)
        
        



def main():
    from ImageD11.rsv import readvol
    from fabio.openimage import openimage
    import numpy as np
    import sys
    #v = readvol( sys.argv[1] )
    v = openimage( sys.argv[1] ).data.astype( np.float32 )
    root = plotvol(v)
    root.mainloop()

if __name__=="__main__":
    main()
