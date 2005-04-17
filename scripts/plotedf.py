
#!

# This is statement is required by the build system to query build info
if __name__ == '__build__':
	raise Exception


## example python/pyopengl script to do tiled texturemapping.
## By david konerding (dek@cgl.ucsf.edu)

import sys
#from Image import *
from OpenGL.GL import *
from OpenGL.Tk import *
try:
	from Numeric import *
except:
	print "This demo requires the Numeric extension, sorry."
	sys.exit()
import math

const = math.pi

class edfFile:
	"""
	Absolute minimum edf file format, UInt16, rows and columns
	"""
	def __init__(self,filename):
		f=open(filename,"rb")
		h=f.read(1024)
		rows=None
		cols=None
		for item in h.split(";"):
			if item.find("Dim_1")>=0:
				rows=int(item.split()[-1])
			if item.find("Dim_2")>=0:
				cols=int(item.split()[-1])
		f.seek(-rows*cols*2,2)
		self.rows=rows
		self.cols=cols
		self.data=array(fromstring(f.read(rows*cols*2),UInt16),savespace=1)
		self.minI=minimum.reduce(self.data)
		self.maxI=maximum.reduce(self.data)
		print "Opened",filename,"max=",self.maxI,"min=",self.minI


class checker:

	def makeImage(self):
		try:
			mi=int(self.minI.get())
			mx=int(self.maxI.get())
		except:
			mi=self.edfFile.minI
			mx=self.edfFile.maxI
		d=clip(self.edfFile.data,mi,mx) # makes a clipped copy
		print mx,mi,maximum.reduce(d),minimum.reduce(d),d.typecode(),
		d=255.*(d-mi)/(mx-mi)
		d=ravel(d.astype(UInt8))
		print maximum.reduce(d),minimum.reduce(d),d.typecode()
		self.image=d.tostring()
		self.imageWidth = self.edfFile.rows
		self.imageHeight = self.edfFile.cols

		
	def display(self, event=None):
		glClearColor( .7, 0.8, 0.9, 0)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		glBegin(GL_QUADS)


		glTexCoord2f(0.0, 0.0);		glVertex3f(-10.0,-10., 0.0)
		glTexCoord2f(0.0, 1.0);		glVertex3f(-10.0, 10., 0.0)
		glTexCoord2f(1.0, 1.0);		glVertex3f( 10. , 10., 0.0)
		glTexCoord2f(1.0, 0.0);		glVertex3f( 10., -10., 0.0)

		glEnd()
		glFlush()

	def change(self,event=None):
		self.SetupTexture()
		self.ogl.tkRedraw()

	def SetupWindow(self):



		self.OglFrame = Frame()
		self.OglFrame.pack(side = 'top', expand=1 ,fill='both')

		self.ogl = Opengl(master=self.OglFrame, width = 500, height = 500, double = 1)


		self.ogl.pack(side = 'top', expand = 1, fill = 'both')
		self.ogl.distance=100.
		self.ogl.fovy=11.5
#		self.ogl.set_centerpoint(0, 0, 0)
                self.ogl.autospin_allowed=1
		self.ogl.redraw = self.display


		# Control buttons for scaling
		self.bf=Frame()
		self.QuitButton = Button(self.bf, {'text':'Quit'})
		self.QuitButton.bind('<ButtonRelease-1>', sys.exit)
		self.QuitButton.pack(side=RIGHT)

		Label(self.bf,text="MIN:").pack(side=LEFT)
		self.minI=StringVar()
		self.minI.set(str(self.edfFile.minI))
		self.minIentry=Entry(self.bf, textvariable=self.minI)
		self.minIentry.bind('<KeyPress-Return>', self.change)
		self.minIentry.pack(side=LEFT)

		Label(self.bf,text="   MAX:").pack(side=LEFT)
		self.maxI=StringVar()
		top=self.edfFile.minI+(self.edfFile.maxI-self.edfFile.minI)/10.
		self.maxI.set(str(top))
		self.maxIentry=Entry(self.bf, textvariable=self.maxI)
		self.maxIentry.bind('<KeyPress-Return>', self.change)
		self.maxIentry.pack(side=LEFT)

		self.Update = Button(self.bf, text="Update", command=self.change).pack(side=LEFT)
		Button(self.bf, text="Reset", command=self.ogl.reset).pack(side=LEFT)
		self.bf.pack(side=TOP,expand=0,fill=X)		
		helpframe=Frame()
		Label(helpframe,text="Left mouse button to translate, Right to zoom").pack(side=BOTTOM)
		helpframe.pack(side=BOTTOM)




	def SetupTexture(self):
		self.makeImage()
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
##		glTexImage2D(GL_TEXTURE_2D, 0, 3, self.imageWidth, self.imageHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE,  self.image)
		glTexImage2D(GL_TEXTURE_2D, 0, 3, self.imageWidth, self.imageHeight, 0, GL_LUMINANCE ,GL_UNSIGNED_BYTE,  self.image)
##		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP)
##		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL)
		glEnable(GL_TEXTURE_2D)
		glShadeModel(GL_FLAT)


	def __init__(self):
		try:
			self.edfFile = edfFile(sys.argv[1])
		except:
			sys.stderr.write("usage: %s edf_file\n"%(sys.argv[0]))
			raise
		self.SetupWindow()
		self.SetupTexture()
		self.ogl.mainloop()

if __name__ == '__main__':
	checker()



