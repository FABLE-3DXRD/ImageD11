#!/usr/bin/env python

from __future__ import print_function
try:
    from Tkinter import *
    from tkMessageBox import showinfo
    import tkSimpleDialog as simpledialog
    import tkFileDialog as filedialog
    import Dialog as dialog
except:
    from tkinter import *
    import tkinter.simpledialog as simpledialog
    import tkinter.filedialog as filedialog
    from tkinter.messagebox import showinfo
    import tkinter.dialog as dialog

from PIL import ImageTk,Image

import os,sys

import numpy as np

"""
TO DO

Scroll with the keyboard arrows
Zoom in and out with Page up and Page down (or others)
Draw the box via the keyboard
Implement a status bar and cursor
"""


def preprocess(image, xxx_todo_changeme,maxi=None,mini=None):
    """ Makes data into 8 bit array for plotting """
    (scalex,scaley) = xxx_todo_changeme
    MAXW=4096
    MAXH=4096
    MINSIZE=1
    DEFAULT_HEIGHT=100
    MAXSIZE=max(MAXW,MAXH)
    assert len(image.shape) in (1, 2) or \
           len(image.shape) == 3 and image.shape[2] == 3, \
           "image not correct format"
    if maxi==None:
        themin = float(np.minimum.reduce(np.ravel(image)))
        themax = float(np.maximum.reduce(np.ravel(image)))
    else:
        themin = float(mini)
        themax = float(maxi)
    if len(image.shape) == 1:
        len_x = image.shape[0]
        ys = ((image - themin)/(themax-themin)*(DEFAULT_HEIGHT-1)).astype('B')
        image = (zeros((DEFAULT_HEIGHT, len_x))+255).astype('B')
        for x in range(len_x):
            image[DEFAULT_HEIGHT-1-ys[x],len_x-x-1] = 0
        image = np.transpose(image)
    elif image.dtype.char != 'b':
        try:
            pass  ## image.savespace(0)
            image = 255 * (image - themin) / (themax-themin)
        except:
            print("Exception",themax,themin)
        image = np.where(image<256,image,255)
        image = np.where(image>0,image,0).astype('B')

    len_x, len_y = image.shape[:2]
#    print "Image dimensions...  x=",len_x,"  y=",len_y
    if scalex < MINSIZE:
        scalex = MINSIZE
    if scaley < MINSIZE:
        scaley = MINSIZE
    if scalex > MAXSIZE:
        scalex = MAXSIZE
    if scaley > MAXSIZE:
        scaley = MAXSIZE
    return image, (int(scalex), int(scaley))


def NumerictoImage( data, xxx_todo_changeme1,maxi=None,mini=None):
    (scalex,scaley) = xxx_todo_changeme1
    image, (scalesx,scalesy) = preprocess(data, (scalex,scaley),maxi,mini)
    width, height = image.shape[:2]
    if len(image.shape) == 3:
        mode = rawmode = 'RGB'
        bits = np.transpose(image, (1,0,2)).tostring()
    else:
        mode = rawmode = 'L'
        bits = np.transpose(image, (1,0)).tostring()
    image2 = Image.frombytes(mode, (width, height),
                                      bits, "raw", rawmode)
    image2=image2.resize((int(scalesx),int(scalesy)))
    return image2,(scalesx,scalesy)






class rubber(Frame):

    def __init__(self, datafile=None, bkgfile=None, master=None):
        self.master=master
        Frame.__init__(self, master)
        if datafile==None:
            datafile=filedialog.askopenfilename(initialdir=os.getcwd())
        if type(datafile)==type("string"):
            self.datafile = datafile
            dataobj=openimage(datafile)
            self.data=dataobj.data.astype(int)
            try:
                self.omega=float(dataobj.header["Omega"])
            except:
                self.omega = None
        else:
            if type(datafile)==type(np.array([10,11,12])):
                self.data=datafile
                self.datafile="unknown0000.edf"
        pass  ## self.data.savespace(0)
        if type(bkgfile)==type("string"):
            self.bkgfile=bkgfile
            bkgobj=openimage(bkgfile)
            self.bkg=bkgobj.data.astype(int)
            self.data=self.data.astype(int)-self.bkg
            print("Got your background from",bkgfile,self.bkg.shape)
        else:
            self.bkg = None
        Pack.config(self,expand=1,fill=BOTH)
        self.b=Frame(self)
        Button(self.b,text="Quit",command=self.close).pack(side=RIGHT)
        Button(self.b,text="Toggle log scaled",command=self.togglelog).pack(side=RIGHT)
        Button(self.b,text="Set Maximum",command=self.setmax).pack(side=RIGHT)
        Button(self.b,text="Set Minimum",command=self.setmin).pack(side=RIGHT)
        Button(self.b,text="Autorange",command=self.autorange).pack(side=RIGHT)
        Button(self.b,text="ZoomIn",command=self.zoomin).pack(side=LEFT)
        Button(self.b,text="ZoomOut",command=self.zoomout).pack(side=LEFT)
        Button(self.b,text="ShowROI",command=self.showroi).pack(side=LEFT)
        self.b.pack(side=BOTTOM)

        self.b1=Frame(self)
        Button(self.b1,text="Background",command=self.rdbkg).pack(side=LEFT)
        Button(self.b1,text="Rock",command=self.makerockingcurve).pack(side=LEFT)
        Button(self.b1,text="Prev",command=self.prev).pack(side=LEFT)
        Button(self.b1,text="Next",command=self.mynext).pack(side=LEFT)
        Label(self.b1,text="File number=").pack(side=LEFT)
        self.filenum=StringVar()
        Entry(self.b1,textvariable=self.filenum).pack(side=LEFT)
        Button(self.b1,text="Jump",command=self.jump).pack(side=LEFT)
        Button(self.b1,text="Open",command=self.file).pack(side=LEFT)
        Button(self.b1,text="Pars",command=self.rdpars).pack(side=LEFT)
        Button(self.b1,text="UBIS",command=self.rdubis).pack(side=LEFT)
        self.b1.pack(side=BOTTOM)
        try:
            self.filenum.set("%d"%(int(self.datafile[-8:-4])))
        except:
            try:
                self.filenum.set("%d"%(int(self.datafile[-9:-5])))
            except:
                try:
                    self.filenum.set("%d"%(int(self.datafile[-4:])))
                except:
                    self.filenum.set("0")

        self.status=Label(self,text=self.datafile)
        self.status.pack(side=BOTTOM)
        self.statistics=Label(self,text=self.getstats())
        self.statistics.pack(side=BOTTOM)

        self.key=Label(self,bg='green',bd=0)
        self.key.pack(side=BOTTOM)
        self.keyt=Label(self,bd=0)
        self.keyt.pack(side=BOTTOM,fill=BOTH)

        self.xscrollbar=Scrollbar(self,orient=HORIZONTAL)
        self.xscrollbar.pack(side=BOTTOM,fill=X)
        self.yscrollbar=Scrollbar(self,orient=VERTICAL)
        self.yscrollbar.pack(side=RIGHT,fill=Y)

#        print image.size
        self.scale=[self.data.shape[0],self.data.shape[1]]
        self.zoom=1.0
        self.log=0
        self.canvasObject = Canvas(self,width=512,height=512,
              xscrollcommand=self.xscrollbar.set,yscrollcommand=self.yscrollbar.set,
              scrollregion=(0, 0, self.data.shape[0], self.data.shape[1]))
        self.canvasObject.pack(side=LEFT,fill="both",expand=1)
        self.xscrollbar.config(command=self.canvasObject.xview)
        self.yscrollbar.config(command=self.canvasObject.yview)

        # this is a "tagOrId" for the rectangle we draw on the canvas
        self.rubberbandBox = None
        self.piccy=None
        self.startx = None
        self.mini=self.maxi=None
        while self.scale[0]*self.zoom < 256:
            self.zoom*=2.
        while self.scale[0]*self.zoom > 512:
            self.zoom/=2.
        self.autorange()
        self.viewinself()
        step=float(self.maxi-self.mini)/self.scale[0]
        l=[]
        x=float(self.mini)
        while x <= self.maxi:
            l.append(x)
            x+=step
        self.legend=np.transpose(np.array([l],float))
        self.legendimage,(scalesx,scalesy)=NumerictoImage(self.legend,(self.scale[0],10))
        self.legendimage=ImageTk.PhotoImage(self.legendimage)
        self.key.configure(image=self.legendimage)
        self.key.update()
        stri="<---%f          %f--->"%(self.mini,self.maxi)
        self.keyt.configure(text=stri)
        self.keyt.update()


        # and the bindings that make it work..
        self.bind_all("<Prior>", lambda e : self.zoomout())
        self.bind_all("<Next>" ,  lambda e : self.zoomin())
        self.bind("q", self.printkey )
        Widget.bind(self.canvasObject, "<Button-1>", self.mouseDown)
        Widget.bind(self.canvasObject, "<Button1-Motion>", self.mouseMotion)
        Widget.bind(self.canvasObject, "<Button1-ButtonRelease>", self.mouseUp)
        Widget.bind(self.canvasObject, "<Button-3>", self.get_hkl)
        self.parameters={}

    def file(self):
        self.datafile=filedialog.askopenfilename(initialdir=os.getcwd())
        self.readdata()


    def rdpars(self,parfile=None):
        if parfile==None:
            parfile=filedialog.askopenfilename(initialdir=os.getcwd())
        from ImageD11 import parameters
        o=parameters.parameters()
        o.loadparameters(parfile)
        self.parameters=o.parameters

    def rdubis(self,ubifile=None):
        from ImageD11 import indexing
        if ubifile==None:
            ubifile=filedialog.askopenfilename(initialdir=os.getcwd())
        self.ubisread = indexing.readubis(ubifile)
        print("Read into self.ubisread",self.ubisread)

    def rdbkg(self, bkgfile=None):
        if bkgfile == None:
            bkgfile=filedialog.askopenfilename(initialdir=os.getcwd())
        self.bkgfile=bkgfile
        bkgobj=openimage(bkgfile)
        self.bkg=bkgobj.data.astype(int)
        self.data=self.data.astype(int)-self.bkg
        print("Got your background from",bkgfile,self.bkg.shape)

    def jump(self):
        number=int(self.filenum.get())
        self.datafile=self.name( number )
        self.readdata()

    def name(self,num):
        try:
            current=int(self.datafile[-8:-4])
        except:
            current=int(self.datafile[-9:-5])
        print("current is",current)
        return "%s%04d%s"%(self.datafile[:-8],current+num,
                self.datafile[-4:])

    def mynext(self):
        self.datafile=self.name(1)
        self.readdata()

    def prev(self):
        self.datafile=self.name(-1)
        self.readdata()

    def readdata(self):
        dataobj=openimage(self.datafile)
        self.status.config(text=self.datafile)
        self.data=dataobj.data.astype(int)
        if self.bkg is not None:
            try:
                self.data=self.data.astype(int)-self.bkg
            except:
                print("Failed to subtract bkg",self.bkg.shape,self.data.shape)
        try:
            self.omega=float(dataobj.header["Omega"])
        except:
            self.omega = None
        self.statistics.configure(text=self.getstats())
        try:
            self.filenum.set("%d"%(int(self.datafile[-8:-4])))
        except:
            self.filenum.set("%d"%(int(self.datafile[-9:-5])))
        self.viewinself()

    def printkey(self,event):
        print("Got a key")
        print(event.char)

    def setmax(self):
        try:
            self.maxi=int(simpledialog.askinteger("Maximum data value?","Maximum data value?"))
            self.viewinself()
            stri="<---%f          %f--->"%(self.mini,self.maxi)
            self.keyt.configure(text=stri)
            self.keyt.update()
        except:
            pass
    def setmin(self):
        try:
            self.mini=int(simpledialog.askinteger("Minimum data value?","Minimum data value?"))
            self.viewinself()
            stri="<---%f          %f--->"%(self.mini,self.maxi)
            self.keyt.configure(text=stri)
            self.keyt.update()
        except:
            pass

    def getstats(self):
        # Convert to float to avoid overflow issues
        t=self.data.astype(float)
        sumi=np.sum(np.ravel(t))
        sumisq=np.sum(np.ravel(t*t))
        npixels=t.shape[0]*t.shape[1]
        self.average=sumi/npixels
        self.variance=np.sqrt(sumisq/npixels - self.average*self.average)
        self.mindata=np.minimum.reduce(np.ravel(self.data))
        self.maxdata=np.maximum.reduce(np.ravel(self.data))
        return "%d Pixels, Average = %f, Variance = %f, Min = %f, Max = %f"%(
            npixels,self.average,self.variance,self.mindata,self.maxdata)

    def autorange(self):
        newmin=max(np.minimum.reduce(np.ravel(self.data)),self.average-self.variance)
        newmax=min(np.maximum.reduce(np.ravel(self.data)),self.average+self.variance)
        if self.mini!=newmin or self.maxi != newmax:
            self.mini=newmin
            self.maxi=newmax
            self.viewinself()

    def close(self):
        if self.master==None:
            self.destroy()
        else:
            self.master.destroy()


    def zoomin(self):
        if self.scale[0]*2 > 2050:
            print("Zoom limit of ",self.zoom,"reached")
        else:
            self.zoom*=2
            print("Zooming in to scale ",self.zoom,self.zoom*self.data.shape[0])
            self.viewinself()

    def zoomout(self):
        if self.zoom < 0.02:
            print("Zoom out limit of",self.zoom, "reached")
        else:
            self.zoom/=2
            print("Zooming out to scale ",self.zoom)
            self.viewinself()

    def viewinself(self):
        print(self.maxi,self.mini)
        if self.log==0:
            self.image,(self.scale[0],self.scale[1])=NumerictoImage(self.data,self.data.shape,self.maxi,self.mini)
        else:
            self.image,(self.scale[0],self.scale[1])=NumerictoImage(np.log(self.data+1),self.data.shape,self.maxi,self.mini)
        self.scale[0]*=self.zoom
        self.scale[1]*=self.zoom
        self.scale = [int(v) for v in self.scale]
        print("Image size",self.scale)
        self.displayimage = ImageTk.PhotoImage(self.image.resize(self.scale))
        if self.piccy==None:
            self.piccy=self.canvasObject.create_image(2,2,image=self.displayimage,
                                             anchor=NW,tags="picture")
        self.canvasObject.itemconfig(self.piccy,image=self.displayimage)
        self.canvasObject.config(scrollregion=(0,0,self.scale[0],self.scale[1]))
        self.canvasObject.delete(self.rubberbandBox)
        stri="<---%f          %f--->"%(self.mini,self.maxi)
        self.keyt.configure(text=stri)
        self.keyt.update()
        if self.rubberbandBox is not None:
            z=self.zoom
            self.rubberbandBox = self.canvasObject.create_rectangle(
              z*self.startx, z*self.starty, z*self.endx, z*self.endy, outline='green')

    def getbox(self):
        print(self.startx*self.zoom,self.endx*self.zoom)
        print(self.starty*self.zoom,self.endy*self.zoom)
        xl=int(min(self.startx,self.endx)) ; xh = int(max(self.startx,self.endx))
        yl=int(min(self.starty,self.endy)) ; yh = int(max(self.starty,self.endy))
        print("[",xl,":",xh,",",yl,":",yh,"]",self.data.shape)
        if xl!=xh and yl!=yh :
            return self.data[xl:xh,yl:yh].copy()
        else:
            return None

    def showroi(self):
        selection=self.getbox()
        if selection is not None:
            t = Toplevel()
            junk=rubber(selection,master=t)


    def togglelog(self):
        self.keyt.update()
        if self.log==0:
            self.log=1
            self.legendimage,(scalesx,scalesy)=NumerictoImage(np.log(self.legend+1),(self.scale[0],10))
            self.legendimage=ImageTk.PhotoImage(self.legendimage)
            self.key.configure(image=self.legendimage)
            self.key.update()
            stri="<---%f          %f--->"%(self.mini,self.maxi)
            self.keyt.configure(text=stri)
            self.maxi=np.log(self.maxi+1)
            self.mini=np.log(self.mini+1)
            print(self.maxi,self.mini)
            self.viewinself()
        else:
            self.log=0
            self.legendimage,(scalesx,scalesy)=NumerictoImage(self.legend,(self.scale[0],10))
            self.legendimage=ImageTk.PhotoImage(self.legendimage)
            self.key.configure(image=self.legendimage)
            self.key.update()
            stri="<---%f          %f--->"%(self.mini,self.maxi)
            self.keyt.configure(text=stri)
            self.maxi=np.exp(self.maxi)-1
            self.mini=np.exp(self.mini)-1
            print(self.maxi,self.mini)
            self.viewinself()

    def get_hkl(self, event):
        # canvas x and y take the screen coords from the event and translate
        # them into the coordinate system of the canvas object
        xpos = self.canvasObject.canvasx(event.x)/self.zoom
        ypos = self.canvasObject.canvasy(event.y)/self.zoom
        print("xpos,ypos",xpos,ypos)
        from ImageD11 import transform
        tth,eta = transform.compute_tth_eta( np.array([[xpos], [ypos ]]) ,
                                             **self.parameters )
        print("omega:",self.omega,type(self.omega))
        om = np.array([float(self.omega)])
        print("tth,eta,om",tth,eta,om)
        self.gv = transform.compute_g_vectors(tth,eta,om,float(self.parameters['wavelength']), self.parameters['wedge'])
        self.gv = np.transpose(self.gv)
        s=""
        i=0
        for ubi in self.ubisread:
            h=np.dot(ubi,self.gv.T)
            print("i=%d"%(i),"hkl= %.2f %.2f %.2f"%tuple(h))
            i+=1
            s+="grain %3d\n h = %.2f k=%.2f l = %.2f\n"%(i,h[0],h[1],h[2])
        showinfo("Right clicked",
                 "You click at %f %f, tth=%f eta=%f omega=%f\n %s"%(xpos,ypos,tth,eta,om,s))


    def mouseDown(self, event):
        # canvas x and y take the screen coords from the event and translate
        # them into the coordinate system of the canvas object
        self.canvasObject.delete(self.rubberbandBox)
        self.startx = self.canvasObject.canvasx(event.x)/self.zoom
        self.starty = self.canvasObject.canvasy(event.y)/self.zoom

    def mouseMotion(self, event):
        # canvas x and y take the screen coords from the event and translate
        # them into the coordinate system of the canvas object
        x = self.canvasObject.canvasx(event.x)/self.zoom
        y = self.canvasObject.canvasy(event.y)/self.zoom

        if (self.startx != event.x)  and (self.starty != event.y) :
            self.canvasObject.delete(self.rubberbandBox)
            z=self.zoom
            self.rubberbandBox = self.canvasObject.create_rectangle(self.startx*z, self.starty*z, x*z, y*z, outline='green')
            # this flushes the output, making sure that
            # the rectangle makes it to the screen
            # before the next event is handled
            self.update_idletasks()

    def mouseUp(self, event):
        self.endx = self.canvasObject.canvasx(event.x)/self.zoom
        self.endy = self.canvasObject.canvasy(event.y)/self.zoom
        print(self.startx,self.endx,self.starty,self.endy)

    def makerockingcurve(self):
        if self.startx==None:
            pass
        print(self.startx*self.zoom,self.endx*self.zoom)
        print(self.starty*self.zoom,self.endy*self.zoom)
        xl=int(min(self.startx,self.endx)) ; xh = int(max(self.startx,self.endx))
        yl=int(min(self.starty,self.endy)) ; yh = int(max(self.starty,self.endy))
        print("[",xl,":",xh,",",yl,":",yh,"]",self.data.shape)
        if xl!=xh and yl!=yh :
            print("ROI = ",xl,xh,yl,yh)
            stem = simpledialog.askstring("Filename stem?","Filename stem?",initialvalue=self.datafile[:-8])
            first = simpledialog.askinteger("First File to treat","First file to treat")
            last = simpledialog.askinteger("Last File to treat", "Last file to treat")
            extn=".edf"
            s='files %s%04d%s to %s%04d%s\nArea of interest %d %d %d %d'%(stem,first,extn,stem,last,extn,xl,xh,yl,yh)
            print("Making rocking curve for:",s)

            d = dialog.Dialog(None, {'title': 'Make a rocking curve',
                      'text':  "Sorry - this is not working yet!!!" ,
                      'bitmap': dialog.DIALOG_ICON,
                      'default': 0,
                      'strings': ('Do it now please!',
                                  """That's not what I meant!""")})
            return
            print(d.num)
            if d.num==0:
                import rocker
                line=rocker.rock(stem,first,last,extn,xl,xh,yl,yh)
                import twodplotml
                twodplotml.twodplot(parent=Toplevel(),data=line)
                print(line)
            else:
                print("What a waste of time that was!")



import sys
from fabio.openimage import openimage
if len(sys.argv)>2:
    test = rubber(sys.argv[1],bkgfile=sys.argv[2])
elif len(sys.argv)==2:
    test = rubber(sys.argv[1])
else:
    test = rubber()

test.mainloop()
