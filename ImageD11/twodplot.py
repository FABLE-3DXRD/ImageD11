# ImageD11_v0.4 Software for beamline ID11
# Copyright (C) 2005  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# This was broadly copied from an example in matplotlib, but then
# modified extensively.
# Should go back to something cleaner later... - JPW 2005


"""
From the matplotlib examples - modified for mouse
"""
import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import Tkinter as Tk
import tkFileDialog, tkMessageBox
import sys,os,time

class data:
    def __init__(self,x,y,d={}):
        self.x=x
        self.y=y
        self.d=d

class twodplot(Tk.Frame):
    def __init__(self,parent=None,data=None,quiet="No"):
        Tk.Frame.__init__(self,parent)
        self.quiet=quiet
        self.f = Figure(figsize=(8,5), dpi=100)
        self.a = self.f.add_subplot(111)
        self.plotitems={}
        self.maxpoints=1000000 # work around slow plotting
        # print data
        if data != None:
            self.plotitems[data[0]]=data[1]
        self.title=None
        self.xr=self.yr=None
        # a tk.DrawingArea
        self.canvas = FigureCanvasTkAgg(self.f, master=self)
        self.canvas.show()
        self.tkc=self.canvas.get_tk_widget()
        self.tkc.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.tkc.bind("<ButtonPress-1>",self.on_down)
        self.tkc.bind("<ButtonRelease-1>",self.on_up)
        self.tkc.bind("<Button1-Motion>",self.on_move)
        self.tkc.bind("<ButtonPress-2>",self.on_2)
        self.tkc.bind("<ButtonPress-3>",self.on_3)
        self.when_down = -1
        self.time_down = 0.1
        self.bindkeys()
        self.rubberbandbox=None
        self.pack_opts={'side':Tk.LEFT, 'padx':'2', 'pady':'2'}
        self.bf1=Tk.Frame(self)
        Tk.Button(master=self.bf1, text='Clear', command=self.clear).pack(self.pack_opts)
        Tk.Button(master=self.bf1, text='Save Plot', command=self.printplot).pack(self.pack_opts)
        Tk.Button(master=self.bf1, text='LogY', command=self.logy).pack(self.pack_opts)
        Tk.Button(master=self.bf1, text='LogX', command=self.logx).pack(self.pack_opts)
# FIXME - buttons for panx/y and zoomx/y
        Tk.Button(master=self.bf1,text='>' ,
                    command=lambda : self.keypress(self.a.panx,1 )  ).pack(self.pack_opts)
        Tk.Button(master=self.bf1,text='<',
                    command=lambda : self.keypress(self.a.panx,-1)  ).pack(self.pack_opts)
        Tk.Button(master=self.bf1,text='^'   ,
                    command=lambda : self.keypress(self.a.pany,-1 ) ).pack(self.pack_opts)
        Tk.Button(master=self.bf1,text='v' ,
                    command=lambda : self.keypress(self.a.pany,1)   ).pack(self.pack_opts)
        self.bf1.pack(side=Tk.BOTTOM)
        self.bf2=Tk.Frame(self)
        Tk.Button(master=self.bf2,text='UnZoomX' ,
                    command=lambda : self.keypress(self.a.zoomx,-1 )  ).pack(self.pack_opts)
        Tk.Button(master=self.bf2,text='ZoomX',
                    command=lambda : self.keypress(self.a.zoomx,1)  ).pack(self.pack_opts)
        Tk.Button(master=self.bf2,text='ZoomY'   ,
                    command=lambda : self.keypress(self.a.zoomy,1 ) ).pack(self.pack_opts)
        Tk.Button(master=self.bf2,text='UnZoomY' ,
                    command=lambda : self.keypress(self.a.zoomy,-1)   ).pack(self.pack_opts)
        Tk.Button(master=self.bf2,text='Autoscale' ,
                    command=lambda  : self.keypress(self.autoscale)   ).pack(self.pack_opts)
        Tk.Button(master=self.bf2,text='Autoscale Y',
                    command=lambda  : self.keypress(self.autoscaley, None )  ).pack(self.pack_opts)
        self.bf2.pack(side=Tk.BOTTOM)
        self.label=Tk.Label(self,text="Plot window messages")
        self.label.pack(side=Tk.BOTTOM,fill=Tk.X, expand=0)
        self.pack(side=Tk.TOP,fill=Tk.BOTH,expand=Tk.YES)
        self.hidden=[]
        self.replot()
        self.xd=None
        self.yd=None

    def printplot(self):
        fn = tkFileDialog.asksaveasfilename(title="File name to print to",
              defaultextension="png")
        print fn
        f,e = os.path.splitext(fn)
        extns=['png','ps','eps','bmp','raw','rgb']
        print e
        if e.lower() in ['.png','.ps','.eps','.bmp','.raw','.rgb']:
            self.update_idletasks() # Try to get screen redrawn
            self.canvas.print_figure(fn, dpi=300, orientation='landscape')
        else:
            tkMessageBox.showinfo("Sorry","I can only make output in these formats"+str(extns))

    def keypress(self,*arg):
        if len(arg)>1:
            arg[0](*arg[1:])
        else:
            arg[0]()
        self.canvas.show()


    def bindkeys(self):
        return
        self.bind_all('<Left>' ,lambda e: self.keypress(self.a.panx,1 )  )
        self.bind_all('<Right>',lambda e: self.keypress(self.a.panx,-1)  )
        self.bind_all('<Up>'   ,lambda e: self.keypress(self.a.pany,-1 ) )
        self.bind_all('<Down>' ,lambda e: self.keypress(self.a.pany,1)   )
        self.bind_all('<Shift-Left>' ,lambda e: self.keypress(self.a.zoomx,-1 )  )
        self.bind_all('<Shift-Right>',lambda e: self.keypress(self.a.zoomx,1)  )
        self.bind_all('<Shift-Up>'   ,lambda e: self.keypress(self.a.zoomy,1 ) )
        self.bind_all('<Shift-Down>' ,lambda e: self.keypress(self.a.zoomy,-1)   )
        self.bind_all('<Next>' , lambda e : self.keypress(self.autoscale)   )
        self.bind_all('<Prior>', lambda e : self.keypress(self.autoscaley, e )  )

    def autoscaley(self,e):
        print dir(self.a.dataLim)
        print self.a.dataLim
        yr=self.a.get_ylim()
        self.a.set_ylim(yr)

    def adddata(self,data):
        """
        Takes a tuple of name, data object
        """
        self.plotitems[data[0]]=data[1]
        if data[0] in self.hidden:
            self.hidden.remove(data[0])
        self.replot()

    def hideall(self):
        self.xr = self.yr = None
        for item in self.plotitems.keys():
            self.hidden.append(item)


    def removedata(self,name):
        try:
            self.plotitems.pop(name)
        except KeyError:
            pass


    def replot(self):
        self.a.clear()
        if self.title!= None:
            self.a.set_title(self.title)
#      b  : blue
#      g  : green
#      r  : red
#      c  : cyan
#      m  : magenta
#      y  : yello
        c = ['g','r','c','m','y','b']
        for name in self.plotitems.keys():
            if name in self.hidden:
                continue
            item=self.plotitems[name]
            #print 'x ', item.x
            #print 'y ', item.y
            #print 'd ', item.d
            #print self.plotitems[name].d
            if item.d.has_key('color'):
                pc=item.d['color']
            else:
                c.append(c[0])
                pc=c.pop(0)
            if item.d.has_key('pointtype'):
                pc=item.d['pointtype']+pc
            else:
                pc="."+pc
            if item.d.has_key("xlabel"):
                self.a.set_xlabel(item.d["xlabel"])
            if item.d.has_key("ylabel"):
                self.a.set_ylabel(item.d["ylabel"])
            if item.d.has_key("title"):
                self.a.set_title(item.d["title"])
            try:
                if  item.d.has_key("plottype"):
                    ret = self.a.hist(item.y,item.x)
                elif item.x.shape[0]>self.maxpoints:
                    if self.quiet=="No":
                        if tkMessageBox.askyesno("Slow plotting workaround","Shall I plot only the first %d points for increased speed?"%(self.maxpoints)):
                            ret = self.a.plot(item.x[:self.maxpoints],item.y[:self.maxpoints],pc)
                        else:
                            ret = self.a.plot(item.x,item.y,pc)
                    else:
                        ret = self.a.plot(item.x[:self.maxpoints],item.y[:self.maxpoints],pc)
                else:
                    ret = self.a.plot(item.x,item.y,pc)
            except:
                print "plotting exception ignored"
        if self.xr!=None:
            self.a.set_xlim(self.xr)
        if self.yr!=None:
            self.a.set_ylim(self.yr)
        self.canvas.show()

    def logy(self):
# FIXME - clip negative values before making logscaled?
        if self.a.yaxis.is_log():
            self.a.set_yscale('linear')
        else:
            self.a.set_yscale('log')
        self.canvas.show()

    def logx(self):
# FIXME - clip negative values before making logscaled?
        if self.a.xaxis.is_log():
            self.a.set_xscale('linear')
        else:
            self.a.set_xscale('log')
        self.canvas.show()

    def on_3(self,event):
        self.autoscale()

    def autoscale(self):
        self.a.cla()
        self.xr = self.yr = None
        self.replot()

    def clear(self):
        self.plotitems={}
        self.replot()

    def on_2(self,event):
        try:
            height = self.f.bbox.height()
            x, y = event.x, height-event.y
            (xd,yd)= self.a.transData.inverse_xy_tup( (x,y) )
        except:
            height = self.f.bbox.height
            x, y = event.x, height-event.y
            (xd,yd)= self.a.transData.inverted().transform((x,y))
        self.label.config(text="Clicked at x=%f y=%f"%(xd,yd))

    # Callback functions for mouse
    def on_down(self,event):
        # get the x and y coords, flip y from top to bottom
        self.when_down = time.time()
        try:
            height = self.f.bbox.height()
            x, y = event.x, height-event.y
            (self.xd,self.yd)= self.a.transData.inverse_xy_tup( (x,y) )
        except:
            height = self.f.bbox.height
            x, y = event.x, height-event.y
            (self.xd,self.yd)= self.a.transData.inverted().transform((x,y))

        # transData transforms data coords to display coords.  Use the
        # inverse method to transform back
        # print "print mouse down at", t, val
        # rubber banding:
        if self.rubberbandbox!=None: self.tkc.delete(self.rubberbandbox)
        self.startx=self.tkc.canvasx(event.x)
        self.starty=self.tkc.canvasx(event.y)

    def on_move(self,event):
        x = self.tkc.canvasx(event.x)
        y = self.tkc.canvasy(event.y)
        if (self.startx != event.x)  and (self.starty != event.y) :
            if self.rubberbandbox!=None:
                self.tkc.delete(self.rubberbandbox)
            self.rubberbandbox = self.tkc.create_rectangle(self.startx, self.starty, x, y, outline='green')
            # this flushes the output, making sure that
            # the rectangle makes it to the screen
            # before the next event is handled

    def on_up(self,event):
        # get the x and y coords, flip y from top to bottom
        if self.xd==None:
            return
        if time.time()-self.when_down < self.time_down and self.when_down>0:
            return
        self.tkc.delete(self.rubberbandbox)
        try:
            height = self.f.bbox.height()
            x, y = event.x, height-event.y
            (self.xu,self.yu) = self.a.transData.inverse_xy_tup( (x,y) )
        except:
            height = self.f.bbox.height
            x, y = event.x, height-event.y
            (self.xu,self.yu)= self.a.transData.inverted().transform((x,y))

        # transData transforms data coords to display coords.  Use the
        # inverse method to transform back
        if self.xu != self.xd and self.yu != self.yd:
            # rescale
            xr=[self.xd,self.xu];xr.sort()
            yr=[self.yd,self.yu];yr.sort()
            self.xr=xr
            self.yr=yr
            self.a.set_xlim(xr)
            self.a.set_ylim(yr)
            self.canvas.show()


if __name__=="__main__":
    import epffile, powbase, mcadata, ciidata
    if len(sys.argv)<3:
        import numpy as np
        from math import pi
        print "Usage: %s filename format"%(sys.argv[0])
        x=np.arange(0.0,3.0,0.01)
        dat=epffile.powderdata(x,
                               np.sin(2*pi*x)+5,
                               np.sqrt(sin(2*pi*x)+5),
                               { "title":"sin x" })
    else:
        try:
            if sys.argv[2]=="powbase":
                dat=powbase.powbase(sys.argv[1])
            if sys.argv[2]=="epf":
                dat=epffile.epffile(sys.argv[1])
            if sys.argv[2]=="mca":
                dat=mcadata.mcadata(sys.argv[1])
        except:
            print "Could not read your file %s" % (sys.argv[1])
            raise

    root = Tk.Tk()
    root.wm_title("Two dimensional plotting")
    p=twodplot(root,data=dat)

    Tk.mainloop()
