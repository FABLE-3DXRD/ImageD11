import ipywidgets, pylab as pl
from ImageD11 import transformer
from tkinter import *
from tkinter import filedialog


class FitGeom(transformer.transformer):
    
    "Gives an IPython notebook UI for the experimental calibrations"
    try:
        if get_ipython().__class__.__name__ == 'ZMQInteractiveShell':
            interactive = True
        else:
            interactive = False
    except NameError:
        interactive = False
    
    def __init__(self):
        transformer.transformer.__init__(self)
        self.vars = "y_center z_center distance tilt_y tilt_z".split()
        self.steps = (1,         1,       100,     0.01, 0.01, 0)
        self.nv = len(self.vars)

    def __parCallBack( self, arg ):
        self.parameterobj.parameters.update( { arg['owner'].description : arg['new'] } )
        self.__drawPlot()

    def __fixVaryCallBack( self, arg ):
        name = arg['owner'].description.split(" ")[1]
        vars = self.getvars()
        if arg.new and not (name in self.vars):
            vars.append( name )
        if name in self.vars and not arg.new:
            vars.remove(name)
        self.parameterobj.set_varylist(vars)

    def __fitCallBack(self, arg):
        """ fit call back - runs fit """
        lo, hi = self.ax1.get_xlim()
        self.fit( lo, hi )
        self.__updateGui()

    def __drawPlot(self):
        tth, eta = self.compute_tth_eta()
        self.pt1.set_data( tth, eta )

    def __loadCallBack(self, arg):
        rootTk = Tk()
        rootTk.withdraw()
        rootTk.call('wm', 'attributes', '.', '-topmost', True)
        filename = filedialog.askopenfilename()
        get_ipython().run_line_magic('gui', 'tk')
        if not filename == '':
            self.loadfileparameters(filename)
            self.__updateGui()

    def __saveCallBack(self, arg):
        rootTk = Tk()
        rootTk.withdraw()
        rootTk.call('wm', 'attributes', '.', '-topmost', True)
        filename = filedialog.asksaveasfilename()
        get_ipython().run_line_magic('gui', 'tk')
        if not filename == '':
            self.parameterobj.saveparameters(filename)

    def __updateGui(self):
        for i, pname in enumerate(self.vars):
            self.layout[i,0].value = self.parameterobj.get(pname)
        self.__drawPlot()
    
    def __drawWidgets(self):
        nv = self.nv
        vars = self.vars
        steps = self.steps
        self.layout = ipywidgets.GridspecLayout(nv+1,3)
        for i,( pname, pstep ) in enumerate( zip( vars, steps ) ) :
            self.layout[i,0] = ipywidgets.FloatText( description=pname, 
                value = self.parameterobj.parameters.get(pname),
                step=pstep)
            self.layout[i,0].observe( self.__parCallBack , names='value' )
            self.layout[i,1] = ipywidgets.ToggleButton( description="Vary "+pname, 
                value = pname in self.getvars() )
            self.layout[i,1].observe( self.__fixVaryCallBack, names='value' )

        self.layout[nv,0] = ipywidgets.FloatText( description='fit_tolerance', 
                value = self.parameterobj.parameters.get("fit_tolerance"), step=0,)
        self.layout[nv,0].observe( self.__parCallBack , names='value' )

        self.layout[nv,1] = ipywidgets.Button(description="Run Fit (blocks)")
        self.layout[nv,1].on_click( self.__fitCallBack )

        self.layout[0,2] = ipywidgets.Button(description="Load parameters")
        self.layout[0,2].on_click( self.__loadCallBack )

        self.layout[1,2] = ipywidgets.Button(description="Save to a file")
        self.layout[1,2].on_click( self.__saveCallBack )

    def fitGui(self):
        if self.__class__.interactive:
            self.fig1 = pl.figure(1, figsize=(9,6))
            self.ax1 = self.fig1.add_subplot()
            tth, eta = self.compute_tth_eta()
            self.addcellpeaks()
            self.pt1, = self.ax1.plot( tth, eta, ",")
            self.ax1.set(xlabel="tth", ylabel="eta")
            self.ax1.plot( self.theorytth, [0,]*len(self.theorytth), "r|", ms=360, alpha=0.2 )
            # Add controls
            self.__drawWidgets()
            display(self.layout)
        else:
            print('Sorry, this Gui works only in IPython notebooks!')



