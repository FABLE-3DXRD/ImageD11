
import ipywidgets, pylab as pl, numpy as np, time
from ImageD11 import transformer
from IPython.display import display


class FitGeom(transformer.transformer):
    
    "Gives an IPython notebook UI for the experimental calibrations"

    
    def __init__(self):
        transformer.transformer.__init__(self)
        names = "cell__a cell__b cell__c cell_alpha cell_beta cell_gamma cell_lattice_[P,A,B,C,I,F,R]".split() +\
        "distance wavelength y_size z_size y_center z_center tilt_x tilt_y tilt_z t_x t_y t_z omegasign".split() +\
        "o11 o12 o21 o22 fit_tolerance".split()
        self.names = names
        self.gui = None

    def parCallBack( self, arg, name ):
        self.parameterobj.parameters.update( { name : arg['new'] } )
        self.drawPlot()

    def fixVaryCallBack( self, arg, name ):
        name = arg['owner'].description.split(" ")[1]
        vary = arg.new
        if arg.new:
            if name not in self.parameterobj.varylist:
                self.parameterobj.varylist.append(name)
        else:
            self.parameterobj.varylist.remove(name)

    def fitCallBack(self, arg):
        """ fit call back - runs fit """
        lo, hi = self.ax1.get_xlim()
        self.fit( lo, hi )
        self.updateGui()

    def drawPlot(self):
        self.compute_tth_eta()
        self.addcellpeaks()
        self.pt1.set_data( self.colfile.tth, self.colfile.eta )
        self.pt2.set_data(  self.theorytth, np.full_like(self.theorytth,0) )
        self.ax1.set(title='fitted')
        self.fig1.canvas.draw()
        self.fig1.canvas.flush_events()

    def updateGui(self):
        for i, name in enumerate( self.parameterobj.varylist ):
            self.valuewidgets[name].value = self.parameterobj.parameters[name]
        self.drawPlot()
    
    def drawWidgets(self):
        values = {}
        fix = {}
        for i,name in enumerate( self.names ):
            value = self.parameterobj.parameters[name]
            values[name] = ipywidgets.FloatText( 
                description="", 
                value = value,
                layout=ipywidgets.Layout(positioning='left'))
            values[name].observe( lambda a, name=name: self.parCallBack(a, name),
                                 names='value' )
            if self.parameterobj.par_objs[name].can_vary:
                fix[name] = ipywidgets.Checkbox( description="Vary "+name, 
                    value = name in self.parameterobj.varylist,
                                                indent=False,
                    layout=ipywidgets.Layout(positioning='left'))
                fix[name].observe( lambda a, name=name: self.fixVaryCallBack(a, name),
                                  names='value' )
            else:
                fix[name] = ipywidgets.Label( 
                    value=name,
                    layout=ipywidgets.Layout(positioning='left') ) 
        self.valuewidgets = values
        self.fixwidgets = fix
        

    def fitGui(self):
        if self.gui is not None:
            return self.gui
        else:
            with pl.ioff():
                self.fig1 = pl.figure(figsize=(9,6), constrained_layout=True)
            self.ax1 = self.fig1.add_subplot()
            self.compute_tth_eta()
            self.addcellpeaks()
            self.pt1, = self.ax1.plot( self.colfile.tth, self.colfile.eta, ".", 
                                          ms=2, zorder=1)
            self.ax1.set(xlabel="tth", ylabel="eta")
            self.pt2, = self.ax1.plot( self.theorytth, np.full_like(self.theorytth,0), 
                           "g|", ms=360, lw= 0.1, zorder=2 )
                
            self.fitbutton = ipywidgets.Button(description="Fit x plot range",
                                               width='auto')
                                              
            self.fitbutton.on_click( self.fitCallBack )
            self.drawWidgets()
            self.gui = ipywidgets.GridBox( 
                [ ipywidgets.VBox( list(self.fixwidgets.values()) ), 
                  ipywidgets.VBox( list(self.valuewidgets.values()) ), 
                  ipywidgets.VBox( [ self.fitbutton, self.fig1.canvas ], ) ],
                layout=ipywidgets.Layout(  width='100%',
                    grid_template_rows='auto auto auto',
                    grid_template_columns='15% 15% 70%',) )
            display(self.gui)
            ping(self.fig1) # gets lost on run all cells otherwise
            
            


def ping(fig):
    # https://github.com/matplotlib/ipympl/issues/290
    canvas = fig.canvas
    #display(canvas)
    canvas._handle_message(canvas, {'type': 'send_image_mode'}, [])
    canvas._handle_message(canvas, {'type':'refresh'}, [])
    canvas._handle_message(canvas,{'type': 'initialized'},[])
    canvas._handle_message(canvas,{'type': 'draw'},[])