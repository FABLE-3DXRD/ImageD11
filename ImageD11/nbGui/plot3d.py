
import io
import numpy as np
import ipywidgets
import scipy.spatial.transform
from ImageD11 import cImageD11
from PIL import Image

def demodata():
    h,k,l = np.mgrid[-3:4,-3:4,-3:4]
    gve = np.dot( np.eye(3)*0.1 , (h.ravel(), k.ravel(), l.ravel() ))
    return gve

class Plot3d:
    
    def __init__(self, xyz=demodata(), w=386, h=256, rx=0., ry=0., rz=0., npx = 1 ):
        self.rgba = np.empty( (h, w, 4), 'B')
        self.xyz = xyz
        self.npx = npx
        self.ipyimg = ipywidgets.Image( )
        self.wrx = ipywidgets.FloatSlider( value=rx,  min=-360, max=360.0, step=1, description='rx:', disabled=False,
                    continuous_update=True,    orientation='vertical',    readout=True, readout_format='.1f' )
        self.wry = ipywidgets.FloatSlider( value=ry,  min=-360, max=360.0, step=1, description='ry:', disabled=False,
                    continuous_update=True,    orientation='vertical',    readout=True, readout_format='.1f' )
        self.wrz = ipywidgets.FloatSlider( value=rz,  min=-360, max=360.0, step=1, description='rz:', disabled=False,
                    continuous_update=True,    orientation='vertical',    readout=True, readout_format='.1f' )
        self.wrx.observe( self.redraw, names='value' )
        self.wry.observe( self.redraw, names='value' )
        self.wrz.observe( self.redraw, names='value' )
        self.redraw(None)
        self.widget = ipywidgets.HBox([ self.ipyimg, self.wrx, self.wry, self.wrz] )
    
    def redraw(self,change):
        u = scipy.spatial.transform.Rotation.from_euler('XYZ',
                                                        (self.wrx.value, self.wry.value, self.wrz.value), 
                                                        degrees=True).as_matrix()
        rotated = u.dot(self.xyz)
        order = np.argsort( (rotated[2]*100).astype(np.int16) )
        cImageD11.splat( self.rgba, rotated[:,order].T, u.ravel(), self.npx )
        img = Image.fromarray(self.rgba)
        with io.BytesIO() as buffer:
            img.save( buffer, format='gif' )
            self.ipyimg.value = buffer.getvalue()

