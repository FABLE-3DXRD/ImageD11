

from __future__ import print_function, division

import os, sys
from timeit import default_timer as timer
import numpy as np
from scipy.spatial.transform import Rotation
from ImageD11 import columnfile, transform, grain, cImageD11
import nbsplat

# This is not using Qt.py
from PyQt5 import QtWidgets, QtCore, QtGui

class dataset:
    def __init__(self,
                 colfilename,
                 parfilename=None,
                 ubifilename=None ):
        self.loadcolf( colfilename )
        if parfilename is not None:
            self.loadpars( parfilename )
        self.grains = [ ]
        self.tol = 1.0
        if ubifilename is not None:
            self.loadgrains( ubifilename )
        self.names ='gx','gy','gz'
        self.colors=np.empty( (self.nrows, 4), np.uint8 )
        self.colors.fill(255)

    def loadcolf(self, colfilename):
        self.colf = columnfile.columnfile( colfilename )
        self.nrows = self.colf.nrows
                
    def loadpars(self, fname):
        self.colf.parameters.loadparameters( fname )
        self.colf.updateGeometry()
        
    def loadgrains(self, fname):
        self.grains = grain.read_grain_file( fname )
        if len(self.grains)>0:
            self.resetlabels()
            self.assignlabels()
        
    def get(self, name):
        if name in self.colf.titles:
            return self.colf.getcolumn(name)
        else:
            raise Exception("name "+name+" not found")
            
    def xyz(self):
        return [self.colf.getcolumn(name) for name in self.names]
            
    def grain_drlv2(self, grainid):
        x,y,z = self.xyz()
        hkl = np.dot( self.grains[grainid].ubi, self.xyz )
        dhkl = np.round(hkl)-hkl
        drlv2 = (dhkl*dhkl).sum(axis=0)
        self.colf.addcolumn("drlv2_%d"%(grainid))
    
    def resetlabels(self):
        self.labels    = np.empty(  self.nrows, 'i')
        self.labels.fill(-1)
        self.gvecs     = np.array( (self.colf.gx,
                                    self.colf.gy,
                                    self.colf.gz ) ).T.copy()
        self.drlv2     = np.ones (  self.nrows, float)
        #print(self.colf.titles)
        self.peaks_xyz = np.array( (self.colf.xl,
                                    self.colf.yl,
                                    self.colf.zl ) ).T.copy()
        self.colf.addcolumn( self.labels, "labels" )
        self.colf.addcolumn( self.drlv2, "drlv2" )
    
    def assignlabels( self ):
        self.resetlabels()
        nr = self.nrows
        pks_xyz = np.empty((nr,3), float)
        pks_xyz[:,0] = self.colf.xl
        pks_xyz[:,1] = self.colf.yl
        pks_xyz[:,2] = self.colf.zl
        #print(self.gvecs.shape, pks_xyz.shape)
        for i, g in enumerate( self.grains ):
            if g.translation is None:
                t = np.array((0,0,0), np.float)
            else:
                t = np.array((g.translation[0],g.translation[1],g.translation[2]),
                             np.float)
            cImageD11.compute_gv( pks_xyz,
                                  self.colf.omega,
                                  self.colf.parameters.get('omegasign'),
                                  self.colf.parameters.get('wavelength'),
                                  self.colf.parameters.get('wedge'),
                                  self.colf.parameters.get('chi'),
                                  t,
                                  self.gvecs )
            cImageD11.score_and_assign( g.ubi,
                                        self.gvecs,
                                        self.tol,
                                        self.drlv2,
                                        self.labels, i )
        self.colf.labels[:]= self.labels
        self.colf.drlv2[:] = self.drlv2

    
class splat3dview:
    def __init__(self, w=2000, h=1500):
        self.resize(w,h)
        #self.w = w                 # width
        #self.h = h                 # height
        self.bg = (0,0,0,0)        # background
        self.fg = (255,255,255,0)  # foreground
        self.ps = 12                # pointsize
        self.origin = (0,0,0)      # origin in reciprocal space
        self.scale = 0.7*min(w,h)    # scale for Angstrom to px
        self.zsort = True
        self.u0 = np.eye(3, dtype=np.float )
        self.u  = np.eye(3, dtype=np.float )
        self.zname = "gz"

    def resize(self, w, h):
        self.w = w
        self.h = h
        self.rgba = np.zeros((h,w,4),np.uint8)

    def rotate(self, dx, dy, axes="XYZ"):
        r = Rotation.from_euler(axes, (dx, dy, 0),
                                degrees=True )
        self.u = r.as_dcm()

    def resetrot(self):
        self.u0 = np.dot(self.u, self.u0)
        self.u = np.eye(3, dtype=np.float)
        
    def matrix(self):
        return self.scale*np.dot(self.u, self.u0)

    def drawPixMapOnWidget(self, datas, target ):
        start = timer()
        #print("self.zsort",self.zsort)
        nbsplat.nbsplat( self.rgba,
                         datas.xyz(),
                         self.matrix(),
                         self.ps,
                         self.zsort,
                         self.fg,
                         self.bg,
                         datas.colors )
        # recompute or cache ? self.colors
        # order should be here too.
        # print("draw",timer()-start)
        # fixme - can we write on i.bits() directly ?
        i = QtGui.QImage( self.rgba, self.w, self.h,
                          QtGui.QImage.Format_RGB32 )
        p = QtGui.QPixmap.fromImage(i)
        # print("draw",timer()-start)#,self.rx,self.ry,self.rz)
        target.setPixmap( p )




def colormap(x, colors):
    # fixme: load save or matplotlib or select
    xmin = x.min()
    xmax = x.max()
    xcen = 0.5*(xmin+xmax)
    print("Colormap",xmin,xmax,xcen)
    # colors = np.empty( (len(x), 4), np.uint8 )
    colors[:,0] = np.interp( x, [ xmin, xcen, xmax ], [ 64,  64, 255])    # R
    colors[:,1] = np.interp( x, [ xmin, xcen, xmax ], [ 64, 255, 64]) # G
    colors[:,2] = np.interp( x, [ xmin, xcen, xmax ], [255,  64,  64]) # B
    colors[:,3] = 0 # alpha
    

    
    

        
class Example( QtWidgets.QWidget ):
    def __init__(self, vu, da):
        QtWidgets.QWidget.__init__(self)
        self.vu = vu
        self.da = da
        self.initUI()
        self.x0 = 0
        self.y0 = 0
        self.ex = 0
        self.ey = 0
        self.doRot = False
        
    def mousePressEvent(self, e):
        #print('down', e.x(),e.y(),e.button() )
        self.doRot = True
        self.resetrot( e.x(), e.y())
        
    def resetrot(self, x0, y0 ):
        self.x0 = x0
        self.y0 = y0
        self.vu.resetrot()
        
    def mouseReleaseEvent(self, e):
        self.mouseMoveEvent(e)
        #print('up', e.x(),e.y())
        self.vu.resetrot( )
        self.doRot = False
        
    def mouseMoveEvent(self,e):
        # print('move',e.x(),e.y())
        axes = "XYZ"
        if e.modifiers() == QtCore.Qt.ShiftModifier:
            axes="ZYX" 
        if e.modifiers() == QtCore.Qt.ControlModifier:
            axes="YZX" 
        self.rotateview( e.x(), e.y(), axes  )

    def rotateview(self, ex, ey, axes="XYZ"):
        self.ex = ex
        self.ey = ey
        dx = (self.x0 - self.ex)/10.  # 0.1 deg/px
        dy = (self.y0 - self.ey)/10.
        self.vu.rotate(-dy, dx, axes )
        self.vu.drawPixMapOnWidget( self.da, self.lbl )

    rotateactions = {
        QtCore.Qt.Key_X : ( 90, 0, "XYZ"),
        QtCore.Qt.Key_Y : ( 90, 0, "YZX"),
        QtCore.Qt.Key_Z : ( 90, 0, "ZXY"),
    }

    coloractions = {
        QtCore.Qt.Key_A : "avg_intensity",
        QtCore.Qt.Key_D : "drlv2",
        QtCore.Qt.Key_E : "eta",
        QtCore.Qt.Key_T : "tth",
        QtCore.Qt.Key_G : "gz",
        QtCore.Qt.Key_O : "omega",
        QtCore.Qt.Key_L : "labels",
        }


    def hlp(self):
        h = """ X,Y,Z = Rotate by 90 in X, Y, Z
Shift/Control = modify mouse rotate
A = color by avg_intensity
D = color by drlv2
L = color by grain label
E = color by eta
O = color by omega
T = color by tth
G = color by gz
S = Sort on z or not
Up Arrow = Zoom in
Down Arrow = Zoom out
Right Arrow = Bigger Points
Left Arrow = Smaller Points
Q = Quit
"""
        print(h)
        qmb= QtWidgets.QMessageBox.information(
            self, "Help", h )
        #qmb.setText(h)
        #qmb.show()

        
    def keyPressEvent(self,e):
        k = e.key()
        if k == QtCore.Qt.Key_H:
            self.hlp()
        if self.doRot and k == QtCore.Qt.Key_Shift:
            self.x0 = self.ex
            self.y0 = self.ey
            self.vu.resetrot()
        if k == QtCore.Qt.Key_Q:
            self.close()
        if k == QtCore.Qt.Key_S:
            self.vu.zsort = not self.vu.zsort
        if k == QtCore.Qt.Key_Up:
            self.vu.scale = self.vu.scale * ( 1 + 1/8)
        if k == QtCore.Qt.Key_Down:
            self.vu.scale = self.vu.scale * ( 1 - 1/8)
        if k == QtCore.Qt.Key_Left:
            self.vu.ps = max(self.vu.ps-1, 1)
        if k == QtCore.Qt.Key_Right:
            self.vu.ps = min(self.vu.ps+1, 100)
        if k in self.rotateactions:
            d1,d2,a = self.rotateactions[k]
            print("Rotate by",d1,d2,a)
            self.vu.rotate( d1, d2, a)
            self.vu.resetrot()
        if k in self.coloractions:
            if self.coloractions[k] in self.da.colf.titles:
                print("color by",self.coloractions[k])
                colormap( self.da.get( self.coloractions[k]),
                          self.da.colors )
        self.vu.drawPixMapOnWidget( self.da, self.lbl )

    def keyReleaseEvent(self,e):
        if self.doRot and e.key() == QtCore.Qt.Key_Shift:
            self.x0 = self.ex
            self.y0 = self.ey
            self.vu.resetrot()
    
    def initUI(self):
        self.lbl = QtWidgets.QLabel(self)
        self.lbl.resize( self.vu.w, self.vu.h )
        self.vu.drawPixMapOnWidget( self.da, self.lbl)
        self.show()
        
        
        
            
        
        


def run( colf, pars, grains):
    # colf       pars        grains
    da = dataset( colf, pars, grains )
    vu = splat3dview( )
    app = QtWidgets.QApplication( sys.argv )
    l = Example( vu, da )
    sys.exit( app.exec_() )



if __name__=="__main__":
    colf   = sys.argv[1]
    pars   = None
    grains = None
    try:
        pars = sys.argv[2]
        grains = sys.argv[3]
    except:
        pass

    run( colf, pars, grains )
    
