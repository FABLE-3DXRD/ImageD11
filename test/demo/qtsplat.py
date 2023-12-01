

from __future__ import print_function, division

import os, sys
from timeit import default_timer as timer
import numpy as np
from scipy.spatial.transform import Rotation
if not hasattr( Rotation, 'as_matrix' ):
    Rotation.as_matrix = Rotation.as_dcm
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
                t = np.array((0,0,0), float)
            else:
                t = np.array((g.translation[0],g.translation[1],g.translation[2]),
                             float)
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
        self.w = self.h = None
        self.resize(w,h)
        self.bg = (0,0,0,0)        # background
        self.fg = (255,255,255,0)  # foreground
        self.ps = 1                # pointsize
        self.origin = (0,0,0)      # origin in reciprocal space
        self.scale = 0.7*min(w,h)    # scale for Angstrom to px
        self.zsort = True
        self.u0 = np.eye(3, dtype=float )
        self.u  = np.eye(3, dtype=float )
        self.zname = "gz"
        self.need_redraw = True

    def resize(self, w, h):
        if w != self.w or h != self.h:
            self.w = w
            self.h = h
            self.rgba = np.zeros((h,w,4),np.uint8)

    def rotate(self, dx, dy, axes="XYZ"):
        self.need_redraw = True
        r = Rotation.from_euler(axes, (dx, dy, 0),
                                degrees=True )
        if hasattr(r, "as_dcm"):
            self.u = r.as_dcm()
        else:
            self.u = r.as_matrix()

    def resetrot(self):
        self.need_redraw = True
        self.u0 = np.dot(self.u, self.u0)
        self.u = np.eye(3, dtype=float)

    def matrix(self):
        return self.scale*np.dot(self.u, self.u0)



    def drawPixMapOnWidget(self, datas, target ):
        if self.need_redraw:
            start = timer()
            w, h = target.width(), target.height()
            self.resize( w, h )
            #print("wh %d %d %d %d"%(
            #    self.w,target.width(),self.h,target.height()))
            DOPROF=False
            if DOPROF:
                import cProfile, pstats
                x = cProfile.Profile()
                x.enable()
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
            p = rgbtopm( self.rgba )
            target.setPixmap( p )
            if DOPROF:
                x.disable()
                p = pstats.Stats( x, stream = sys.stdout).sort_stats("tottime")
                p.reverse_order()
                p.print_stats()
            print("drawPixMapOnWidget %.3f ms %s"%(
                1e3*(timer()-start), self.need_redraw))#, end="\r")
            sys.stdout.flush()
            self.need_redraw = False
            return(timer()-start)
        return 0

def rgbtopm( rgb ):
    i = QtGui.QImage(rgb, rgb.shape[1], rgb.shape[0],
                     QtGui.QImage.Format_RGB32 )
    return QtGui.QPixmap.fromImage(i)


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
        self.animate()
        self.show()

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
            axes="XZY"
        self.rotateview( e.x(), e.y(), axes  )

    def rotateview(self, ex, ey, axes="XYZ"):
        self.ex = ex
        self.ey = ey
        dx = (self.x0 - self.ex)/10.  # 0.1 deg/px
        dy = (self.y0 - self.ey)/10.
        self.vu.rotate( dy, -dx, axes )

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
        QtCore.Qt.Key_V : "gz",
        QtCore.Qt.Key_O : "omega",
        QtCore.Qt.Key_L : "labels",
        }


    def hlp(self):
        h = """ X,Y,Z = Rotate by 90 in X, Y, Z
Shift/Control = modify mouse rotate
+ = Zoom in
- = Zoom out
[0-9] = Pointsize
A : "avg_intensity",
D : "drlv2",
E : "eta",
T : "tth",
V : "gz",
O : "omega",
L : "labels",
S = Toggle zsort
Q = Quit
"""
        print(h)
        qmb= QtWidgets.QMessageBox.information(
            self, "Help", h )
        #qmb.setText(h)
        #qmb.show()


    def keyPressEvent(self,e):
        k = e.key()
        print("Got key",k)
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
            self.status.setText( "zsort %s"%(str(self.vu.zsort)))
        if k == QtCore.Qt.Key_Plus:
            self.vu.scale = self.vu.scale * ( 1 + 1/8)
            self.status.setText( "Scale %f"%(self.vu.scale))
        if k == QtCore.Qt.Key_Minus:
            self.vu.scale = self.vu.scale * ( 1 - 1/8)
            self.status.setText( "Scale %f"%(self.vu.scale))
        if k >= QtCore.Qt.Key_0 and k <= QtCore.Qt.Key_9:
            self.vu.ps = int(k - QtCore.Qt.Key_0)
            self.status.setText( "Pointsize %d"%(self.vu.ps))
        if k in self.rotateactions:
            d1,d2,a = self.rotateactions[k]
            self.status.setText("Rotated by %d %d %s"%(d1,d2,a))
            self.vu.rotate( d1, d2, a)
            self.vu.resetrot()
        if k in self.coloractions:
            if self.coloractions[k] in self.da.colf.titles:
                self.status.setText("color by %s"%(self.coloractions[k]))
                colormap( self.da.get( self.coloractions[k]),
                          self.da.colors )
            else:
                self.status.setText("No column %s"%(self.coloractions[k]))
        self.vu.need_redraw=True

    def keyReleaseEvent(self,e):
        if self.doRot and e.key() == QtCore.Qt.Key_Shift:
            self.x0 = self.ex
            self.y0 = self.ey
            self.vu.resetrot()


    def initUI(self):

        self.layout = QtWidgets.QVBoxLayout()

        self.lbl = QtWidgets.QLabel(self)
        self.lbl.resize( self.vu.w, self.vu.h )
        self.lbl.setMinimumSize( 128, 128 )
        e = QtWidgets.QSizePolicy.Expanding
        self.lbl.setSizePolicy( e, e )
        self.layout.addWidget( self.lbl )

        self.status = QtWidgets.QLabel(self)
        self.status.setText("Hello")
        self.layout.addWidget( self.status )

        self.buttons = QtWidgets.QHBoxLayout()
        hb= QtWidgets.QPushButton("Help", self)
        hb.clicked.connect( self.hlp )
        self.buttons.addWidget(hb)
        qb= QtWidgets.QPushButton("Quit", self)
        qb.clicked.connect( self.close )
        self.buttons.addWidget(qb)

        self.layout.addLayout( self.buttons )
        self.setLayout( self.layout )


    def animate(self):
        self.vu.drawPixMapOnWidget( self.da, self.lbl)
        QtCore.QTimer.singleShot(1000//24, self.animate )


    def resizeEvent(self, evt):
        self.vu.need_redraw=True




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

