

from __future__ import print_function, division

import os, sys
from timeit import default_timer as timer
import numpy as np
from scipy.spatial.transform import Rotation
from ImageD11 import columnfile, transform
import nbsplat

# This is not using Qt.py
from PyQt5 import QtWidgets, QtCore, QtGui

class dataset:
    def __init__(self,
                 colfilename,
                 parfilename=None,
                 ubifilename=None ):
        self.colf = columnfile.columnfile( colfilename )
        self.nrows = self.colf.nrows
        if parfilename is not None:
            self.loadpars( parfilename )
        self.grains = [ ]
        if ubifilename is not None:
            self.grains = self.loadgrains( ubifilename )
        self.tol = 0.05
        self.update()
        self.names ='gx','gy','gz'
        
    def loadpars(self, fname):
        self.colf.parameters.loadparameters( fname )
        
    def loadgrains(self, fname):
        self.grains = read_grain_file( fname )

    def update(self):
        self.colf.updateGeometry()
        if len(self.grains)>0:
            self.resetlabels()

    def xyz(self):
        return [self.colf.getcolumn(name) for name in self.names]
            
    def grain_drlv2(self, grainid ):
        x,y,z = self.xyz()
        hkl = np.dot( self.grains[grainid].ubi, self.xyz )
        dhkl = np.round(hkl)-hkl
        drlv2 = (dhkl*dhkl).sum(axis=0)
        return drlv2

    def resetlabels(self):
        self.labels    = np.empty(  nr    , 'i')
        self.labels.fill(-1)
        self.gvecs     = np.array( (self.colf.gx,
                                    self.colf.gy,
                                    self.colf.gz ) ).T.copy()
        self.drlv2     = np.ones (  nr    , float)
        self.peaks_xyz = np.array( (self.colf.XL,
                                    self.colf.YL,
                                    self.colf.ZL ) ).T.copy()
        self.colf.addcolumn( "labels", labels )
        self.colf.addcolumn( "drlv2", drlv2 )
    
    def assignlabels( self, tol=None ):
        if tol is not None:
            self.tol = tol
        self.resetlabels()
        nr = self.colf.nrows
        for i, g in enumerate( self.grains ):
            cImageD11.compute_gv( self.pks_xyz,
                                  self.colf.omega,
                                  self.colf.parameters.get('omegasign'),
                                  self.colf.parameters.get('wavelength'),
                                  self.colf.parameters.get('wedge'),
                                  self.colf.parameters.get('chi'),
                                  g.translation,
                                  self.gvecs )
            cImageD11.score_and_assign( g.ubi,
                                        self.gvecs,
                                        self.tol,
                                        self.drlv2,
                                        self.labels, i )
        self.colf.labels[:]= self.labels
        self.colf.drlv2[:] = self.drlv2

    def peaksize(self):
        return self.colf.Number_of_pixels[:]
    
class splat3dview:
    def __init__(self, w=2000, h=1500):
        self.resize(w,h)
        #self.w = w                 # width
        #self.h = h                 # height
        self.bg = (0,0,0,0)        # background
        self.fg = (255,255,255,0)  # foreground
        self.ps = 2                # pointsize
        self.origin = (0,0,0)      # origin in reciprocal space
        self.scale = 0.4*min(w,h)    # scale for Angstrom to px
        self.zsort = True
        self.u0 = np.eye(3, dtype=np.float )
        self.u  = np.eye(3, dtype=np.float )

    def score(self, datas):
        return datas.peaksize()

    def resize(self, w, h):
        self.w = w
        self.h = h
        self.rgba = np.zeros((h,w,4),np.uint8)

    def rotate(self, dx, dy):
        r = Rotation.from_euler("XYZ", (dx, dy, 0),
                                degrees=True )
        self.u = r.as_dcm()

    def resetrot(self):
        self.u0 = np.dot(self.u, self.u0)
        self.u = np.eye(3, dtype=np.float)
        
    def matrix(self):
        return self.scale*np.dot(self.u, self.u0)
    
    def drawPixMapOnWidget(self, datas, target):
        start = timer()
        nbsplat.nbsplat( self.rgba,
                         datas.xyz(),
                         self.matrix(),
                         self.ps,
                         self.zsort,
                         self.fg,
                         self.bg )
        # recompute or cache ? self.colors
        # order should be here too.
        # print("draw",timer()-start)
        # fixme - can we write on i.bits() directly ?
        i = QtGui.QImage( self.rgba, self.w, self.h,
                          QtGui.QImage.Format_RGB32 )
        p = QtGui.QPixmap.fromImage(i)
        # print("draw",timer()-start)#,self.rx,self.ry,self.rz)
        target.setPixmap( p )


        
class grainScoreView(splat3dview):
    def __init__(self, grain_id, atol=0.01):
        self.tol2 = atol*atol
        self.grain_id = grain_id
        super().__init__(self)

    def score(self ):
        """ log( drlv2 + a ) with a small """
        return -np.log( self.model.grain_drlv2( self.grainid ) + self.tol2 )

    def setGrain(self, grain_id):
        self.grain_id = grain_id

    def setTol(self, tol):
        self.tol2 = tol*tol
        
               
class grainIdView(splat3dview):
    def __init__(self, tol=0.05):
        self.tol = tol
        super().__init__(self)

    def score(self, grainid, atol=0.01 ):
        """ log( drlv2 + a ) with a small """
        return -np.log( self.model.grain_drlv2( grainid ) + self.tol*self.tol )
    
    def set_tol(self, tol):
        self.tol = tol
    

        
class Example( QtWidgets.QWidget ):
    def __init__(self, vu, da):
        QtWidgets.QWidget.__init__(self)
        self.vu = vu
        self.da = da
        self.initUI()
        self.x0 = 0
        self.y0 = 0
        
    def mousePressEvent(self, e):
        #print('down', e.x(),e.y(),e.button() )
        self.x0 = e.x()
        self.y0 = e.y()
        
    def mouseReleaseEvent(self, e):
        self.mouseMoveEvent(e)
        #print('up', e.x(),e.y())
        self.vu.resetrot( )
        
    def mouseMoveEvent(self,e):
        # print('move',e.x(),e.y())
        dx = (self.x0 - e.x())/10.  # 0.1 deg/px
        dy = (self.y0 - e.y())/10.
        self.vu.rotate( dy, dx )
        self.vu.drawPixMapOnWidget( self.da, self.lbl )
    
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
    
