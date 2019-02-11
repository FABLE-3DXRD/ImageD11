

import sys, os
import numpy as np, fabio
import silx.gui.qt , silx.gui.plot

from ImageD11.peakmerge import peakmerger


class sptview(object):

    def __init__(self, fname=None):
        self.pm = peakmerger()
        self.fname = fname
        self.pm.readpeaks( fname )
        self.select_image(0)
        self.drawUI()

    def select_image(self, i):
        im = self.pm.images[i]
        self.currentnum = i
        j = im.imagenumber
        self.pm.harvestpeaks(  numlim=(j-.1,j+0.1) )
        self.frame = fabio.open(im.name)       
        self.pkx = [ p.x for p in self.pm.allpeaks ]
        self.pky = [ p.y for p in self.pm.allpeaks ]
        self.pkI = [ p.avg for p in self.pm.allpeaks ]

    def next(self):
        self.select_image( self.currentnum + 1)
        self.changeframe()

    def prev(self):
        self.select_image( self.currentnum - 1)
        self.changeframe()

    def drawUI(self):
        w1 = silx.gui.qt.QWidget()
        l = silx.gui.qt.QLabel( self.fname )
        self.framelabel = silx.gui.qt.QLabel( self.frame.filename )
        n = silx.gui.qt.QPushButton( "next" )
        p = silx.gui.qt.QPushButton( "prev" )
        n.clicked.connect( self.next )
        p.clicked.connect( self.prev )
        ui = silx.gui.qt.QVBoxLayout(w1)
        ui.addWidget(l)
        ui.addWidget(self.framelabel)
        ui.addWidget(n)
        ui.addWidget(p)
        self.w = w1
        self.win = silx.gui.plot.Plot2D( parent=None, backend='gl' )
        self.imglabel = self.win.addImage( self.frame.data, z=0, origin=(-0.5,-0.5),
            colormap=silx.gui.colors.Colormap(name='magma'))
        self.plotlabel = self.win.addScatter( self.pky, self.pkx, self.pkI, z=1, symbol='o',
            colormap=silx.gui.colors.Colormap(name='viridis'))
        self.win.getScatter( self.plotlabel ).setSymbolSize(10)
        self.win.show()
        self.w.show()
        
    def changeframe(self):
        d = self.frame.data
        im = self.axes.getImage( self.imglabel )
        im.setData(d)
        pl = self.axes.getScatter( self.plotlabel )
        pl.setData( self.pkx, self.pky, self.pkI )
        self.framelabel.setText( self.frame.filename )


if __name__=="__main__":
    qapp = silx.gui.qt.QApplication([])
    sptview( sys.argv[1] )
    qapp.exec_()