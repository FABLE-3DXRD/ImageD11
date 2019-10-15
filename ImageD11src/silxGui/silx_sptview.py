

import sys, os
import numpy as np, fabio
import silx.gui.qt , silx.gui.plot

from ImageD11.peakmerge import peakmerger
qapp=None

class sptview(object):

    def __init__(self, fname=None):
        self.pm = peakmerger()
        self.fname = fname
        self.pm.readpeaks( fname )
        self.select_image(0)
        self.drawUI()

    def select_image(self, i):
        if i < 0 or  i > len(self.pm.images)-1:
            return False
        im = self.pm.images[i]
        self.currentnum = i
        j = im.imagenumber
        self.pm.harvestpeaks(  numlim=(j-.1,j+0.1) )
        self.frame = fabio.open(im.name)       
        self.pkx = [ p.x for p in self.pm.allpeaks ]
        self.pky = [ p.y for p in self.pm.allpeaks ]
        self.pkI = [ p.avg for p in self.pm.allpeaks ]
        return True

    def next(self):
        if self.select_image( self.currentnum + 1):
            self.changeframe()
        else:
            self.warn("No next image")

    def prev(self):
        if self.select_image( self.currentnum - 1):
            self.changeframe()
        else:
            self.warn("No previous image")

    def drawUI(self):
        self.widget = silx.gui.qt.QWidget()
        self.widget.setWindowTitle("spot viewer")
        self.layout = silx.gui.qt.QGridLayout()
        self.widget.setLayout( self.layout )
        self.spot_label = silx.gui.qt.QLabel( self.fname )
        self.framelabel = silx.gui.qt.QLabel( self.frame.filename )
        self.framelabel.setAlignment( silx.gui.qt.Qt.AlignCenter )
        self.layout.addWidget(self.spot_label, 0, 0, 1, 3)
        self.layout.addWidget(self.framelabel, 1, 1)
        self.n = silx.gui.qt.QPushButton( "next" )
        self.p = silx.gui.qt.QPushButton( "prev" )
        self.n.clicked.connect( self.next )
        self.p.clicked.connect( self.prev )
        self.layout.addWidget(self.p, 1, 0 )
        self.layout.addWidget(self.n, 1, 2 )
        self.plot = silx.gui.plot.Plot2D(  backend='gl')
        self.imglabel = self.plot.addImage( self.frame.data, z=0, origin=(-0.5,-0.5),
            colormap=silx.gui.colors.Colormap(name='magma'))
        self.plotlabel = self.plot.addScatter( self.pky, self.pkx, self.pkI, z=1, symbol='+',
            colormap=silx.gui.colors.Colormap(name='viridis'))
        self.plot.getScatter( self.plotlabel ).setSymbolSize(10)
        self.plot.setKeepDataAspectRatio(True)
        self.layout.addWidget(self.plot, 2, 0, 1, 3)
        self.widget.show()
        self.plot.show()

    def changeframe(self):
        d = self.frame.data
        im = self.plot.getImage( self.imglabel )
        im.setData(d)
        pl = self.plot.getScatter( self.plotlabel )
        pl.setData( self.pky, self.pkx, self.pkI )
        self.framelabel.setText( self.frame.filename )

    def warn(self, message):
        silx.gui.qt.QMessageBox.warning(None, "Warning", message)

if __name__=="__main__":
    qapp = silx.gui.qt.QApplication([])
    s=sptview( sys.argv[1] )
    s.widget.show()
    sys.exit(qapp.exec_())
