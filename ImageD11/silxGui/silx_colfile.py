
from ImageD11.columnfile import columnfile 

import sys, os
import numpy as np
import silx.gui.qt , silx.gui.plot

class silxqtcolfile( object ):

    def __init__(self, filename=None, xlabel=None, ylabel=None, zlabel=None ):
        """
        Reads an ImageD11 peak search output for applying masks on 2D scatter plots
        """
        self.xlabel = None
        self.ylabel = None
        self.zlabel = None
        self.drawUI()
        self.loadColfile(filename)
        
    def drawUI( self ):
        """
        Sets up the UI
        | cfname        |  load | save
        | parname       |  load |
        | xaxis | yaxis | color | apply_mask |
        """
        self.win = silx.gui.qt.QWidget()
        self.scatter_widget = silx.gui.plot.ScatterView( backend='gl' )
        # import pdb; pdb.set_trace()
        if not self.scatter_widget._plot()._backend.isValid():
            self.scatter_widget = silx.gui.plot.ScatterView( backend='matplotlib' )
        self.scatter_widget.setColormap( silx.gui.colors.Colormap(name='viridis') )
        self.scatter_widget.getMaskToolsWidget().parent().setFloating(True)
        self.scatter_widget.getMaskToolsWidget().show()
        # glayout.addWidget( self.scatter_widget, 3, 0, 1, 4 )
        
        glayout = silx.gui.qt.QGridLayout(self.win)
        #   0     1    2     3
        #  fname       load  save             0
        #  pname       load                   1
        #  x      y    z     domask           2
        #  plots

        self.fnamelabel = silx.gui.qt.QLabel( 'Columnfile: ' )
        blc=silx.gui.qt.QPushButton( "Load colfile" )
        glayout.addWidget( self.fnamelabel, 0, 0, 1, 2 )
        blc.clicked.connect( self.loadColfile )
        glayout.addWidget(blc, 0, 2)
        bsc=silx.gui.qt.QPushButton( "Save colfile" )
        bsc.clicked.connect( self.savecolfile )
        glayout.addWidget(bsc, 0, 3)

        self.pnamelabel = silx.gui.qt.QLabel( 'Parameter file: ' )
        blp=silx.gui.qt.QPushButton( "Load parfile" )
        glayout.addWidget( blp, 1, 3)
        blp.clicked.connect(self.loadparameters)
        glayout.addWidget( self.pnamelabel, 1, 0, 1, 2 )

        bxa=silx.gui.qt.QComboBox( None )
        bya=silx.gui.qt.QComboBox( None )
        bza=silx.gui.qt.QComboBox( None )
        self.axisboxes = [bxa,bya,bza]
        self.ignoreselect= True
        for b in self.axisboxes:
            b.currentIndexChanged.connect( self.select )
        glayout.addWidget(silx.gui.qt.QLabel("x-plot"),2,0)
        glayout.addWidget(silx.gui.qt.QLabel("y-plot"),2,1)
        glayout.addWidget(silx.gui.qt.QLabel("color"),2,2)
        glayout.addWidget(bxa,3,0)
        glayout.addWidget(bya,3,1)
        glayout.addWidget(bza,3,2)

        bdel=silx.gui.qt.QPushButton( "Delete selected" )
        bdel.clicked.connect( self.applymask )
        glayout.addWidget(bdel,2,3)

        self.scatter_widget.show()
        self.win.show()

    def setTitles(self):
        """ Read the colfile titles into the dropdowns """
        self.ignoreselect=True
        for i,b in enumerate(self.axisboxes):
            t = b.currentText()
            b.clear()
            b.addItems( self.colfile.titles )
            if t in self.colfile.titles:
                b.setCurrentIndex( self.colfile.titles.index( t ) )
            else:
                b.setCurrentIndex(i)
        self.ignoreselect=False

    def select(self,col):
        """
        Choose the x,y,z axes for plotting
        """
        if self.ignoreselect:
            return
        self.update()
        self.scatter_widget.resetZoom()

    def savecolfile(self):
        """ Write after editing """
        filename = silx.gui.qt.QFileDialog(None,"Columnfile").getSaveFileName()
        print(filename)
        try:
            self.colfile.writefile( filename[0] )
        except:
            m = silx.gui.qt.QMessageBox.about(self.scatter_widget, 
            "Fail", "Saving failed - bad filename?"  )

    def loadColfile(self, filename=None):
        """ Read in a new file """
        if filename is None or filename is False:
            filename = silx.gui.qt.QFileDialog(self.win,"Columnfile").getOpenFileName()[0]
        try:
            c = columnfile( filename )
            self.colfile = c
        except:
            print("problem opening file",filename)
            return
        print("loaded file",filename)
        self.fnamelabel.setText( "Columnfile: " + self.colfile.filename )
        self.setTitles()
        self.update()
        self.scatter_widget.resetZoom()
        self.scatter_widget.activateWindow()
        self.win.activateWindow()

    def update(self):
        """ Refreshes the plot """
        self.x = self.colfile.getcolumn(self.axisboxes[0].currentText())
        self.y = self.colfile.getcolumn(self.axisboxes[1].currentText())
        self.z = self.colfile.getcolumn(self.axisboxes[2].currentText())
        self.scatter_widget.getXAxis().setLabel(self.axisboxes[0].currentText())
        self.scatter_widget.getYAxis().setLabel(self.axisboxes[1].currentText())
        self.scatter_widget.setGraphTitle('color: ' + self.axisboxes[2].currentText())
        self.scatter_widget.setData( self.x, self.y, self.z )

    def applymask(self):
        """ Applies the scatterplot mask to the colfile """
        m = self.scatter_widget.getMaskToolsWidget().getSelectionMask()
        print("Masking",(m!=0).sum())
        self.colfile.filter( m == 0 )
        self.update()

    def loadparameters(self):
        fname =  silx.gui.qt.QFileDialog(self.win,"Parfile").getOpenFileName()[0]
        if not os.path.exists(fname):
            m = silx.gui.qt.QMessageBox.about(self.scatter_widget, 
            "Fail", "Setting parameters failed - bad filename?"  )
        self.colfile.parameters.loadparameters( fname )
        self.colfile.updateGeometry()
        self.setTitles()
        self.pnamelabel.setText("Parameters: "+fname)
        self.update()

if __name__=="__main__":
    app = silx.gui.qt.QApplication( sys.argv )
    if len(sys.argv) > 1:
        qp = silxqtcolfile( sys.argv[1] )
    else:
        qp = silxqtcolfile( )
    # Avoid segfault on exceptions
    sys.excepthook = silx.gui.qt.exceptionHandler
    app.exec_()
