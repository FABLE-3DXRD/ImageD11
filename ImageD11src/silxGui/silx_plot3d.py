


from ImageD11.columnfile import columnfile 
import sys
import numpy as np
import silx.gui.qt , silx.gui.plot3d.SceneWindow
from silx.gui.plot3d.tools.PositionInfoWidget import PositionInfoWidget
from silx.gui.widgets.BoxLayoutDockWidget import BoxLayoutDockWidget




if __name__ == "__main__":
    colf = columnfile( sys.argv[1] )
    colf.parameters.loadparameters( sys.argv[2] )
    colf.updateGeometry()

    qapp = silx.gui.qt.QApplication([])

    # Create a SceneWindow widget
    window = silx.gui.plot3d.SceneWindow.SceneWindow()

    sceneWidget = window.getSceneWidget()
    sceneWidget.setBackgroundColor((0.1, 0.18, 0.08, 1.))
    sceneWidget.setForegroundColor((1., 1., 1., 1.))
    sceneWidget.setTextColor((0.5, 0.7, 0.7, 1.))
    sceneWidget.setProjection('orthographic')

    positionInfo = PositionInfoWidget()
    positionInfo.setSceneWidget(sceneWidget)
    dock = BoxLayoutDockWidget()
    dock.setWindowTitle("Selection Info")
    dock.setWidget(positionInfo)
    window.addDockWidget(silx.gui.qt.Qt.BottomDockWidgetArea, dock)

    x,y,z,values = colf.gx, colf.gy, colf.gz, colf.avg_intensity

    scatter3d = silx.gui.plot3d.items.Scatter3D()
    scatter3d.setData(x, y, z, values)
    # Set scatter3d properties
    scatter3d.getColormap().setName('viridis')  # Use 'magma' colormap
    scatter3d.setSymbol('.') # Use point markers
    scatter3d.setSymbolSize(11)  # Set the size of the markers

    # Add scatter3d to the scene
    sceneWidget.addItem( scatter3d )

    # Set scatter3d transform
    SIZE=1
    scatter3d.setScale(SIZE, SIZE, SIZE)

    window.show()
    # Avoid segfault on exceptions
    sys.excepthook = silx.gui.qt.exceptionHandler
    qapp.exec_()