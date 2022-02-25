#!/usr/bin/env python
import sys, time
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from ImageD11.columnfile import columnfile 
 
class PeaksDisplay(QMainWindow):
    def __init__(self,parent=None):
        super(PeaksDisplay,self).__init__(parent)
        self.setWindowTitle('Peaks file')
        view = QTableView()
        tableData = PeaksTableModel()
        view.setModel(tableData)
        self.setCentralWidget(view)
 
 
class PeaksTableModel(QAbstractTableModel):
    """Model class that drives the population of tabular display"""
    def __init__(self):
        super(PeaksTableModel,self).__init__()
        start = time.time()
        self.colfile = columnfile( sys.argv[1] )
        print( time.time()-start)
 
    def rowCount(self,index=QModelIndex()):
        return self.colfile.nrows
 
    def readfile(self,fname):
        self.beginResetModel()
#	self.colfile = columnfile(fname)
        self.endResetModel()
 
    def columnCount(self,index=QModelIndex()):
        return self.colfile.ncols
 
    def data(self,index,role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            return QVariant( float(self.colfile.bigarray[ index.column(), index.row()] ))
        else:
            return QVariant()
 
    def headerData(self,section,orientation,role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return QVariant() 
        if orientation == Qt.Horizontal:
            return QVariant(self.colfile.titles[section])
        return QVariant(int(section + 1))
 
def start():
    app  = QApplication(sys.argv)
    appWin = PeaksDisplay()
    appWin.show()
    app.exec_()
 
if __name__ == "__main__":
    start()
