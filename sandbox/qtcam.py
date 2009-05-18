
from PyQt4 import Qt
import numpy, sys, time, os

"""
Read an image from a Tango Device
Display it in an opencv window, tracking features

Assumes an RGB picture (=> dim_x is 3*image size)
"""

if not os.environ.has_key('TANGO_HOST'):
    os.environ['TANGO_HOST'] = "crunch:20000"
print "Using TANGO_HOST = ", os.environ['TANGO_HOST']

import PyTango





def normalise( im ):
    """ Set mean to 0, std to 1 """
    avg = numpy.ravel(im).mean()
    std = numpy.ravel(im).std()
    if std > 0:
        return ( im - avg ) / std
    else:
        return im - avg

def cross_correlate(a2, b2, fftshape):
    """ Make 2D cross correlation image """
    # FFT:
    a2 = numpy.fft.rfft2(a2, s = fftshape)
    b2 = numpy.fft.rfft2(b2, s = fftshape)
    # Cross correlate
    c = numpy.fft.irfft2( a2 * b2.conj() )
    return numpy.fft.fftshift( c )

def find_offset( c ):
    """ Co-ordinates of max in a 2D cross correlation """
    flatmax = c.argmax()
    # print flatmax, c.flat[flatmax]
    dim0 = int(flatmax/c.shape[1])
    dim1 = flatmax - c.shape[1] * dim0
    # print dim0,dim1,c[dim0,dim1]
    roi = c[ dim0-6:dim0+7, dim1-6:dim1+7 ]
    troi = roi.ravel().sum()
    x  = numpy.dot( roi.sum(axis=1), numpy.arange(dim0-6,dim0+7) ) / troi
    y  = numpy.dot( roi.sum(axis=0), numpy.arange(dim1-6,dim1+7) ) / troi
    # Average the local area ?
    return x, y


def register_image( im1, im2 ):
    """ We take the current and place it onto reference """

    ref = normalise(numpy.asarray(im1).sum(axis=2, dtype=numpy.float))
    cur = normalise(numpy.asarray(im2).sum(axis=2, dtype=numpy.float))

    fftshape = ( ref.shape[0] + cur.shape[0] ,
                 ref.shape[1] + cur.shape[1]  )
    cor = cross_correlate( ref, cur, fftshape)
    x, y = find_offset( cor )
    #print x, y, ref.shape, cur.shape
    return x, y







class tangocam:
    def __init__(self, CAMERA_DEVICE = "id11/ccd1394/1"):
        self.cam = PyTango.DeviceProxy(CAMERA_DEVICE)

    def getimage_numpy(self):
        """ return the image as a numeric array """
        res = self.cam.read_attribute_as_str("Image")
        ar = numpy.reshape(numpy.fromstring(res.value,
                                            numpy.uint8) ,
                                   ( res.dim_y , res.dim_x/3 , 3 ))
        return ar


class camreader( Qt.QThread ):
    def __init__(self, lock, cam, im, parent = None):
        super(camreader,self).__init__(parent)
        self.lock = lock
        self.cam = cam
        self.im = im
    def run(self):
        while 1:
            try:
                self.lock.lockForWrite()
                self.im[:,:] = self.cam.getimage_numpy()
            except:
                pass
            self.lock.unlock()
            time.sleep(0.5)


class ImageQt(Qt.QImage):

    def __init__(self, npyar):

        data = None
        colortable = None
        assert type(npyar) == numpy.ndarray
        assert len(npyar.shape) == 3
        assert npyar.shape[2] == 3
        assert npyar.dtype == numpy.uint8

        h , w = npyar.shape[0:2]
        im = numpy.zeros( (h, w, 4),
                          numpy.uint8)
        im[:,:,0] = npyar[:,:,2]
        im[:,:,1] = npyar[:,:,1]
        im[:,:,2] = npyar[:,:,0]  
        data = im.tostring()
        self.__data = data or im.tostring()
        Qt.QImage.__init__(self, self.__data, w, h, Qt.QImage.Format_RGB32 )
        


class qtcam( Qt.QApplication ):
    def __init__(self, args):

        self.cam = tangocam()
        self.numpyarray = self.cam.getimage_numpy()
        # reference is div/2 + 128 for later subtracting data/2
        self.refdata = self.numpyarray.copy()
        self.refimage = self.numpyarray/2 + 128
        self.last = self.numpyarray.copy()
        self.scalfac = 0.5
        self.iframe = 0

        Qt.QApplication.__init__(self, args )
    
        self.lock = Qt.QReadWriteLock()
        self.camreader = camreader( self.lock, self.cam, self.numpyarray )
        self.camreader.start()

        self.addwidgets()

        self.connect( self.loadbutton,
                Qt.SIGNAL("clicked()"),
                self.loadreference )

        self.connect( self.savebutton,
                Qt.SIGNAL("clicked()"),
                self.savereference )
        
        self.connect( self.zoombutton,
                Qt.SIGNAL("clicked()"),
                self.setscalfac )

        self.connect( self.calcbutton,
                Qt.SIGNAL("clicked()"),
                self.compute )


        self.timer = Qt.QTimer()
        self.timer.setInterval( 300 )
        self.connect( self.timer,
                Qt.SIGNAL("timeout()"),
                self.update )
        self.update()
        self.timer.start()

        self.exec_()


    def addwidgets(self):
        self.form = Qt.QDialog()
        self.form.setWindowTitle( "ID11 Camera" )
        self.vl = Qt.QGridLayout()

        self.vl.addWidget(  Qt.QLabel ("<h2>Live Feed</h2>" ),0,0 )
        self.zoombutton = Qt.QPushButton("Set zoom")
        self.vl.addWidget( self.zoombutton,0,1 )

        self.savebutton = Qt.QPushButton( "Save Current")
        self.vl.addWidget( self.savebutton, 0, 2)

        self.vl.addWidget(  Qt.QLabel ("<h2>Difference</h2>" ),0,3 )

        self.loadbutton = Qt.QPushButton( "Load Reference")
        self.vl.addWidget( self.loadbutton, 0,4)

        self.calcbutton = Qt.QPushButton("Compute")
        self.vl.addWidget( self.calcbutton,0,5 )

        self.camlabel = Qt.QLabel()
        self.vl.addWidget( self.camlabel,1,0,1,3 )

        self.difflabel = Qt.QLabel()
        self.vl.addWidget( self.difflabel,1,3,1,3 )

        self.form.setLayout(self.vl)
        self.form.show()


    def compute(self):
        self.lock.lockForRead()
        if self.numpyarray[0,0,0] != self.last[0,0,0]:
            self.last = self.numpyarray.copy()
        self.lock.unlock()
        x, y = register_image( self.refdata, self.last )
        print x, y
        Qt.QMessageBox.information (self.form, 
                "Computed Offset", 
                "Vertical  = %.5f  \nHorizontal = %.5f"%(
                    self.refdata.shape[0] - x ,
                    self.refdata.shape[1] - y), 
                Qt.QMessageBox.Ok)


    def setscalfac(self):
        val, ok = Qt.QInputDialog.getDouble(self.form, 
                "Scale factor", 
                "Please enter the image scale",
                self.scalfac, #value = 
                0.01, #minValue = 
                10 , #maxValue 
                 3)       #         int decimals 
        if ok:
            self.scalfac = val
        self.vl.invalidate()
        self.update(Forced=True)

    def loadreference(self):
        fn = Qt.QFileDialog.getOpenFileName(
                self.form,
                "Choose a file to open", 
                os.getcwd()
                )   
        if fn == "":
            return
        ar = numpy.fromstring( open(fn,"rb").read(),
                               numpy.uint8)
        self.refdata = ar.reshape( self.refimage.shape )
        self.refimage = self.refdata/2 + 128
        self.update(Forced=True)

    def savereference(self):
        fn = Qt.QFileDialog.getSaveFileName(
                self.form,
                "Give the filename to save", 
                os.getcwd()
                )   
        if fn == "":
            return
        open( fn, "wb" ).write( self.last.tostring() )
        


    def array_to_display(self, ar):
        i = ImageQt( ar )
        h, w, junk = ar.shape
        j =  Qt.QImage.scaled( i,
                               int(w*self.scalfac),
                               int(h*self.scalfac) )
        pm = Qt.QPixmap.fromImage(  j )
        return pm

    

    def update(self, Forced=False):
        import time
        start = time.clock()
        try:
            self.lock.lockForRead()
            if self.numpyarray[0,0,0] != self.last[0,0,0]:
                self.last = self.numpyarray.copy()
                if not Forced:
                    self.lock.unlock()
                    return
        except:
            pass
        self.lock.unlock()
        diff = numpy.subtract( self.refimage, self.last/2 ) 
        self.camlabel.setPixmap( self.array_to_display(self.last) )
        self.difflabel.setPixmap(self.array_to_display(diff) )
        self.camlabel.show()
        self.difflabel.show()
        self.iframe += 1
        print "Nframes %5d  last %.5f/s\r"%(self.iframe, time.clock()-start),


if __name__=="__main__":
    qtcam(sys.argv)
