
from fabio import tifimage
from ImageD11 import columnfile
import numpy as np
import os, unittest


class test_tifs(unittest.TestCase):

    FIRST=6
    LAST =16
    NPK  =0 

    def setUp(self):
        blank = np.zeros((24,24),np.uint16)
        peak   = blank.copy()
        peak[ 12:14, 15:17 ] = 42
        
        # put a peak on every other image
        for i in range(self.FIRST, self.LAST+1):
            if i%2 == 0:
                obj = tifimage.tifimage( peak, { } )
                self.NPK += 1
                print i, self.NPK
            else:
                obj = tifimage.tifimage( blank, { } )
            obj.write( 'tiftest%04d.tif'%(i))
        
        if os.path.exists('peaks_t5.spt'):
            os.remove('peaks_t5.spt')
        if os.path.exists('peaks_t5.flt'):
            os.remove('peaks_t5.flt')

    def tearDown(self):
        return
        # we'll leave the files around for now

        for i in range(self.FIRST, self.LAST+1):
            os.remove( 'tiftest%04d.tif'%(i) )
        


    def testpeaksearch(self):
        os.system("peaksearch.py -n tiftest " + \
                  " -t 5 -F .tif -S 0.3 " + \
                  " -T 11 -p Y --OmegaOverRide " +\
                  " -f %d -l %d"%(self.FIRST, self.LAST))
        results = columnfile.columnfile( "peaks_t5.flt" )
        self.assertEqual( results.nrows, self.NPK)


        
        
if __name__ == '__main__':
    unittest.main()
    
   
    
    

