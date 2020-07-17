from __future__ import print_function
from fabio import tifimage
from ImageD11 import columnfile
import numpy as np
import os, unittest, sys


class test_tifs(unittest.TestCase):

    FIRST=6
    LAST =16
    NPK  =0 

    CMD  = sys.executable + " " + os.path.join( 
            "..","..", "scripts", "peaksearch.py " )
    
    def setUp(self):
        blank = np.zeros((24,24),np.uint16)
        peak   = blank.copy()
        peak[ 12:14, 15:17 ] = 42
        
        # put a peak on every other image
        for i in range(self.FIRST, self.LAST+1):
            if i%2 == 0:
                obj = tifimage.tifimage( peak, { } )
                self.NPK += 1
                print( i, self.NPK)
            else:
                obj = tifimage.tifimage( blank, { } )
            obj.write( 'tiftest%04d.tif'%(i))
        
       
    def tearDown(self):
        return
        # we'll leave the files around for now

        for i in range(self.FIRST, self.LAST+1):
            os.remove( 'tiftest%04d.tif'%(i) )
        


    def testpeaksearch(self):
        os.system(self.CMD + " -n tiftest " + \
                  " -t 5 -F .tif -S 0.3 " + \
                  " -T 11 -p Y --OmegaOverRide " +\
                  " -f %d -l %d"%(self.FIRST, self.LAST))
        results = columnfile.columnfile( "peaks_t5.flt" )
        self.assertEqual( results.nrows, self.NPK)
        if os.path.exists( 'peaks_t6.flt' ):
            self.assertEqual( open( 'peaks_t5.flt').read().rstrip(),
                              open( 'peaks_t6.flt').read().rstrip() )

    def testpeaksearch_nooverride(self):
        os.system(self.CMD + " -n tiftest " + \
                  " -t 6 -F .tif -S 0.3 " + \
                  " -T 11 -p Y " +\
                  " -f %d -l %d"%(self.FIRST, self.LAST))
        results = columnfile.columnfile( "peaks_t6.flt" )
        self.assertEqual( results.nrows, self.NPK)
        if os.path.exists( 'peaks_t5.flt' ):
            self.assertEqual( open( 'peaks_t5.flt').read().rstrip(),
                              open( 'peaks_t6.flt').read().rstrip() )
        
    def testpeaksearch_singlethread(self):
        os.system(self.CMD + " -n tiftest " + \
                  " -t 7 -F .tif -S 0.3 " + \
                  " -T 11 -p Y --singleThread " +\
                  " -f %d -l %d"%(self.FIRST, self.LAST))
        results = columnfile.columnfile( "peaks_t7.flt" )
        self.assertEqual( results.nrows, self.NPK)
        if os.path.exists( 'peaks_t6.flt' ):
            self.assertEqual( open( 'peaks_t7.flt').read().rstrip(),
                              open( 'peaks_t6.flt').read().rstrip() )
        
        
if __name__ == '__main__':
    unittest.main()
    
   
    
    

