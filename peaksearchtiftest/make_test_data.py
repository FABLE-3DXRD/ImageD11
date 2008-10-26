

from fabio import tifimage
import numpy as np
import os

FIRST=6
LAST =16

def setUp():
    blank = np.zeros((24,24),np.uint16)
    peak   = blank.copy()
    peak[ 12:14, 15:17 ] = 42

    # put a peak on every other image
    for i in range(FIRST, LAST+1):
        if i%2 == 0:
            obj = tifimage.tifimage( peak, { } )
        else:
            obj = tifimage.tifimage( blank, { } )
        obj.write( 'tiftest%04d.tif'%(i))


def tearDown():

    for i in range(FIRST, LAST+1):
        os.remove( 'tiftest%04d.tif'%(i) ) 


def testpeaksearch():
    if os.path.exists('peaks_t5.spt'):
        os.remove('peaks_t5.spt')
    if os.path.exists('peaks_t5.flt'):
        os.remove('peaks_t5.flt')
    os.system('peaksearch.py -n tiftest -f %d -l %d -t 5 -F .tif -S 0.3 -T 11 -p Y --OmegaOverRide'%(FIRST, LAST))

if __name__ == '__main__':
    setUp()
    testpeaksearch()
    
   
    
    

