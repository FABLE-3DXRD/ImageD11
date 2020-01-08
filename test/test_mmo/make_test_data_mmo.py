
from __future__ import print_function, division
################################################# make_test_data.py
#
#  save as a .py file and run it from a command line via:
#  python make_test_data.py
"""
Generate a series of test images
"""
import os
from ImageD11.columnfile import columnfile
import unittest
from PIL import Image
import platform, sys

def runtest(start, step):
    """ check tif range """
    pksh = sys.executable + ' ' + os.path.join("..","..","scripts","peaksearch.py")
    os.system(pksh + 
          " -n test -F .tif -f 0 -l 9 -t 1.0 -D 0 " + 
          "-T %f -S %f > testpksearch.log"%(start,step))
    c = columnfile("peaks_t1.flt")
    print ("start",start,"step", step,"finds:", end=" " )
    print ("Min_o",c.Min_o,"Max_o",c.Max_o,"omega",c.omega)
    if len(c.Min_o) == 3:
        print("I think it might be OK")
    return len(c.Min_o)==3


class test_mmo(unittest.TestCase):
    def setUp(self):
        # make a test image
        im = Image.new( "L", ( 264, 299) ) 
        for i in [0, 1,            8, 9]: 
            im.save("test%04d.tif"%(i))
        # Put a 2x3 peak at (23, 45) for images 3-7 in a 0-9 series
        # peak 1 - simple blob
        for p in [  (23,45),   (23,46), (23,47), 
                    (24,45),   (24,46), (24,47)]:
            im.putpixel( p , 128)
        # peak 2 - two will merge to one
        for p in [  (123,45),          (123,47), 
                    (124,45),          (124,47)]:
            im.putpixel( p , 129)
        # peak 3 - one becomes two
        for p in [  (23,145), (23,146), (23,147), 
                    (24,145), (24,146), (24,147)]:
            im.putpixel( p , 127)
        for i in [2, 3, 4]: 
            im.save("test%04d.tif"%(i))
        for p in  [  (123,46),  
                     (124,46)   ]:
            im.putpixel( p , 129)
        # peak 3 - one becomes two
        for p in [   (23,146),  
                     (24,146)]:
            im.putpixel( p , 0)
        for i in [5, 6, 7]: 
            im.save("test%04d.tif"%(i))
            
    def test1(self):
        runtest( -4, +1) 
        runtest( -40, +1) 
        runtest( 40, +1) 

if __name__=="__main__":
    unittest.main()
    
