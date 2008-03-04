

"""
Generate a series of test images
"""
import Image
# Make an image
im = Image.new( "L", ( 64, 99) ) 
for i in [0, 1, 2, 8, 9]:  # 0 , 1 , 2
    im.save("test%04d.tif"%(i))
# Put a 2x3 peak at (23, 45) for images 3-7 in a 0-9 series
im.putpixel( (23,45) , 128)
im.putpixel( (23,46) , 128)
im.putpixel( (24,45) , 128)
im.putpixel( (24,46) , 128)
im.putpixel( (23,47) , 128)
im.putpixel( (24,47) , 128)
for i in range(3, 8):  # 3, 4, 5, 6, 7
    im.save("test%04d.tif"%(i))
    
import os
from ImageD11.columnfile import columnfile

def test(start, step):
    os.remove("peaks_t1.flt")
    os.system("c:\python25\python c:\python25\scripts\peaksearch.py " +
          "-n test -F .tif -f 0 -l 9 -t 1.0 -D 0 " + 
          "-T %f -S %f >> testpksearch.log"%(start,step))
    c = columnfile("peaks_t1.flt")
    print "start",start,"step", step,"finds:",
    print "Min_o",c.Min_o[0],"Max_o",c.Max_o[0],"omega",c.omega[0]

test( -4, +1) 
test( -40, +1) 
test( 40, +1) 


