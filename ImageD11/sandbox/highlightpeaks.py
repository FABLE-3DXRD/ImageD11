
from __future__ import print_function



# Kindly contributed by Martijn Van Hulzen, Delft, 2016.


# script that increases signal to background ratio by 
# looking for the minimum and maximum value of a pixel in a set of edfs 
# (usually rotated along the z-axis over an omega range) and then subtracting 
# the min from max it also produces the average edf per cycle



import fabio, glob, numpy as np, sys, os

stem, f, l, n = sys.argv[1:] # define base filename (stem), first framenr (f), 
         # last framenr (l) and nr of frames in a cycle (n) as arguments
f=int(f) # turn first framenr into integer
l=int(l) # turn last framenr  into integer
n=int(n) # turn number of frames in a cycle into integer
m = 0    # correction added to the framenr in case a framenr is skipped

# loop over the range of edfs from first to last frame with step n
for i in range(f,l,n):
    fn = "%s%04d.edf"%(stem,i+m)     # set the filename
    # loop until a valid filename is found
    while not os.path.isfile(fn) and m < 10:  
        # print the filename that does not exist
        print("%s does not exist" % fn)        
        # increase correction by one because the frame does not exist 
        m = m + 1                            
        # set the new filename now that m = m + 1 
        fn = "%s%04d.edf"%(stem,i+m)         
    if m > 9:
        print ("Stopping, as too many files do not exist")
        break              # break when too many filenames are missing
    im = fabio.open(fn)            # read the frame data
    s = im.data.astype(np.float32) # assign floating point data to s
                                   # used to determine the average
    lo = s.copy()                  # copy for determining minimum
    hi = s.copy()                  # copy for determining maximum
    for j in range(n):                     # loop over a cycle of frames
        fn = "%s%04d.edf" % (stem,i+j+m)         # set frame filename
        while not os.path.isfile(fn) and m < 10: # check whether filename exists
            # file does not exist, increase m and try again
            m = m + 1                                
            # print filename that does not exist
            print("%s does not exist" % fn)      
            # set new filename      
            fn = "%s%04d.edf" % (stem,i+j+m)    
        if m > 9: 
            print ("Stopping, as too many files do not exist")
            break               # break when too many filenames are missing
        print("%s" % fn)         # print frame to be processed
        f = fabio.open(fn).data # retrieve frame data
        s = s + f               # add new frame to previously added frames
        # determine the min between the current minimum and the current frame
        lo = np.minimum( lo, f )
        # determine the max between the current maximum and the current frame
        hi = np.maximum( hi, f )                
    if m > 9: 
        break                   # break out if files are not found
    s = s / n   # determine the average by dividing by the number of files
    # assign back to original structure that was created at the beginning 
    # of the process (reuse)
    im.data = s  
    im.write( "avg%04d.edf"%(i/n)) # write the average to file
    print("avg%04d.edf"%(i/n))      # print the average filename
    im.data = hi-lo                # assign the max - min to the data structure
    im.write( "pks%04d.edf"%(i/n)) # write the peaks to file
    print("pks%04d.edf"%(i/n))      # print the peaks filename
