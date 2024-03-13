
import numpy as np
from ImageD11 import labelimage
from ImageD11 import columnfile


def testcase():
        d0 = np.array(
           [[0,0,0,0,0],
            [0,1,0,1,0],
            [0,0,0,0,0],
            [0,1,0,1,0],
            [0,0,0,0,0]],   # 4 peaks
           np.float32)
        d1 = np.array(
           [[0,0,0,0,0],
            [0,1,1,1,0],
            [0,0,0,0,0],
            [0,1,0,1,0],
            [0,0,0,0,0]],   # 3 peaks, merges top
           np.float32)
        d2 = np.array(
           [[0,0,0,0,0],
            [0,1,0,1,0],
            [0,1,0,0,0],
            [0,1,0,1,0],
            [0,0,0,0,0]],   # 3 peaks, merges left
           np.float32)
        d3 = np.array(
           [[0,0,0,0,0],
            [0,1,0,1,0],
            [0,0,0,1,0],
            [0,1,0,1,0],
            [0,0,0,0,0]],   # 3 peaks, merges right
           np.float32 )
        d4 = np.array(
           [[0,0,0,0,0],
            [0,1,0,1,0],
            [0,0,0,0,0],
            [0,1,1,1,0],
            [0,0,0,0,0]],   # 3 peaks, merges bottom
           np.float32)
        z = np.zeros_like(d0)
        frame_series = [ z, d0, d1, 
                         z, d0, d1, d0,
                         z, d0, d2, 
                         z, d0, d2, d0,
                         z, d0, d3, 
                         z, d0, d3, d0,
                         z, d0, d4, 
                         z, d0, d4, d0,
                         z, d0, d1, d2, d3, d4, z ]
        fltfile = "test3d.flt"
        with open("test3d.spt","w") as sptfile:
            lio = labelimage.labelimage(d0.shape , fltfile )
            for i, frame in enumerate(frame_series):
                sptfile.write("\n# Frame %d\n"%(i))
                lio.peaksearch(frame, 0.1, i)
                if lio.npk > 0:
                    lio.output2dpeaks( sptfile )
                lio.mergelast()
            lio.finalise()
    

if __name__=="__main__":
    testcase()
