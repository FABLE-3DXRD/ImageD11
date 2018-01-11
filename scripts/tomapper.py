#!/usr/bin/env python


"""
Script to convert and ImageD11 grain file (ubi file)
into the format needed by the GrainMapper program
"""


if __name__=="__main__":
    from ImageD11.grain import *
    import sys
    gl = read_grain_file( sys.argv[1] )
    f = open(sys.argv[2],"w")
    f.write("%d\n"%(len(gl)))
    for g in gl:
        u = g.ubi
        t = g.translation
        f.write("%f %f %f "%(u[0,0],u[0,1],u[0,2]))
        f.write("%f %f %f "%(u[1,0],u[1,1],u[1,2]))
        f.write("%f %f %f "%(u[2,0],u[2,1],u[2,2]))
        f.write("%f %f \n" %(t[0]/1000, t[1]/1000))
    
