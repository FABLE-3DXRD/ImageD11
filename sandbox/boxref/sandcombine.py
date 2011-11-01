 
import glob, os
def sand_special_case(fname):
    """
    Create an ascii file holding fltname, x, y, z in microns
    one file per line
    """
    rootdir = "/data/id11/inhouse/jon/sand0910/newpeaksearches/peaksearches"
    stem = "sand_097_0910"
    scans = "a_","b_"
    xr = range(7)
    zr = range(1,2)
    fltnames = [ ("%s%s%d_%d__all.flt"%(stem, s, xp, zp),xp,zp)
                 for s in scans
                 for zp in zr
                 for xp in xr ]
    #for f,x,z in fltnames:
    #    print f, x, z, os.path.exists(os.path.join(rootdir, f))
    assert len(fltnames) == len(glob.glob(os.path.join(rootdir,"*_all.flt")))
    f = open(fname,"w")
    y = 0
    for flt, x,z in fltnames:
        f.write("%s %f %f %f\n"%(
            os.path.join(rootdir, flt),
            x * 1500 - 4500,
            y,
            z * 500 - 750
            ))
    f.close()
    # print open(fname).read()
    print "Hardwired sand special case!"
    print "using default file with fltname samtx samty samtz"

if __name__=="__main__":
    from combine_flt import *
    sand_special_case( "sand.xyz")
    combine_flt_from_file( "sand.xyz", "sand.hdf", "sand0910" )

